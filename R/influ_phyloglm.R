#' Influential species detection - Phylogenetic Logistic Regression
#'
#' Performs leave-one-out deletion analyis for phylogenetic logistic regression,
#' and detects influential species.
#'
#' @param formula The model formula
#' @param data Data frame containing species traits with row names matching tips
#' in \code{phy}.
#' @param phy A phylogeny (class 'phylo') matching \code{data}.
#' @param btol Bound on searching space. For details see \code{phyloglm}
#' @param cutoff The cutoff value used to identify for influential species
#' (see Details)
#' @param track Print a report tracking function progress (default = TRUE)
#' @param ... Further arguments to be passed to \code{phylolm}
#' @details
#' This function sequentially removes one species at a time, fits a phylogenetic
#' logistic regression model using \code{\link[phylolm]{phyloglm}}, stores the
#' results and detects influential species.
#'
#' Currently only logistic regression using the "logistic_MPLE"-method from
#' \code{phyloglm} is implemented.
#'
#' \code{influ_phyglm} detects influential species based on the standardised
#' difference in intercept and/or slope when removing a given species compared
#' to the full model including all species. Species with a standardised difference
#' above the value of \code{cutoff} are identified as influential. The default
#' value for the cutoff is 2 standardised differences change.
#'
#' Currently, this function can only implement simple models (i.e. 
#' \eqn{y = a + bx}). In the future we will implement more complex models.
#'
#' Output can be visualised using \code{sensi_plot}.
#'
#' @return The function \code{influ_phyglm} returns a list with the following
#' components:
#' @return \code{cutoff}: The value selected for \code{cutoff}
#' @return \code{formula}: The formula
#' @return \code{full.model.estimates}: Coefficients, aic and the optimised
#' value of the phylogenetic parameter (i.e. \code{alpha}) for the full model
#' without deleted species.
#' @return \code{influential_species}: List of influential species, both
#' based on standardised difference in interecept and in the slope of the
#' regression. Species are ordered from most influential to less influential and
#' only include species with a standardised difference > \code{cutoff}.
#' @return \code{influ.model.estimates}: A data frame with all simulation
#' estimates. Each row represents a deleted clade. #' Columns report the calculated
#' regression intercept (\code{intercept}), difference between simulation
#' intercept and full model intercept (\code{DFintercept}), the standardised
#' difference (\code{sDFintercept}), the percentage of change in intercept compared
#' to the full model (\code{intercept.perc}) and intercept p-value
#' (\code{pval.intercept}). All these parameters are also reported for the regression
#' slope (\code{DFslope} etc.). Additionally, model aic value (\code{AIC}) and
#' the optimised value (\code{optpar}) of the phylogenetic parameter
#' (i.e. \code{alpha}) are reported.
#' @return \code{data}: Original full dataset.
#' @return \code{errors}: Species where deletion resulted in errors.
#' @author Gustavo Paterno & Gijsbert D.A. Werner
#' @seealso \code{\link[phylolm]{phyloglm}}, \code{\link{samp_phyglm}},
#' \code{\link{influ_phylm}}, \code{\link{sensi_plot}}
#' @references Here still: reference to phylolm paper + our own?
#' @examples 
#'\dontrun{
#'# Simulate Data:
#'set.seed(6987)
#'phy = rtree(150)
#'x = rTrait(n=1,phy=phy)
#'X = cbind(rep(1,150),x)
#'y = rbinTrait(n=1,phy=phy, beta=c(-1,0.5), alpha=.7 ,X=X)
#'dat = data.frame(y, x)
#'# Run sensitivity analysis:
#'influ <- influ_phyglm(y ~ x, data = dat, phy = phy) 
#'# To check summary results and most influential species:
#'summary(influ)
#'# Visual diagnostics for clade removal:
#'sensi_plot(influ)
#'}
#' @export

influ_phyglm <- function(formula,data,phy,btol=50,cutoff=2,track=TRUE,...){
        if(class(formula)!="formula") stop("formula must be class 'formula'")
        if(class(data)!="data.frame") stop("data must be class 'data.frame'")
        if(class(phy)!="phylo") stop("phy must be class 'phylo'")
        else
        
        # Check match between data and phy 
        data_phy <- match_dataphy(formula, data, phy)
        #Calculates the full model, extracts model parameters
        full.data <- data_phy$data
        phy <- data_phy$phy
        #Calculates the full model, extracts model parameters
        N               <- nrow(full.data)
        mod.0           <- phylolm::phyloglm(formula, data=full.data,
                                   phy=phy,method="logistic_MPLE",btol=btol,...)
        intercept.0      <- mod.0$coefficients[[1]]
        slope.0          <- mod.0$coefficients[[2]]
        pval.intercept.0 <- phylolm::summary.phyloglm(mod.0)$coefficients[[1,4]]
        pval.slope.0     <- phylolm::summary.phyloglm(mod.0)$coefficients[[2,4]]
        optpar.0         <- mod.0$alpha
        if(isTRUE(mod.0$convergence!=0)) stop("Full model failed to converge,
                                        consider changing btol. See ?phyloglm")
        else

        #Creates empty data frame to store model outputs
        influ.model.estimates<-
                data.frame("species" =numeric(), "intercept"=numeric(),
                           "DFintercept"=numeric(),"intercept.perc"=numeric(),
                           "pval.intercept"=numeric(),"slope"=numeric(),
                           "DFslope"=numeric(),"slope.perc"=numeric(),
                           "pval.slope"=numeric(),"AIC"=numeric(),
                           "optpar" = numeric())

        #Loops over all species, and removes each one individually
        counter <- 1
        errors <- NULL
        pb <- txtProgressBar(min = 0, max = N, style = 1)
        for (i in 1:N){
                crop.data <- full.data[c(1:N)[-i],]
                crop.phy <-  ape::drop.tip(phy,phy$tip.label[i])
                mod=try(phylolm::phyloglm(formula, data=crop.data,phy=crop.phy,
                                          method="logistic_MPLE",btol=btol,...),
                        TRUE)
                if(isTRUE(class(mod)=="try-error")) {
                        error <- i
                        names(error) <- rownames(full.data$data)[i]
                        errors <- c(errors,error)
                        next }
                else {  sp                   <- phy$tip.label[i]
                        intercept            <- mod$coefficients[[1]]
                        slope                <- mod$coefficients[[2]]
                        DFintercept          <- intercept - intercept.0
                        DFslope              <- slope - slope.0
                        intercept.perc       <- round((abs(DFintercept/intercept.0))*100,digits=1)
                        slope.perc           <- round((abs(DFslope/slope.0))*100,digits=1)
                        pval.intercept       <- phylolm::summary.phyloglm(mod)$coefficients[[1,4]]
                        pval.slope           <- phylolm::summary.phyloglm(mod)$coefficients[[2,4]]
                        aic.mod              <- mod$aic
                        optpar               <- mod$alpha

                        if(track==TRUE) setTxtProgressBar(pb, i)

                        #Stores values for eacht simulation
                        estim.simu <- data.frame(sp, intercept, DFintercept, intercept.perc,
                                                 pval.intercept, slope, DFslope, slope.perc,
                                                 pval.slope, aic.mod, optpar,
                                                 stringsAsFactors = F)
                        influ.model.estimates[counter, ]  <- estim.simu
                        counter=counter+1
                }
        }
        on.exit(close(pb))
        #Calculates standardized DFbeta and DFintercept
        sDFintercept <- influ.model.estimates$DFintercept/
                stats::sd(influ.model.estimates$DFintercept)
        sDFslope     <- influ.model.estimates$DFslope/
                stats::sd(influ.model.estimates$DFslope)

        influ.model.estimates$sDFslope     <- sDFslope
        influ.model.estimates$sDFintercept <- sDFintercept

        #Creates a list with full model estimates:
        param0 <- list(coef=phylolm::summary.phyloglm(mod.0)$coefficients,
                       aic=phylolm::summary.phyloglm(mod.0)$aic,
                       optpar=phylolm::summary.phyloglm(mod.0)$alpha)

        #Identifies influencital species (sDF > cutoff) and orders by influence
        reorder.on.slope         <-influ.model.estimates[order(abs(
                influ.model.estimates$sDFslope),decreasing=T),c("species","sDFslope")]
        influ.sp.slope           <-as.character(reorder.on.slope$species[abs(
                reorder.on.slope$sDFslope)>cutoff])
        reorder.on.intercept     <-influ.model.estimates[order(abs(
                influ.model.estimates$sDFintercept),decreasing=T),c("species","sDFintercept")]
        influ.sp.intercept       <-as.character(reorder.on.intercept$species[abs(
                reorder.on.intercept$sDFintercept)>cutoff])

        #Generates output:
        res <- list(cutoff=cutoff,
                    formula=formula,
                    full.model.estimates=param0,
                    influential.species= list(influ.sp.slope=influ.sp.slope,
                                              influ.sp.intercept=influ.sp.intercept),
                    influ.model.estimates=influ.model.estimates,
                    data=full.data,errors=errors)
        class(res) <- "sensiInflu"
        ### Warnings:
        if (length(res$errors) >0){
                warning("Some species deletion presented errors, please check: output$errors")}
        else {
                message("No errors found. All single deletions were performed and stored successfully. Please, check outpu$influ.model.estimates.")
                res$errors <- "No errors found."
        }

        return(res)

}



