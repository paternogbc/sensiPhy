#' Leave-one-out-deletion analysis for gls phylogenetic regression.
#'
#' \code{influ_pgls} performs leave-one-out-deletion analysis for
#' \code{\link{gls}} phylogenetic linear regression. It removes one species at a time, fits
#' a gls model without the species and store the results of the model
#' estimates. The function repeats this procedure for all species in the dataset
#' and store all the results in a data.frame.
#' @aliases influ_pgls
#' @param formula A model formula
#' @param data A data frame containing variables that can be attributed to the taxa
#' at the tips of a phylogeny
#' @param phy A phylogeny (class 'phylo') to be matched to the data above
#' @details This functions only works for simple linear regression \eqn{y = bx +
#'   a}. Future implementation will deal with more complex models. This functions only
#'   works with Pagel correlation (see \code{corPagel}). If any error occur during simulation,
#'   check 'output$errors' for details.
#' @return The function \code{influ_gls} returns a list with the following
#'   components:
#'
#' @return \code{formula} The model formula
#' @return \code{model_estimates} Full model estimates
#' @return \code{influential_species} Most influential species for beta
#' @return \code{results} A data frame with all simulation estimates. DFbeta and
#'   DFintercept represent absolute difference between full model and simulation.
#'   sDFintercept and sDFbeta represent the standardized differences.
#' @return \code{data} Original dataset
#' @return \code{errors} Species that showed erros during pgls fit
#' @section Warning: This code is note fully checked. Please be aware.
#' @seealso \code{\link[caper]{pgls}}, \code{\link{samp_pgls}}
#' @export

influ_phylolm <- function(formula,data,phy,model="lambda",cutoff=2,...)
{
        # Basic error checking:
        if(class(formula)!="formula") stop("Please formula must be class
                                           'forumla'")
        if(class(data)!="data.frame") stop("Please data must be class
                                           'data.frame'")
        if(class(phy)!="phylo") stop("Please phy must be class
                                           'phylo'")
        if (sum(rownames(data) != phy$tip.label) > 0) stop("Species must be at the same order
                                                      in data and phy")
        if ((model == "trend") & (is.ultrametric(phy)))
                stop("the trend is unidentifiable for ultrametric trees.")
        else

        # FULL MODEL calculations:
        full.data <- data
        N         <- nrow(full.data)
        mod.0     <- phylolm::phylolm(formula, data=full.data,model=model,phy=phy)

        intercept.0      <- mod.0$coefficients[[1]]             # Intercept (full model)
        slope.0          <- mod.0$coefficients[[2]]             # Beta (full model)
        pval.intercept.0 <- phylolm::summary.phylolm(mod.0)$coefficients[[1,4]] # p.value (intercept)
        pval.slope.0     <-     phylolm::summary.phylolm(mod.0)$coefficients[[2,4]] # p.value (slope)
        optpar.0         <- mod.0$optpar

        #Create the influ.model.estimates data.frame
        influ.model.estimates<-data.frame("species" =numeric(), "intercept"=numeric(),
                            "DFintercept"=numeric(),"intercept.perc"=numeric(),"pval.intercept"=numeric(),
                            "slope"=numeric(),"DFslope"=numeric(),"slope.perc"=numeric(),
                            "pval.slope"=numeric(),"AIC"=numeric(),
                             "optpar" = numeric())
        #Loop:
        counter <- 1
        errors <- NULL

        for (i in 1:N){
                crop.data <- full.data[c(1:N)[-i],]
                crop.phy <-  ape::drop.tip(phy,phy$tip.label[i])

                mod=try(phylolm::phylolm(formula, data=crop.data,model=model,phy=crop.phy),TRUE)

                if(isTRUE(class(mod)=="try-error")) {
                        error <- i
                        names(error) <- rownames(full.data$data)[i]
                        errors <- c(errors,error)
                        next }

                else {
                        ### Calculating model estimates:

                        sp                   <- phy$tip.label[i]      # species removed
                        intercept            <- mod$coefficients[[1]] # Intercept (crop model)
                        slope                <- mod$coefficients[[2]] # Beta (crop model)
                        DFintercept          <- intercept - intercept.0 # DF intercept
                        DFslope              <- slope - slope.0 # DF beta
                        intercept.perc       <- round((abs(DFintercept/intercept.0))*100,digits=1)  # Percentage of intercept change
                        slope.perc           <- round((abs(DFslope/slope.0))*100,digits=1)  # Percentage of beta change
                        pval.intercept       <- phylolm::summary.phylolm(mod)$coefficients[[1,4]] # p.value (intercept)
                        pval.slope           <- phylolm::summary.phylolm(mod)$coefficients[[2,4]] # p.value
                        aic.mod              <- mod$aic # Model AIC
                        optpar               <- mod$optpar# Estimated lambda
                        print(i)

                        ### Storing values for each simulation
                        influ.model.estimates[counter,1]  <- sp
                        influ.model.estimates[counter,2]  <- intercept
                        influ.model.estimates[counter,3]  <- DFintercept
                        influ.model.estimates[counter,4]  <- intercept.perc
                        influ.model.estimates[counter,5]  <- pval.intercept
                        influ.model.estimates[counter,6]  <- slope
                        influ.model.estimates[counter,7]  <- DFslope
                        influ.model.estimates[counter,8]  <- slope.perc
                        influ.model.estimates[counter,9]  <- pval.slope
                        influ.model.estimates[counter,10] <- aic.mod
                        influ.model.estimates[counter,11] <- optpar
                        counter=counter+1
                }
        }

        ### Calculating Standardized DFbeta and DFintercept
        sDFintercept <- influ.model.estimates$DFintercept/sd(influ.model.estimates$DFintercept)
        sDFslope     <- influ.model.estimates$DFslope/sd(influ.model.estimates$DFslope)


        influ.model.estimates$sDFslope     <- sDFslope;
        influ.model.estimates$sDFintercept <- sDFintercept

        ### Original model estimates:
        param0 <- list(coef=phylolm::summary.phylolm(mod.0)$coefficients,
                       aic=phylolm::summary.phylolm(mod.0)$aic,
                       optpar=phylolm::summary.phylolm(mod.0)$optpar)


        ### Statistically Influential species for Beta (sDFbetas > 2)
        influ.sp.b <- as.character(influ.model.estimates$species[which(abs(influ.model.estimates$sDFslope) > cutoff)])

        ### Output:
        res <- list(analysis.type="influ_phylolm",
             formula=formula,
             full.model.estimates=param0,
             influential.species= influ.sp.b,
             influ.model.estimates=influ.model.estimates,
             data=full.data,errors=errors)

        ### Warnings:
        if (length(res$errors) >0){
                warning("Some species deletion presented errors, please check: output$errors")}
        else {
                message("No errors found. All single deletions were performed and stored successfully. Please, check outpu$influ.model.estimates.")
                res$errors <- "No errors found."
        }

        return(res)

}
