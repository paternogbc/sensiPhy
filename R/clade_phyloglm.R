#' Influential clade detection - Phylogenetic Logistic Regression
#'
#' Estimate the impact on model estimates of phylogenetic logistic regression after 
#' removing clades from the analysis. 
#'
#' @param formula The model formula
#' @param data Data frame containing species traits with row names matching tips
#' in \code{phy}.
#' @param phy A phylogeny (class 'phylo') matching \code{data}.
#' @param btol Bound on searching space. For details see \code{phyloglm}
#' @param track Print a report tracking function progress (default = TRUE)
#' @param clade.col The name of a column in the provided data frame with clades 
#' specification (a character vector with clade names).
#' @param n.species Minimum number of species in the clade in order to include
#' this clade in the leave-one-out deletion analyis. Default is \code{5}.
#' @param ... Further arguments to be passed to \code{phyloglm}
#' 
#' @details
#' Currently only logistic regression using the "logistic_MPLE"-method from
#' \code{phyloglm} is implemented.
#'
#' \code{clade_phyglm} detects influential clades based on
#' difference in intercept and/or slope when removing a given clade compared
#' to the full model including all species.
#' 
#' Currently, this function can only implement simple linear models (i.e. 
#' \eqn{y = a + bx}). In the future we will implement more complex models.
#'
#' Output can be visualised using \code{clade_plot}.
#'
#' @return The function \code{clade_phyglm} returns a list with the following
#' components:
#' @return \code{formula}: The formula
#' @return \code{full.model.estimates}: Coefficients, aic and the optimised
#' value of the phylogenetic parameter (e.g. \code{alpha}) for the full model
#' without deleted species.
#' @return \code{clade.model.estimates}: A data frame with all simulation
#' estimates. Each row represents a deleted clade. Columns report the calculated
#' regression intercept (\code{intercept}), difference between simulation
#' intercept and full model intercept (\code{DFintercept}), the percentage of change
#' in intercept compared to the full model (\code{intercept.perc}) and intercept
#' p-value (\code{pval.intercept}). All these parameters are also reported for the regression
#' slope (\code{DFslope} etc.). Additionally, model aic value (\code{AIC}) and
#' the optimised value (\code{optpar}) of the phylogenetic parameter are
#' reported.
#' @return \code{data}: Original full dataset.
#' @return \code{errors}: Clades where deletion resulted in errors.
#' @author Gustavo Paterno & Gijsbert Werner
#' @seealso \code{\link[sensiPhy]{clade_phylm}}, \code{\link{influ_phyglm}},
#' \code{\link{sensi_plot}}
#' @references Here still: reference to phylolm paper + our own?
#' @examples 
#' \dontrun{
#'# Simulate Data:
#'set.seed(6987)
#'phy = rtree(150)
#'x = rTrait(n=1,phy=phy)
#'X = cbind(rep(1,150),x)
#'y = rbinTrait(n=1,phy=phy, beta=c(-1,0.5), alpha=.7 ,X=X)
#'cla <- rep(c("A","B","C","D","E"), each = 30)
#'dat = data.frame(y, x, cla)
#'# Run sensitivity analysis:
#'clade <- clade_phyglm(y ~ x, phy = phy, data = dat, clade.col = "cla")
#'# To check summary results and most influential clades:
#'summary(clade)
#'# Visual diagnostics for clade removal:
#'sensi_plot(clade)
#'# Specify which clade removal to plot:
#'sensi_plot(clade, "B")
#'sensi_plot(clade, "C")
#'}
#' @export

clade_phyglm <- function(formula, data, phy, btol=50, track = TRUE,
                           clade.col, n.species = 5, ...){
    
    if(class(formula)!="formula") stop("formula must be class 'formula'")
    if(!is.data.frame(data)) stop("data must be class 'data.frame'")
    if(class(phy)!="phylo") stop("phy must be class 'phylo'")
    if(missing(clade.col)) stop("clade.col not defined. Please, define the",
                                " column with clade names.")
    
    # Check match between data and phy 
    data_phy <- match_dataphy(formula, data, phy)
    #Calculates the full model, extracts model parameters
    full.data <- data_phy$data
    phy <- data_phy$phy
    namesInd <- match(clade.col, names(full.data))
    if (is.na(namesInd)) {
        stop("Names column '", clade.col, "' not found in data frame'")
    }
    
    N               <- nrow(full.data)
    mod.0           <- phylolm::phyloglm(formula, data = full.data, 
                                         phy = phy, method = "logistic_MPLE",
                                         btol = btol)
    intercept.0      <- mod.0$coefficients[[1]]
    slope.0          <- mod.0$coefficients[[2]]
    pval.intercept.0 <- phylolm::summary.phyloglm(mod.0)$coefficients[[1,4]]
    pval.slope.0     <- phylolm::summary.phyloglm(mod.0)$coefficients[[2,4]]
    optpar.0 <- mod.0$alpha
    
    
    if(isTRUE(mod.0$convergence!=0)) stop("Full model failed to converge,
                                          consider changing btol. See ?phyloglm")
    else
        
        #Creates empty data frame to store model outputs
        clade.model.estimates<-
        data.frame("clade" = I(as.character()), "intercept" = numeric(),
                   "DFintercept" = numeric(), "intercept.perc" = numeric(),
                   "pval.intercept" = numeric(), "slope" = numeric(),
                   "DFslope" = numeric(), "slope.perc" = numeric(),
                   "pval.slope" = numeric(), "AIC" = numeric(),
                   "optpar" = numeric())
    
    #Loops over all clades, and removes each one individually
    counter <- 1
    errors <- NULL
    
    all.clades <- levels(full.data[ ,clade.col])
    k <- names(which(table(full.data[ ,clade.col]) > n.species ))
    if (length(k) == 0) stop(paste("There is no clade with more than ",
                                   n.species," species. Change 'n.species' 
                                   to fix this problem",sep = ""))
    # Loop:
    pb <- txtProgressBar(min = 0, max = length(k), style = 1)
    for (i in k){
        if (length(k) > 1) {
            crop.data <- full.data[full.data[ ,clade.col] %in% setdiff(all.clades,i), ]
            crop.sp <-   which(!full.data[ ,clade.col] %in% setdiff(all.clades,i))
        }
        if (length(k) == 1) {
            crop.data <- full.data[!full.data[ ,clade.col] %in% k, ]
            crop.sp <-   which(full.data[ ,clade.col] %in% k)
        }
        
        crop.phy <-  ape::drop.tip(phy,phy$tip.label[crop.sp])
        mod=try(phylolm::phyloglm(formula, data = crop.data, 
                                  phy = crop.phy, method = "logistic_MPLE",
                                  btol = btol), TRUE)
        if(isTRUE(class(mod) == "try-error")) {
            
            error <- i
            errors <- c(errors,error)
            next }
        else {  
            
            intercept   <- mod$coefficients[[1]]
            slope       <- mod$coefficients[[2]]
            DFintercept <- intercept - intercept.0
            DFslope     <- slope - slope.0
            intercept.perc  <- round((abs(DFintercept / intercept.0)) * 100, digits = 1)
            slope.perc      <- round((abs(DFslope / slope.0)) * 100, digits = 1)
            pval.intercept  <- phylolm::summary.phyloglm(mod)$coefficients[[1,4]]
            pval.slope      <- phylolm::summary.phyloglm(mod)$coefficients[[2,4]]
            aic.mod         <- mod$aic
            optpar          <- mod$alpha
            
            if(track==TRUE) (setTxtProgressBar(pb, counter))
            
            # Stores values for each simulation
            estim.simu <- data.frame(i, intercept, DFintercept, intercept.perc,
                                     pval.intercept, slope, DFslope, slope.perc,
                                     pval.slope, aic.mod, optpar,
                                     stringsAsFactors = F)
            clade.model.estimates[counter, ]  <- estim.simu
            counter=counter+1
        }
    }
    on.exit(close(pb))
    #Creates a list with full model estimates:
    param0 <- list(coef=phylolm::summary.phylolm(mod.0)$coefficients,
                   aic=phylolm::summary.phylolm(mod.0)$aic,
                   optpar = optpar.0)
    
    #Generates output:
    res <- list(formula = formula,
                full.model.estimates = param0,
                clade.model.estimates = clade.model.estimates,
                data=full.data,
                errors = errors,
                clade.col = clade.col)
    class(res) <- c("sensiClade","sensiCladeL")
    ### Warnings:
    if (length(res$errors) > 0){
        warning("Some clades deletion presented errors, please check: output$errors")}
    else {
        message("No errors found. All deletions were performed and stored successfully. Please, check outpu$clade.model.estimates.")
        res$errors <- "No errors found."
    }
    return(res)
    
}
