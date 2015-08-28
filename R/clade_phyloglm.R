#' Influential clade detection - Phylogenetic Logistic Regression
#'
#' Estimate the impact on model estimates phylogenetic Logist regression after 
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
#' \code{clade_phyloglm} detects influential clades based on
#' difference in intercept and/or slope when removing a given clade compared
#' to the full model including all species.
#' 
#' Currently, this function can only implement simple linear models (i.e. 
#' \eqn{y = a + bx}). In the future we will implement more complex models.
#'
#' Output can be visualised using \code{clade_plot}.
#'
#' @return The function \code{clade_phyloglm} returns a list with the following
#' components:
#' @return \code{formula}: The formula
#' @return \code{full.model.estimates}: Coefficients, aic and the optimised
#' value of the phylogenetic parameter (e.g. \code{alpha}) for the full model
#' without deleted species.
#' @return \code{clade.model.estimates}: A data frame with all simulation
#' estimates. Each row represents a deleted clade. Reported are the calculated
#' regression intercept (\code{intercept}), difference between simulation
#' intercept and full model intercept (\code{DFintercept}), the percentage change
#'  in intercept compared to the full model (\code{intercept.perc}) and intercept 
#'  p-value (\code{pval.intercept}). All of these are also reported for the regression
#' slope (\code{DFslope} etc.). Additonally, model aic value (\code{AIC}) and
#' the optimised value (\code{optpar}) of the phylogenetic parameter are
#' reported.
#' @return \code{data}: Original full dataset.
#' @return \code{errors}: Clades where deletion resulted in errors.
#' @author Gustavo Paterno & Gijsbert Werner
#' @seealso \code{\link[sensiPhy]{clade_phylolm}}, \code{\link{influ_phyloglm}},
#' \code{\link{sensi_plot}}
#' @references Here still: reference to phylolm paper + our own?
#' @export

clade_phyloglm <- function(formula, data, phy, btol=50, track = TRUE,
                          clade.col, n.species = 5, ...){
    
    if(class(formula)!="formula") stop("formula must be class 'formula'")
    if(!is.data.frame(data)) stop("data must be class 'data.frame'")
    if(class(phy)!="phylo") stop("phy must be class 'phylo'")
    if(missing(clade.col)) stop("clade.col not defined. Please, define the",
                                " column with clade names.")
    
    #Calculates the full model, extracts model parameters
    data_phy <- match_dataphy(formula, data, phy)
    #Calculates the full model, extracts model parameters
    full.data <- data_phy$data
    phy <- data_phy$phy
    full.data <- data
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
    
    k <- names(which(table(full.data[ ,clade.col]) > n.species ))
    if (length(k) == 0) stop(paste("There is no clade with more than ",
                                   n.species," species. Change 'n.species' 
                                   to fix this problem",sep = ""))
    # Loop:
    for (i in k){
        if (length(k) > 1) {
            crop.data <- full.data[full.data[ ,clade.col] %in% setdiff(k,i), ]
            crop.sp <-   which(!full.data[ ,clade.col] %in% setdiff(k,i))
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
            
            if(track==TRUE) (print(i))
            
            # Stores values for each simulation
            clade.model.estimates[counter,1]  <- i
            clade.model.estimates[counter,2]  <- intercept
            clade.model.estimates[counter,3]  <- DFintercept
            clade.model.estimates[counter,4]  <- intercept.perc
            clade.model.estimates[counter,5]  <- pval.intercept
            clade.model.estimates[counter,6]  <- slope
            clade.model.estimates[counter,7]  <- DFslope
            clade.model.estimates[counter,8]  <- slope.perc
            clade.model.estimates[counter,9]  <- pval.slope
            clade.model.estimates[counter,10] <- aic.mod
            clade.model.estimates[counter,11] <- optpar
            counter=counter+1
        }
    }
    
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
