#' Influential clade detection - Phylogenetic Linear Regression
#'
#' Estimate the impact on model estimates phylogenetic linear regression after 
#' removing clades from the analysis. 
#'
#' @param formula The model formula
#' @param data Data frame containing species traits with row names matching tips
#' in \code{phy}.
#' @param phy A phylogeny (class 'phylo') matching \code{data}.
#' @param model The phylogenetic model to use (see Details). Default is \code{lambda}.
#' @param track Print a report tracking function progress (default = TRUE)
#' @param clade.col The name of a column in the provided data frame with clades 
#' specification (a character vector with clade names).
#' @param n.species Minimum number of species in the clade in order to include
#' this clade in the leave-one-out deletion analyis. Default is \code{5}.
#' @param ... Further arguments to be passed to \code{phylolm}
#' @details
#' This function sequentially removes one clade at a time, fits a phylogenetic
#' linear regression model using \code{\link[phylolm]{phylolm}} and stores the
#' results. The impact of of a specific clade on model estimates is calculated by the
#'  comparison between the full model (with all species) and the model without 
#'  the species belonging to a clade.
#'
#' All phylogenetic models from \code{phylolm} can be used, i.e. \code{BM},
#' \code{OUfixedRoot}, \code{OUrandomRoot}, \code{lambda}, \code{kappa},
#' \code{delta}, \code{EB} and \code{trend}. See ?\code{phylolm} for details.
#'
#' \code{clade_phylolm} detects influential clades based on
#' difference in intercept and/or slope when removing a given clade compared
#' to the full model including all species.
#' 
#' Currently, this function can only implement simple linear models (i.e. 
#' \eqn{y = a + bx}). In the future we will implement more complex models.
#'
#' Output can be visualised using \code{sensi_plot}.
#'
#' @return The function \code{clade_phylolm} returns a list with the following
#' components:
#' @return \code{formula}: The formula
#' @return \code{full.model.estimates}: Coefficients, aic and the optimised
#' value of the phylogenetic parameter (e.g. \code{lambda}) for the full model
#' without deleted species.
#' @return \code{clade.model.estimates}: A data frame with all simulation
#' estimates. Each row represents a deleted clade. Reported are the calculated
#' regression intercept (\code{intercept}), difference between simulation
#' intercept and full model intercept (\code{DFintercept}), the percentage change
#'  in intercept compared to the full model (\code{intercept.perc}) and intercept 
#'  p-value (\code{pval.intercept}). All of these are also reported for the regression
#' slope (\code{DFslope} etc.). Additonally, model aic value (\code{AIC}) and
#' the optimised value (\code{optpar}) of the phylogenetic parameter
#' (e.g. \code{kappa} or \code{lambda}, depends on phylogeneticmodel used) are
#' reported.
#' @return \code{data}: Original full dataset.
#' @return \code{errors}: Clades where deletion resulted in errors.
#' @examples
#' \dontrun{
#' library(sensiPhy);library(phylolm)
#'
#' #Generate a random tree
#' set.seed(2468)
#' tree <- rtree(100)
#'
#' #Generate random predictor variable (pred), evolving according to a BM model.
#' pred<- rTraitCont(tree,root.value=0,sigma=1,model="BM")
#'
#' #Generate one continous traits
#' cont_trait <- pred + rTraitCont(tree,model="BM",sigma=2.5)
#' fam <- rep(c("fam1","fam2","fam3","fam4","fam5"),each=20)
#' dat<-data.frame(pred,cont_trait,fam)
#'
#' #Determine influential clades:
#' clade.test <- clade_phylolm(cont_trait~pred,data=dat,phy=tree,clade.col="fam")
#' }
#' @author Gustavo Paterno
#' @seealso \code{\link[phylolm]{phylolm}}, \code{\link{influ_phylolm}},
#' \code{\link{sensi_plot}}
#' @references Here still: reference to phylolm paper + our own?
#' @export

clade_phylolm <- function(formula, data, phy, model = "lambda", track = TRUE,
                        clade.col, n.species = 5, ...){
    if(!is.data.frame(data)) stop("data must be class 'data.frame'")
    if(missing(clade.col)) stop("clade.col not defined. Please, define the",
                                " column with clade names.")
    data_phy <- match_dataphy(formula, data, phy)
    #Calculates the full model, extracts model parameters
    full.data <- data_phy$data
    phy <- data_phy$phy
    clade.col <- clade.col
    namesInd <- match(clade.col, names(full.data))
    if (is.na(namesInd)) {
        stop("Names column '", clade.col, "' not found in data frame'")
    }
    
    N               <- nrow(full.data)
    mod.0           <- phylolm::phylolm(formula, data=full.data,
                                        model=model,phy=phy)
    intercept.0      <- mod.0$coefficients[[1]]
    slope.0          <- mod.0$coefficients[[2]]
    pval.intercept.0 <- phylolm::summary.phylolm(mod.0)$coefficients[[1,4]]
    pval.slope.0     <- phylolm::summary.phylolm(mod.0)$coefficients[[2,4]]
    optpar.0 <- mod.0$optpar

    
    #Creates empty data frame to store model outputs
    clade.model.estimates<-
        data.frame("clade" =I(as.character()), "intercept"=numeric(),
                   "DFintercept"=numeric(),"intercept.perc"=numeric(),
                   "pval.intercept"=numeric(),"slope"=numeric(),
                   "DFslope"=numeric(),"slope.perc"=numeric(),
                   "pval.slope"=numeric(),"AIC"=numeric(),
                   "optpar" = numeric())
    
    #Loops over all clades, and removes each one individually
    counter <- 1
    errors <- NULL
    
    k <- names(which(table(full.data[,clade.col]) > n.species ))
    if (length(k) == 0) stop(paste("There is no clade with more than ",
                          n.species," species. Change 'n.species' to fix this
                          problem",sep=""))
    # Loop:

    for (i in k){
        if (length(k) > 1) {
            crop.data <- full.data[full.data[ ,clade.col] %in% setdiff(k,i),]
            crop.sp <-   which(!full.data[ ,clade.col] %in% setdiff(k,i))
        }
        if (length(k) == 1) {
            crop.data <- full.data[!full.data[ ,clade.col] %in% k,]
            crop.sp <-   which(full.data[ ,clade.col] %in% k)
        }
        
        crop.phy <-  ape::drop.tip(phy,phy$tip.label[crop.sp])
        mod=try(phylolm::phylolm(formula, data=crop.data,model=model,
                                 phy=crop.phy),TRUE)
        if(isTRUE(class(mod)=="try-error")) {
            
            error <- i
            errors <- c(errors,error)
            next }
        else {  
            
            intercept            <- mod$coefficients[[1]]
            slope                <- mod$coefficients[[2]]
            DFintercept          <- intercept - intercept.0
            DFslope              <- slope - slope.0
            intercept.perc       <- round((abs(DFintercept / intercept.0)) * 100,
                                          digits = 1)
            slope.perc           <- round((abs(DFslope / slope.0)) * 100,
                                          digits = 1)
            pval.intercept       <- phylolm::summary.phylolm(mod)$coefficients[[1,4]]
            pval.slope           <- phylolm::summary.phylolm(mod)$coefficients[[2,4]]
            aic.mod              <- mod$aic
            if (model == "BM" | model == "trend"){
                optpar <- NA
            }
            if (model != "BM" & model != "trend" ){
                optpar               <- mod$optpar
            }
            
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
                   optpar=mod.0$optpar)

    #Generates output:
    res <- list(model = model,
                formula = formula,
                full.model.estimates = param0,
                clade.model.estimates = clade.model.estimates,
                data = full.data,
                errors = errors,
                clade.col = clade.col)
    class(res) <- "sensiClade"
    ### Warnings:
    if (length(res$errors) >0){
        warning("Some clades deletion presented errors, please check: output$errors")}
    else {
        message("No errors found. All deletions were performed and stored successfully. Please, check outpu$clade.model.estimates.")
        res$errors <- "No errors found."
    }
    return(res)
    
}


