#' Influential species detection - Phylogenetic Linear Regression
#'
#' Performs leave-one-out deletion analyis for phylogenetic linear regression,
#' and detects influential clades.
#'
#' @param formula The model formula
#' @param data Data frame containing species traits with row names matching tips
#' in \code{phy}.
#' @param phy A phylogeny (class 'phylo') matching \code{data}.
#' @param model The phylogenetic model to use (see Details). Default is \code{lambda}.
#' @param cutoff The cutoff value used to identify for influential clades
#' (see Details)
#' @param track Print a report tracking function progress (default = TRUE)
#' @param clade.col The name of a column in the provided data frame with clades 
#' specification.
#' @n.species Minimum required number of species in the clade in order to include
#' this clade in the leave-one-out deletion analyis. Default is \code{10}.
#' @param ... Further arguments to be passed to \code{phylolm}
#' @details
#' This function sequentially removes one clade at a time, fits a phylogenetic
#' linear regression model using \code{\link[phylolm]{phylolm}}, stores the
#' results and detects influential clades.
#'
#' All phylogenetic models from \code{phylolm} can be used, i.e. \code{BM},
#' \code{OUfixedRoot}, \code{OUrandomRoot}, \code{lambda}, \code{kappa},
#' \code{delta}, \code{EB} and \code{trend}. See ?\code{phylolm} for details.
#'
#' \code{clade_phylolm} detects influential clades based on the standardised
#' difference in intercept and/or slope when removing a given clade compared
#' to the full model including all species. Clades with a standardised difference
#' above the value of \code{cutoff} are identified as influential. The default
#' value for the cutoff is 2 standardised differences change.
#'
#' Currently, this function can only implement simple linear models (i.e. \eqn{trait~
#' predictor}). In the future we will implement more complex models.
#'
#' Output can be visualised using \code{sensi_plot}.
#'
#' @return The function \code{clade_phylolm} returns a list with the following
#' components:
#' @return \code{cutoff}: The value selected for \code{cutoff}
#' @return \code{formula}: The formula
#' @return \code{full.model.estimates}: Coefficients, aic and the optimised
#' value of the phylogenetic parameter (e.g. \code{lambda}) for the full model
#' without deleted species.
#' @return \code{influential_clades}: List of influential clades, both
#' based on standardised difference in interecept and in the slope of the
#' regression. Clades are ordered from most influential to less influential and
#' only include clades with a standardised difference > \code{cutoff}.
#' @return \code{clade.model.estimates}: A data frame with all simulation
#' estimates. Each row represents a deleted clade. Reported are the calculated
#' regression intercept (\code{intercept}), difference between simulation
#' intercept and full model intercept (\code{DFintercept}), the standardised
#' difference (\code{sDFintercept}), the percentage change in intercept compared
#' to the full model (\code{intercept.perc}) and intercept p-value
#' (\code{pval.intercept}). All of these are also reported for the regression
#' slope (\code{DFslope} etc.). Additonally, model aic value (\code{AIC}) and
#' the optimised value (\code{optpar}) of the phylogenetic parameter
#' (e.g. \code{kappa} or \code{lambda}, depends on phylogeneticmodel used) are
#' reported.
#' @return \code{data}: Original full dataset.
#' @return \code{errors}: Clades where deletion resulted in errors.
#' @examples
#' library(sensiPhy);library(phylolm)
#'
#' #Generate a random tree
#' set.seed(2468)
#' tree <- rtree(100)
#'
#' #Generate random predictor variable (pred), evolving according to a BM model.
#' pred<- rTraitCont(tree,root.value=0,sigma=1,model="BM")
#'
#' #Generate two continous traits, one evolving highly correlated with the
#' #predictor (trait 1), and one evolving more randomly (trait 2)
#' cont_trait1 <- pred + rTraitCont(tree,model="BM",sigma=0.1)
#' cont_trait2 <- pred + rTraitCont(tree,model="BM",sigma=10)
#'
#' #Generate two binary traits, one highly correlated to pred (trait 1), the other less.
#' bin_trait1<-rbinTrait(n=1,tree,beta=c(-1,0.5),alpha=0.1,
#'                      X=cbind(rep(1,length(tree$tip.label)),pred))
#' bin_trait2<-rbinTrait(n=1,tree,beta=c(-1,0.5),alpha=5,
#'                       X=cbind(rep(1,length(tree$tip.label)),pred))
#' dat<-data.frame(pred,cont_trait1,cont_trait2,bin_trait1,bin_trait2)
#'
#' #Determine influential species for both regressions.
#' fit1<-influ_phylolm(cont_trait1~pred,data = dat,phy = tree)
#' fit2<-influ_phylolm(cont_trait2~pred,data = dat,phy = tree)
#'
#' #For purposes of comparison the full model output from phylolm:
#' summary(phylolm(cont_trait1~pred,data = dat,phy = tree,model = "lambda"))
#' summary(phylolm(cont_trait2~pred,data = dat,phy = tree,model = "lambda"))
#' @author Gustavo Paterno & Gijsbert D.A. Werner
#' @seealso \code{\link[phylolm]{phylolm}}, \code{\link{samp_phylolm}},
#' \code{\link{influ_phyloglm}},\code{\link{sensi_plot}}
#' @references Here still: reference to phylolm paper + our own?
#' @export

clade_phylolm <- function(formula,data,phy,model="lambda",cutoff=2,track=TRUE,
                        clade, n.species = 10, ...){
    if(class(formula)!="formula") stop("formula must be class 'formula'")
    if(class(data)!="data.frame") stop("data must be class 'data.frame'")
    if(class(phy)!="phylo") stop("phy must be class 'phylo'")
    if(class(clade)!="character") stop("clade must be class 'character'")
    
    #Calculates the full model, extracts model parameters
    full.data <- data
    clade <- clade
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
    
    k <- names(which(table(full.data[,clade]) > n.species ))
    if (length(k) == 0) stop(paste("There is no clade with more than ",
                          n.species," species. Change 'n' to fix this problem",sep=""))
    # Loop:
    for (i in k){
        if (length(k) > 1) {
            crop.data <- full.data[full.data[ ,clade] %in% setdiff(k,i),]
            crop.sp <-   which(!full.data[ ,clade] %in% setdiff(k,i))
        }
        if (length(k) == 1) {
            crop.data <- full.data[full.data[ ,clade] %in% k,]
            crop.sp <-   which(!full.data[ ,clade] %in% k)
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
            intercept.perc       <- round((abs(DFintercept/intercept.0))*100,digits=1)
            slope.perc           <- round((abs(DFslope/slope.0))*100,digits=1)
            pval.intercept       <- phylolm::summary.phylolm(mod)$coefficients[[1,4]]
            pval.slope           <- phylolm::summary.phylolm(mod)$coefficients[[2,4]]
            aic.mod              <- mod$aic
            if (model == "BM"){
                optpar <- NA
            }
            if (model != "BM"){
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
    
    #Calculates Standardized DFbeta and DFintercept
    sDFintercept <- clade.model.estimates$DFintercept/
        sd(clade.model.estimates$DFintercept)
    sDFslope     <- clade.model.estimates$DFslope/
        sd(clade.model.estimates$DFslope)
    
    clade.model.estimates$sDFslope     <- sDFslope
    clade.model.estimates$sDFintercept <- sDFintercept
    
    #Creates a list with full model estimates:
    param0 <- list(coef=phylolm::summary.phylolm(mod.0)$coefficients,
                   aic=phylolm::summary.phylolm(mod.0)$aic,
                   optpar=mod.0$optpar)
    
    #Identifies influencital clade (sDF > cutoff) and orders by influence
    reorder.on.slope         <-clade.model.estimates[order(abs(
        clade.model.estimates$sDFslope),decreasing=T),c("clade","sDFslope")]
    influ.clade.slope           <-as.character(reorder.on.slope$clade[abs(
        reorder.on.slope$sDFslope)>cutoff])
    reorder.on.intercept     <-clade.model.estimates[order(abs(
        clade.model.estimates$sDFintercept),decreasing=T),c("clade","sDFintercept")]
    influ.clade.intercept       <-as.character(reorder.on.intercept$clade[abs(
        reorder.on.intercept$sDFintercept)>cutoff])
    
    #Generates output:
    res <- list(analysis.type="clade_phylolm",
                cutoff=cutoff,
                formula=formula,
                full.model.estimates=param0,
                influential.clades= list(influ.clade.slope=influ.clade.slope,
                                          influ.clade.intercept=influ.clade.intercept),
                clade.model.estimates=clade.model.estimates,
                data=full.data,errors=errors)
    
    ### Warnings:
    if (length(res$errors) >0){
        warning("Some clades deletion presented errors, please check: output$errors")}
    else {
        message("No errors found. All deletions were performed and stored successfully. Please, check outpu$clade.model.estimates.")
        res$errors <- "No errors found."
    }
    
    return(res)
    
}


