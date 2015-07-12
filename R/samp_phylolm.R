#' Sensitivity Analysis Species Sampling  - Phylogenetic Linear Regression
#'
#' Performs analyses of sensitivity to species sampling by randomly removing
#' species and detecting the effects on parameter estimates in a phylogenetic
#' linear regression.
#'
#' @param formula The model formula
#' @param data Data frame containing species traits with row names matching tips
#' in \code{phy}.
#' @param phy A phylogeny (class 'phylo') matching \code{data}.
#' @param times The number of times species are randomly deleted for each
#' \code{break}.
#' @param breaks A vector containing the percentages of species to remove.
#' @param model The phylogenetic model to use (see Details). Default is \code{lambda}.
#' @param track Print a report tracking function progress (default = TRUE)
#' @param ... Further arguments to be passed to \code{phylolm}
#' @details
#'
#' This function randomly removes a given percentage of species (controlled by
#' \code{breaks}) from the full phylogenetic linear regression, fits a phylogenetic
#' linear regression model without these species using \code{\link[phylolm]{phylolm}},
#' repeates this many times (controlled by \code{times}), stores the results and
#' calculates the effects on model parameters.
#'
#' All phylogenetic models from \code{phylolm} can be used, i.e. \code{BM},
#' \code{OUfixedRoot}, \code{OUrandomRoot}, \code{lambda}, \code{kappa},
#' \code{delta}, \code{EB} and \code{trend}. See ?\code{phylolm} for details.
#'
#' Currently, this function can only implement simple linear models (i.e. \eqn{trait~
#' predictor}). In the future we will implement more complex models.
#'
#' Output can be visualised using \code{sensi_plot}.
#' @return The function \code{samp_phylolm} returns a list with the following
#' components:
#' @return \code{formula}: The formula
#' @return \code{full.model.estimates}: Coefficients, aic and the optimised
#' value of the phylogenetic parameter (e.g. \code{lambda} or \code{kappa}) for
#' the full model without deleted species.
#' @return \code{samp.model.estimates}: A data frame with all simulation
#' estimates. Each row represents a model rerun with a given number of species
#' \code{n.remov} removed, representing \code{n.percent} of the full dataset.
#' Reported are the calculated regression intercept (\code{intercept}),
#' difference between simulation intercept and full model intercept (\code{DFintercept}),
#' the percentage change in intercept compared to the full model (\code{intercept.perc})
#' and intercept p-value (\code{pval.intercept}). All of these are also reported
#' for the regression slope (\code{DFslope} etc.). Additonally, model aic value
#' (\code{AIC}) and the optimised value (\code{optpar}) of the phylogenetic
#' parameter (e.g. \code{kappa} or \code{lambda}, depends on phylogeneticmodel
#' used) are reported.
#' @return \code{sign.analysis} For each break (i.e. each percentage of species
#' removed) this reports the percentage of statistically signficant (at p<0.05)
#' intercepts (\code{perc.sign.intercept}) over all repititions as well as the
#' percentage of statisticaly significant (at p<0.05) slopes (\code{perc.sign.slope}).
#' @return \code{data}: Original full dataset.
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
#' #For both regressions, determine sensitivity of results to species sampled.
#' fit1<-samp_phylolm(cont_trait1~pred,data = dat,phy = tree)
#' fit2<-samp_phylolm(cont_trait2~pred,data = dat,phy = tree)
#'
#' #It is possible to change the species removal percentages and number of repeats
#' fit3<-samp_phylolm(cont_trait2~pred,data = dat,phy = tree,
#'      breaks = c(0.25,0.5,0.75),times=100)
#'
#' @author Gustavo Paterno & Gijsbert D.A. Werner
#' @seealso \code{\link[phylolm]{phylolm}}, \code{\link{samp_phyloglm}},
#' \code{\link{influ_phylolm}},\code{\link{sensi_plot}}
#' @references Ho, L. S. T. and Ané, C. 2014. "A linear-time algorithm for 
#' Gaussian and non-Gaussian trait evolution models". Systematic Biology 63(3):397-408.
#' @export

samp_phylolm <- function(formula,data,phy,times=20,
                         breaks=seq(.1,.7,.1),model="lambda",track=TRUE,...)
        {if(class(formula)!="formula") stop("formula must be class 'formula'")
         if(class(data)!="data.frame") stop("data must be class 'data.frame'")
         if(class(phy)!="phylo") stop("phy must be class 'phylo'")
         if(length(breaks)<2) stop("Please include more than one break,
                                  e.g. breaks=c(.3,.5)")
         if ((model == "trend") & (is.ultrametric(phy)))
                 stop("Trend is unidentifiable for ultrametric trees.,
                      see ?phylolm for details")
         else

        #Calculates the full model, extracts model parameters
        full.data       <- data
        N               <- nrow(full.data)
        mod.0           <- phylolm::phylolm(formula, data=full.data,
                                             model=model,phy=phy)
        intercept.0      <- mod.0$coefficients[[1]]
        slope.0          <- mod.0$coefficients[[2]]
        pval.intercept.0 <- phylolm::summary.phylolm(mod.0)$coefficients[[1,4]]
        pval.slope.0     <- phylolm::summary.phylolm(mod.0)$coefficients[[2,4]]
        optpar.0         <- mod.0$optpar
        aic.0                   <- mod.0$aic

        #Creates empty data frame to store model outputs
        samp.model.estimates<-
                data.frame("n.remov" =numeric(), "n.percent"=numeric(),
                           "intercept"=numeric(),"DFintercept"=numeric(),
                           "intercept.perc"=numeric(),"pval.intercept"=numeric(),
                           "slope"=numeric(),"DFslope"=numeric(),
                           "slope.perc"=numeric(),"pval.slope"=numeric(),
                           "AIC"=numeric(),"optpar" = numeric())

        #Loops over breaks, remove percentage of species determined by 'breaks
        #and repeat determined by 'times'.
        counter=1
        limit <- sort(round((breaks)*nrow(full.data),digits=0))
        for (i in limit){
                for (j in 1:times){
                        exclude <- sample(1:N,i)
                        crop.data <- full.data[-exclude,]
                        crop.phy <-  ape::drop.tip(phy,phy$tip.label[exclude])
                        mod=try(phylolm::phylolm(formula, data=crop.data,
                                                 model=model,phy=crop.phy),TRUE)
                        if(isTRUE(class(mod)=="try-error")) { next }
                        else {  intercept             <- mod$coefficients[[1]]
                                slope                 <- mod$coefficients[[2]]
                                optpar                <- mod$optpar
                                pval.intercept        <-
                                        phylolm::summary.phylolm(mod)$coefficients[[1,4]]
                                pval.slope            <-
                                        phylolm::summary.phylolm(mod)$coefficients[[2,4]]
                                aic                   <- mod$aic
                                DFintercept           <- intercept - intercept.0
                                DFslope               <- slope - slope.0
                                intercept.perc        <- round((abs(
                                        DFintercept/intercept.0))*100,digits=1)
                                slope.perc            <- round((abs(
                                        DFslope/slope.0))*100,digits=1)
                                aic                   <- mod$aic
                                optpar                <- mod$optpar
                                n.remov <- i
                                n.percent <- round((n.remov/N)*100,digits=0)
                                rep <- j

                                if(track==TRUE) (
                                print(paste("Break = ",n.percent,". Repetition = ",j,sep="")))

                                # Stores values for each simulation
                                samp.model.estimates[counter,1]<- n.remov
                                samp.model.estimates[counter,2]<- n.percent
                                samp.model.estimates[counter,3]<- intercept
                                samp.model.estimates[counter,4]<- DFintercept
                                samp.model.estimates[counter,5]<- intercept.perc
                                samp.model.estimates[counter,6]<- pval.intercept
                                samp.model.estimates[counter,7]<- slope
                                samp.model.estimates[counter,8]<- DFslope
                                samp.model.estimates[counter,9]<- slope.perc
                                samp.model.estimates[counter,10]<- pval.slope
                                samp.model.estimates[counter,11]<- aic
                                samp.model.estimates[counter,12]<- optpar
                                counter=counter+1
                        }
                }
        }

        #Calculates percentages of signficant intercepts & slopes within breaks.
        res                     <- samp.model.estimates
        times                   <- table(res$n.remov)
        breaks                  <- unique(res$n.percent)
        sign.intercept          <- res$pval.intercept > .05
        sign.slope              <- res$pval.slope > .05
        res$sign.intercept      <- sign.intercept
        res$sign.slope          <- sign.slope
        perc.sign.intercept     <- 1-(with(res,tapply(sign.intercept,n.remov,sum)))/times
        perc.sign.slope         <- 1-(with(res,tapply(sign.slope,n.remov,sum)))/times
        perc.sign.tab       <- data.frame(percent_sp_removed=breaks,
                                          perc.sign.intercept=as.numeric(perc.sign.intercept),
                                          perc.sign.slope=as.numeric(perc.sign.slope))

        #Creates a list with full model estimates:
        param0<-list(coef=phylolm::summary.phylolm(mod.0)$coefficients,
                     aic=phylolm::summary.phylolm(mod.0)$aic,
                     optpar=mod.0$optpar)

        #Generates output:
        return(list(analyis.type = "samp_phylolm",
                    formula=formula,
                    full.model.estimates=param0,
                    samp.model.estimates=samp.model.estimates,
                    sign.analysis=perc.sign.tab,
                    data=full.data))

}
