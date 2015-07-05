#' Sampling effort analysis for pgls linear regression.
#'
#' \code{samp_pgls} performs sample size sensitive analysis for \code{pgls}
#' linear regressions. It removes species at random, fits a pgls model without the
#' species and store the results of the model estimates. The percentage of
#' species removed is specified with \code{breaks} and the number of simulations
#' per break is defined by \code{times}.
#' @aliases samp_pgls
#' @inheritParams influ_pgls
#' @param breaks Percentage intervals to remove species. For example:
#'   \code{breaks = c(.1,.2,.3)},removes 10,20 and 30 percentage of species at
#'   random in each simulation.
#' @param times The number of times to repeat each simulation (per
#'   \code{breaks}) interval.
#' @details This functions only works for simple linear regression \eqn{y = bx +
#'   a}. Future implementation will deal with more complex models.
#' @return The function \code{samp_gls} returns a list with the following
#'   components:
#' @return \code{model_estimates} Full model estimates
#' @return \code{results} A data frame with all simulation estimates.
#' @return \code{power_analysis} A data frame with power analysis
#' @return \code{data} Original dataset
#' @seealso \code{\link[caper]{pgls}}, \code{\link{influ_pgls}}
#' @examples
#' \dontrun{
#' library(caper);library(sensiC)
#' data(shorebird)
#' comp.data <- comparative.data(shorebird.tree, shorebird.data, Species, vcv=TRUE, vcv.dim=3)
# # First we need to match tip.labels with rownames in data:
#' sp.ord <- match(shorebird.tree$tip.label, rownames(shorebird.data))
#' shorebird.data <- shorebird.data[sp.ord,]
#' samp1 <- samp_pgls(log(Egg.Mass) ~ log(M.Mass),data=shorebird.data,phy=shorebird.tree)
#' }
#' @export

samp_phylolm <- function(formula,data,phy,times=20,breaks=seq(.1,.7,.1),model="lambda",
                      ...)
{
        ### Basic error checking:
        if(class(formula)!="formula") stop("Please formula must be
                                           class 'formula'")
        if(class(data)!="data.frame") stop("data data must be of class
                                                 'data.frame'. See function `comparative.data`.")
        if(length(breaks)<2) stop("please include more then one break
                                  (eg. breaks=c(.3,.5)")
        if(class(phy) != "phylo") stop("Please phy must be of class 'phylo'")

        if (sum(rownames(data) != phy$tip.label) > 0) stop("Species must be at the same order
                                                      in data and phy")
        if ((model == "trend") & (is.ultrametric(phy)))
                stop("the trend is unidentifiable for ultrametric trees.")
        else

        # FULL MODEL calculations:
        full.data <- data
        N <- nrow(full.data)
        mod.0 <- phylolm::phylolm(formula, data=full.data,model=model,phy=phy)

        intercept.0             <- mod.0$coefficients[[1]] #Intercept (full model)
        slope.0                 <- mod.0$coefficients[[2]] #Slope (full model)
        optpar.0                <- mod.0$alpha             #The optimisation paratemer alpha (phylogenetic correlation parameter)
        pval.intercept.0        <- phylolm::summary.phylolm(mod.0)$coefficients[[1,4]] #P-value intercept (full model)
        pval.slope.0            <- phylolm::summary.phylolm(mod.0)$coefficients[[2,4]] #P-value slope (full model)
        aic.0                   <- mod.0$aic

        #Create the results data.frame
        samp.model.estimates<-data.frame("n.remov" =numeric(), "n.percent"=numeric(),
                                         "intercept"=numeric(),"DFintercept"=numeric(),"intercept.perc"=numeric(),"pval.intercept"=numeric(),
                                         "slope"=numeric(),"DFslope"=numeric(),"slope.perc"=numeric(),"pval.slope"=numeric(),
                                         "AIC"=numeric(),"optpar" = numeric())

        # Loop:
        counter=1
        limit <- sort(round((breaks)*nrow(full.data),digits=0))
        for (i in limit){
                for (j in 1:times){
                        exclude <- sample(1:N,i)
                        crop.data <- full.data[-exclude,]
                        crop.phy <-  ape::drop.tip(phy,phy$tip.label[exclude])

                        mod=try(phylolm::phylolm(formula, data=crop.data,model=model,phy=crop.phy),TRUE)
                        if(isTRUE(class(mod)=="try-error")) { next }
                        else {
                                intercept             <- mod$coefficients[[1]] #Intercept (crop model)
                                slope                 <- mod$coefficients[[2]] #Slope (crop model)
                                optpar                <- mod$optpar            #The optimisation paratemer alpha (phylogenetic correlation parameter)
                                pval.intercept        <- phylolm::summary.phylolm(mod)$coefficients[[1,4]] #P-value intercept (crop model)
                                pval.slope            <- phylolm::summary.phylolm(mod)$coefficients[[2,4]] #P-value slope (crop model)
                                aic                   <- mod$aic
                                DFintercept           <- intercept - intercept.0   # DF intercept
                                DFslope               <- slope - slope.0           # DF slope
                                intercept.perc        <- round((abs(DFintercept/intercept.0))*100,digits=1)  # Percentage of intercept change
                                slope.perc            <- round((abs(DFslope/slope.0))*100,digits=1)  # Percentage of slope change
                                pval.intercept        <- phylolm::summary.phylolm(mod)$coefficients[[1,4]] #P-value intercept (full model)
                                pval.slope            <- phylolm::summary.phylolm(mod)$coefficients[[2,4]] #P-value slope (full model)
                                aic                   <- mod$aic
                                optpar                <- mod$optpar             #The optimisation paratemer alpha (phylogenetic correlation parameter)
                                n.remov <- i
                                n.percent <- round((n.remov/N)*100,digits=0)
                                rep <- j
                                print(paste("Break = ",n.percent,sep=""));
                                print(paste("Repetition= ",j,""));

                                ### Storing ### Storing values for each simulation
                                #write in a table
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

        #Percentage Significant:
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

        #Summary of all full model details for output
        param0<-list(coef=phylolm::summary.phylolm(mod.0)$coefficients,
                     aic=phylolm::summary.phylolm(mod.0)$aic,
                     optpar=summary(mod.0)$optpar)

        # Function output:
        return(list(analyis.type = "samp_phylolm",
                    formula=formula,
                    full.model.estimates=param0,
                    samp.model.estimates=samp.model.estimates,
                    sign.analysis=perc.sign.tab,
                    data=full.data))

}
