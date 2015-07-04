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


samp_pgls <- function(formula,data,phy,times=20,breaks=seq(.1,.7,.1),model="lambda",
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
        c.data <- data
        N <- nrow(c.data)
        mod.0 <- phylolm::phylolm(formula, data=c.data,model=model,phy=phy)
        a.0 <- mod.0$coefficients[[1]]             # Intercept (full model)
        b.0 <- mod.0$coefficients[[2]]             # Beta (full model)
        p.val.a0 <-    phylolm::summary.phylolm(mod.0)$coefficients[[1,4]] # p.value (intercept)
        p.val.b0 <-    phylolm::summary.phylolm(mod.0)$coefficients[[2,4]] # p.value (slope)
        optpar.0 <- mod.0$optpar
        aic.0    <-mod.0$aic

        #Create the results data.frame
        results<-data.frame("n.removs" =numeric(), "n.percents"=numeric(),
                            "intercept"=numeric(),"DFintercept"=numeric(),
                            "intercept.change"=numeric(),
                            "beta"=numeric(),"DFbeta"=numeric(),
                            "beta.change"=numeric(),"pval.intercept"=numeric(),
                            "pval.beta"=numeric(),"AIC"=numeric(),
                            "optpar"=numeric())

        # Loop:
        counter=1
        limit <- sort(round((breaks)*nrow(c.data),digits=0))
        for (i in limit){
                for (j in 1:times){
                        exclude <- sample(1:N,i)
                        crop.data <- c.data[-exclude,]
                        crop.phy <-  ape::drop.tip(phy,phy$tip.label[exclude])

                        mod=try(phylolm::phylolm(formula, data=crop.data,model=model,phy=crop.phy),TRUE)
                        if(isTRUE(class(mod)=="try-error")) { next }
                        else {
                                ### Calculating model estimates
                                a <-    mod$coefficients[[1]]          # Intercept (crop model)
                                b <-    mod$coefficients[[2]]       # Beta (crop model)
                                DFa <- a - a.0                 # DF intercept
                                DFb <- b - b.0                 # DF beta
                                a.change <- round((abs(DFa/a.0))*100,digits=1)  # Percentage of intercept change
                                b.change <- round((abs(DFb/b.0))*100,digits=1)  # Percentage of beta change
                                pval.a <-    phylolm::summary.phylolm(mod.0)$coefficients[[1,4]] # p.value (intercept)
                                pval.b <-    phylolm::summary.phylolm(mod.0)$coefficients[[2,4]] # p.value (slope)
                                aic.mod <-   mod$aic           # Model AIC
                                optpar <-    mod$optpar# Estimated lambda
                                n.remov <- i
                                n.percent <- round((n.remov/N)*100,digits=0)
                                rep <- j
                                print(paste("Species removed= ",i,sep=""));
                                print(paste("Repetition= ",j,""));


                                ### Storing values for each simulation
                                #write in a table
                                results[counter,1]<- n.remov
                                results[counter,2]<- n.percent
                                results[counter,3]<- a
                                results[counter,4]<- DFa
                                results[counter,5]<- a.change
                                results[counter,6]<- b
                                results[counter,7]<- DFb
                                results[counter,8]<- b.change
                                results[counter,9]<- pval.a
                                results[counter,10]<- pval.b
                                results[counter,11]<- aic.mod
                                results[counter,12]<- optpar


                                counter=counter+1



                        }


                }
        }

        ## Power Analysis:
        res <- results
        times <- table(res$n.removs)
        breaks <- unique(res$n.percents)
        sig.a <- res$pval.intercept > .05
        sig.b <- res$pval.beta > .05
        res$sig.a <- sig.a
        res$sig.b <- sig.b
        power.intercept <- 1-(with(res,tapply(sig.a,n.removs,sum)))/times
        power.beta <- 1-(with(res,tapply(sig.b,n.removs,sum)))/times
        power.tab <- data.frame(percent_sp_removed=breaks,
                                power.intercept=as.numeric(power.intercept),
                                power.beta=as.numeric(power.beta))

        # Function output:
        return(list(output = "samp_pgls",
                    model_estimates=data.frame(intercept = a.0,beta=b.0),
                    results=results,power.analysis=power.tab,data=c.data))

}
