#' Sampling effort analysis for gls phylogenetic regression.
#'
#' \code{samp_gls} performs sample size sensitive analysis for \code{pgls}
#' regressions. It removes species at random, fits a pgls model without the
#' species and store the results of the model estimates. The percentage of
#' species removed is specified with \code{breaks} and the number of simulations
#' per break is defined by \code{times}.
#' @aliases samp_gls
#' @inheritParams influ_gls
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
#' @return \code{beta95_IC} Full model beta 95 confidence interval
#' @return \code{results} A data frame with all simulation estimates.
#' @return \code{power_analysis} A data frame with power analysis for each
#' @return \code{data} Original dataset
#' @section Warning: This code is note fully checked. Please be aware.
#' @seealso \code{\link{pgls}}, \code{\link{influ_pgls}}, \code{\link{samp_pgls}}
#' @examples
#' library(caper);library(ggplot2);library(gridExtra)
#' data(shorebird)
#' comp.data <- comparative.data(shorebird.tree, shorebird.data, Species, vcv=TRUE, vcv.dim=3)
# # First we need to match tip.labels with rownames in data:
#' sp.ord <- match(shorebird.tree$tip.label, rownames(shorebird.data))
#' shorebird.data <- shorebird.data[sp.ord,]
#' samp1 <- samp_gls(log(Egg.Mass) ~ log(M.Mass),data=shorebird.data,phy=shorebird.tree)
#' sensi_plot(samp1)
#' @export


samp_gls <- function(formula,data,phy,times=20,breaks=seq(.1,.7,.1))
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
        else

        # FULL MODEL calculations:

        c.data <- data
        N <- nrow(c.data)
        cor.0 <- ape::corPagel(1,phy=phy,fixed=F)

        mod.0 <- nlme::gls(formula, data=c.data,method="ML",correlation=cor.0)
        sumMod <- as.data.frame(summary(mod.0)$tTable)

        intercept.0 <-    sumMod[1,1]       # Intercept (full model)
        beta.0 <-    sumMod[2,1]            # Beta (full model)
        pval.0 <-    sumMod[2,4]            # p.value (full model)
        sd.beta.0 <- sumMod[2,2]            # Standart Error (full model)
        df.0 <- N-2                         # Degrees of Freedon (full model))
        if (df.0 >= 30) {
                beta.IC <- 1.96*sd.beta.0
        }
        else {
        beta.IC <- qt(0.975,df.0)*sd.beta.0 # Beta CI (full model)
        }
        beta.0.low <- beta.0 - beta.IC      # Low limit of beta CI (full model)
        beta.0.up <-  beta.0 + beta.IC      # Up limit of beta CI (full model)

        # Sampling effort analysis:
        intercepts <- as.numeric()
        betas <- as.numeric()
        DFbetas <- as.numeric()
        beta.change <- NULL
        p.values <- as.numeric()
        n.removs <- as.numeric()
        n.percents <- as.numeric()

        # Loop:
        limit <- sort(round((breaks)*nrow(c.data),digits=0))
        for (i in limit){
                for (j in 1:times){
                        exclude <- sample(1:N,i)
                        crop.data <- c.data[-exclude,]
                        crop.phy <-  ape::drop.tip(phy,phy$tip.label[exclude])
                        crop.cor <- ape::corPagel(1,phy=crop.phy,fixed=F)

                        mod=try(nlme::gls(formula, data=crop.data,
                                          correlation=crop.cor,method="ML"),TRUE)
                        if(isTRUE(class(mod)=="try-error")) { next }
                        else {
                                ### Calculating model estimates:
                                sum.Mod <- as.data.frame(summary(mod)$tTable)
                                intercept <-    sum.Mod[1,1]       # Intercept (crop model)
                                beta <-    sum.Mod[2,1]            # Beta (crop model)
                                DFbeta <- beta - beta.0
                                if (abs(DFbeta) < 0.05*beta.0)
                                        b.change = "within 5%"

                                if (abs(DFbeta) > 0.05*beta.0)
                                        b.change = "higher than 5%"
                                if (abs(DFbeta) > 0.1*beta.0){
                                        b.change = "higher than 10%"
                                }
                                if (abs(DFbeta) > 0.3*beta.0){
                                        b.change = "higher than 30%"
                                }



                                pval <-    sum.Mod[2,4]            # p.value (crop model)
                                n.remov <- i
                                n.percent <- round((n.remov/N)*100,digits=0)
                                rep <- j

                                ### Storing values for each simulation
                                intercepts <- c(intercepts,intercept)
                                betas <- c(betas,beta)
                                DFbetas <- c(DFbetas,DFbeta)
                                beta.change <- c(beta.change,b.change)
                                p.values <- c( p.values,pval)
                                n.removs <- c(n.removs,n.remov)
                                n.percents <- c(n.percents,n.percent)
                        }


                }
        }

        # Data frame with results:
        estimates <- data.frame(intercepts,betas,DFbetas,beta.change,p.values,n.removs,n.percents)

        ## Power Analysis:
        times <- table(estimates$n.removs)[1]
        breaks <- unique(estimates$n.percents)
        simu.sig <- estimates$p.values > .05
        estimates$simu.sig <- simu.sig
        power <- 1-(with(estimates,tapply(simu.sig,n.removs,sum)))/times
        power.tab <- data.frame(percent_sp_removed=breaks,power)
        estimates <- estimates[,-ncol(estimates)]

        param0 <- data.frame(beta.0,intercept.0)
        beta_IC <- data.frame(beta.low=beta.0.low,beta.up=beta.0.up)
        return(list(model_estimates=param0,beta_95_IC=beta_IC,
                    results=estimates,power_analysis=power.tab,data=c.data))


}
