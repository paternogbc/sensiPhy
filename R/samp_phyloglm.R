#' Sampling effort analysis for phylogenetic generalized linear models (logistic regression).
#'
#' \code{samp_phyloglm} performs sample size sensitivity analysis for logistic regresssions using \code{phyloglm}.
#' It removes species at random, fits a phyloglm model without the
#' species and stores the results of the model estimates. The percentage of
#' species removed is specified with \code{breaks} and the number of simulations
#' per break is defined by \code{times}.
#' @aliases samp_phyloglm
#' @inheritParams samp_pgls
#' @param breaks Percentage intervals to remove species. For example:
#'   \code{breaks = c(.1,.2,.3)},removes 10, 20 and 30 percent of species at
#'   random in each simulation.
#' @param times The number of times to repeat each simulation (per
#'   \code{breaks}) interval.
#' @param btol Bound on the linear predictor to bound the searching space (see ?phyloglm)
#' @details This functions currently only works for logistic regressions, Poisson regressions
#' will be implemented later (see R-package phylolm, Ives and Garland 2002 and Ives and Garland 2010).
#' Also it currently only uses logistic_MPLE as a method, and can only deal with simple logistic regression models.
#' @return The function \code{samp_phyloglm} returns a list with the following
#'   components:
#' @return \code{model_estimates} Full model estimates. This contains slope and intercept for the full model and associated p-values, as well as the phylogenetic correlation parameter alpha from \code{phyloglm}.
#' @return \code{results} A data frame with all simulation estimates.
#' @return \code{power_analysis} A data frame with power analysis for each
#' @return \code{data} Original dataset
#' @section Warning: This code is not fully checked. Please be aware.
#' @seealso \code{\link{pgls}}, \code{\link{influ_pgls}}, \code{\link{samp_pgls}}, \code{\link{samp_phyloglm}}
#' @examples
#' library(caper);library(ggplot2);library(gridExtra);library(phylolm)
#' data(shorebird)
#  #First, we need to match tip.labels with rownames in data:
#' sp.ord <- match(shorebird.tree$tip.label, rownames(shorebird.data))
#' shorebird.data <- shorebird.data[sp.ord,]
#' #Create a binary variable (large egg / small egg), for illustration purposes.
#' mean(shorebird.data$Egg.Mass)
#' shorebird.data$Egg.Mass.binary<-ifelse(shorebird.data$Egg.Mass>30,1,0) #Turn egg mass into a binary variable
#' table(shorebird.data$Egg.Mass.binary) #Mostly small eggs.
#' #Run samp_phyloglm
#' samp_phyloglm1<-samp_phyloglm(formula = Egg.Mass.binary~M.Mass,data = shorebird.data,phy=shorebird.tree,btol = 50)
#' sensi_plot(samp_phyloglm1)
#' @export


samp_phyloglm <- function(formula,data,phy,times=20,breaks=seq(.1,.7,.1),btol=50,...)
{### Basic error checking:
        if(class(formula)!="formula") stop("Please formula must be
                                           class 'formula'")
        if(class(data)!="data.frame") stop("data data must be of class
                                                 'data.frame'. See ?phyloglm")
        if(length(breaks)<2) stop("Please include more than one break
                                  (eg. breaks=c(.3,.5)")
        if(class(phy) != "phylo") stop("phy must be of class 'phylo'")

        if (sum(rownames(data) != phy$tip.label) > 0) stop("Species must be in the same order
                                                      in data and phy")
        else

        # FULL MODEL calculations:
        full.data <- data
        N <- nrow(full.data)

        mod.0 <- phylolm::phyloglm(formula, data=full.data,
                                   phy=phy,method="logistic_MPLE",btol=btol,...)
        if(isTRUE(mod.0$convergence!=0)) stop("Null model failed to converge, consider changing btol. See ?phyloglm")
        #The above line checks if the null model converges, and if not terminates with the suggestion to change btol.
        else

        intercept.0 <-    mod.0$coefficients[[1]]       # Intercept (full model)
        slope.0 <-    mod.0$coefficients[[2]]            # slope (full model)
        alpha.0 <-    mod.0$alpha                #Alpha (phylogenetic correlation parameter)
        pval.intercept.0 <- phylolm::summary.phyloglm(mod.0)$coefficients[[1,4]] #P-value intercept (full model)
        pval.slope.0 <- phylolm::summary.phyloglm(mod.0)$coefficients[[2,4]]  #P-value slope (full model)

        # Sampling effort analysis:
        intercepts <- as.numeric()
        slopes <- as.numeric()
        DFslopes <- as.numeric()
        slope.change <- NULL
        p.values.slope <- as.numeric()
        p.values.intercept <- as.numeric()
        alphas <- as.numeric()
        n.removs <- as.numeric()
        n.percents <- as.numeric()

        # Loop:
        limit <- sort(round((breaks)*nrow(full.data),digits=0))
        for (i in limit){
                for (j in 1:times){
                        exclude <- sample(1:N,i)
                        crop.data <- full.data[-exclude,]
                        crop.phy <-  ape::drop.tip(phy,phy$tip.label[exclude])

                        mod=try(phylolm::phyloglm(formula, data=crop.data,
                                                  phy=crop.phy,method="logistic_MPLE",btol=btol,...),TRUE)
                        if(isTRUE(class(mod)=="try-error")) { next }
                        else {
                                ### Calculating model estimates:
                                intercept <-    mod$coefficients[[1]]       # Intercept (crop model)
                                slope <-    mod$coefficients[[2]]            # slope (crop model)
                                alpha <-    mod$alpha                #Alpha (phylogenetic correlation parameter)
                                DFslope <- slope - slope.0

                                if (abs(DFslope) < 0.05*slope.0)
                                        b.change = "within 5%"

                                if (abs(DFslope) > 0.05*slope.0)
                                        b.change = "higher than 5%"
                                if (abs(DFslope) > 0.1*slope.0){
                                        b.change = "higher than 10%"
                                }
                                else

                                pval.slope <- phylolm::summary.phyloglm(mod)$coefficients[[2,4]]
                                pval.intercept <- phylolm::summary.phyloglm(mod)$coefficients[[1,4]]
                                alpha<-mod$alpha
                                n.remov <- i
                                n.percent <- round((n.remov/N)*100,digits=0)
                                rep <- j

                                ### Storing values for each simulation
                                intercepts <- c(intercepts,intercept)
                                slopes <- c(slopes,slope)
                                DFslopes <- c(DFslopes,DFslope)
                                slope.change <- c(slope.change,b.change)
                                p.values.slope <- c( p.values.slope,pval.slope)
                                p.values.intercept <- c(p.values.intercept,pval.intercept)
                                alphas <- c(alphas,alpha)
                                n.removs <- c(n.removs,n.remov)
                                n.percents <- c(n.percents,n.percent)
                        }


                }
        }

        # Data frame with results:
        estimates <- data.frame(intercepts,slopes,DFslopes,slope.change,p.values.slope,
                                p.values.intercept,alphas,n.removs,n.percents)

        ## Power Analysis:
        times <- table(estimates$n.removs)[1]
        breaks <- unique(estimates$n.percents)
        simu.sig <- estimates$p.values.slope > .05
        estimates$simu.sig <- simu.sig
        power <- 1-(with(estimates,tapply(simu.sig,n.removs,sum)))/times
        power.tab <- data.frame(percent_sp_removed=breaks,power)
        estimates <- estimates[,-ncol(estimates)]

        param0 <- data.frame(slope.0,pval.slope.0,intercept.0,pval.intercept.0,alpha.0)
        return(list(model_estimates=param0,
                    results=estimates,power_analysis=power.tab,data=full.data))
#Consider: print also the fitted formula to the output, as a reminder for people.

}


