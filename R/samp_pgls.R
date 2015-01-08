#' Sampling effort analysis for pgls regression.
#'
#' \code{samp_pgls} performs sampling effort diagnostics for \code{pgls}
#' regressions. It removes species at random, fits a pgls model without the
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
#' @return The function \code{samp_pgls} returns a list with the following
#'   components:
#' @return \code{model_estimates} Full model estimates
#' @return \code{beta95_IC} Full model beta 95 confidence interval
#' @return \code{results} A data frame with all simulation estimates.
#' @return \code{power_analysis} A data frame with power analysis for each
#' @section Warning: This code is note fully checked. Please be aware.
#' @seealso \code{\link{pgls}}, \code{\link{influ_pgls}}
#' @examples
#' library(caper);library(ggplot2);library(gridExtra)
#' data(shorebird)
#' comp.data <- comparative.data(shorebird.tree, shorebird.data, Species, vcv=TRUE, vcv.dim=3)
#' samp1 <- samp_pgls(log(Egg.Mass) ~ log(M.Mass),data=comp.data)
#' # You can specify the number of replicates and break intervals:
#' samp2 <- samp_pgls(log(Egg.Mass) ~ log(M.Mass),data=comp.data,times=20,breaks=c(.1,.3,.5))
#' # Plot results
#' sensi_plot(samp1)
#' sensi_plot(samp2)
#' @export

samp_pgls <- function(formula,data,times=20,breaks=seq(.1,.7,.1),lambda="ML")
{
          ### Basic error checking:
          if(class(formula)!="formula") stop("Please formula must be
                                             class 'forumla'")
          if(class(data)!="comparative.data") stop("data data must be of class
                        'comparative.data'. See function `comparative.data`.")
          if(length(breaks)<2) stop("please include more then one break
                                    (eg. breaks=c(.3,.5)")
          else

          # FULL MODEL calculations:
          c.data <- data
          N <- nrow(c.data$data)             # Sample size
          mod.0 <- caper::pgls(formula, data=c.data,lambda=lambda)
          sumMod <- summary(mod.0)
          intercept.0 <-    sumMod[[c(5,1)]] # Intercept (full model)
          beta.0 <-    sumMod[[c(5,2)]]      # Beta (full model)
          pval.0 <-    sumMod[[c(5,8)]]      # p.value (full model)
          sd.beta.0 <- sumMod[[c(5,4)]]      # Standart Error (full model)
          df.0 <- sumMod[[2]][2] # Degree if Freedon (full model))
          beta.IC <- qt(0.975,df.0)*sd.beta.0 # Beta CI (full model)
          beta.0.low <- beta.0 - beta.IC  # Low limit of beta CI (full model)
          beta.0.up <- beta.0 + beta.IC   # Up limit of beta CI (full model)

          # Sampling effort analysis:
          intercepts <- as.numeric()
          betas <- as.numeric()
          p.values <- as.numeric()
          n.removs <- as.numeric()
          n.percents <- as.numeric()

          # Loop:
          limit <- sort(round((breaks)*nrow(c.data$data),digits=0))
          for (i in limit){
                    for (j in 1:times){
                              exclude <- sample(1:N,i)
                              crop.data <- c.data[-exclude,]
                              mod=try(caper::pgls(formula, data=crop.data,lambda),TRUE)
                              if(isTRUE(class(mod)=="try-error")) { next }
                              else {
                                        ### Calculating model estimates:
                                        sum.Mod <- summary(mod)
                                        beta <-    sum.Mod[[c(5,2)]]     # Beta
                                        intercept <-    sum.Mod[[c(5,1)]]# Intercept
                                        pval <-    sum.Mod[[c(5,8)]] # p.value
                                        n.remov <- i
                                        n.percent <- round(n.remov/N,digits=1)*10
                                        rep <- j

                                        ### Storing values for each simulation
                                        intercepts <- c(intercepts,intercept)
                                        betas <- c(betas,beta)
                                        p.values <- c( p.values,pval)
                                        n.removs <- c(n.removs,n.remov)
                                        n.percents <- c(n.percents,n.percent)
                              }


                    }
          }


          # Data frame with results:
          estimates <- data.frame(intercepts,betas,p.values,n.removs,n.percents)

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
                      results=estimates,power_analysis=power.tab,data=c.data$data))


}
