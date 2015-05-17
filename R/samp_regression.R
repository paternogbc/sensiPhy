#' Sampling effort analysis for gls phylogenetic regression.
#'
#' \code{samp_regression} performs sample size sensitive analysis for \code{pgls}
#' regressions with continuous variables. It removes species at random, fits a pgls model without the
#' species and store the results of the model estimates. The percentage of
#' species removed is specified with \code{breaks} and the number of simulations
#' per break is defined by \code{times}.
#' @aliases samp_regression
#' @inheritParams influ_pgls
#' @param breaks Percentage intervals to remove species. For example:
#'   \code{breaks = c(.1,.2,.3)},removes 10,20 and 30 percentage of species at
#'   random in each simulation.
#' @param times The number of times to repeat each simulation (per
#'   \code{breaks}) interval.
#' @details This functions only works with continuous variables on Y and X. More
#' then one explanatory variable and interactions are allowed in the model.
#' Future implementation will deal categorical explanatory variables.
#' @return The function \code{samp_gls} returns a list with the following
#'   components:
#' @return \code{model_estimates} Full model estimates
#' @return \code{intervals} Full model 95% confidence interval
#' @return \code{results} A data frame with all simulation estimates.
#' @return \code{power_analysis} A data frame with power analysis
#' @return \code{data} Original dataset
#' @section Warning: This code is note fully checked. Please be aware.
#' @seealso \code{\link{pgls}}, \code{\link{influ_pgls}}, \code{\link{samp_pgls}}
#' @export


samp_regression <- function(formula,data,phy,times=50,breaks=seq(.1,.5,.1))
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

### Full model estimates:
        c.data <- data
        N <- nrow(c.data)
        cor.0 <- ape::corPagel(1,phy=phy,fixed=F)

        ### Full model:
        mod.0 <- nlme::gls(formula, data=c.data,method="ML",correlation=cor.0)
        sumMod <- as.data.frame(summary(mod.0)$tTable)
        intervals.0 <- nlme::intervals(mod.0)
        params <- rownames(sumMod)
        n.params <- length(rownames(sumMod))


# Sampling effort analysis:
results <- NULL
n.removs <- as.numeric()
n.percents <- as.numeric()

limit <- sort(round((breaks)*nrow(c.data),digits=0))
for (i in limit){
        for (j in 1:times){
                exclude <- sample(1:N,i)
                crop.data <- c.data[-exclude,]
                crop.phy <-  ape::drop.tip(phy,phy$tip.label[exclude])
                crop.cor <- ape::corPagel(0,phy=crop.phy,fixed=F)

                mod=try(nlme::gls(formula, data=crop.data,
                                  correlation=crop.cor,method="ML"),TRUE)
                if(isTRUE(class(mod)=="try-error")) { next }
                else {
                        ### Calculating model estimates:
                        sum.Mod <- as.data.frame(summary(mod)$tTable)
                        n.remov <- i
                        n.percent <- round((n.remov/N)*100,digits=0)


                        estim <- data.frame(n.remov,
                                            rm.percent = rep(n.percent,n.params),
                                            parameter=rownames(sum.Mod),
                                            estimate=sum.Mod[,1],
                                            sdError = sum.Mod[,2],
                                            DFestimate = sum.Mod[,1]-sumMod[,1],
                                            p.val=sum.Mod[,4])
                        results <- rbind(results,estim)
                        rep <- j

                }


        }
}

## Power Analysis:
reps <- as.numeric(table(results$rm.percent))/n.params
non.sig <- with(results,
                tapply(p.val,list(parameter,rm.percent),
                       function(x)sum(x > 0.05)))
power <- 1- (non.sig/reps)
power.values <- NULL
for(k in 1:n.params){
        value <- as.numeric(power[k,])
        power.values <- c(power.values,value)
}
power.tab <- data.frame(percent = rep(breaks,n.params),
                        parameter = rep(params,each=length(breaks)),
                        power = power.values)

## Estimate differences:
estimate.tab <- dplyr::summarise(dplyr::group_by(results,rm.percent,parameter),
                                 Estimate = mean(estimate),
                                 sdEstimate = sd(estimate),
                                 mDFestimate = mean(abs(DFestimate)))

#### Output:
print(sumMod);
print(paste("See $results; $summary_results and $power_analysis for details"))

return(list(model_estimates=sumMod, intervals=intervals.0,
            results = results,
            summary_results = estimate.tab,
            power_analysis=power.tab,
            data= c.data)
)
}
