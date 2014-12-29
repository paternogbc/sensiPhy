#' Plot results from package `sensiC`
#'
#' \code{sensi_plot} Performs Sampling effort analysis for \code{pgls}
#' regression. It removes species at random, fits a pgls model without the
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
#' @details This functions only works for simple linear regression (y = bx + a).
#'   Future implementation will deal with more complex models.
#' @section Warning: This code is note fully checked. Please be aware.
#' @seealso \code{\link{[caper]pgls}}
#' @return The function \code{influpgls} returns a list with the following
#'   components:
#' @return \code{errors} Species removal that showed erros during pgls fit
#' @return \code{formula} The model formula
#' @return \code{model_estimates} Full model estimates
#' @return \code{beta95_IC} Full model beta 95 confidence interval
#' @return \code{influential_species} Most influential species for beta and
#'   intercept
#' @return \code{results} A data frame with all simulation estimates. DFbeta and
#'   DFintercept represent absolute difference between full model and simulation
#'   estimate, beta and intercept respectively.
#' @export

### Start:
sensi_plot <- function(x, method="sampling"){
          if (method == "sampling"){
                    result <- x[[3]]
                    beta.0 <- as.numeric(x[[1]][1])
                    beta.0.low <- as.numeric(x[[2]][1])
                    beta.0.up <- as.numeric(x[[2]][2])
                    .e <- environment()

                    ## Graphs: Estimated betas ~ % species removed
                    p2 <- ggplot(result,aes(y=betas,x=n.percents))+
                              geom_point(size=3,alpha=.7)+
                              scale_x_continuous(breaks=result$n.percents)+
                              ylab("Estimated Betas")+
                              xlab("% of Species Removed ")+
                              geom_hline(yintercept=beta.0.low,linetype=2,color="red")+
                              geom_hline(yintercept=beta.0.up,linetype=2,color="red")+
                              geom_hline(yintercept=beta.0,linetype=2,color="red",size=1.1)

                    ## Mean estimated Betas:
                    med <- with(result,tapply(betas,n.removs,mean))
                    Sdev <- with(result,tapply(betas,n.removs,sd))
                    n.sp <- nrow(x[[5]]$data)-as.numeric(rownames(med))
                    result.med <- data.frame(med,Sdev,n.sp,beta.0)
                    p1 <- ggplot(result.med,aes(y=med,x=n.sp),environment=.e)+
                              geom_point(size=3,alpha=.7)+
                              scale_x_continuous(breaks=c(n.sp,nrow(x[[5]]$data)))+
                              geom_errorbar(aes(ymin=med-Sdev, ymax=med+Sdev), width=.1)+
                              geom_hline(yintercept= beta.0,linetype=2,color="red")+
                              geom_hline(yintercept=beta.0.low,linetype=2,color="red")+
                              geom_hline(yintercept=beta.0.up,linetype=2,color="red")+
                              xlab("Number of Species") + ylab("Mean Beta (+-SD)")+
                              geom_point(aes(x=nrow(x[[5]]$data),y=beta.0,size=3,colour="red"))+
                              theme(legend.position="none")+
                              xlab("Number of species")

                    ## Power Analysis: p.value
                    times <- table(result$n.removs)[1]
                    breaks <- unique(result$n.percents)
                    simu.sig <- result$p.values > .05
                    result$simu.sig <- simu.sig
                    power <- 1-(with(result,tapply(simu.sig,n.removs,sum)))/times
                    power.tab <- data.frame(breaks,power)
                    p3 <-ggplot(power.tab,aes(y=power,x=breaks))+
                              scale_y_continuous(limits=c(0,1))+
                              scale_x_continuous(breaks=breaks)+
                              xlab("% Species removed")+
                              geom_point(size=5,colour="red")+
                              geom_line(colour="red")+
                              ylab("Power  [p-value]")

                    ## Power Analysis: beta (percentage of betas > or < then CI)
                    beta.high <- result$betas > beta.0.up
                    beta.low <- result$betas < beta.0.low
                    result$beta.out.CI <- beta.high+beta.low
                    power <- as.numeric(1-(with(result,tapply(beta.out.CI,n.removs,sum)))/times)
                    power.tab <- data.frame(breaks,power)
                    p4 <- ggplot(power.tab,aes(y=power,x=breaks))+
                              scale_y_continuous(limits=c(0,1))+
                              scale_x_continuous(breaks=breaks)+
                              xlab("% Species removed")+
                              geom_point(size=5,colour="red")+
                              geom_line(colour="red")+
                              ylab("Power  [Beta]")
                    grid.arrange(p1,p2,p3,p4,ncol=2,nrow=2)
          }
          if (method == "influence"){
                    .e <- environment()
                    result <- x[[6]]
                    vars <- all.vars(x[[2]])
                    intercept.0 <-  as.numeric(x[[3]][1])
                    beta.0 <-  as.numeric(x[[3]][2])
                    p1 <- ggplot(result,aes(x=betas,y=..density..),environment=.e)+
                              geom_histogram(fill="lightyellow", alpha=.9,colour="grey60", size=.2) +
                              geom_density(size=.2) +
                              geom_vline(xintercept = beta.0,color="red",linetype=2,size=.7)+
                              xlab("Estimated Betas")
                    p2 <- ggplot(result,aes(x=intercepts,y=..density..),environment=.e)+
                              geom_histogram(fill="lightyellow", alpha=.9,colour="grey60", size=.2) +
                              geom_density(size=.2) +
                              geom_vline(xintercept = intercept.0,color="red",linetype=2,size=.7)+
                              xlab("Estimated Intercepts")
                    if (isTRUE(class(x[[1]]) != "character" )){
                              x$data <- x$data[-as.numeric(x[[1]]),]
                    }
                    result.tab <- data.frame(x$results,x$data[vars])
                    p3<-ggplot(result.tab,aes(y=get(vars[1]),
                                              x=get(vars[2]),
                                              colour=abs(DFbetas)),environment=.e,)+
                              geom_point(size=3,alpha=.8)+
                              scale_colour_gradient( low="black",high="red",name="")+
                              theme(legend.key.width = unit(.2,"cm"),
                                    panel.background=element_rect(fill="white",colour="black"),
                                    panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank())+
                              ylab(vars[1])+
                              xlab(vars[2])+
                              ggtitle("Change in Beta estimate")
                    p4<-ggplot(result.tab,aes(y=get(vars[1]),
                                              x=get(vars[2]),
                                              colour=abs(DFintercepts)),environment=.e,)+
                              geom_point(size=3,alpha=.8)+
                              scale_colour_gradient( low="black", high="red",name="")  +
                              theme(legend.key.width = unit(.2,"cm"),
                                    panel.background=element_rect(fill="white",colour="black"),
                                    panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank())+
                              ylab(all.vars(x$formula)[1])+
                              xlab(all.vars(x$formula)[2])+
                              ggtitle("Change in Intercept estimate")

                    suppressWarnings(grid.arrange(p1,p2,p3,p4,nrow=2,ncol=2))
          }



}

