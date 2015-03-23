#' Plot results from package `sensiC`
#'
#' \code{sensi_plot} Plot results from \code{samp_pgls}, \code{samp_gls},
#' \code{influ_pgls},\code{influ_gls}
#' @aliases sensi_plot
#' @param x output from \code{samp_pgls}, \code{influ_pgls}, \code{samp_gls} and
#' code{influ_gls}
#' @export

### Start:
sensi_plot <- function(x){
          if (length(x) == 5){
                    result <- x[[3]]
                    beta.0 <- as.numeric(x[[1]][1])
                    beta.0.low <- as.numeric(x[[2]][1])
                    beta.0.up <- as.numeric(x[[2]][2])
                    .e <- environment()

                    ## Graphs: Estimated betas ~ % species removed
                    p2 <- ggplot2::ggplot(result,aes(y=betas,x=n.percents))+
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
                    n.sp <- nrow(x[[5]])-as.numeric(rownames(med))
                    result.med <- data.frame(med,Sdev,n.sp,beta.0)
                    p1 <- ggplot2::ggplot(result.med,aes(y=med,x=n.sp),environment=.e)+
                              geom_point(size=3,alpha=.7)+
                              scale_x_continuous(breaks=c(n.sp,nrow(x[[5]]$data)))+
                              geom_errorbar(aes(ymin=med-Sdev, ymax=med+Sdev), width=.1)+
                              geom_hline(yintercept= beta.0,linetype=2,color="red")+
                              geom_hline(yintercept=beta.0.low,linetype=2,color="red")+
                              geom_hline(yintercept=beta.0.up,linetype=2,color="red")+
                              xlab("Number of Species") + ylab("Mean Beta (+-SD)")+
                              geom_point(aes(x=nrow(x[[5]]),y=beta.0,size=3,colour="red"))+
                              theme(legend.position="none")+
                              xlab("Number of species")

                    ## Power Analysis: p.value
                    times <- table(result$n.removs)[1]
                    breaks <- unique(result$n.percents)
                    simu.sig <- result$p.values > .05
                    result$simu.sig <- simu.sig
                    power <- 1-(with(result,tapply(simu.sig,n.removs,sum)))/times
                    power.tab <- data.frame(breaks,power)
                    p3 <-ggplot2::ggplot(power.tab,aes(y=power,x=breaks))+
                              scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.05))+
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
                    p4 <- ggplot2::ggplot(power.tab,aes(y=power,x=breaks))+
                              scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.05))+
                              scale_x_continuous(breaks=breaks)+
                              xlab("% Species removed")+
                              geom_point(size=5,colour="red")+
                              geom_line(colour="red")+
                              ylab("Power  [Beta]")
                    suppressWarnings(gridExtra::grid.arrange(p1,p2,p3,p4,ncol=2,nrow=2))
          }
          else      {

                    .e <- environment()
                    result <- x[[6]]
                    vars <- all.vars(x[[2]])
                    intercept.0 <-  as.numeric(x[[3]][1])
                    beta.0 <-  as.numeric(x[[3]][2])

                    # Distribution of estimated betas:
                    p1 <- ggplot2::ggplot(result,aes(x=betas,y=..density..),environment=.e)+
                              geom_histogram(fill="lightyellow", alpha=.9,colour="grey60", size=.2) +
                              geom_density(size=.2) +
                              geom_vline(xintercept = beta.0,color="red",linetype=2,size=.7)+
                              xlab("Estimated Betas")
                    # Distribution of estimated intercepts:
                    p2 <- ggplot2::ggplot(result,aes(x=intercepts,y=..density..),environment=.e)+
                              geom_histogram(fill="lightyellow", alpha=.9,colour="grey60", size=.2) +
                              geom_density(size=.2) +
                              geom_vline(xintercept = intercept.0,color="red",linetype=2,size=.7)+
                              xlab("Estimated Intercepts")
                    if (isTRUE(class(x[[1]]) != "character" )){
                              x$data <- x$data[-as.numeric(x[[1]]),]
                    }
                    result.tab <- data.frame(x$results,x$data[vars])

                    # Influential points for beta estimation:
                    p3<-ggplot2::ggplot(result.tab,aes(y=get(vars[1]),
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

                    # Influential points for intercept estimation
                    p4<-ggplot2::ggplot(result.tab,aes(y=get(vars[1]),
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

                    suppressWarnings(gridExtra::grid.arrange(p1,p2,p3,p4,nrow=2,ncol=2))
          }



}

