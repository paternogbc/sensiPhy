#' Graphical sensitive analysis for comparative methods
#'
#' \code{sensi_plot} Plot results from \code{samp_gls},
#' \code{influ_gls}
#' @aliases sensi_plot
#' @param x output from \code{samp_gls}, \code{influ_gls}
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
                              geom_hline(yintercept=beta.0,linetype=2,color="red",size=1.1)+
                            theme(axis.text=element_text(size=14),
                                  axis.title=element_text(size=16))

                    ## Power Analysis: p.value
                    times <- table(result$n.removs)
                    breaks <- unique(result$n.percents)
                    simu.sig <- result$p.values > .05
                    result$simu.sig <- simu.sig
                    p.out <- (with(result,tapply(simu.sig,n.removs,sum))/times)
                    power <- as.numeric(1-p.out)
                    power.tab <- data.frame(breaks,power)
                    p3 <-ggplot2::ggplot(power.tab,aes(y=power,x=breaks))+
                              scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.1))+
                              scale_x_continuous(breaks=breaks)+
                              xlab("% Species removed")+
                              geom_point(size=5,colour="red")+
                              geom_line(colour="red")+
                              ylab("Power  [p-value]")+
                            theme(axis.text=element_text(size=14),
                                  axis.title=element_text(size=16))

                    ## Power Analysis: beta (percentage of betas > or < then CI)
                    beta.high <- result$betas > beta.0.up
                    beta.low <- result$betas < beta.0.low
                    result$beta.out.CI <- beta.high+beta.low
                    b.out <-(with(result,tapply(beta.out.CI,n.removs,sum))/times)
                    p.b.out <- as.numeric(b.out)
                    p.b.in <- 1 -p.b.out
                    proportion <- c(p.b.in,p.b.out)
                    b.class <- rep(c("Within 95% CI" ,"Out of 95% CI"),each=length(breaks))
                    beta.tab <- data.frame(breaks,b.class,proportion)

                    p4 <- ggplot(beta.tab,aes(y=proportion,x=as.factor(breaks),fill=b.class))+
                            geom_bar(stat="identity")+
                            scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.1))+
                            theme( legend.position = "top",
                                   legend.direction = "horizontal",
                                   legend.text=element_text(size=14),
                                   legend.title = element_blank(),
                                   axis.text=element_text(size=14),
                                   axis.title=element_text(size=16))+
                            xlab("% of Species Removed")+
                            ylab("Proportion of estimated betas")
                    suppressWarnings(gridExtra::grid.arrange(p2,p4,p3,ncol=2,nrow=2))
          }
          else      {

                    .e <- environment()
                    result <- x[[6]]
                    vars <- all.vars(x[[2]])
                    vars2 <- gsub("list","",attr(terms(x[[2]]),"variables"))[-1]
                    intercept.0 <-  as.numeric(x[[3]][1])
                    beta.0 <-  as.numeric(x[[3]][2])

                    # Distribution of estimated betas:
                    p1 <- ggplot2::ggplot(result,aes(x=betas,y=..density..),environment=.e)+
                              geom_histogram(fill="lightyellow", alpha=.9,colour="grey60", size=.2) +
                              geom_density(size=.2) +
                              geom_vline(xintercept = beta.0,color="red",linetype=2,size=.7)+
                              xlab("Estimated Betas")+
                            theme(axis.text = element_text(size=14),
                                  axis.title = element_text(size=16))
                    # Distribution of estimated intercepts:
                    p2 <- ggplot2::ggplot(result,aes(x=intercepts,y=..density..),environment=.e)+
                              geom_histogram(fill="lightyellow", alpha=.9,colour="grey60", size=.2) +
                              geom_density(size=.2) +
                              geom_vline(xintercept = intercept.0,color="red",linetype=2,size=.7)+
                              xlab("Estimated Intercepts")+
                            theme(axis.text = element_text(size=14),
                                  axis.title = element_text(size=16))
                    if (isTRUE(class(x[[1]]) != "character" )){
                              x$data <- x$data[-as.numeric(x[[1]]),]
                    }
                    result.tab <- data.frame(x$results,x$data[vars])

                    # Influential points for beta estimation:
                    p3<-ggplot2::ggplot(result.tab,aes(eval(parse(text=vars2[2])),
                                                       eval(parse(text=vars2[1])),
                                                       colour=abs(sDFbetas)),environment=.e)+
                            geom_point(size=3,alpha=.8)+
                            scale_colour_gradient( low="black",high="red",name="")+
                            theme(legend.key.width = unit(.2,"cm"),
                                  panel.background=element_rect(fill="white",colour="black"),
                                  legend.text = element_text(size=14),
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank())+
                            ylab(vars2[1])+
                            xlab(vars2[2])+
                            ggtitle("Standardized Difference in Beta")+
                            theme(axis.text = element_text(size=14),
                                  axis.title = element_text(size=16));p3

                    # Statistically influential points for Beta estimate
                    p4 <- ggplot2::ggplot(result,aes(x=sDFbetas),environment=.e)+
                            geom_histogram(fill="red",color="black",binwidth=.5) +
                            xlab("Standardized Difference in Beta")+
                            geom_histogram(data=subset(result,sDFbetas<2&sDFbetas>-2),
                                           colour="black", fill="white",binwidth=.5)+
                            theme(axis.text = element_text(size=14),
                                  axis.title = element_text(size=16))+
                            geom_vline(xintercept = -2,color="red",linetype=2,size=.7)+
                            geom_vline(xintercept = 2,color="red",linetype=2,size=.7)



                    suppressWarnings(gridExtra::grid.arrange(p1,p2,p3,p4,nrow=2,ncol=2))
          }



}

