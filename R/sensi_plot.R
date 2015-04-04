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
                    result    <- x[[3]]
                    beta.0    <- as.numeric(x[[1]][1])
                    beta.5    <- .05*beta.0
                    beta.10   <- .1*beta.0
                    result$beta.change = with(result, factor(beta.change,
                                        levels = rev(levels(beta.change))))
                    .e <- environment()

                    ## Graphs: Estimated betas ~ % species removed

                    p1 <- ggplot2::ggplot(result,aes(y=betas,x=n.percents,colour=beta.change))+
                            geom_hline(yintercept=beta.0,linetype=1,color="red",
                                       size=1,alpha=.6,name="Original Beta")+
                            geom_point(size=2,position = "jitter")+
                            scale_colour_brewer(palette="Reds",
                                                name="Deviation from original beta")+
                              scale_x_continuous(breaks=result$n.percents)+
                              ylab("Estimated Betas")+
                              xlab("% of Species Removed ")+

                              geom_hline(yintercept=beta.0+beta.5,linetype=2,
                                         alpha=.6)+
                              geom_hline(yintercept=beta.0-beta.5,linetype=2,
                                         alpha=.6)+
                              geom_hline(yintercept=beta.0+beta.10,linetype=2,
                                         alpha=.6)+
                              geom_hline(yintercept=beta.0-beta.10,linetype=2,
                                         alpha=.6)+
                            theme( legend.position = "top",
                                   legend.direction = "horizontal",
                                   legend.text=element_text(size=14),
                                   legend.title=element_text(size=14),
                                   axis.text=element_text(size=14),
                                   axis.title=element_text(size=16),
                                   legend.key.width=unit(.5,"line"),
                                   legend.key.size = unit(.5,"cm"),
                                   panel.background = element_rect(fill="white",
                                                                   colour="black"))

                    ## Power Analysis: p.value
                    times <- table(result$n.percents)
                    breaks <- unique(result$n.percents)
                    simu.sig <- result$p.values > .05
                    result$simu.sig <- simu.sig
                    p.out <- (with(result,tapply(simu.sig,n.removs,sum))/times)
                    power <- as.numeric(1-p.out)
                    power.tab <- data.frame(breaks,power)
                    p4 <-ggplot2::ggplot(power.tab,aes(y=power,x=breaks))+
                              scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.1))+
                              scale_x_continuous(breaks=result$n.percents)+
                              xlab("% Species removed")+
                              geom_point(size=5,colour="red")+
                              geom_line(colour="red")+
                              ylab("Power  [p-value]")+
                            theme(axis.text=element_text(size=14),
                                  axis.title=element_text(size=16),
                                  panel.background = element_rect(fill="white",
                                                                  colour="black"))

                    ## Power Analysis: beta (percentage of betas > or < then 5% or 10%)

                    beta.tab <- dplyr::summarise(dplyr::group_by(result,beta.change,n.percents),
                                     proportion=n())
                    ## Correcting for the number of replications per n.percent interval:
                    attach(beta.tab)
                    for (jj in 1:length(times)){
                        a <- beta.tab[n.percents==unique(beta.tab$n.percent)[jj],]$proportion
                        a <- a/times[jj]
                        beta.tab[n.percents==unique(beta.tab$n.percent)[jj],]$proportion <- a

                    }
                    detach(beta.tab)
                    p2 <- ggplot(beta.tab,
                                 aes(y=proportion,x=n.percents,fill=factor(beta.change)))+
                            geom_bar(stat="identity",alph=.5)+
                            scale_fill_brewer(palette="Reds",
                                              name="")+
                            scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.1))+
                            scale_x_continuous(breaks=result$n.percents)+
                            theme( legend.position = "top",
                                   legend.direction = "horizontal",
                                   legend.text=element_text(size=14),
                                   legend.title = element_text(size=12),
                                   axis.text=element_text(size=14),
                                   axis.title=element_text(size=16),
                                   legend.key.width=unit(.5,"line"),
                                   legend.key.size = unit(.5,"cm"),
                                   panel.background = element_rect(fill="white",
                                                                   colour="black"))+
                            xlab("% of Species Removed")+
                            ylab("Proportion of estimated betas")

                    ### standardized Beta across % of species removed:
                    m.DFbetas <- dplyr::summarise(dplyr::group_by(result,n.percents),
                                            mDFbetas = mean(abs(DFbetas)),
                                            sd = as.numeric(sd(abs(DFbetas))))
                    p3 <- ggplot(m.DFbetas,aes(y=mDFbetas,x=n.percents))+
                            geom_point(size=5,colour="red",alpha=.6)+
                            geom_errorbar(aes(ymin=mDFbetas-sd, ymax=mDFbetas+sd),
                                          colour="red", width=.8)+
                            scale_x_continuous(breaks=result$n.percents)+
                            theme(axis.title=element_text(size=16),
                                  axis.text = element_text(size=14),
                                  panel.background = element_rect(fill="white",
                                                                  colour="black"))+
                            xlab("% of Species Removed")+
                            ylab("Mean DFbetas (+- SD)")
                    suppressWarnings(gridExtra::grid.arrange(p1,p2,p3,p4,ncol=2,nrow=2))
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
                            theme(axis.text = element_text(size=14,colour="black"),
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

