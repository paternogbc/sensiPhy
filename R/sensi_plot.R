#' Graphical sensitive analysis for comparative methods
#'
#' \code{sensi_plot} Plot results from \code{samp_pgls},
#' \code{influ_pgls}
#' @aliases sensi_plot
#' @param x output from \code{samp_gls}, \code{influ_gls}
#' @param graphs choose which graphs should be printed on the output ("all", 1,2,3 or 4)
#' @export

### Start:
sensi_plot <- function(x,graphs="all"){

########### samp_pgls graphs ##################

          if (x[[1]] == "samp_pgls"){
                    result    <- x[[3]]
                    power.tab <- x[[4]]

                    # classes of beta.change:
                    result$beta.class <- "class"
                    result[result$beta.change <= 5,]$beta.class <- "within 5%"
                    result[result$beta.change > 5
                           & result$beta.change <= 10 ,]$beta.class <- "higher than 5%"
                    result[result$beta.change > 10,]$beta.class <- "higher than 10%"
                    result$beta.class <- as.factor(result$beta.class)
                    beta.0    <- as.numeric(x[[2]][2])
                    beta.5    <- .05*beta.0
                    beta.10   <- .1*beta.0

                    # reverting the order of the levels
                    result$beta.class =
                            with(result, factor(beta.class,
                                        levels = rev(levels(result$beta.class))))
                    .e <- environment()

                    ## Graphs: Estimated betas ~ % species removed
                    if(length(levels(result$beta.class)) == 3){
                            col = c("skyblue","orange","red2")
                    }
                    if(length(levels(result$beta.class)) == 2){
                            col = c("skyblue","orange")
                    }
                    if(length(levels(result$beta.class)) == 1){
                            col = c("skyblue")
                    }

                    p1 <- ggplot2::ggplot(result,aes(y=beta,x=n.percents,colour=beta.class))+

                            geom_point(size=4,position = "jitter",alpha=.5)+
                            scale_x_continuous(breaks=result$n.percents)+
                              ylab("Estimated Betas")+
                              xlab("% of Species Removed ")+
                            scale_colour_manual(values=col)+
                            geom_hline(yintercept=beta.0,linetype=1,color="red",
                                       size=1,alpha=.6,name="Original Beta")+

                              geom_hline(yintercept=beta.0+beta.5,linetype=2,
                                         alpha=.6)+
                              geom_hline(yintercept=beta.0-beta.5,linetype=2,
                                         alpha=.6)+
                              geom_hline(yintercept=beta.0+beta.10,linetype=2,
                                         alpha=.6)+
                              geom_hline(yintercept=beta.0-beta.10,linetype=2,
                                         alpha=.6)+
                            theme( legend.position = "none",
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

                    p4 <-ggplot2::ggplot(power.tab,
                                         aes(y=power.beta,x=percent_sp_removed))+
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

                    beta.tab <- dplyr::summarise(dplyr::group_by(result,beta.class,n.percents),
                                     proportion=n())
                    ## Correcting for the number of replications per n.percent interval:

                    times <- table(result$n.percents)
                    for (jj in 1:length(times)){
                        a <- beta.tab[beta.tab$n.percents==unique(beta.tab$n.percent)[jj],]$proportion
                        a <- a/times[jj]
                        beta.tab[beta.tab$n.percents==unique(beta.tab$n.percent)[jj],]$proportion <- a

                    }

                    p2 <- ggplot(beta.tab,
                                 aes(y=proportion,x=n.percents,fill=factor(beta.class)))+
                            geom_bar(stat="identity",alpha=.5)+
                            scale_fill_manual(values=col,name="Change in beta")+
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
                                            mDFbetas = mean(abs(DFbeta)),
                                            sd = as.numeric(sd(abs(DFbeta))))
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

                    if (graphs == "all"){
                        gridExtra::grid.arrange(p1,p2,p3,p4,ncol=2,nrow=2)
                    }
                    if (graphs == 1){
                        print(p1)
                    }
                    if (graphs == 2){
                        print(p2)
                    }
                    if (graphs == 3){
                        print(p3)
                    }
                    if (graphs == 4){
                        print(p4)
                    }



          }

########### influ_pgls graphs ##################

          if  (x[[1]] == "influ_pgls")    {

                    .e <- environment()
                    result <- x[[5]]
                    vars <- all.vars(x[[2]])
                    vars2 <- gsub("list","",attr(terms(x[[2]]),"variables"))[-1]
                    intercept.0 <-  as.numeric(x[[3]][1])
                    beta.0 <-  as.numeric(x[[3]][2])

                    # Distribution of estimated betas:
                    p1 <- ggplot2::ggplot(result,aes(x=beta,y=..density..),environment=.e)+
                              geom_histogram(fill="lightyellow", alpha=.9,colour="grey60", size=.2) +
                              geom_density(size=.2) +
                              geom_vline(xintercept = beta.0,color="red",linetype=2,size=.7)+
                              xlab("Estimated Betas")+
                            theme(axis.text = element_text(size=14),
                                  axis.title = element_text(size=16))
                    # Distribution of estimated intercepts:
                    p2 <- ggplot2::ggplot(result,aes(x=intercept,y=..density..),environment=.e)+
                              geom_histogram(fill="lightyellow", alpha=.9,colour="grey60", size=.2) +
                              geom_density(size=.2) +
                              geom_vline(xintercept = intercept.0,color="red",linetype=2,size=.7)+
                              xlab("Estimated Intercepts")+
                            theme(axis.text = element_text(size=14),
                                  axis.title = element_text(size=16))

                    ### Removing species with error:
                    if (isTRUE(class(x[[7]]) != "character" )){
                              x$data <- x$data[-as.numeric(x[[7]]),]
                    }
                    result.tab <- data.frame(x$results,x$data[vars])

                    # Influential points for beta estimation:
                    p3<-ggplot2::ggplot(result.tab,aes(eval(parse(text=vars2[2])),
                                                       eval(parse(text=vars2[1])),
                                                       colour=abs(sDFbeta)),environment=.e)+
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
                                  axis.title = element_text(size=16))

                    # Statistically influential points for beta estimate
                    p4 <- ggplot2::ggplot(result,aes(x=sDFbeta),environment=.e)+
                            geom_histogram(fill="red",color="black",binwidth=.5) +
                            xlab("Standardized Difference in Beta")+
                            geom_histogram(data=subset(result,sDFbeta<2&sDFbeta>-2),
                                           colour="black", fill="white",binwidth=.5)+
                            theme(axis.text = element_text(size=14),
                                  axis.title = element_text(size=16))+
                            geom_vline(xintercept = -2,color="red",linetype=2,size=.7)+
                            geom_vline(xintercept = 2,color="red",linetype=2,size=.7)

                    if (graphs == "all"){
                            suppressMessages(gridExtra::grid.arrange(p1,p2,p3,p4,ncol=2,nrow=2))
                    }
                    if (graphs == 1){
                            suppressMessages(print(p1))
                    }
                    if (graphs == 2){
                            suppressMessages(print(p2))
                    }
                    if (graphs == 3){
                            suppressMessages(print(p3))
                    }
                    if (graphs == 4){
                            suppressMessages(print(p4))
                    }
          }



}

