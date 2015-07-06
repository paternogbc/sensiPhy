#' Graphical diagnostics for \code{samp_phylolm}
#'
#' \code{plot_samp_phylolm} Plot results from \code{samp_phylolm},
#' @param x output from \code{samp_phylolm}
#' @param graphs choose which graph should be printed on the output ("all", 1,2,3 or 4)
#' @param param choose which parameter ("intercept" or "slope") should be printed

### Start:
plot_samp_phylolm <- function(x,graphs="all",param="slope"){

        ########### samp_pgls graphs ##################
        if (x[[1]] != "samp_phylolm" & x[[1]] != "samp_phyloglm")
                stop("x must be an output from samp_phylolm or samp_phyloglm!")
        else

        x <- test.samp.glm
        result    <- x$samp.model.estimates
        sig.tab <- x$sign.analysis

        # classes of slope.perc:
        result$slope.class <- "class"
        result[result$slope.perc <= 5,]$slope.class <- "within 5%"
        result[result$slope.perc > 5
               & result$slope.perc <= 10 ,]$slope.class <- "higher than 5%"
        result[result$slope.perc > 10,]$slope.class <- "higher than 10%"
        result$slope.class <- as.factor(result$slope.class)
        slope.0    <- as.numeric(x$full.model.estimates$coef[2])
        slope.5    <- .05*slope.0
        slope.10   <- .1*slope.0

        # classes of intercept.perc:
        result$intercept.class <- "class"
        result[result$intercept.perc <= 5,]$intercept.class <- "within 5%"
        result[result$intercept.perc > 5
               & result$intercept.perc <= 10 ,]$intercept.class <- "higher than 5%"
        result[result$intercept.perc > 10,]$intercept.class <- "higher than 10%"
        result$intercept.class <- as.factor(result$intercept.class)
        intercept.0    <- as.numeric(x$full.model.estimates$coef[1])
        intercept.5    <- .05*intercept.0
        intercept.10   <- .1*intercept.0

        # reverting the order of the levels
        result$slope.class =
                with(result, factor(slope.class,
                                    levels = rev(levels(result$slope.class))))
        result$intercept.class =
                with(result, factor(intercept.class,
                                    levels = rev(levels(result$intercept.class))))

        ## Organizing colours:
        if(length(levels(result$slope.class)) == 3){
                col = c("skyblue","orange","red2")
        }
        if(length(levels(result$slope.class)) == 2){
                col = c("skyblue","orange")
        }
        if(length(levels(result$slope.class)) == 1){
                col = c("skyblue")
        }

        ### Graphs--------------------------------------------------------------

        ### Estimated slopes across n.percent:
        s1 <- ggplot2::ggplot(result,aes(y=slope,x=n.percent,
                                         colour=slope.class),
                              environment = environment())+

                geom_point(size=4,position = "jitter",alpha=.5)+
                scale_x_continuous(breaks=result$n.percent)+
                ylab("Estimated slopes")+
                xlab("% of Species Removed ")+
                scale_colour_manual(values=col)+
                geom_hline(yintercept=slope.0,linetype=1,color="red",
                           size=1,alpha=.6,name="Original slope")+

                geom_hline(yintercept=slope.0+slope.5,linetype=2,
                           alpha=.6)+
                geom_hline(yintercept=slope.0-slope.5,linetype=2,
                           alpha=.6)+
                geom_hline(yintercept=slope.0+slope.10,linetype=2,
                           alpha=.6)+
                geom_hline(yintercept=slope.0-slope.10,linetype=2,
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
        ### Estimated intercept across n.percent
        i1<- ggplot2::ggplot(result,aes(y=intercept,x=n.percent,
                                   colour=intercept.class),
                        environment = environment())+

                geom_point(size=4,position = "jitter",alpha=.5)+
                scale_x_continuous(breaks=result$n.percent)+
                ylab("Estimated intercepts")+
                xlab("% of Species Removed ")+
                scale_colour_manual(values=col)+
                geom_hline(yintercept=intercept.0,linetype=1,color="red",
                           size=1,alpha=.6,name="Original intercept")+

                geom_hline(yintercept=intercept.0+intercept.5,linetype=2,
                           alpha=.6)+
                geom_hline(yintercept=intercept.0-intercept.5,linetype=2,
                           alpha=.6)+
                geom_hline(yintercept=intercept.0+intercept.10,linetype=2,
                           alpha=.6)+
                geom_hline(yintercept=intercept.0-intercept.10,linetype=2,
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

        ### Proportion of change.classes across n.percent
        n.perc.times <- as.numeric(table(result$n.percent))
        slope.perc <- with(result,aggregate(data=result,slope ~ slope.class*n.percent,FUN=length))
        a <- colSums(table(slope.perc$slope.class,slope.perc$n.percent))
        slope.perc$slope <- (slope.perc$slope/rep(n.perc.times,
                                                  times=a))*100
        intercept.perc <- with(result,aggregate(data=result,intercept ~ intercept.class*n.percent,FUN=length))
        b <- colSums(table(intercept.perc$intercept.class,intercept.perc$n.percent))
        intercept.perc$intercept <- (intercept.perc$intercept/rep(n.perc.times,
                                                  times=b))*100
        ### Graph: Slope
        s2 <- ggplot(slope.perc,
                     aes(y=slope,x=n.percent,
                         fill=factor(slope.class)),
                     environment = environment())+
                geom_bar(stat="identity",alpha=.5)+
                scale_fill_manual(values=col,name="Change in beta")+
                scale_y_continuous(limits=c(0,100),breaks=seq(0,100,10))+
                scale_x_continuous(breaks=result$n.percent)+
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
                ylab("Proportion of estimated slope")

        ### Graph: Intercept
        i2 <- ggplot(intercept.perc,
                     aes(y=intercept,x=n.percent,
                         fill=factor(intercept.class)),
                     environment = environment())+
                geom_bar(stat="identity",alpha=.5)+
                scale_fill_manual(values=col,name="Change in beta")+
                scale_y_continuous(limits=c(0,100),breaks=seq(0,100,10))+
                scale_x_continuous(breaks=result$n.percent)+
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
                ylab("Proportion of estimated intercept")

        ### Optpar acros % removed species:
        opt <- ggplot(result,aes(y=optpar,x=n.percent,group=as.factor(n.percent)))+
                geom_point()+
                geom_boxplot(fill="red",alpha=.5)+
                scale_x_continuous(breaks=result$n.percent)+
                theme(axis.title=element_text(size=16),
                      axis.text = element_text(size=14),
                      panel.background = element_rect(fill="white",
                                                      colour="black"))+
                xlab("% of Species Removed")+
                ylab("Phylogenetic model parameter")

        ## Significance Analysis : p.value of slope

        s4 <-ggplot2::ggplot(sig.tab,
                             aes(y=perc.sign.slope*100,x=percent_sp_removed),
                             environment = environment())+
                scale_y_continuous(limits=c(0,100),breaks=seq(0,100,10))+
                scale_x_continuous(breaks=result$n.percent)+
                xlab("% Species removed")+
                geom_point(size=5,colour="red")+
                geom_line(colour="red")+
                ylab("% of significant Slopes")+
                theme(axis.text=element_text(size=14),
                      axis.title=element_text(size=16),
                      panel.background = element_rect(fill="white",
                                                      colour="black"))

        i4 <-ggplot2::ggplot(sig.tab,
                             aes(y=perc.sign.intercept*100,x=percent_sp_removed),
                             environment = environment())+
                scale_y_continuous(limits=c(0,100),breaks=seq(0,100,10))+
                scale_x_continuous(breaks=result$n.percent)+
                xlab("% Species removed")+
                geom_point(size=5,colour="red")+
                geom_line(colour="red")+
                ylab("% of significant Intercept")+
                theme(axis.text=element_text(size=14),
                      axis.title=element_text(size=16),
                      panel.background = element_rect(fill="white",
                                                      colour="black"))

        ### Ploting:
        if (param == "slope" & graphs=="all")
                print(grid.arrange(s1,s2,opt,s4,ncol=2))
        if (param == "slope" & graphs==1)
                print(s1)
        if (param == "slope" & graphs==2)
                print(s2)
        if (param == "slope" & graphs==3)
                print(opt)
        if (param == "slope" & graphs==4)
                print(s4)
        if (param == "intercept" & graphs=="all")
                print(grid.arrange(i1,i2,opt,i4,ncol=2))
        if (param == "intercept" & graphs==1)
                print(i1)
        if (param == "intercept" & graphs==2)
                print(i2)
        if (param == "intercept" & graphs==3)
                print(opt)
        if (param == "intercept" & graphs==4)
                print(i4)
}
