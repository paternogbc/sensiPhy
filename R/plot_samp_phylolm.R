if (x[[1]] == "samp_pgls"){

        x <- samp2
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
        ## Graphs: Estimated slopes ~ % species removed
        if(length(levels(result$slope.class)) == 3){
                col = c("skyblue","orange","red2")
        }
        if(length(levels(result$slope.class)) == 2){
                col = c("skyblue","orange")
        }
        if(length(levels(result$slope.class)) == 1){
                col = c("skyblue")
        }

        p1 <- ggplot2::ggplot(result,aes(y=slope,x=n.percent,
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

        ## Significance Analysis : p.value of slope

        p4 <-ggplot2::ggplot(sig.tab,
                             aes(y=perc.sign.slope,x=percent_sp_removed),
                             environment = environment())+
                scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.1))+
                scale_x_continuous(breaks=result$n.percent)+
                xlab("% Species removed")+
                geom_point(size=5,colour="red")+
                geom_line(colour="red")+
                ylab("Power  [p-value]")+
                theme(axis.text=element_text(size=14),
                      axis.title=element_text(size=16),
                      panel.background = element_rect(fill="white",
                                                      colour="black"))


        p2 <- ggplot(beta.tab,
                     aes(y=proportion,x=n.percents,
                         fill=factor(beta.class)),
                     environment = environment())+
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
        p3 <- ggplot(m.DFbetas,aes(y=mDFbetas,x=n.percents),
                     environment = environment())+
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
