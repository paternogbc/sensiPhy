#' Graphical diagnostics for class 'sensiSamp'
#'
#' \code{plot_samp_phylm} Plot results from \code{samp_phylm} and
#' \code{influ_phyloglm}
#' @param x output from \code{samp_phylm}
#' @param graphs choose which graph should be printed on the output ("all", 1,2,3 or 4)
#' @param param choose which model parameter should be ploted  ("intercept" or "estimate")
#' @param ... further arguments to methods
#' @importFrom ggplot2 scale_x_continuous scale_colour_manual geom_hline 
#' geom_bar scale_fill_manual scale_y_continuous geom_boxplot geom_line 
#' @author Gustavo Paterno
#' @seealso \code{\link[ggplot2]{ggplot}}, \code{\link[sensiPhy]{samp_phylm}}
#' \code{\link[sensiPhy]{samp_phyglm}}
#' @details For 'x' from samp_phylm or samp_phyglm:
#' 
#' \strong{Graph 1:} Estimated slopes or intercepts for each simulation across  
#' percentages of species removed. Colours represent percentage 
#' of change in comparison with the full model (blue = lower than 5, orange = 
#' between 5 and 10 and red = higher than 10).
#' The red horizontal line represents the original slope or 
#' intercept from the full model (with all species). 
#' 
#' \strong{Graph 2:} The proportion of estimated slopes and intercepts in each category 
#' across the percentage of species removed.
#' 
#' \strong{Graph 3:} Estimated phylogenetic model parameter for each simulation across
#' the percentage of species removed.
#' 
#' \strong{Graph 4:} The percentage of significant slopes or intercepts across the 
#' percentage of species removed.  
#' 
#' @note If model = "BM", only plots 1, 2 and 4 are printed. Plot 3, phylogenetic
#'  model parameter is not available for model = "BM"
#' @export

sensi_plot.sensiSamp <- function(x, graphs = "all", param = "estimate", ...)
{

   # x <- samp
# nulling variables:
estimate <- n.percent <- estimate.class <- intercept <- model <- intercept.class <- NULL
optpar <- perc.sign.estimate <- percent_sp_removed <- perc.sign.intercept <- NULL

        result    <- x$sensi.estimates
        sig.tab <- x$sign.analysis
        
        # classes of slope.perc:
        result$estimate.class <- "class"
        ### Within 5%:
        if (length(result[result$estimate.perc <= 5 ,]$estimate.class) >= 1){
                result[result$estimate.perc <= 5,]$estimate.class <- "within 5%"
        }
        ### Higher than 5%
        if (length(result[result$estimate.perc > 5
                & result$estimate.perc <= 10 ,]$estimate.class) >= 1){
                result[result$estimate.perc > 5
                       & result$estimate.perc <= 10 ,]$estimate.class <- "higher than 5%"
        }
        ### Higher than 10%
        if (length(result[result$estimate.perc > 10,]$estimate.class) >= 1){
                result[result$estimate.perc > 10,]$estimate.class <- "higher than 10%"
        }

        result$estimate.class <- as.factor(result$estimate.class)
        estimate.0    <- as.numeric(x$full.model.estimates$coef[2])
        estimate.5    <- .05*estimate.0
        estimate.10   <- .1*estimate.0

        # classes of intercept.perc:
        result$intercept.class <- "class"
        ### Within 5%:
        if (length(result[result$intercept.perc <= 5 ,]$intercept.class) >= 1){
                result[result$intercept.perc <= 5,]$intercept.class <- "within 5%"
        }
        ### Higher than 5%
        if (length(result[result$intercept.perc > 5
                          & result$intercept.perc <= 10 ,]$intercept.class) >= 1){
                result[result$intercept.perc > 5
                       & result$intercept.perc <= 10 ,]$intercept.class <- "higher than 5%"
        }
        ### Higher than 10%
        if (length(result[result$intercept.perc > 10,]$intercept.class) >= 1){
                result[result$intercept.perc > 10,]$intercept.class <- "higher than 10%"
        }

        result$intercept.class <- as.factor(result$intercept.class)
        intercept.0    <- as.numeric(x$full.model.estimates$coef[1])
        intercept.5    <- .05*intercept.0
        intercept.10   <- .1*intercept.0

        # reverting the order of the levels
        result$estimate.class =
                with(result, factor(estimate.class,
                                    levels = rev(levels(result$estimate.class))))
        result$intercept.class =
                with(result, factor(intercept.class,
                                    levels = rev(levels(result$intercept.class))))

        ## Organizing colours: slope
        if(length(levels(result$estimate.class)) == 3){
                colS = c("skyblue","orange","red2")
        }
        if(length(levels(result$estimate.class)) == 2){
                colS = c("skyblue","orange")
        }
        if(length(levels(result$estimate.class)) == 1){
                colS = c("skyblue")
        }
        ## Organizing colours: intercept
        if(length(levels(result$intercept.class)) == 3){
                colI = c("skyblue","orange","red2")
        }
        if(length(levels(result$intercept.class)) == 2){
                colI = c("skyblue","orange")
        }
        if(length(levels(result$intercept.class)) == 1){
                colI = c("skyblue")
        }

        ### Graphs--------------------------------------------------------------

        ### Estimated slopes across n.percent:
        s1 <- ggplot2::ggplot(result,aes(y=estimate,x=n.percent,
                                         colour=estimate.class),
                              environment = parent.frame())+

                geom_point(size=4,position = "jitter",alpha=.5)+
                scale_x_continuous(breaks=result$n.percent)+
                ylab("Estimates")+
                xlab("% of Species Removed ")+
                scale_colour_manual(values=colS)+
                geom_hline(yintercept=estimate.0,linetype=1,color="red",
                           size=1, alpha = .6)+

                geom_hline(yintercept=estimate.0+estimate.5,linetype=2,
                           alpha=.6)+
                geom_hline(yintercept=estimate.0-estimate.5,linetype=2,
                           alpha=.6)+
                geom_hline(yintercept=estimate.0+estimate.10,linetype=2,
                           alpha=.6)+
                geom_hline(yintercept=estimate.0-estimate.10,linetype=2,
                           alpha=.6)+
                theme( legend.position = "none",
                       legend.direction = "horizontal",
                       legend.text=element_text(size=12),
                       legend.title=element_text(size=12),
                       axis.text=element_text(size=12),
                       axis.title=element_text(size=12),
                       legend.key.width=unit(.5,"line"),
                       legend.key.size = unit(.5,"cm"),
                       panel.background = element_rect(fill="white",
                                                       colour="black"))
        ### Estimated intercept across n.percent
        i1<- ggplot2::ggplot(result,aes(y=intercept,x=n.percent,
                                   colour=intercept.class),
                             environment = parent.frame())+

                geom_point(size=4,position = "jitter",alpha=.5)+
                scale_x_continuous(breaks=result$n.percent)+
                ylab("Intercepts")+
                xlab("% of Species Removed ")+
                scale_colour_manual(values=colI)+
                geom_hline(yintercept=intercept.0,linetype=1,color="red",
                           size=1,alpha=.6)+

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
                       legend.text=element_text(size=12),
                       legend.title=element_text(size=12),
                       axis.text=element_text(size=12),
                       axis.title=element_text(size=12),
                       legend.key.width=unit(.5,"line"),
                       legend.key.size = unit(.5,"cm"),
                       panel.background = element_rect(fill="white",
                                                       colour="black"))

        ### Proportion of change.classes across n.percent
        n.perc.times <- as.numeric(table(result$n.percent))
        estimate.perc <- with(result,aggregate(data=result,estimate ~ estimate.class*n.percent,FUN=length))
        a <- colSums(table(estimate.perc$estimate.class,estimate.perc$n.percent))
        estimate.perc$estimate <- (estimate.perc$estimate/rep(n.perc.times,
                                                  times=a))*100
        intercept.perc <- with(result,aggregate(data=result,intercept ~ intercept.class*n.percent,FUN=length))
        b <- colSums(table(intercept.perc$intercept.class,intercept.perc$n.percent))
        intercept.perc$intercept <- (intercept.perc$intercept/rep(n.perc.times,
                                                  times=b))*100
        ### Graph: Slope
        s2 <- ggplot2::ggplot(estimate.perc,
                     aes(y=estimate,x=n.percent,
                         fill=factor(estimate.class)),
                     environment = parent.frame())+
                geom_bar(stat="identity",alpha=.5)+
                scale_fill_manual(values=colS,name="Change in beta")+
                scale_y_continuous(breaks=seq(0,100,10))+
                scale_x_continuous(breaks=result$n.percent)+
                theme( legend.position = "top",
                       legend.direction = "horizontal",
                       legend.text=element_text(size=12),
                       legend.title = element_text(size=12),
                       axis.text=element_text(size=12),
                       axis.title=element_text(size=12),
                       legend.key.width=unit(.5,"line"),
                       legend.key.size = unit(.5,"cm"),
                       panel.background = element_rect(fill="white",
                                                       colour="black"))+
                xlab("% of Species Removed")+
                ylab("% of Estimates")

        ### Graph: Intercept
        i2 <- ggplot2::ggplot(intercept.perc,
                     aes(y=intercept,x=n.percent,
                         fill=factor(intercept.class)),
                     environment = parent.frame())+
                geom_bar(stat="identity",alpha=.5)+
                scale_fill_manual(values=colI,name="Change in beta")+
                scale_y_continuous(breaks=seq(0,100,10))+
                scale_x_continuous(breaks=result$n.percent)+
                theme( legend.position = "top",
                       legend.direction = "horizontal",
                       legend.text=element_text(size=12),
                       legend.title = element_text(size=12),
                       axis.text=element_text(size=12),
                       axis.title=element_text(size=12),
                       legend.key.width=unit(.5,"line"),
                       legend.key.size = unit(.5,"cm"),
                       panel.background = element_rect(fill="white",
                                                       colour="black"))+
                xlab("% of Species Removed")+
                ylab("% of Intercepts")

        ### Optpar acros % removed species:
        opt <- ggplot2::ggplot(result,aes(y=optpar,x=n.percent,group=as.factor(n.percent)))+
                geom_point()+
                geom_boxplot(fill="red",alpha=.5)+
                scale_x_continuous(breaks=result$n.percent)+
                theme(axis.title=element_text(size=12),
                      axis.text = element_text(size=12),
                      panel.background = element_rect(fill="white",
                                                      colour="black"))+
                xlab("% of Species Removed")+
                ylab("Phylogenetic model parameter")

        ## Significance Analysis : p.value of slope

        s4 <-ggplot2::ggplot(sig.tab,
                             aes(y=perc.sign.estimate*100,x=percent_sp_removed),
                             environment = parent.frame())+
                scale_y_continuous(limits=c(0,100),breaks=seq(0,100,10))+
                scale_x_continuous(breaks=result$n.percent)+
                xlab("% Species removed")+
                geom_point(size=5,colour="red")+
                geom_line(colour="red")+
                ylab("% of significant estimates")+
                theme(axis.text=element_text(size=12),
                      axis.title=element_text(size=12),
                      panel.background = element_rect(fill="white",
                                                      colour="black"))

        i4 <-ggplot2::ggplot(sig.tab,
                             aes(y=perc.sign.intercept*100,x=percent_sp_removed),
                             environment = parent.frame())+
                scale_y_continuous(limits=c(0,100),breaks=seq(0,100,10))+
                scale_x_continuous(breaks=result$n.percent)+
                xlab("% Species removed")+
                geom_point(size=5,colour="red")+
                geom_line(colour="red")+
                ylab("% of significant intercepts")+
                theme(axis.text=element_text(size=12),
                      axis.title=element_text(size=12),
                      panel.background = element_rect(fill="white",
                                                      colour="black"))
        
        which_plot(param = param, graphs = graphs,
                   s1 = s1, s2 = s2, s4 = s4, opt = opt,
                   i1 = i1, i2 = i2, i4 = i4, model = x$model)
}


#####Plotting method for sensi_plot.sensiSamp.TraitEvol
#' @export

sensi_plot.sensiSamp.TraitEvol <- function(x, graphs = "all", ...){
  if(as.character(x$call[[1]])=="samp_continuous"){ #Check what type of TraitEvolution is evaluated
    ### Nulling variables:
    estimate.sigsq <- n.percent <- perc.sign <- percent_sp_removed <- NULL
    result.sigsq    <- x$sensi.estimates
    
    # classes of perc:
    result.sigsq$class <- "class"
    ### Within 5%:
    if (length(result.sigsq[result.sigsq$sigsq.perc <= 5 ,]$class) >= 1){
      result.sigsq[result.sigsq$sigsq.perc <= 5,]$class <- "within 5%"
    }
    ### Higher than 5%
    if (length(result.sigsq[result.sigsq$sigsq.perc> 5
                            & result.sigsq$sigsq.perc <= 10 ,]$class) >= 1){
      result.sigsq[result.sigsq$sigsq.perc > 5
                   & result.sigsq$sigsq.perc <= 10 ,]$class <- "higher than 5%"
    }
    ### Higher than 10%
    if (length(result.sigsq[result.sigsq$sigsq.perc > 10,]$class) >= 1){
      result.sigsq[result.sigsq$sigsq.perc > 10,]$class <- "higher than 10%"
    }
    
    result.sigsq$class <- as.factor(result.sigsq$class)
    e.0    <- as.numeric(x$full.model.estimates$sigsq)
    e.5    <- .05*e.0
    e.10   <- .1*e.0
    
    # reverting the order of the levels
    result.sigsq$class =
      with(result.sigsq, factor(class,
                                levels = rev(levels(result.sigsq$class))))
    
    ## Organizing colours
    if(length(levels(result.sigsq$class)) == 3){
      colS.sigsq = c("skyblue","orange","red2")
    }
    if(length(levels(result.sigsq$class)) == 2){
      colS.sigsq = c("skyblue","orange")
    }
    if(length(levels(result.sigsq$class)) == 1){
      colS.sigsq = c("skyblue")
    }
    
    ### Graphs--------------------------------------------------------------
    
    ### Estimated across n.percent:
    s1 <- ggplot2::ggplot(result.sigsq,aes(y=sigsq,x=n.percent,
                                           colour=class),
                          environment = parent.frame())+
      
      geom_point(size=4,position = "jitter",alpha=.5)+
      scale_x_continuous(breaks=result.sigsq$n.percent)+
      ylab("Estimated sigsq")+
      xlab("% of Species Removed ")+
      scale_colour_manual(values=colS.sigsq)+
      geom_hline(yintercept=e.0,linetype=1,color="red",
                 size=1, alpha = .6)+
      
      geom_hline(yintercept=e.0+e.5,linetype=2,
                 alpha=.6)+
      geom_hline(yintercept=e.0-e.5,linetype=2,
                 alpha=.6)+
      geom_hline(yintercept=e.0+e.10,linetype=2,
                 alpha=.6)+
      geom_hline(yintercept=e.0-e.10,linetype=2,
                 alpha=.6)+
      theme( legend.position = "none",
             legend.direction = "horizontal",
             legend.text=element_text(size=12),
             legend.title=element_text(size=12),
             axis.text=element_text(size=12),
             axis.title=element_text(size=12),
             legend.key.width=unit(.5,"line"),
             legend.key.size = unit(.5,"cm"),
             panel.background = element_rect(fill="white",
                                             colour="black"))
    
    
    ### Graph2
    ### Proportion of change.classes across n.percent
    n.perc.times <- as.numeric(table(result.sigsq$n.percent))
    perc <- with(result.sigsq,aggregate(data=result.sigsq,sigsq ~ class*n.percent,FUN=length))
    a <- colSums(table(perc$class,perc$n.percent))
    perc$sigsq <- (perc$sigsq/rep(n.perc.times,
                                  times=a))*100
    s2 <- ggplot2::ggplot(perc,
                          aes(y=sigsq,x=n.percent,
                              fill=factor(class)),
                          environment = parent.frame())+
      geom_bar(stat="identity",alpha=.5)+
      scale_fill_manual(values=colS.sigsq,name="Change in sigsq")+
      scale_y_continuous(breaks=seq(0,100,10))+
      scale_x_continuous(breaks=result.sigsq$n.percent)+
      theme( legend.position = "top",
             legend.direction = "horizontal",
             legend.text=element_text(size=12),
             legend.title = element_text(size=12),
             axis.text=element_text(size=12),
             axis.title=element_text(size=12),
             legend.key.width=unit(.5,"line"),
             legend.key.size = unit(.5,"cm"),
             panel.background = element_rect(fill="white",
                                             colour="black"))+
      xlab("% of Species Removed")+
      ylab(paste("% of estimated sigsq change"))
    
    ### Export two graphs:
    if (graphs == 1) 
      suppressMessages(return(s1))
    if (graphs == 2) 
      suppressMessages(return(s2))
    if (graphs == "all")
      suppressMessages(return(multiplot(s1,s2, cols = 2)))
  }
  ##When a samp_discrete object has been called
  if(as.character(x$call[[1]])=="samp_discrete"){ #Check what type of TraitEvolution is evaluated
    ### Nulling variables:
    estimate.q12 <- n.percent <- perc.sign <- percent_sp_removed <- NULL
    result.q12    <- x$sensi.estimates
    
    # classes of perc:
    result.q12$class <- "class"
    ### Within 5%:
    if (length(result.q12[result.q12$q12.perc <= 5 ,]$class) >= 1){
      result.q12[result.q12$q12.perc <= 5,]$class <- "within 5%"
    }
    ### Higher than 5%
    if (length(result.q12[result.q12$q12.perc> 5
                          & result.q12$q12.perc <= 10 ,]$class) >= 1){
      result.q12[result.q12$q12.perc > 5
                 & result.q12$q12.perc <= 10 ,]$class <- "higher than 5%"
    }
    ### Higher than 10%
    if (length(result.q12[result.q12$q12.perc > 10,]$class) >= 1){
      result.q12[result.q12$q12.perc > 10,]$class <- "higher than 10%"
    }
    
    result.q12$class <- as.factor(result.q12$class)
    e.0    <- as.numeric(x$full.model.estimates$q12)
    e.5    <- .05*e.0
    e.10   <- .1*e.0
    
    # reverting the order of the levels
    result.q12$class =
      with(result.q12, factor(class,
                              levels = rev(levels(result.q12$class))))
    
    ## Organizing colours
    if(length(levels(result.q12$class)) == 3){
      colS.q12 = c("skyblue","orange","red2")
    }
    if(length(levels(result.q12$class)) == 2){
      colS.q12 = c("skyblue","orange")
    }
    if(length(levels(result.q12$class)) == 1){
      colS.q12 = c("skyblue")
    }
    
    ### Graphs--------------------------------------------------------------
    
    ### Estimated across n.percent:
    s1 <- ggplot2::ggplot(result.q12,aes(y=q12,x=n.percent,
                                         colour=class),
                          environment = parent.frame())+
      
      geom_point(size=4,position = "jitter",alpha=.5)+
      scale_x_continuous(breaks=result.q12$n.percent)+
      ylab("Estimated q12")+
      xlab("% of Species Removed ")+
      scale_colour_manual(values=colS.q12)+
      geom_hline(yintercept=e.0,linetype=1,color="red",
                 size=1, alpha = .6)+
      
      geom_hline(yintercept=e.0+e.5,linetype=2,
                 alpha=.6)+
      geom_hline(yintercept=e.0-e.5,linetype=2,
                 alpha=.6)+
      geom_hline(yintercept=e.0+e.10,linetype=2,
                 alpha=.6)+
      geom_hline(yintercept=e.0-e.10,linetype=2,
                 alpha=.6)+
      theme( legend.position = "none",
             legend.direction = "horizontal",
             legend.text=element_text(size=12),
             legend.title=element_text(size=12),
             axis.text=element_text(size=12),
             axis.title=element_text(size=12),
             legend.key.width=unit(.5,"line"),
             legend.key.size = unit(.5,"cm"),
             panel.background = element_rect(fill="white",
                                             colour="black"))
    
    
    ### Graph2
    ### Proportion of change.classes across n.percent
    n.perc.times <- as.numeric(table(result.q12$n.percent))
    perc <- with(result.q12,aggregate(data=result.q12,q12 ~ class*n.percent,FUN=length))
    a <- colSums(table(perc$class,perc$n.percent))
    perc$q12 <- (perc$q12/rep(n.perc.times,
                              times=a))*100
    s2 <- ggplot2::ggplot(perc,
                          aes(y=q12,x=n.percent,
                              fill=factor(class)),
                          environment = parent.frame())+
      geom_bar(stat="identity",alpha=.5)+
      scale_fill_manual(values=colS.q12,name="Change in q12")+
      scale_y_continuous(breaks=seq(0,100,10))+
      scale_x_continuous(breaks=result.q12$n.percent)+
      theme( legend.position = "top",
             legend.direction = "horizontal",
             legend.text=element_text(size=12),
             legend.title = element_text(size=12),
             axis.text=element_text(size=12),
             axis.title=element_text(size=12),
             legend.key.width=unit(.5,"line"),
             legend.key.size = unit(.5,"cm"),
             panel.background = element_rect(fill="white",
                                             colour="black"))+
      xlab("% of Species Removed")+
      ylab(paste("% of estimated q12 change"))
    
    ####Redo the exact same thing for the q21
    ### Nulling variables:
    estimate.q21 <- n.percent <- perc.sign <- percent_sp_removed <- NULL
    result.q21    <- x$sensi.estimates
    
    # classes of perc:
    result.q21$class <- "class"
    ### Within 5%:
    if (length(result.q21[result.q21$q21.perc <= 5 ,]$class) >= 1){
      result.q21[result.q21$q21.perc <= 5,]$class <- "within 5%"
    }
    ### Higher than 5%
    if (length(result.q21[result.q21$q21.perc> 5
                          & result.q21$q21.perc <= 10 ,]$class) >= 1){
      result.q21[result.q21$q21.perc > 5
                 & result.q21$q21.perc <= 10 ,]$class <- "higher than 5%"
    }
    ### Higher than 10%
    if (length(result.q21[result.q21$q21.perc > 10,]$class) >= 1){
      result.q21[result.q21$q21.perc > 10,]$class <- "higher than 10%"
    }
    
    result.q21$class <- as.factor(result.q21$class)
    e.0    <- as.numeric(x$full.model.estimates$q21)
    e.5    <- .05*e.0
    e.10   <- .1*e.0
    
    # reverting the order of the levels
    result.q21$class =
      with(result.q21, factor(class,
                              levels = rev(levels(result.q21$class))))
    
    ## Organizing colours
    if(length(levels(result.q21$class)) == 3){
      colS.q21 = c("skyblue","orange","red2")
    }
    if(length(levels(result.q21$class)) == 2){
      colS.q21 = c("skyblue","orange")
    }
    if(length(levels(result.q21$class)) == 1){
      colS.q21 = c("skyblue")
    }
    
    ### Graphs--------------------------------------------------------------
    
    ### Estimated across n.percent:
    s3 <- ggplot2::ggplot(result.q21,aes(y=q21,x=n.percent,
                                         colour=class),
                          environment = parent.frame())+
      
      geom_point(size=4,position = "jitter",alpha=.5)+
      scale_x_continuous(breaks=result.q21$n.percent)+
      ylab("Estimated q21")+
      xlab("% of Species Removed ")+
      scale_colour_manual(values=colS.q21)+
      geom_hline(yintercept=e.0,linetype=1,color="red",
                 size=1, alpha = .6)+
      
      geom_hline(yintercept=e.0+e.5,linetype=2,
                 alpha=.6)+
      geom_hline(yintercept=e.0-e.5,linetype=2,
                 alpha=.6)+
      geom_hline(yintercept=e.0+e.10,linetype=2,
                 alpha=.6)+
      geom_hline(yintercept=e.0-e.10,linetype=2,
                 alpha=.6)+
      theme( legend.position = "none",
             legend.direction = "horizontal",
             legend.text=element_text(size=12),
             legend.title=element_text(size=12),
             axis.text=element_text(size=12),
             axis.title=element_text(size=12),
             legend.key.width=unit(.5,"line"),
             legend.key.size = unit(.5,"cm"),
             panel.background = element_rect(fill="white",
                                             colour="black"))
    
    
    ### Graph2
    ### Proportion of change.classes across n.percent
    n.perc.times <- as.numeric(table(result.q21$n.percent))
    perc <- with(result.q21,aggregate(data=result.q21,q21 ~ class*n.percent,FUN=length))
    a <- colSums(table(perc$class,perc$n.percent))
    perc$q21 <- (perc$q21/rep(n.perc.times,
                              times=a))*100
    s4 <- ggplot2::ggplot(perc,
                          aes(y=q21,x=n.percent,
                              fill=factor(class)),
                          environment = parent.frame())+
      geom_bar(stat="identity",alpha=.5)+
      scale_fill_manual(values=colS.q21,name="Change in q21")+
      scale_y_continuous(breaks=seq(0,100,10))+
      scale_x_continuous(breaks=result.q21$n.percent)+
      theme( legend.position = "top",
             legend.direction = "horizontal",
             legend.text=element_text(size=12),
             legend.title = element_text(size=12),
             axis.text=element_text(size=12),
             axis.title=element_text(size=12),
             legend.key.width=unit(.5,"line"),
             legend.key.size = unit(.5,"cm"),
             panel.background = element_rect(fill="white",
                                             colour="black"))+
      xlab("% of Species Removed")+
      ylab(paste("% of estimated q21 change"))
    
    ### Export four graphs:
    if (graphs == 1) 
      suppressMessages(return(s1))
    if (graphs == 2) 
      suppressMessages(return(s2))
    if (graphs == 3) 
      suppressMessages(return(s3))
    if (graphs == 4) 
      suppressMessages(return(s4))
    if (graphs == "all"){
      if(x$optpar !="BM"){
        suppressMessages(return(multiplot(s1,s3,s2,s4, cols = 2)))
      } else
        suppressMessages(return(multiplot(s1,s2, cols = 2)))
    }
  }
}

