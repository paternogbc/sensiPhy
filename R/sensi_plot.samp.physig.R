#' Graphical diagnostics for class 'samp.physig'
#'
#' \code{plot_samp_phylm} Plot results from \code{samp_physig}
#' @param x output from \code{samp_physig}
#' @param graphs choose which graph should be printed on the output ("all", 1,2 and 3 )
#' @param ... further arguments to methods
#' @importFrom ggplot2 scale_x_continuous scale_colour_manual geom_hline 
#' geom_bar scale_fill_manual scale_y_continuous geom_boxplot geom_line 
#' @author Gustavo Paterno
#' @seealso \code{\link[ggplot2]{ggplot}}, \code{\link[sensiPhy]{samp_phylm}}
#' \code{\link[sensiPhy]{samp_physig}}
#' @details For 'x' from samp_physig:
#' 
#' * Graph 1: Estimated phylogenetic signal for each simution across  
#' percentages of species removed. Colours represent percentage 
#' of change in comparison with the full data estimate (blue = lower than 5, orange = 
#' between 5 and 10 and red = higher than 10).
#' The red horizontal line represents the original phylogenetic signal 
#' estimated from the full model (with all species). 
#' 
#' * Graph 2: The proportion of estimated signal in each category 
#' across the percentage of species removed.
#' 
#' * Graph 3: The percentage of significant signal estimates across the 
#' percentage of species removed. 
#' @md 
#' @export

sensi_plot.samp.physig <- function(x, graphs = "all", ...){
  
result    <- x$sensi.estimates
sig.tab <- x$sign.analysis

# classes of perc:
result$class <- "class"
### Within 5%:
if (length(result[result$perc <= 5 ,]$class) >= 1){
  result[result$perc <= 5,]$class <- "within 5%"
}
### Higher than 5%
if (length(result[result$perc > 5
                  & result$perc <= 10 ,]$class) >= 1){
  result[result$perc > 5
         & result$perc <= 10 ,]$class <- "higher than 5%"
}
### Higher than 10%
if (length(result[result$perc > 10,]$class) >= 1){
  result[result$perc > 10,]$class <- "higher than 10%"
}

result$class <- as.factor(result$class)
e.0    <- as.numeric(x$full.data.estimates[[1]])
e.5    <- .05*e.0
e.10   <- .1*e.0

# reverting the order of the levels
result$class =
  with(result, factor(class,
                      levels = rev(levels(result$class))))

## Organizing colours
if(length(levels(result$class)) == 3){
  colS = c("skyblue","orange","red2")
}
if(length(levels(result$class)) == 2){
  colS = c("skyblue","orange")
}
if(length(levels(result$class)) == 1){
  colS = c("skyblue")
}

### Graphs--------------------------------------------------------------

### Estimated across n.percent:
method <- x$call$method
if(is.null(x$call$method)) method <- "K"

s1 <- ggplot2::ggplot(result,aes(y=estimate,x=n.percent,
                                 colour=class),
                      environment = parent.frame())+
  
  geom_point(size=4,position = "jitter",alpha=.5)+
  scale_x_continuous(breaks=result$n.percent)+
  ylab(paste("Estimated", method, sep = " "))+
  xlab("% of Species Removed ")+
  scale_colour_manual(values=colS)+
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
n.perc.times <- as.numeric(table(result$n.percent))
perc <- with(result,aggregate(data=result,estimate ~ class*n.percent,FUN=length))
a <- colSums(table(perc$class,perc$n.percent))
perc$estimate <- (perc$estimate/rep(n.perc.times,
                                    times=a))*100
s2 <- ggplot2::ggplot(perc,
                      aes(y=estimate,x=n.percent,
                          fill=factor(class)),
                      environment = parent.frame())+
  geom_bar(stat="identity",alpha=.5)+
  scale_fill_manual(values=colS,name=paste("Change in ", method, sep = ""))+
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
  ylab(paste("% of estimated ", method))

#### Graph significance
s3 <-ggplot2::ggplot(sig.tab,
                     aes(y=perc.sign*100,x=percent_sp_removed),
                     environment = parent.frame())+
  scale_y_continuous(limits=c(0,100),breaks=seq(0,100,10))+
  scale_x_continuous(breaks=result$n.percent)+
  xlab("% Species removed")+
  geom_point(size=5,colour="red")+
  geom_line(colour="red")+
  ylab(paste("% of significant", method))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        panel.background = element_rect(fill="white",
                                        colour="black"))

### Export two graphs:
if (graphs == 1) 
  suppressMessages(return(s1))
if (graphs == 2) 
  suppressMessages(return(s2))
if (graphs == 3) 
  suppressMessages(return(s3))
if (graphs == "all")
  suppressMessages(return(multiplot(s1,s3,s2, cols = 2)))
}