#' Graphical diagnostics for class 'sensiIntra'
#'
#' \code{plot_tree.intra_phylm} Plot results from \code{tree_phylm},
#' \code{intra_phylm} and \code{intra_phyglm}
#' @param x output from \code{tree_phylm}, \code{tree_phyglm},
#' \code{intra_phylm} or \code{intra_phyglm}
#' @param graphs choose which graph should be printed in the output ("all", 1, 2 or 3)
#' @param ... further arguments to methods
#' @importFrom ggplot2 scale_color_manual geom_histogram geom_abline geom_density 
#' geom_vline xlab geom_point theme
#' @author Caterina Penone and Gustavo Paterno
#' @seealso \code{\link[ggplot2]{ggplot}}, \code{\link[sensiPhy]{tree_phylm}}
#' \code{\link[sensiPhy]{intra_phylm}}
#' @details For 'x' from \code{tree_phylm}, \code{tree_phyglm},
#' \code{intra_phylm} or \code{intra_phyglm}:
#' 
#' \strong{Graphs 1 and 2:} Distribution of estimated slopes (estimates) and intercepts for each tree (for \code{tree_phylm}) or 
#' value generated within a given interval (\code{intra_phylm})
#' Red vertical line represents the mean estimate or intercept for all models. 
#'
#' \strong{Graph 3:} Scatterplot with mean regression (black line) and standard deviation of the regression (red dotted lines).
#' @importFrom grid unit 
#' @importFrom stats plogis
#' @export

sensi_plot.sensiIntra <- function(x, graphs="all", ...){
  
    # nulling variables
    formula <- estimate <- ..density.. <- intercept <- NULL
    s3 <- yy <- linety <- pval.estimate <- NULL

    mappx <- x$formula[[3]]
    mappy <- x$formula[[2]]
    full.data <- x$data[all.vars(x$formula)]
    
    if(!is.null(x$y.transf))
    {full.data[,1] <- x$y.transf(full.data[,1])}
    
    if(!is.null(x$x.transf))
    {full.data[,2] <- x$x.transf(full.data[,2])}
    
    result <- x$sensi.estimates
    statm<- x$all.stats
    estimate.0 <-  as.numeric(statm[4,3])
    intercept.0 <-  as.numeric(statm[1,3])
    sensi.estimates<-x$sensi.estimates
    
    xf <- model.frame(formula = x$formula, data = full.data)[,2]
    yf <- plogis(intercept.0 + estimate.0 * xf)
    yp <- plogis((statm[1,6]) + (statm[4,6] * xf)) 
    ym <- plogis((statm[1,5]) + (statm[4,5] * xf))
    
    plot_data <- data.frame("xf" = c(xf,xf,xf),
                            "yy" = c(yf, yp, ym),
                            linety = rep(c("Mean","High","Low"),each = length(yf)))

    #Distribution of estimated estimates:
    s1 <- ggplot2::ggplot(sensi.estimates,aes(x=estimate),
                          environment = parent.frame())+
      geom_histogram(fill="yellow", colour="black", size=.2, alpha = .3) +
      geom_vline(xintercept = estimate.0,color="red",linetype=2,size=.7)+
      xlab("Estimates")+
      ylab("Frequency")+
        theme(axis.title=element_text(size=12),
              axis.text = element_text(size=12),
              panel.background = element_rect(fill="white",
                                              colour="black"))

    #Distribution of estimated intercepts:
    i1 <- ggplot2::ggplot(sensi.estimates,aes(x=intercept),
                          environment = parent.frame())+
      geom_histogram(fill="yellow", colour="black", size=.2, alpha = .3) +
      geom_vline(xintercept = intercept.0,color="red",linetype=2,size=.7)+
      xlab("Estimated Intercepts")+
      ylab("Frequency")+
        theme(axis.title=element_text(size=12),
              axis.text = element_text(size=12),
              panel.background = element_rect(fill="white",
                                              colour="black"))
 
    #third plot: data visualisation
    s2 <- ggplot2::ggplot(full.data, aes_string(y = mappy, x = mappx),
                          environment = parent.frame())+
      geom_point(size=3,alpha=.8)+
      theme(panel.background=element_rect(fill="white",colour="black"),
            axis.title=element_text(size=12),
            axis.text = element_text(size=12),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position="none")

    if(length(class(x)) == 1){
      s2.out <- s2 +
        geom_abline(intercept = intercept.0, slope = estimate.0)+
        geom_abline(intercept = statm[1,5], 
                    slope     = statm[4,5], colour = "red",
                    linetype  = 2, size = 1, show.legend = T)+
        geom_abline(intercept = statm[1,6],
                    slope     = statm[4,6], colour = "red",
                    linetype  = 2, size = 1, show.legend = T)
    }
    
    
    if(length(class(x)) == 2){
      s2.out <- s2 + geom_line(data = plot_data, aes(x = xf, y = yy, linetype=linety, col=linety),size=1)+
        scale_linetype_manual(values = c(2,2,1)) +
        scale_color_manual("",values = c("red","red","black"),guide=F)
                   
    }
    
    #Distribution of p-values (estimate)
    p1 <- ggplot2::ggplot(sensi.estimates,aes(x=pval.estimate),
                          environment = parent.frame())+
        geom_histogram(fill="yellow", colour="black", size=.2, alpha = .3) +
        xlab("Distribution of P-values")+
        ylab("Frequency")+
        theme(axis.title=element_text(size=12),
              axis.text = element_text(size=12),
              panel.background = element_rect(fill="white",
                                              colour="black"))
    
    
    
    ### Plotting:
    if (graphs=="all")
      suppressMessages(return(multiplot(s1,s2.out,i1,p1, cols=2)))
    if (graphs==1)
      suppressMessages(return(s1))
    if (graphs==2)
      suppressMessages(return(i1))
    if (graphs==3)
      suppressMessages(return(s2.out))
    if (graphs==4)
        suppressMessages(return(p1))

}

#' Graphical diagnostics for class 'sensiTree'
#'
#' \code{plot_tree.intra_phylm} Plot results from \code{tree_phylm},
#' \code{intra_phylm} and \code{intra_phyglm}
#' @param x output from \code{tree_phylm}, \code{tree_phyglm},
#' \code{intra_phylm} or \code{intra_phyglm}
#' @param graphs choose which graph should be printed in the output ("all", 1, 2 or 3)
#' @param ... further arguments to methods
#' @importFrom ggplot2 scale_color_manual geom_histogram geom_abline geom_density 
#' geom_vline xlab geom_point theme
#' @author Caterina Penone and Gustavo Paterno
#' @seealso \code{\link[ggplot2]{ggplot}}, \code{\link[sensiPhy]{tree_phylm}}
#' \code{\link[sensiPhy]{intra_phylm}}
#' @details For 'x' from \code{tree_phylm}, \code{tree_phyglm},
#' \code{intra_phylm} or \code{intra_phyglm}:
#' 
#' \strong{Graphs 1 and 2:} Distribution of estimated estimates and intercepts for each tree (for \code{tree_phylm}) or 
#' value generated within a given interval (\code{intra_phylm})
#' Red vertical line represents the mean estimate or intercept for all models. 
#'
#' \strong{Graph 3:} Scatterplot with mean regression (black line) and standard deviation of the regression (blue dotted lines).
#' 
#' @importFrom grid unit 
#' @export
sensi_plot.sensiTree <- function(x, graphs = "all", ...){
    sensi_plot.sensiIntra(x, graphs, ...)
}


#' Graphical diagnostics for class 'sensiIntra_Tree'
#'
#' \code{plot_tree.intra_phylm} Plot results from \code{tree_intra_phylm},
#' \code{tree_intra_phylm} and \code{tree_intra_phyglm}
#' @param x output from \code{tree_intra_phylm}, \code{tree_intra_phyglm}
#' @param graphs choose which graph should be printed in the output ("all", 1, 2 or 3)
#' @param uncer.type chosse which uncertainty type should be printed ("all", "intra", "tree")
#' @param ... further arguments to methods
#' @importFrom ggplot2 scale_color_manual geom_histogram geom_abline geom_density 
#' geom_vline xlab geom_point theme
#' @author Caterina Penone and Gustavo Paterno
#' @seealso \code{\link[ggplot2]{ggplot}}, \code{\link[sensiPhy]{tree_phylm}}
#' \code{\link[sensiPhy]{intra_phylm}}
#' @details For 'x' from \code{tree_intra_phylm} or \code{tree_intra_phyglm}:
#' 
#' \strong{Graphs 1 and 2:} Distribution of estimated estimates and intercepts for each tree (for \code{tree_phylm}) or 
#' value generated within a given interval (\code{interaction_intra_tree_phylm})
#' Red vertical line represents the mean estimate or intercept for all models. 
#'
#' \strong{Graph 3:} Scatterplot with mean regression (black line) and standard deviation of the regression (blue dotted lines).
#' @importFrom grid unit 
#' @importFrom stats plogis
#' @export

sensi_plot.sensiTree_Intra <- function(x, graphs="all", uncer.type = "all",...){
  
  # nulling variables
  formula <- estimate <- ..density.. <- intercept <- NULL
  s3 <- yy <- linety <- pval.estimate <- NULL
  
  mappx <- x$formula[[3]]
  mappy <- x$formula[[2]]
  full.data <- x$data[all.vars(x$formula)]
  
  if(!is.null(x$y.transf))
  {full.data[,1] <- x$y.transf(full.data[,1])}
  
  if(!is.null(x$x.transf))
  {full.data[,2] <- x$x.transf(full.data[,2])}
  
  result <- x$sensi.estimates
  if(uncer.type == "all") {statm <- x$all.stats[,1:6]}
  if(uncer.type == "intra") {statm <- x$all.stats[,7:12]}
  if(uncer.type == "tree") {statm <- x$all.stats[,13:18]}
  
  estimate.0 <-  as.numeric(statm[4,3])
  intercept.0 <-  as.numeric(statm[1,3])
  
  if(uncer.type == "all") {sensi.estimates <- x$sensi.estimates}
  if(uncer.type == "intra") {sensi.estimates <- stats::aggregate(.~n.intra, data = x$sensi.estimates,mean)}
  if(uncer.type == "tree") {sensi.estimates <- stats::aggregate(.~n.tree, data = x$sensi.estimates,mean)}

  xf <- model.frame(formula = x$formula, data = full.data)[,2]
  yf <- plogis(intercept.0 + estimate.0 * xf)
  yp <- plogis((statm[1,6]) + (statm[4,6] * xf)) 
  ym <- plogis((statm[1,5]) + (statm[4,5] * xf))
  
  plot_data <- data.frame("xf" = c(xf,xf,xf),
                          "yy" = c(yf, yp, ym),
                          linety = rep(c("Mean","High","Low"),each = length(yf)))
  
  #Distribution of estimated estimates:
  s1 <- ggplot2::ggplot(sensi.estimates,aes(x=estimate),
                        environment = parent.frame())+
    geom_histogram(fill="yellow", colour="black", size=.2, alpha = .3) +
    geom_vline(xintercept = estimate.0,color="red",linetype=2,size=.7)+
    xlab("Estimated estimates")+
    ylab("Frequency")+
    theme(axis.title=element_text(size=12),
          axis.text = element_text(size=12),
          panel.background = element_rect(fill="white",
                                          colour="black"))
  
  #Distribution of estimated intercepts:
  i1 <- ggplot2::ggplot(sensi.estimates,aes(x=intercept),
                        environment = parent.frame())+
    geom_histogram(fill="yellow", colour="black", size=.2, alpha = .3) +
    geom_vline(xintercept = intercept.0,color="red",linetype=2,size=.7)+
    xlab("Estimated Intercepts")+
    ylab("Frequency")+
    theme(axis.title=element_text(size=12),
          axis.text = element_text(size=12),
          panel.background = element_rect(fill="white",
                                          colour="black"))
  
  #third plot: data visualisation
  s2 <- ggplot2::ggplot(full.data, aes_string(y = mappy, x = mappx),
                        environment = parent.frame())+
    geom_point(size=3,alpha=.8)+
    theme(panel.background=element_rect(fill="white",colour="black"),
          axis.title=element_text(size=12),
          axis.text = element_text(size=12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position="none")
  
  if(length(class(x)) == 1){
    s2.out <- s2 +
      geom_abline(intercept = intercept.0, slope = estimate.0)+
      geom_abline(intercept = statm[1,5], 
                  slope     = statm[4,5], colour = "red",
                  linetype  = 2, size = 1, show.legend = T)+
      geom_abline(intercept = statm[1,6],
                  slope     = statm[4,6], colour = "red",
                  linetype  = 2, size = 1, show.legend = T)
  }
  
  
  if(length(class(x)) == 2){
    s2.out <- s2 + geom_line(data = plot_data, aes(x = xf, y = yy, linetype=linety, col=linety),size=1)+
      scale_linetype_manual(values = c(2,2,1)) +
      scale_color_manual("",values = c("red","red","black"),guide=F)
    
  }
  
  #Distribution of p-values (estimate)
  p1 <- ggplot2::ggplot(sensi.estimates,aes(x=pval.estimate),
                        environment = parent.frame())+
    geom_histogram(fill="yellow", colour="black", size=.2, alpha = .3) +
    xlab("Distribution of P-values")+
    ylab("Frequency")+
    theme(axis.title=element_text(size=12),
          axis.text = element_text(size=12),
          panel.background = element_rect(fill="white",
                                          colour="black"))
  
  
  
  ### Plotting:
  if (graphs=="all")
    suppressMessages(return(multiplot(s1,s2.out,i1,p1, cols=2)))
  if (graphs==1)
    suppressMessages(return(s1))
  if (graphs==2)
    suppressMessages(return(i1))
  if (graphs==3)
    suppressMessages(return(s2.out))
  if (graphs==4)
    suppressMessages(return(p1))
  
}




#########
###Still to fully write Check sensi.Intra for general strategy

sensi_plot.sensiTree.TraitEvol <- function(x, graphs="all", ...){
  q12_fig<-ggplot()+
    geom_histogram(aes(x$sensi.estimates$q12))+
    geom_vline(aes(xintercept=mean(x$sensi.estimates$q12)),colour="red")+
    geom_vline(aes(xintercept=median(x$sensi.estimates$q12)),colour="blue")+
    theme_minimal()
  q21_fig<-ggplot()+
    geom_histogram(aes(x$sensi.estimates$q21))+
    geom_vline(aes(xintercept=mean(x$sensi.estimates$q21)),colour="red")+
    geom_vline(aes(xintercept=median(x$sensi.estimates$q21)),colour="blue")+
    theme_minimal()
  aicc_fig<-ggplot()+
    geom_histogram(aes(x$sensi.estimates$aicc))+
    geom_vline(aes(xintercept=mean(x$sensi.estimates$aicc)),colour="red")+
    geom_vline(aes(xintercept=mediean(x$sensi.estimates$aicc)),colour="blue")+
    theme_minimal()
  optpar_fig<-ggplot()+
    geom_histogram(aes(x$sensi.estimates$optpar))+
    geom_vline(aes(xintercept=mean(x$sensi.estimates$optpar)),colour="red")+
    geom_vline(aes(xintercept=mediean(x$sensi.estimates$optpar)),colour="blue")+
    theme_minimal()
  
  if (graphs=="all")
    suppressMessages(return(multiplot(q12_fig,q21_fig,aicc_fig,optpar_fig, cols=2)))
}
