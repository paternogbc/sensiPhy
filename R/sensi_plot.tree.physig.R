#' Graphical diagnostics for class 'tree.physig'
#'
#' \code{sensi_plot_tree.physig} Plot results from \code{tree_physig}.
#' @param x output from \code{tree_physig}
#' @param graphs choose which graph should be printed in the output ("all", 1 or 2)
#' @param ... further arguments to methods
#' @importFrom ggplot2 scale_color_manual geom_histogram geom_abline geom_density 
#' geom_vline xlab geom_point theme
#' @author Caterina Penone and Gustavo Paterno
#' @seealso \code{\link[ggplot2]{ggplot}}, \code{\link[sensiPhy]{tree_phylm}}
#' \code{\link[sensiPhy]{intra_phylm}}
#' @details For 'x' from \code{tree_physig}
#' 
#' \strong{Graphs 1:} Distribution of estimated phylosgentic signal (lambda or K) for each tree
#' Red vertical line represents the average signal among all estimates. 
#'
#' \strong{Graph 2:} Distribution of p-values for the phylogenetic signal (K or lambda) 
#' for each tree. Red vertical line represents the alpha significance level = 0.05.
#' @importFrom grid unit 
#' @importFrom stats plogis
#' @export
sensi_plot.tree.physig <- function(x, graphs="all", ...){
  
  ### Basic checking:
  method <- x$call$method
  if(is.null(x$call$method)) method <- "K"
  ### Distribution of K values estimated
  e1 <- ggplot2::ggplot(x$tree.physig.estimates, aes(x=estimate))+
    geom_histogram(fill="yellow",colour="black", size=.2,
                   alpha = .3) +
    geom_vline(xintercept = x$stats$mean[1], color="red",linetype=2,size=.7)+
    xlab(paste("Estimated", method, "values", sep = " "))+
    ylab("Frequency")+
    theme(axis.title=element_text(size=12),
          axis.text = element_text(size=12),
          panel.background = element_rect(fill="white",
                                          colour="black"))
  
  ### Distribution of Values estimated
  e2 <- ggplot2::ggplot(x$tree.physig.estimates, aes(x = pval))+
    geom_histogram(fill="yellow",colour="black", size=.2,
                   alpha = .3) +
    geom_vline(xintercept = 0.05,color="red",linetype=1,size=.7)+
    #scale_x_continuous(limits = c(0,1), breaks = seq(0,1,.1))+
    xlab("Estimated P values")+
    ylab("Frequency")+
    theme(axis.title=element_text(size=12),
          axis.text = element_text(size=12),
          panel.background = element_rect(fill="white",
                                          colour="black"))
  ### Plotting:
  if (graphs=="all")
    suppressMessages(return(multiplot(e1,e2, cols=2)))
  if (graphs==1)
    suppressMessages(return(e1))
  if (graphs==2)
    suppressMessages(return(e2))
} 

#' Graphical diagnostics for class 'intra.physig'
#'
#' \code{sensi_plot_intra.physig} Plot results from \code{intra_physig}.
#' @param x output from \code{intra_physig}
#' @param graphs choose which graph should be printed in the output ("all", 1 or 2)
#' @param ... further arguments to methods
#' @importFrom ggplot2 scale_color_manual geom_histogram geom_abline geom_density 
#' geom_vline xlab geom_point theme
#' @author Caterina Penone and Gustavo Paterno
#' @seealso \code{\link[ggplot2]{ggplot}}, \code{\link[sensiPhy]{intra_phylm}}
#' \code{\link[sensiPhy]{intra_phylm}}
#' @details For 'x' from \code{intra_physig}
#' 
#' Graphs 1: Distribution of estimated phylosgentic signal (lambda or K) for each simulation
#' Red vertical line represents the average signal among all estimates. 
#'
#' Graph 2: Distribution of p-values for the phylogenetic signal (K or lambda) 
#' for each simulation. Red vertical line represents the alpha significance level = 0.05.
#' @importFrom grid unit 
#' @importFrom stats plogis
#' @export
sensi_plot.intra.physig <- function(x, graphs="all", ...){
  
  ### Basic checking:
  method <- x$call$method
  if (is.null(x$call$method)) method <- "K"
  ### Distribution of K values estimated
  e1 <- ggplot2::ggplot(x$intra.physig.estimates, aes(x = estimate)) +
    geom_histogram(fill = "yellow",colour = "black", size = .2,
                   alpha = .3) +
    geom_vline(xintercept = x$stats$mean[1], color = "red", linetype = 2, size = .7) +
    xlab(paste("Estimated", method, "values", sep = " ")) +
    ylab("Frequency") +
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          panel.background = element_rect(fill = "white",
                                          colour = "black")); 
  
  ### Distribution of Values estimated
  e2 <- ggplot2::ggplot(x$intra.physig.estimates, aes(x = pval)) +
    geom_histogram(fill = "yellow", colour = "black", size = .2,
                   alpha = .3) +
    geom_vline(xintercept = 0.05, color = "red", linetype = 1, size = .7) +
    xlab("Estimated P values") +
    ylab("Frequency") +
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          panel.background = element_rect(fill = "white",
                                          colour = "black"));
  ### Plotting:
  if (graphs == "all")
    suppressMessages(return(multiplot(e1,e2, cols = 2)))
  if (graphs == 1)
    suppressMessages(return(e1))
  if (graphs == 2)
    suppressMessages(return(e2))
} 
