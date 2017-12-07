#' Graphical diagnostics for class 'tree.physig'
#'
#' \code{sensi_plot_tree.bd} Plot results from \code{tree_bd}.
#' @param x output from \code{tree_bd}
#' @param graphs choose which graph should be printed in the output ("all", 1 or 2)
#' @param ... further arguments to methods
#' @importFrom ggplot2 scale_color_manual geom_histogram geom_abline geom_density 
#' geom_vline xlab geom_point theme
#' @author Caterina Penone and Gustavo Paterno
#' @seealso \code{\link[ggplot2]{ggplot}}, \code{\link[sensiPhy]{tree_phylm}}
#' \code{\link[sensiPhy]{intra_phylm}}
#' @details For 'x' from \code{tree_bd}
#' 
#' \strong{Graphs 1:} Distribution of estimated diversification or speciation rate
#'  for each tree. Red vertical line represents the average signal among all estimates. 
#'  
#' \strong{Graphs 2:} Estimates across each phylogenetic tree.
#'
#' @importFrom grid unit 
#' @importFrom stats plogis
#' @importFrom stats reorder
#' @export
sensi_plot.tree.bd <- function(x, graphs="all", ...){
  
  ### Nulling variables
  estimate <- n.tree <- model <- NULL
  ### Basic checking:
  method <- x$call$method
  if(is.null(x$call$method)) method <- "ms"
  if (method == "ms") label = "estimated diversification rate [ms]"
  if (method == "km") label = "estimated speciation rate [km]"
  ### Distribution of ms values estimated
  e1 <- ggplot2::ggplot(x$tree.bd.estimates, aes(x=estimate))+
    geom_histogram(fill="yellow",colour="black", size=.2,
                   alpha = .3) +
    geom_vline(xintercept = x$stats$mean, color="red",linetype=2,size=.7)+
    xlab(label)+
    ylab("Frequency")+
    theme(axis.title=element_text(size=12),
          axis.text = element_text(size=12),
          panel.background = element_rect(fill="white",
                                          colour="black"))
  
  ### Distribution of Values estimated
  e2 <- ggplot2::ggplot(x$tree.bd.estimates, 
                        aes(y = estimate, x = reorder(n.tree, estimate)))+
    geom_point(size = 3, color = "tomato") +
    xlab("tree")+
    ylab("estimate")+
    theme(axis.title=element_text(size=12),
          axis.text.x = element_text(size=8),
          axis.text.y = element_text(size=12),
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