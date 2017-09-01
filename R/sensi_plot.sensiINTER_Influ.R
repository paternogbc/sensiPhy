#' Graphical diagnostics for class 'sensiIntra_Influ'
#'
#' \code{sensi_plot.intra_influ} Plot results from \code{intra_influ_phylm} and 
#' \code{intra_influ_phyglm}
#' @param x output from \code{influ_phylm} or \code{intra_influ_phyglm}
#' @param graphs choose which graph should be printed on the output ("all", 1,2,3 or 4)
#' @param param choose which parameter ("intercept" or "estimate" should be printed)
#' @param ... further arguments to methods
#' @importFrom ggplot2 aes geom_histogram geom_density geom_vline 
#' xlab theme element_text geom_point scale_colour_gradient element_rect ylab xlab
#' ggtitle element_blank
#' @author Gustavo Paterno, Caterina Penone
#' @seealso \code{\link[sensiPhy]{intra_influ_phylm}}; 
#' \code{\link[sensiPhy]{intra_influ_phyglm}}
#' @details For 'x' from intra_influ_phylm or intra_influ_phyglm:
#' 
#' \strong{Graph 1:} Most common influential species on model estimates.
#'  Percentage of iterations (n.tree or n.intra) where the removal of individual
#'   species caused significant change in model estimate (sDIFestimate > cutoff). 
#' 
#' \strong{Graph 2:} Shift on model estimate after species removal for most influential
#'  species. Small red dots represent individual reruns between different trees or 
#'  simulations while large red dots represent the average DIFestimate among all 
#'  iterations.
#' @importFrom grid unit 
#' @export

### Start:
sensi_plot.sensiIntra_Influ <- function(x, graphs="all", ...){
  ### Nulling variables:
  Significant <- Species.removed <- iteration <- DIFestimate <- species <- NULL
  species <- 
  ### Graph one
  n.tree <- length(unique(x$sensi.estimates$iteration))
  sp.estimate <- unlist(as.list(x$influential.species$influ.sp.estimate$influ.sp.estimate))
  sp.estimate.tab <- table(sp.estimate)
  sp.estimate <- sp.estimate.tab[order(sp.estimate.tab,decreasing=T)] 
  sp.estimate.tab <- data.frame("Species removed" = names(sp.estimate), 
                                "Significant" = (as.numeric(sp.estimate)/n.tree)*100)
  g1 <- 
    ggplot2::ggplot(sp.estimate.tab, aes(y = (Significant), 
                                x = reorder(Species.removed, - Significant))) +
    geom_bar(stat = "identity", fill = "lightblue") + 
    xlab("Species removed") + ylab("(%) iterations") +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=12),
          panel.background = element_rect(fill="white",
                                          colour="black"),
          legend.position = "none",
          axis.text.x=element_text(angle=65,hjust=1)) +
    scale_y_continuous(limits = c(0,100), breaks = seq(0,100,10)) +
    ggplot2::annotate(geom = "text",x=Inf, y = Inf, vjust= 1.5, hjust= 1.2,
             label = paste("Number of iterations =", n.tree))
  
  ### Graph two
  es <- x$sensi.estimates
  es <- es[es$species %in% summary(x)[[1]][[1]], ]
  es$species <- as.factor(es$species)
  
  ### ord species in the same order as graph 1
  es$species <- factor(es$species, levels = sp.estimate.tab$Species.removed)
  g2 <-
    ggplot2::ggplot(es, aes(y = DIFestimate, x = species)) + 
    geom_point(color = "red") +
    #scale_x_discrete(labels = c(strtrim(es$species, 4))) +
    ggplot2::stat_summary(fun.y = "mean", size = 5, color = "red", geom = "point") +
    geom_hline(yintercept = 0, color = "black") +
    theme(axis.text=element_text(size= 12),
          axis.title=element_text(size=12),
          panel.background = element_rect(fill="white",
                                          colour="black"),
          legend.position = "none",
          axis.text.x=element_text(angle=65,hjust=1)) +
    xlab("Species removed") + ylab("DIFestimate")
  
  ### Output-------------------------------------------------------------------------------
  ### Plotting:
  if (graphs=="all")
    return(suppressMessages(multiplot(g1,g2, cols = 2)))
  if (graphs==1)
    return(suppressMessages(g1))
  if (graphs==2)
    return(suppressMessages(g2))
}


#' Graphical diagnostics for class 'sensiTree_Influ'
#'
#' \code{sensi_plot.sensiTree_Influ} Plot results from \code{tree_influ_phylm} and 
#' \code{tree_influ_phyglm}
#' @param x output from \code{tree_influ_phylm} or \code{tree_influ_phyglm}
#' @param graphs choose which graph should be printed on the output ("all", 1,2,3 or 4)
#' @param param choose which parameter ("intercept" or "estimate" should be printed)
#' @param ... further arguments to methods
#' @importFrom ggplot2 aes geom_histogram geom_density geom_vline 
#' xlab theme element_text geom_point scale_colour_gradient element_rect ylab xlab
#' ggtitle element_blank
#' @author Gustavo Paterno, Caterina Penone
#' @seealso \code{\link[ggplot2]{ggplot}}
#' @details For 'x' from sensiTree_Influ_phylm or sensiTree_Influ_phyglm:
#' 
#' \strong{Graph 1:} Most common influential species on model estimates.
#'  Percentage of iterations (n.tree or n.intra) where the removal of individual
#'   species caused significant change in model estimate (sDIFestimate > cutoff). 
#' 
#' \strong{Graph 2:} Shift on model estimate after species removal for most influential
#'  species. Small red dots represent individual reruns between different trees or 
#'  simulations while large red dots represent the average DIFestimate among all 
#'  iterations.
#' @importFrom grid unit 
#' @export
#' @export
sensi_plot.sensiTree_Influ <- function(x, graphs="all", ...){
  sensi_plot.sensiIntra_Influ(x, graphs, ...)
}
