#' Graphical diagnostics for class 'sensiTree_Intra'
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
#' @seealso \code{\link[ggplot2]{ggplot}}
#' @details For 'x' from intra_influ_phylm or intra_influ_phyglm:
#' 
#' Graph 1: Distribution of estimated slopes or intercepts for each 
#' simulation (leave-one-out deletion). Red vertical line represents the original
#' slope or intercept from the full model (with all species). 
#' 
#' Graph 2: Original regression plot (\eqn{trait~predictor}). Standardized 
#' difference in slope or intercept is represented by a continous size scale. 
#' The names of the most influential species (sDF > cutoff) are ploted in the
#' graph. 
#' 
#' Graph 3: Distribution of standardized difference in slope or intercept. Red 
#' colour indicates inbfluential species (with a standardised difference above 
#' the value of \code{cutoff}).
#' 
#' Graph 4: Distribution of the percentage of change in slope or intercept.
#' @importFrom grid unit 
#' @export

### Start:
sensi_plot.sensiTree_Influ <- function(x, graphs="all", ...){
  ### Graph one
  n.tree <- length(unique(x$sensi.estimates$iteration))
  sp.estimate <- unlist(as.list(x$influential.species$influ.sp.estimate$influ.sp.estimate))
  sp.estimate.tab <- table(sp.estimate)
  sp.estimate <- sp.estimate.tab[order(sp.estimate.tab,decreasing=T)] 
  sp.estimate.tab <- data.frame("Species removed" = names(sp.estimate), 
                                "Significant" = (as.numeric(sp.estimate)/n.tree)*100)
  g1 <- 
    ggplot(sp.estimate.tab, aes(y = (Significant), 
                                x = reorder(Species.removed, - Significant))) +
    geom_bar(stat = "identity", fill = "lightblue") + 
    xlab("Species removed") + ylab("(%) iterations") +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=12),
          panel.background = element_rect(fill="white",
                                          colour="black"),
          legend.position = "none",
          axis.text.x=element_text(angle=65,hjust=1)) +
    annotate(geom = "text",x=Inf, y = Inf, vjust= 1.5, hjust= 1.2,
             label = paste("Number of iterations =", n.tree))
  
  ### Graph two
  es <- x$sensi.estimates
  es <- es[es$species %in% summary(x)[[1]][[1]], ]
  es$species <- as.factor(es$species)
  
  ### ord species in the same order as graph 1
  es$species <- factor(es$species, levels = sp.estimate.tab$Species.removed)
  g2 <-
    ggplot(es, aes(y = DIFestimate, x = species)) + 
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


#' Graphical diagnostics for class 'sensiIntraInflu'
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
#' @details For 'x' from sensiTree_Intra_phylm or sensiTree_Intra_phyglm:
#' 
#' Graph 1: Distribution of estimated estimates (slopes) or intercepts for each 
#' simulation (leave-one-out deletion). Red vertical line represents the original
#' slope or intercept from the full model (with all species). 
#' 
#' Graph 2: Original regression plot (\eqn{trait~predictor}). Standardized 
#' difference in slope or intercept is represented by a continous size scale. 
#' The names of the most influential species (sDF > cutoff) are ploted in the
#' graph. 
#' 
#' Graph 3: Distribution of standardized difference in slope or intercept. Red 
#' colour indicates inbfluential species (with a standardised difference above 
#' the value of \code{cutoff}).
#' 
#' Graph 4: Distribution of the percentage of change in slope or intercept.
#' @importFrom grid unit 
#' @export



#' @export
sensi_plot.sensiIntra_Influ <- function(x, graphs="all", ...){
  sensi_plot.sensiTree_Influ(x, graphs, ...)
}
