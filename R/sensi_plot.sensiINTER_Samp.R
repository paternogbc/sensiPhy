#' Graphical diagnostics for class 'sensiIntra_Samp' and 'sensiTree_samp'
#'
#' \code{sensi_plot.sensiIntra_Samp} Plot results from \code{sensiIntra_Samp_phylm} and
#' \code{sensiIntra_Samp_phyglm}
#' @param x output from \code{sensiIntra_Samp_phylm} or \code{sensiIntra_Samp_phyglm}
#' @param graphs choose which graph should be printed on the output ("all", 1,2,3 or 4)
#' @param param choose which model parameter should be ploted  ("intercept" or "estimate")
#' @param ... further arguments to methods
#' @importFrom ggplot2 scale_x_continuous scale_colour_manual geom_hline 
#' geom_bar scale_fill_manual scale_y_continuous geom_boxplot geom_line 
#' @author Gustavo Paterno
#' @seealso \code{\link[ggplot2]{ggplot}}, \code{\link[sensiPhy]{samp_phylm}}
#' \code{\link[sensiPhy]{samp_phyglm}}
#' @details For 'x' from sensiIntra_Samp_phylm or sensiIntra_Samp_phyglm:
#' 
#' Graph 1: Estimated estimates (slopes) or intercepts for each simution across  
#' percentages of species removed and intra simulations. Colours represent percentage 
#' of change in comparison with the full model (blue = lower than 5, orange = 
#' between 5 and 10 and red = higher than 10).
#' The red horizontal line represents the original slope or 
#' intercept from the full model (with all species). 
#' 
#' Graph 2: The proportion of estimated slopes and intercepts in each category 
#' across the percentage of species removed.
#' 
#' Graph 3: Estimated phylogenetic model parameter for each simulation across
#' the percentage of species removed.
#' 
#' Graph 4: The percentage of significant slopes or intercepts across the 
#' percentage of species removed.  
#' 
#' @note If model = "BM", only plots 1, 2 and 4 are printed. Plot 3, phylogenetic
#'  model parameter is not available for model = "BM"
#' @export

sensi_plot.sensiIntra_Samp <- function(x, graphs = "all", ...)
{
  
  ### prepare data
  fes <- x$full.model.estimates
  es <- x$sensi.estimates
  n.tree <- nrow(fes)/2
  est.or <- seq(2, nrow(fes), by = 2)
  e.0.mean <- mean(fes[est.or, ]$Estimate)
  
  ### simulated estimates
  es$iteration <- as.factor(es$iteration)
  es$n.percent <- as.factor(es$n.percent)
  
  ### sig table
  sig <- x$sign.analysis
  
  ### plots ---------------------------------------------------------------------------------
  # 1 Estimates across trees and % of removal
  g1 <- ggplot(es, aes(y = estimate, x = n.percent, group = reorder(iteration, estimate))) + 
    geom_jitter(size = .8, position= position_dodge(width = .8), color = "red") +
    ggplot2::stat_summary(fun.y = "mean", size = 3, geom = "point", position= position_dodge(width = .8), color = "red") +
    geom_hline(yintercept = e.0.mean, size = 1) +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=12),
          panel.background = element_rect(fill="white",
                                          colour="black"),
          legend.position = "none") +
    xlab("% Species removed") + ylab("estimate")
  
  
  # 2 Sigificance table across trees and % of removal
  g2 <- ggplot(sig, aes(y = perc.sign.estimate*100, x = percent_sp_removed)) +
    geom_point(size = 1, position= position_dodge(width = .8), color = "red") +
    ggplot2::stat_summary(fun.y = "mean", size = 5, geom = "point",
                          position= position_dodge(width = .8), color = "red", alpha = .5) +
    ggplot2::stat_summary(fun.y = "mean", size = 1, geom = "line", 
                          position= position_dodge(width = .8), color = "red") +
    scale_y_continuous(limits=c(0,100),breaks=seq(0,100,10))+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=12),
          panel.background = element_rect(fill="white",
                                          colour="black"),
          legend.position = "none") +
    xlab("% Species removed") + ylab("% Significant estimate")
    
  ### Output-------------------------------------------------------------------------------
  ### Plotting:
  if (graphs=="all")
    return(suppressMessages(multiplot(g1,g2, cols = 2)))
  if (graphs==1)
    return(suppressMessages(g1))
  if (graphs==2)
    return(suppressMessages(g2))
}



#' Graphical diagnostics for class 'sensiTree_samp'
#'
#' \code{sensi_plot.sensiTree_Samp} Plot results from \code{sensiTree_Samp_phylm} and
#' \code{sensiTree_Samp_phyglm}
#' @param x output from \code{sensiTree_Samp_phylm} or \code{sensiTree_Samp_phyglm}
#' @param graphs choose which graph should be printed on the output ("all", 1,2,3 or 4)
#' @param param choose which model parameter should be ploted  ("intercept" or "estimate")
#' @param ... further arguments to methods
#' @importFrom ggplot2 scale_x_continuous scale_colour_manual geom_hline 
#' geom_bar scale_fill_manual scale_y_continuous geom_boxplot geom_line 
#' @author Gustavo Paterno
#' @seealso \code{\link[ggplot2]{ggplot}}, \code{\link[sensiPhy]{samp_phylm}}
#' \code{\link[sensiPhy]{samp_phyglm}}
#' @details For 'x' from sensiTree_Samp_phylm or sensiTree_Samp_phyglm:
#' 
#' Graph 1: Estimated slopes (slopes) or intercepts for each simution across  
#' percentages of species removed and intra simulations. Colours represent percentage 
#' of change in comparison with the full model (blue = lower than 5, orange = 
#' between 5 and 10 and red = higher than 10).
#' The red horizontal line represents the original slope or 
#' intercept from the full model (with all species). 
#' 
#' Graph 2: The proportion of estimated slopes and intercepts in each category 
#' across the percentage of species removed.
#' 
#' Graph 3: Estimated phylogenetic model parameter for each simulation across
#' the percentage of species removed.
#' 
#' Graph 4: The percentage of significant slopes or intercepts across the 
#' percentage of species removed.  
#' 
#' @note If model = "BM", only plots 1, 2 and 4 are printed. Plot 3, phylogenetic
#'  model parameter is not available for model = "BM"
#' @export

sensi_plot.sensiTree_Samp <- function(x, graphs = "all", param = "estimate", ...){
  sensi_plot.sensiIntra_Samp(x, graphs, ...)
}