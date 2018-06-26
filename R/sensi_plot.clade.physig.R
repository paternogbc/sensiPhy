#' Graphical diagnostics for class 'clade.physig'
#'
#' \code{plot_clade_physig} Plot results from \code{clade_physig}
#' @param x output from \code{influ_physig}
#' @param clade The name of the clade to be evaluated (see details)
#' @param ... further arguments to methods.
#' @importFrom ggplot2 aes geom_histogram geom_density geom_vline 
#' xlab theme element_text geom_point scale_colour_gradient element_rect ylab xlab
#' ggtitle element_blank
#' @author Gustavo Paterno
#' @seealso \code{\link[ggplot2]{ggplot}}, \code{\link[phytools]{phylosig}}
#' @details For 'x' from clade_physig:
#' 
#' \strong{Graph 1:} Distribution of the simulated signal estimates (Null distribution
#' for a given clade sample size).
#' The red dashed line represents the estimated signal for the reduced data 
#' (without the focal clade) and the black line represents the signal estimate
#'  for the full data.
#' 
#' @export

### Start:
sensi_plot.clade.physig <- function(x, clade = NULL, ...) {
  ### nulling variables:
  estimate <- model <- NULL
  
  # check clade
  clades.names <- x$sensi.estimates$clade
  if (is.null(clade) == TRUE){
    clade <- clades.names[1]
    message("Clade argument was not defined. Plotting results for clade: ",
            clade,"
                Use clade = 'clade name' to plot results for other clades")
  }
  
  method <- x$method
  times <- x$call$n.sim
  if(is.null(x$call$n.sim)) times <- 1000
    
  ces       <- x$sensi.estimates 
  N.species <- ces[ces$clade==clade, ]$N.species
  p.values <- summary(x)[[2]]
  P <- p.values[p.values$clade.removed == clade, ]$Pval.randomization   
    
  wcf <- x$sensi.estimates$clade %in% clade
  wcn <- x$null.dist$clade %in% clade
  if(sum(wcf) == 0) stop(paste(clade, "is not a valid clade name", sep = " "))
    
  # Full estimate
  e.0 <- x$full.data.estimates$estimate 
  # Withou clade estimate
  c.e <- x$sensi.estimates[wcf, ]$estimate
  nd <- x$null.dist[wcn, ] ### CLADE NULL DIST
    
  ## Estimates dataframe
  vl <- data.frame(model = c("Full data", "Without clade"), 
                     estimate = c(e.0,c.e))
  leg.title <- paste("Estimated", method)
    
  ### Graph 1
  p1 <- ggplot2::ggplot(nd, aes(x = nd$estimate)) + 
      geom_histogram(fill="yellow",colour="black", size=.2, alpha = .3) +
      geom_vline(data = vl, aes(xintercept = estimate,
                                colour = model, linetype = model),
                 size = 1.2) + 
      scale_color_manual(leg.title, values = c("black","red")) +
      scale_linetype_manual(leg.title, values = c(1,2)) +
      theme(axis.title=element_text(size=12),
            axis.text = element_text(size=12),
            panel.background = element_rect(fill="white",
                                            colour="black"),
            legend.key = element_rect(fill = "transparent"),
            legend.background = element_rect(fill = "transparent"),
            legend.position = c(.9,.9))+
      ylab("Frequency")+
      xlab(paste("Simulated", method, "| N.species = ", 
                 N.species, "| N.sim = ", times)) +
      ggtitle(paste("Randomization test for", clade, " | P = ", 
                    sprintf("%.3f", P)))
  suppressMessages(print(p1)) 
}
  

 
