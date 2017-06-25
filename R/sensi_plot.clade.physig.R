#' Graphical diagnostics for class 'clade.physig'
#'
#' \code{plot_influ_physig} Plot results from \code{influ_physig}
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
#' Graph 2: Distribution of the simulated slopes (Null distribution
#' for a given clade sample size).
#' The red dashed line represents the estimated slope for the reduced model 
#' (without the focal clade) and the black line represents the slope for the 
#' full model.
#' 
#' @export

### Start:
sensi_plot.clade.physig <- function(x, clade = NULL, ...) {
  clade = NULL
  # check clade
  clades.names <- x$sensi.estimates$sensi.clade$clade
  if (is.null(clade) == TRUE){
    clade <- clades.names[1]
    message("Clade argument was not defined. Plotting results for clade: ",
            clade,"
                Use clade = 'clade name' to plot results for other clades")
  }
  
  method <- x$method
  times <- x$call$times
  if(is.null(x$call$times)) times <- 1000
    
  ces       <- x$sensi.estimates$sensi.clade 
  N.species <- ces[ces$clade==clade, ]$N.species
  p.values <- summary(x)[[2]]
  P <- p.values[p.values$clade.removed == clade, ]$Pval.randomization   
    
  wcf <- x$sensi.estimates$sensi.clade$clade %in% clade
  wcn <- x$sensi.estimates$null.dist$clade %in% clade
  if(sum(wcf) == 0) stop(paste(clade, "is not a valid clade name", sep = " "))
    
  # Full estimate
  e.0 <- x$full.data.estimates$estimate 
  # Withou clade estimate
  c.e <- x$sensi.estimates$sensi.clade[wcf, ]$estimate
  nd <- x$sensi.estimates$null.dist[wcn, ] ### CLADE NULL DIST
    
  ## Estimates dataframe
  vl <- data.frame(model = c("Full data", "Without clade"), 
                     estimate = c(e.0,c.e))
  leg.title <- paste("Estimated", method)
    
  ### Graph 1
  p1 <- ggplot(nd, aes(x = nd$estimate)) + 
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
  

 
