#' Graphical diagnostics for class 'influ.physig'
#'
#' \code{sensi_plot.influ_physig} Plot results from \code{influ_physig}
#' @param x output from \code{influ_physig}
#' @param graphs choose which graph should be printed on the output ("all", 1, 2)
#' @param ... further arguments to methods.
#' @importFrom ggplot2 aes geom_histogram geom_density geom_vline 
#' xlab theme element_text geom_point scale_colour_gradient element_rect ylab xlab
#' ggtitle element_blank
#' @author Gustavo Paterno
#' @seealso \code{\link[ggplot2]{ggplot}}, \code{\link[phytools]{phylosig}}
#' @details For 'x' from influ_physig:
#' 
#' \strong{Graph 1:} Distribution of estimated phylogenetic signal (K or lambda) for each 
#' simulation (leave-one-out deletion). Dashed red vertical line represents the original
#'  estimate of phylogenetic signal with the full data (with all species). 
#' 
#' \strong{Graph 2:} Distribution of P values for the phylogenetic signal (K or lambda) 
#' for each simulation (leave-one-out deletion). Red vertical line represents 
#' the alpha significance level = 0.05.
#' 
#' @export

### Start:
sensi_plot.influ.physig <- function(x, graphs = "all", ...){
  ### Nulling variables:
  estimate <- pval <- NULL
  ### Organizing values:
  result <- x$influ.physig.estimates
  method <- x$call$method
  if(is.null(x$call$method)) method <- "K"
  est.0  <- as.numeric(x$full.data.estimates$estimate)
  cutoff <-  x$cutoff
  
  # Plot 1: Distribution of estimated phylogenetic signal metrics
  i1 <- ggplot2::ggplot(result, aes(x = estimate))+
    geom_histogram(fill="yellow",colour="black", size=.2,
                   alpha = .3) +
    geom_vline(aes(xintercept = est.0, color = "red"), linetype = 2, size =.7)+
    scale_color_manual(values = "red", name = "", labels = c("Full data")) +
    xlab(paste("Estimated", method, "values", sep = " "))+
    ylab("Frequency")+
    theme(axis.title=element_text(size=12),
          axis.text = element_text(size=12),
          panel.background = element_rect(fill="white",
                                          colour="black"),
          legend.key = element_rect(fill = "transparent"),
          legend.background = element_rect(fill = "transparent"),
          legend.position = c(.9,.9)) 
  ### Distribution of P-values
  i2 <- ggplot2::ggplot(result,aes(x = pval))+
    geom_histogram(fill="yellow",colour="black", size=.2,
                   alpha = .3) +
    geom_vline(xintercept = 0.05,color="red",linetype=1,size=.7)+
    xlab(paste("P values for", method, sep = " "))+
    ylab("Frequency")+
    theme(axis.title=element_text(size=12),
          axis.text = element_text(size=12),
          panel.background = element_rect(fill="white",
                                          colour="black"))
  ### Expor two graphs:
  if (graphs == 1) 
    suppressMessages(return(i1))
  if (graphs == 2) 
    suppressMessages(return(i2))
  if (graphs == "all")
    suppressMessages(return(multiplot(i1, i2, cols = 2)))
}

