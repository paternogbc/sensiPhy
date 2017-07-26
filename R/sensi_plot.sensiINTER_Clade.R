#' Graphical diagnostics for class 'sensiIntra_Clade'
#'
#' Plot results from \code{intra_clade_phylm} and \code{intra_clade_phyglm}
#' @param x output from \code{intra_clade_phylm} or \code{intra_clade_phyglm}
#' @param clade The name of the clade to be evaluated (see details)
#' @param ... further arguments to methods.
#' @importFrom ggplot2 aes theme element_text geom_point element_rect ylab xlab
#' ggtitle element_blank geom_abline scale_shape_manual scale_linetype_manual 
#' guide_legend element_rect
#' guides
#' 
#' @author Gustavo Paterno, Caterina Penone
#' @seealso \code{\link[sensiPhy]{intra_clade_phylm}} 
#' @details For 'x' from intra_clade_phylm or intra_clade_phyglm:
#' 
#' Graph 1: The original scatterplot \eqn{y = a + bx} (with the 
#' full dataset) and a comparison between the regression lines of the full model
#' and the model without the selected clade (set by \code{clade}). For further
#' details about this method, see \code{\link[sensiPhy]{clade_phylm}}.
#' 
#' Species from the selected clade are represented in red (removed species), black
#' solid line represents the regression with the full model and red dashed line represents
#' the regression of the model without the species from the selected clade.
#' To check the available clades to plot, see \code{x$clade.model.estimates$clade} 
#' in the object returned from \code{intra_clade_phylm} or \code{intra_clade_phyglm}. 
#' 
#' Graph 2: Distribution of the simulated slopes (Null distribution
#' for a given clade sample size).
#' The red dashed line represents the estimated slope for the reduced model 
#' (without the focal clade) and the black line represents the slope for the 
#' full model.
#'  
#' @importFrom ggplot2 aes_string
#' @importFrom stats model.frame qt plogis 
#' @export

sensi_plot.sensiIntra_Clade <- function(x, clade = NULL, ...){
  yy <- NULL; slope <- NULL
  
  #aggregate results (needs to be improved in future!)
  x$full.model.estimates$type <- c("inter","slope")
  x$full.model.estimates <- stats::aggregate(.~type,x$full.model.estimates,FUN=mean)
  x$clade.model.estimates <- stats::aggregate(.~clade,x$clade.model.estimates,FUN=mean)
  x$null.dist$simul <- 1:x$call$n.sim
  x$null.dist <- stats::aggregate(.~clade+simul,x$null.dist,FUN=mean)
  
  # start:
  full.data <- x$data
  mappx <- x$formula[[3]]
  mappy <- x$formula[[2]]
  vars <- all.vars(x$formula)
  clade.col <- x$clade.col
  
  clades.names <- unique(x$clade.model.estimates$clade)
  if (is.null(clade) == T){
    clade <- clades.names[1]
    warning("Clade argument was not defined. Plotting results for clade: ",
            clade,"
            Use clade = 'clade name' to plot results for other clades")
  }
  clade.n <- which(clade == clades.names)
  if (length(clade.n) == 0) stop(paste(clade,"is not a valid clade name"))
  
  ### Organizing values:
  result <- x$clade.model.estimates
  intercept.0 <-  as.numeric(x$full.model.estimates$Estimate[1])
  slope.0     <-  as.numeric(x$full.model.estimates$Estimate[2])
  
  inter <- c(x$clade.model.estimates$intercept[clade.n ],
             intercept.0)
  slo <-  c(x$clade.model.estimates$slope[clade.n ],
            slope.0)
  model <- NULL
  estimates <- data.frame(inter,slo, model=c("Without clade", "Full model"))
  
  xf <- model.frame(formula = x$formula, data = full.data)[,2]
  yf <- plogis(estimates[2,1] + estimates[2,2] * xf)
  yw <- plogis(estimates[1,1] + estimates[1,2] * xf)
  plot_data <- data.frame("xf" = c(xf,xf),
                          "yy" = c(yw, yf),
                          model = rep(c("Without clade","Full model"),
                                      each = length(yf)))
  
  match.y <- which(full.data[, clade.col] == clade)
  match.n <- which(full.data[, clade.col] != clade)
  
  g1 <- ggplot2::ggplot(full.data, aes_string(y = mappy, x = mappx),
                        environment = parent.frame())+
    geom_point(data = full.data[match.n, ], alpha = .7,
               size = 4)+
    geom_point(data = full.data[match.y, ],alpha = .5,
               size = 4, aes(shape = "Removed species"), colour = "red")+
    
    scale_shape_manual(name = "", values = c("Removed species" = 16))+
    guides(shape = guide_legend(override.aes = list(linetype = 0)))+
    scale_linetype_manual(name = "Model", values = c("Full model" = "solid",
                                                     "Without clade" = "dashed"))+
    scale_color_manual(name = "Model", values = c("Full model" = "black",
                                                  "Without clade" = "red"))+
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = 12),
          panel.background = element_rect(fill = "white", colour = "black"),
          legend.position="bottom", legend.box = "horizontal")+
    ggtitle(paste("Clade removed: ", clade, sep = ""))
  
  ### Permuation Test plot:
  nd <- x$null.dist
  ces <- x$clade.model.estimates
  
  nes <- nd[nd$clade == clade, ]
  slob <- ces[ces$clade == clade ,]$slope
  slfu <- x$full.model.estimates$Estimate[[2]]
  
  ### P.value permutation test:
  p.values <- summary(x)[[1]]
  P <- p.values[p.values$clade.removed == clade, ]$Pval.randomization
  
  g2 <- ggplot2::ggplot(nes ,aes(x=slope))+
    geom_histogram(fill="yellow",colour="black", size=.2,
                   alpha = .3) +
    geom_vline(xintercept = slob, color="red",linetype=2,size=.7)+
    geom_vline(xintercept = slfu, color="black",linetype=1,size=.7)+
    xlab(paste("Simulated slopes | N.species = ", 
               ces[ces$clade==clade, ]$N.species, "| N.sim = ", 
               nrow(nes))) +
    ylab("Frequency")+
    theme(axis.title=element_text(size=12),
          axis.text = element_text(size=12),
          panel.background = element_rect(fill="white",
                                          colour="black"))+
    ggtitle(paste("Randomization test for", clade, " | P = ", 
                  sprintf("%.3f", P)))
  
  ### plot lines: linear or logistic depending on output class
  if(length(class(x)) == 1){
    g.out <- g1 + geom_abline(data = estimates, aes(intercept = inter, slope = slo,
                                                    linetype = factor(model),color=factor(model)),
                              size=.8)
  }
  if(length(class(x)) == 2){
    g.out <- g1 + geom_line(data = plot_data, aes(x = xf, y = yy, linetype = factor(model),color=factor(model)))
  }
  return(suppressMessages(multiplot(g.out, g2, cols=2)))
  }



#' Graphical diagnostics for class 'sensiTree_Clade'
#'
#' Plot results from \code{tree_clade_phylm} and \code{tree_clade_phyglm}
#' @param x output from \code{tree_clade_phylm} or \code{tree_clade_phyglm}
#' @param clade The name of the clade to be evaluated (see details)
#' @param ... further arguments to methods.
#' @importFrom ggplot2 aes theme element_text geom_point element_rect ylab xlab
#' ggtitle element_blank geom_abline scale_shape_manual scale_linetype_manual 
#' guide_legend element_rect
#' guides
#' 
#' @author Gustavo Paterno, Caterina Penone
#' @seealso \code{\link[sensiPhy]{tree_clade_phylm}} 
#' @details For 'x' from tree_clade_phylm or tree_clade_phyglm:
#' 
#' Graph 1: The original scatterplot \eqn{y = a + bx} (with the 
#' full dataset) and a comparison between the regression lines of the full model
#' and the model without the selected clade (set by \code{clade}). For further
#' details about this method, see \code{\link[sensiPhy]{clade_phylm}}.
#' 
#' Species from the selected clade are represented in red (removed species), black
#' solid line represents the regression with the full model and red dashed line represents
#' the regression of the model without the species from the selected clade.
#' To check the available clades to plot, see \code{x$clade.model.estimates$clade} 
#' in the object returned from \code{tree_clade_phylm} or \code{tree_clade_phyglm}. 
#' 
#' Graph 2: Distribution of the simulated slopes (Null distribution
#' for a given clade sample size).
#' The red dashed line represents the estimated slope for the reduced model 
#' (without the focal clade) and the black line represents the slope for the 
#' full model.
#'  
#' @importFrom ggplot2 aes_string
#' @importFrom stats model.frame qt plogis 
#' @export
sensi_plot.sensiTree_Clade <- function(x, graphs = "all", ...){
  sensi_plot.sensiIntra_Clade(x, clade = NULL, ...)
}