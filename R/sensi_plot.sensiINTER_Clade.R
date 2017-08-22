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
#' To check the available clades to plot, see \code{x$sensi.estimates$clade} 
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

sensi_plot.sensiTree_Clade <- function(x, clade = NULL, graphs = "all", ...){
    ### Get clade name:
  clades.names <- unique(x$sensi.estimates$clade)
  if (is.null(clade) == T){
    clade <- clades.names[1]
    warning("Clade argument was not defined. Plotting results for clade: ",
            clade,"
            Use clade = 'clade name' to plot results for other clades")
  }
  clade.n <- which(clade == clades.names)
  if (length(clade.n) == 0) stop(paste(clade,"is not a valid clade name"))

### organize data:
fes   <- x$full.model.estimates
f     <- data.frame(fes)[seq(2,nrow(fes), 2), ]
f.ord <- match(sort(f$iteration), f$iteration)
ces   <- x$sensi.estimates
ces.c <- ces[ces$clade == clade, ]
nd    <-  x$null.dist
nd.c  <- nd[nd$clade == clade, ]
s.est <- summary(x)[[1]]

### null > | < then clade removed
if (s.est[s.est$clade.removed == clade, ]$DIFestimate > 0)            
  out <- which(nd.c$estimate >= rep(ces.c$estimate, each = x$call$n.sim))
if (s.est[s.est$clade.removed == clade, ]$DIFestimate < 0)            
  out <- which(nd.c$estimate <= rep(ces.c$estimate, each = x$call$n.sim))
nd.c.out <- nd.c[out, ]

### Transform iteration intofactor:
fes$iteration   <- as.factor(fes$iteration)
ces.c$iteration <- as.factor(ces.c$iteration)
nd.c$iteration  <- as.factor(nd.c$iteration)


### Plot 1: Estimated slopes for each clade between trees ----------------------
g1 <-
  ggplot(ces, aes(y = estimate, x = reorder(as.factor(strtrim(clade, 4)), estimate))) + 
  geom_point(color = "red") +
  ggplot2::stat_summary(fun.data = "mean_cl_boot", size = 1, color = "red") +
  geom_hline(yintercept = mean(f$Estimate), size = 1.5) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        panel.background = element_rect(fill="white",
                                        colour="black"),
        legend.position = "none",
        axis.text.x=element_text(angle=45,hjust=1)) +
  xlab("Clade removed") + ylab("estimate")

### Plot 2: Estimated slopes against for null distribution/removed clade ~ trees 
cols <- c("Without clade" = "red", "Full data" = "black", "Null distribution" = "lightblue")
g2 <- 
  ggplot(nd.c) +       
  geom_jitter(width = .35, aes(y = estimate, x = as.factor(iteration), color = "Null distribution")) +
  geom_line(data = f[f.ord, ], aes(y = Estimate, x = as.factor(iteration), color = "Full data",  group = 1), size = 1) +
  geom_point(data = f[f.ord, ] , aes(y = Estimate, x = as.factor(iteration),  color = "Full data", group = 1), size = 1) +
  geom_line(data = ces.c, aes(y = estimate, x = as.factor(iteration), color = "Without clade", group = 1), size = 1) + 
  geom_point(data = ces.c, aes(y = estimate, x = as.factor(iteration), color = "Without clade", group = 1), size = 1) + 
 # geom_point(data = nd.c.out , aes(y = estimate, x = as.factor(iteration)), color = "blue",  size = .5) +
  
  theme_bw() + 
  scale_color_manual(name = "", values = cols) +
  theme(legend.position = c(.15, .95),
        legend.direction = "vertical",
        legend.text=element_text(size=12),
        axis.text=element_text(size=10),
        axis.title=element_text(size=12),
        legend.key.width=unit(.5,"line"),
        legend.key.size = unit(.5,"cm"),
        legend.background = element_blank(),
        panel.background = element_rect(fill="white",
                                        colour="black"),
        axis.line = element_line(colour = "black"),
        panel.grid = element_blank()) +
  xlab("Tree") + ylab("estimate")+
  ggtitle(paste("Clade = ", clade, "| ", "n.sim = ", x$call$n.sim, " | ",
                " n.tree = ", x$call$n.tree,
                "| Sig. iterations =", s.est[s.est$clade.removed == clade, ]$`Significant (%)`, "%")) 
if(length(levels(ces.c$iteration)) > 30) g2 <- g2 +  theme(axis.text.x = element_blank()) 

### Output-------------------------------------------------------------------------------
### Plotting:
if (graphs=="all")
  return(suppressMessages(multiplot(g1,g2, cols = 2)))
if (graphs==1)
  return(suppressMessages(g1))
if (graphs==2)
  return(suppressMessages(g2))
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
#' To check the available clades to plot, see \code{x$sensi.estimates$clade} 
#' in the object returned from \code{tree_clade_phylm} or \code{tree_clade_phyglm}. 
#' 
#' Graph 2: Distribution of the simulated slopes (Null distribution
#' for a given clade sample size).
#' The red dashed line represents the estimated estimate for the reduced model 
#' (without the focal clade) and the black line represents the estimate for the 
#' full model.
#'  
#' @importFrom ggplot2 aes_string
#' @importFrom stats model.frame qt plogis 
#' @export
sensi_plot.sensiIntra_Clade <- function(x, clade = NULL, graphs = "all", ...){
  sensi_plot.sensiIntra_Clade(x, clade = NULL, graphs, ...)
}