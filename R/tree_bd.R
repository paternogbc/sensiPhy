#' Phylogenetic uncertainty - net diversification rate
#'
#' Performs estimates of diversification rate evaluating
#' uncertainty in trees topology.
#'
#' @param phy A phylogeny (class 'multiPhylo', see ?\code{ape}).
#' @param n.tree Number of times to repeat the analysis with n different trees picked 
#' randomly in the multiPhylo file. (If \code{n.tree} = "all", diversification will be estimated
#' among the set of trees provided in \code{phy})
#' @param method the method for estimating diversification rate ("ms" or "km") (see Details). 
#' @param track Print a report tracking function progress (default = TRUE)
#' @param ... Further arguments to be passed to \code{phylosig}
#' @details
#' This function estimates net diversification rate using \code{\link[geiger]{bd.ms}} 
#' (Magallon and Sanderson (2000) method) or speciation rate using \code{\link[geiger]{bd.km}}
#' (Kendall-Moran method) for n trees, randomly picked from a multiPhylo file.
#'
#' Output can be visualised using \code{sensi_plot}.
#'
#' @return The function \code{tree_bd} returns a list with the following
#' components:
#' @return \code{tree.bd.estimates}: Three number, diversification/speciation rate estimate 
#' ("Magallon and Sanderson" or "Kendall-Moran") for each run with a different phylogenetic tree.
#' @return \code{stats}: Main statistics for estimates across trees.\code{CI_low} and \code{CI_high} are the lower 
#' and upper limits of the 95% confidence interval.
#' @author Gustavo Paterno
#' @seealso \code{\link[geiger]{bd.ms}},
#' \code{\link{tree_phylm}},\code{\link{sensi_plot}}
#' @references 
#' Magallon S and MJ Sanderson. 2000. Absolute diversification rates in 
#' angiosperm clades. Evolution 55:1762-1780.
#' 
#' @examples 
#'data("primates")
#'# To estimate diversification rate with Magallon and Sanderson method:
#'fit <- tree_bd(phy = primates.phy, n.tree = 30, method = "ms")
#'fit$stats
#'# To estimate speciation rate Kendall-Moran method
#'fit <- tree_bd(phy = primates.phy, n.tree = 30, method = "km")
#'fit$stats
#' @export
#'@importFrom geiger bd.ms bd.km
tree_bd <- function(phy, n.tree = "all", method = "ms", track = F, ...){
  
  #Error check
  if(!method %in% c("ms", "km")) stop("method must be 'ms' or 'km'")
  if (n.tree == "all") n.tree <- length(phy)
  if(class(phy)!="multiPhylo") stop("phy must be class 'multiPhylo'")
  if(length(phy)<n.tree) stop("'n.tree' must be smaller (or equal) than the number
                              of trees in the 'multiPhylo' object")
  
  # Pick n=n.tree random trees or all
  trees <- sample(length(phy), n.tree, replace = F)
  
  #Create the results data.frame
  tree.bd.estimates <- data.frame("n.tree" = numeric(),
                                  "estimate" = numeric())
  
  #Model calculation
  counter = 1
  if(track==TRUE) pb <- utils::txtProgressBar(min = 0, max = n.tree, style = 3)
  
  for (j in trees){
    
    if (method == "ms") { mod.s <- geiger::bd.ms(phy = phy[[j]]) }
    if (method == "km") { mod.s <- geiger::bd.km(phy = phy[[j]]) }
    
    
    estimate <- mod.s
    
    if(track==TRUE) (utils::setTxtProgressBar(pb, counter))
    
    #write in a table
    estim.simu <- data.frame(j, estimate)
    tree.bd.estimates[counter, ]  <- estim.simu
    counter = counter + 1
    
  }
  
  if(track==TRUE) on.exit(close(pb))
  #calculate mean and sd for each parameter
  #mean_by_tree <- stats::aggregate(. ~ n.tree, data = tree.physig.estimates, mean)
  
  statresults <- data.frame(min = apply(tree.bd.estimates, 2, min),
                            max = apply(tree.bd.estimates, 2, max),
                            mean = apply(tree.bd.estimates, 2, mean),
                            sd_tree = apply(tree.bd.estimates, 2, stats::sd))[-1, ]
  
  statresults$CI_low  <- statresults$mean - qt(0.975, df = n.tree-1) * statresults$sd_tree / sqrt(n.tree)
  statresults$CI_high <- statresults$mean + qt(0.975, df = n.tree-1) * statresults$sd_tree / sqrt(n.tree)
  
  stats <- round(statresults[c(1),c(3,5,6,1,2)],digits=7)
  
  cl <- match.call()
  res <- list(   call = cl,
                 tree.bd.estimates = tree.bd.estimates,
                 stats = stats)
  class(res) <- "tree.bd"
  return(res)
}

