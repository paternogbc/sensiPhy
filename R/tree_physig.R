#' Phylogenetic uncertainty - Phylogenetic signal
#'
#' Performs phylogenetic signal estimates evaluating
#' uncertainty in trees topology.
#'
#' @param trait.col The name of a column in the provided data frame with trait 
#'  to be analyzed  (e.g. "Body_mass").
#' @param data Data frame containing species traits with row names matching tips
#' in \code{phy}.
#' @param phy A phylogeny (class 'phylo') matching \code{data}.
#' @param method Method to compute signal: can be "K" or "lambda".
#' @param n.tree Number of times to repeat the analysis with n different trees picked 
#' randomly in the multiPhylo file. (If \code{n.tree} = "all", phylosgentic signal will be estimated
#' among the all set of trees provided in \code{phy})
#' @param track Print a report tracking function progress (default = TRUE)
#' @param ... Further arguments to be passed to \code{phylosig}
#' @details
#' This function estimates phylogenetic signal using \code{\link[phytools]{phylosig}}
#' to n trees, randomly picked in a multiPhylo file.
#'
#' Output can be visualised using \code{sensi_plot}.
#'
#' @return The function \code{tree_physig} returns a list with the following
#' components:
#' @return \code{Trait}: Column name of the trait analysed
#' @return \code{data}: Original full dataset
#' @return \code{tree.physig.estimates}: Three number, phylogenetic signal estimate 
#' (lambda or K) and the p-value for each run with a different phylogenetic tree.
#' @return \code{N.obs}: Size of the dataset after matching it with tree tips and removing NA's.
#' @return \code{stats}: Main statistics for phylogenetic estimates.\code{CI_low} and \code{CI_high} are the lower 
#' and upper limits of the 95% confidence interval.
#' @author Caterina Penone & Gustavo Paterno
#' @seealso \code{\link[phytools]{phylosig}}, 
#' \code{\link{tree_phylm}},\code{\link{sensi_plot}}
#' @references 
#' Donoghue, M.J. & Ackerly, D.D. (1996). Phylogenetic Uncertainties and 
#' Sensitivity Analyses in Comparative Biology. Philosophical Transactions:
#'  Biological Sciences, pp. 1241-1249.
#'  
#' Blomberg, S. P., T. Garland Jr., A. R. Ives (2003) 
#' Testing for phylogenetic signal in comparative data: 
#' Behavioral traits are more labile. Evolution, 57, 717-745.
#' 
#' Pagel, M. (1999) Inferring the historical patterns of biological evolution. 
#' Nature, 401, 877-884.
#' 
#' Kamilar, J. M., & Cooper, N. (2013). Phylogenetic signal in primate behaviour,
#'  ecology and life history. Philosophical Transactions of the Royal Society 
#'  B: Biological Sciences, 368: 20120341.
#'  
#' @examples 
#'# Load data:
#'data(alien)
#'# Logtransform data
#'alien.data$logMass <- log(alien.data$adultMass) 
#'# Run sensitivity analysis:
#'tree <- tree_physig(trait.col = "logMass", data = alien.data, 
#'phy = alien.phy, n.tree = 10)
#'summary(tree)
#'sensi_plot(tree)
#'sensi_plot(tree, graphs = 1)
#'sensi_plot(tree, graphs = 2)
#' @export
tree_physig <- function(trait.col, data, phy, n.tree = "all", method = "K", track = TRUE, ...){

  #Error check
  if (n.tree == "all") n.tree <- length(phy)
  if(class(data)!="data.frame") stop("data must be class 'data.frame'")
  if(class(phy)!="multiPhylo") stop("phy must be class 'multiPhylo'")
  if(length(phy)<n.tree) stop("'n.tree' must be smaller (or equal) than the number of trees in the 'multiPhylo' object")
  
  # Check match between data and phy 
  datphy <- match_dataphy(get(trait.col) ~ 1, data, phy)
  full.data <- datphy$data
  phy <- datphy$phy
  trait     <- full.data[[trait.col]]
  names(trait)  <- phy[[1]]$tip.label
  N <- nrow(full.data)
  
  
    # Pick n=n.tree random trees or all
    trees <- sample(length(phy), n.tree, replace = F)
    #Create the results data.frame
    tree.physig.estimates <- data.frame("n.tree" = numeric(), "estimate" = numeric(),
                                        "pval" = numeric())
    
    #Model calculation
    counter = 1
    if(track==TRUE) pb <- utils::txtProgressBar(min = 0, max = n.tree, style = 3)
    for (j in trees){
      
      mod.s    <- phytools::phylosig(tree = phy[[j]], x = trait, method = method, test = TRUE, ...)
      estimate <- mod.s[[1]]
      pval     <- mod.s$P
      
      if(track==TRUE) (utils::setTxtProgressBar(pb, counter))
      #write in a table
      estim.simu <- data.frame(j, estimate, pval)
      tree.physig.estimates[counter, ]  <- estim.simu
      counter = counter + 1
      
    }
    
    if(track==TRUE) on.exit(close(pb))
    #calculate mean and sd for each parameter
    #mean_by_tree <- stats::aggregate(. ~ n.tree, data = tree.physig.estimates, mean)
    
    statresults <- data.frame(min = apply(tree.physig.estimates, 2, min),
                              max = apply(tree.physig.estimates, 2, max),
                              mean = apply(tree.physig.estimates, 2, mean),
                              sd_tree = apply(tree.physig.estimates, 2, stats::sd))[-1, ]
    
    statresults$CI_low  <- statresults$mean - qt(0.975, df = n.tree-1) * statresults$sd_tree / sqrt(n.tree)
    statresults$CI_high <- statresults$mean + qt(0.975, df = n.tree-1) * statresults$sd_tree / sqrt(n.tree)
    
    stats <- round(statresults[c(1:2),c(3,5,6,1,2)],digits=5)
    
    cl <- match.call()
    res <- list(   call = cl,
                   Trait = trait.col,
                   tree.physig.estimates = tree.physig.estimates,
                   N.obs = N,
                   stats = stats,
                   data = full.data)
    class(res) <- "tree.physig"
    return(res)
  }

