#' Interaction between phylogenetic uncertainty and sensitivity to species sampling  - Phylogenetic Logistic Regression
#' 
#' Performs analyses of sensitivity to species sampling by randomly removing
#' species and detecting the effects on parameter estimates in phylogenetic
#' logistic regression, while evaluating uncertainty in trees topology.
#'
#' @param formula The model formula
#' @param data Data frame containing species traits with row names matching tips
#' in \code{phy}.
#' @param phy A phylogeny (class 'phylo') matching \code{data}.
#' @param n.sim The number of times species are randomly deleted for each
#' \code{break}.
#' @param n.tree Number of times to repeat the analysis with n different trees picked 
#' randomly in the multiPhylo file.
#' @param breaks A vector containing the percentages of species to remove.
#' @param btol Bound on searching space. For details see \code{phyloglm}
#' @param track Print a report tracking function progress (default = TRUE)
#' @param ... Further arguments to be passed to \code{phyloglm}
#' @details
#'
#' This function randomly removes a given percentage of species (controlled by
#' \code{breaks}) from the full phylogenetic logistic regression, fits a phylogenetic
#' logistic regression model without these species using \code{\link[phylolm]{phyloglm}},
#' repeats this many times (controlled by \code{times}), stores the results and
#' calculates the effects on model parameters. It repeats this operation using n trees, 
#' randomly picked in a multiPhylo file.
#'
#' Only logistic regression using the "logistic_MPLE"-method from
#' \code{phyloglm} is implemented.
#'
#' Currently, this function can only implement simple logistic models (i.e. \eqn{trait~
#' predictor}). In the future we will implement more complex models.
#'
#' Output can be visualised using \code{sensi_plot}.
#' @return The function \code{samp_phylm} returns a list with the following
#' components:
#' @return \code{formula}: The formula
#' @return \code{full.model.estimates}: Coefficients, aic and the optimised
#' value of the phylogenetic parameter (e.g. \code{lambda} or \code{kappa}) for
#' the full model without deleted species.
#' @return \code{sensi.estimates}: A data frame with all simulation
#' estimates. Each row represents a model rerun with a given number of species
#' \code{n.remov} removed, representing \code{n.percent} of the full dataset.
#' Columns report the calculated regression intercept (\code{intercept}),
#' difference between simulation intercept and full model intercept (\code{DIFintercept}),
#' the percentage of change in intercept compared to the full model (\code{intercept.perc})
#' and intercept p-value (\code{pval.intercept}). All these parameters are also reported
#' for the regression slope (\code{DIFestimate} etc.). Additionally, model aic value
#' (\code{AIC}) and the optimised value (\code{optpar}) of the phylogenetic
#' parameter (e.g. \code{kappa} or \code{lambda}, depending on the phylogenetic model
#' used) are reported. Lastly we reported the standardised difference in intercept 
#' (\code{sDIFintercept}) and slope (\code{sDIFestimate}). 
#' @return \code{sign.analysis} For each break (i.e. each percentage of species
#' removed) this reports the percentage of statistically signficant (at p<0.05)
#' intercepts (\code{perc.sign.intercept}) over all repititions as well as the
#' percentage of statisticaly significant (at p<0.05) slopes (\code{perc.sign.estimate}).
#' @return \code{data}: Original full dataset.
#' @note Please be aware that dropping species may reduce power to detect 
#' significant slopes/intercepts and may partially be responsible for a potential 
#' effect of species removal on p-values. Please also consult standardised differences
#' in the (summary) output. 
#' @author Gustavo Paterno, Gijsbert D.A. Werner & Caterina Penone
#' @seealso \code{\link[phylolm]{phyloglm}}, \code{\link{samp_phylm}},
#' \code{\link{influ_phyglm}}, \code{\link{sensi_plot}}
#' @references 
#' Werner, G.D.A., Cornwell, W.K., Sprent, J.I., Kattge, J. & Kiers, E.T. (2014).
#' A single evolutionary innovation drives the deep evolution of symbiotic N2-fixation
#' in angiosperms. Nature Communications, 5, 4087.
#'   
#' #' Ho, L. S. T. and Ane, C. 2014. "A linear-time algorithm for 
#' Gaussian and non-Gaussian trait evolution models". Systematic Biology 63(3):397-408.
#' @examples
#'# Simulate Data:
#'set.seed(6987)
#'mphy = rmtree(100, N = 30)
#'x = rTrait(n=1,phy=mphy[[1]])
#'X = cbind(rep(1,100),x)
#'y = rbinTrait(n=1,phy=mphy[[1]], beta=c(-1,0.5), alpha=.7 ,X=X)
#'dat = data.frame(y, x)
#'# Run sensitivity analysis:
#'tree_samp <- tree_samp_phyglm(y ~ x, data = dat, phy = mphy, n.tree = 3, n.sim=10) 
#'summary(tree_samp)
#'sensi_plot(tree_samp)
#'sensi_plot(tree_samp, graphs = 1)
#'sensi_plot(tree_samp, graphs = 2)
#' @export

tree_samp_phyglm <- function(formula, data, phy, n.sim = 30, n.tree = 2, 
                             breaks=seq(.1, .5, .1), btol = 50, track = TRUE,...) {
  
  # Error checking:
  if(!is.data.frame(data)) stop("data must be class 'data.frame'")
  if(class(formula)!="formula") stop("formula must be class 'formula'")
  if(class(phy)!="multiPhylo") stop("phy must be class 'multiPhylo'")
  if(length(phy)<n.tree) stop("'times' must be smaller (or equal) than the number of trees in the 'multiPhylo' object")
  if(length(breaks) < 2) stop("Please include more than one break, e.g. breaks=c(.3,.5)")
  
  #Match data and phy
  data_phy <- match_dataphy(formula, data, phy, ...)
  phy <- data_phy$phy
  full.data <- data_phy$data

  # If the class of tree is multiphylo pick n=n.tree random trees
  trees<-sample(length(phy),n.tree,replace=F)
  
  
  #List to store information
  tree.influ <- list ()
  
  #Start tree loop here
  errors <- NULL
  if(track==TRUE) pb <- utils::txtProgressBar(min = 0, max = n.tree, style = 3)
  counter = 1
  
  for (j in trees){
    #Match data order to tip order
    full.data <- full.data[phy[[j]]$tip.label,]
    
    #Select tree
    tree <- phy[[j]]
    
    tree.influ[[counter]] <- samp_phyglm(formula, data = full.data, phy=tree, n.sim = n.sim,
                                    breaks=breaks, btol = btol, track = FALSE, verbose = FALSE, ...)
    
    if(track==TRUE) utils::setTxtProgressBar(pb, counter)
    counter = counter + 1
  }
  
  if(track==TRUE) close(pb)
  names(tree.influ) <- trees
  
  # Merge lists into data.frames between iterations:
  full.estimates  <- suppressWarnings(recombine(tree.influ, slot1 = 4, slot2 = 1))
  influ.estimates <- recombine(tree.influ, slot1 = 5)
  influ.estimates$info <- NULL
  perc.sign <- recombine(tree.influ, slot1 = 6)
  perc.sign$info <- NULL
  
  
  #Generates output:
  res <- list(call = match.call(),
              model = "logistic_MPLE",
              formula = formula,
              full.model.estimates = full.estimates,
              sensi.estimates = influ.estimates,
              sign.analysis = perc.sign,
              data = full.data)
  
  
  class(res) <- c("sensiTree_Samp","sensiTree_SampL")
  ### Warnings:
  if (length(res$errors) >0){
    warning("Some species deletion presented errors, please check: output$errors")}
  else {
    res$errors <- "No errors found."
  }
  return(res)
}