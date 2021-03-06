#' Interaction between phylogenetic uncertainty and sensitivity to species sampling - Phylogenetic Linear Regression
#'
#' Performs analyses of sensitivity to species sampling by randomly removing
#' species and detecting the effects on parameter estimates in a phylogenetic
#' linear regression, while evaluating uncertainty in trees topology.
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
#' @param model The phylogenetic model to use (see Details). Default is \code{lambda}.
#' @param track Print a report tracking function progress (default = TRUE)
#' @param ... Further arguments to be passed to \code{phylolm}
#' @details
#'
#' This function randomly removes a given percentage of species (controlled by
#' \code{breaks}) from the full phylogenetic linear regression, fits a phylogenetic
#' linear regression model without these species using \code{\link[phylolm]{phylolm}},
#' repeats this many times (controlled by \code{n.sim}), stores the results and
#' calculates the effects on model parameters. It repeats this operation using n trees, 
#' randomly picked in a multiPhylo file.
#'
#' All phylogenetic models from \code{phylolm} can be used, i.e. \code{BM},
#' \code{OUfixedRoot}, \code{OUrandomRoot}, \code{lambda}, \code{kappa},
#' \code{delta}, \code{EB} and \code{trend}. See ?\code{phylolm} for details.
#'
#' Currently, this function can only implement simple linear models (i.e. \eqn{trait~
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
#' removed) this reports the percentage of statistically significant (at p<0.05)
#' intercepts (\code{perc.sign.intercept}) over all repetitions as well as the
#' percentage of statisticaly significant (at p<0.05) slopes (\code{perc.sign.estimate}).
#' @return \code{data}: Original full dataset.
#' #' @note Please be aware that dropping species may reduce power to detect 
#' significant slopes/intercepts and may partially be responsible for a potential 
#' effect of species removal on p-values. Please also consult standardised differences
#' in the (summary) output.
#' @author Gustavo Paterno, Gijsbert D.A. Werner & Caterina Penone
#' @seealso \code{\link[phylolm]{phylolm}}, \code{\link{samp_phylm}},
#' \code{\link{tree_phylm}},\code{\link{tree_samp_phyglm}},\code{\link{sensi_plot}}
#' @references 
#' 
#' Paterno, G. B., Penone, C. Werner, G. D. A. 
#' \href{http://doi.wiley.com/10.1111/2041-210X.12990}{sensiPhy: 
#' An r-package for sensitivity analysis in phylogenetic 
#' comparative methods.} Methods in Ecology and Evolution 
#' 2018, 9(6):1461-1467
#' 
#' Werner, G.D.A., Cornwell, W.K., Sprent, J.I., Kattge, J. & Kiers, E.T. (2014).
#' A single evolutionary innovation drives the deep evolution of symbiotic N2-fixation
#' in angiosperms. Nature Communications, 5, 4087.
#'   
#' Ho, L. S. T. and Ane, C. 2014. "A linear-time algorithm for 
#' Gaussian and non-Gaussian trait evolution models". Systematic Biology 63(3):397-408.
#' 
#' @import ape phylolm
#' 
#' @examples 
#' \dontrun{
#' # Load data:
#' data(alien)
#' # Run analysis:
#' samp <- tree_samp_phylm(log(gestaLen) ~ log(adultMass), phy = alien$phy,
#'                                     data = alien$data, n.tree = 5, n.sim=10)
#' summary(samp)
#' head(samp$sensi.estimates)
#' # Visual diagnostics
#' sensi_plot(samp)
#' sensi_plot(samp, graphs = 1)
#' sensi_plot(samp, graphs = 2)
#' }
#' \dontshow{
#' # Load data:
#' data(alien)
#' # Run analysis:
#' samp <- tree_samp_phylm(log(gestaLen) ~ log(adultMass), phy = alien$phy,
#'                                     data = alien$data, n.tree = 1, n.sim=1)
#' summary(samp)
#' # Visual diagnostics
#' sensi_plot(samp)
#' }
#' @export


tree_samp_phylm <- function(formula, data, phy, n.sim = 30, n.tree = 2, 
                            breaks=seq(.1,.5,.1), model = "lambda", track = TRUE,...) {
  

  # Error checking:
  if (!inherits(data, "data.frame"))
    stop("data must be class 'data.frame'")
  if (!inherits(formula, "formula"))
    stop("formula must be class 'formula'")
  if (!inherits(phy, "multiPhylo"))
    stop("phy must be class 'multiPhylo'")
  if (length(phy) < n.tree)
    stop("'times' must be smaller (or equal) than the number of trees in the 'multiPhylo' object")
  if (length(breaks) < 2)
    stop("Please include more than one break, e.g. breaks=c(.3,.5)")
  if((model == "trend") && (sum(ape::is.ultrametric(phy))>1)) 
    stop("Trend is unidentifiable for ultrametric trees., see ?phylolm for details")
  else
  
  #Match data and phy
  data_phy <- match_dataphy(formula, data, phy, ...)
  phy <- data_phy$phy
  full.data <- data_phy$data
  
  # If the class of tree is multiphylo pick n=n.tree random trees
  trees<-sample(length(phy),n.tree,replace=F)
  
  
  #List to store information
  tree.samp <- list ()
  
  #Start tree loop here
  errors <- NULL
  if(track==TRUE) pb <- utils::txtProgressBar(min = 0, max = n.tree, style = 3)
  counter = 1
  
  for (j in trees){
    #Match data order to tip order
    full.data <- full.data[phy[[j]]$tip.label,]
    
    #Select tree
    tree <- phy[[j]]
    
    tree.samp[[counter]] <- samp_phylm(formula, data = full.data, phy=tree, n.sim = n.sim,
                                   model, breaks=breaks, track = FALSE, verbose = FALSE, ...)
    
    if(track==TRUE) utils::setTxtProgressBar(pb, counter)
    counter = counter + 1
  }
  
  if(track==TRUE) close(pb)
  names(tree.samp) <- trees
  
  # Merge lists into data.frames between iterations:
  full.estimates  <- suppressWarnings(recombine(tree.samp, slot1 = 4, slot2 = 1))
  samp.estimates <- recombine(tree.samp, slot1 = 5)
  samp.estimates$info <- NULL
  perc.sign <- recombine(tree.samp, slot1 = 6)
  perc.sign$info <- NULL

  
  #Generates output:
  res <- list(call = match.call(),
              model = model,
              formula = formula,
              full.model.estimates = full.estimates,
              sensi.estimates = samp.estimates,
              sign.analysis = perc.sign,
              data = full.data)
  

  class(res) <- "sensiTree_Samp"
  ### Warnings:
  if (length(res$errors) >0){
    warning("Some species deletion presented errors, please check: output$errors")}
  else {
    res$errors <- "No errors found."
  }
  return(res)
}