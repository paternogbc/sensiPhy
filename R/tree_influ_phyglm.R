#' Influential species detection and phylogenetic uncertainty - Phylogenetic Logistic Regression
#'
#' Performs leave-one-out deletion analyis for phylogenetic logistic regression,
#' and detects influential species while evaluating uncertainty in trees topology.
#'
#' @param formula The model formula
#' @param data Data frame containing species traits with row names matching tips
#' in \code{phy}.
#' @param phy A phylogeny (class 'phylo') matching \code{data}.
#' @param btol Bound on searching space. For details see \code{phyloglm}
#' @param cutoff The cutoff value used to identify for influential species
#' (see Details)
#' @param n.tree Number of times to repeat the analysis with n different trees picked 
#' randomly in the multiPhylo file.
#' @param track Print a report tracking function progress (default = TRUE)
#' @param ... Further arguments to be passed to \code{phyloglm}
#' @details
#' This function sequentially removes one species at a time, fits a phylogenetic
#' logistic regression model using \code{\link[phylolm]{phyloglm}}, stores the
#' results and detects influential species. It repeats this operation using n trees, 
#' randomly picked in a multiPhylo file.
#'
#' Currently only logistic regression using the "logistic_MPLE"-method from
#' \code{phyloglm} is implemented.
#'
#' \code{influ_phyglm} detects influential species based on the standardised
#' difference in intercept and/or slope when removing a given species compared
#' to the full model including all species. Species with a standardised difference
#' above the value of \code{cutoff} are identified as influential. The default
#' value for the cutoff is 2 standardised differences change.
#'
#' Currently, this function can only implement simple logistic models (i.e. \eqn{trait~
#' predictor}). In the future we will implement more complex models.
#'
#' Output can be visualised using \code{sensi_plot}.
#'
#' @return The function \code{influ_phyglm} returns a list with the following
#' components:
#' @return \code{cutoff}: The value selected for \code{cutoff}
#' @return \code{formula}: The formula
#' @return \code{full.model.estimates}: Coefficients, aic and the optimised
#' value of the phylogenetic parameter (i.e. \code{alpha}) for the full model
#' without deleted species.
#' @return \code{influential_species}: List of influential species, both
#' based on standardised difference in interecept and in the slope of the
#' regression. Species are ordered from most influential to less influential and
#' only include species with a standardised difference > \code{cutoff}.
#' @return \code{sensi.estimates}: A data frame with all simulation
#' estimates. Each row represents a deleted clade. #' Columns report the calculated
#' regression intercept (\code{intercept}), difference between simulation
#' intercept and full model intercept (\code{DIFintercept}), the standardised
#' difference (\code{sDIFintercept}), the percentage of change in intercept compared
#' to the full model (\code{intercept.perc}) and intercept p-value
#' (\code{pval.intercept}). All these parameters are also reported for the regression
#' slope (\code{DIFestimate} etc.). Additionally, model aic value (\code{AIC}) and
#' the optimised value (\code{optpar}) of the phylogenetic parameter
#' (i.e. \code{alpha}) are reported.
#' @return \code{data}: Original full dataset.
#' @return \code{errors}: Species where deletion resulted in errors.
#' @author Gustavo Paterno, Caterina Penone & Gijsbert D.A. Werner
#' @seealso \code{\link[phylolm]{phyloglm}}, \code{\link{samp_phyglm}},
#' \code{\link{influ_phylm}}, \code{\link{sensi_plot}}
#' @references Ho, L. S. T. and Ane, C. 2014. "A linear-time algorithm for 
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
#'influ <- tree_influ_phyglm(y ~ x, data = dat, phy = mphy, n.tree = 20)
#'summary(influ)
#'sensi_plot(influ)
#' @export

tree_influ_phyglm <- function(formula, data, phy, n.tree = 2, 
                                         cutoff = 2, btol = 50, track = TRUE,...) {
  # Error checking:
  if(!is.data.frame(data)) stop("data must be class 'data.frame'")
  if(class(formula)!="formula") stop("formula must be class 'formula'")
  if(class(phy)!="multiPhylo") stop("phy must be class 'multiPhylo'")
  if(length(phy)<n.tree) stop("'times' must be smaller (or equal) than the number of trees in the 'multiPhylo' object")
  else
  
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
  counter = 1
  
  for (j in trees){
    #Match data order to tip order
    full.data <- full.data[phy[[j]]$tip.label,]
    
    #Select tree
    tree <- phy[[j]]
    
    tree.influ[[counter]] <- influ_phyglm(formula, data = full.data, phy=tree,
                                   verbose = FALSE, track = FALSE,...)
    
    counter = counter + 1
  }
  

  names(tree.influ) <- trees
  
  # Merge lists into data.frames between iterations:
  full.estimates  <- suppressWarnings(recombine(tree.influ, slot1 = 3, slot2 = 1))
  
  #influ species slope
  influ.sp.estimate <- (lapply(tree.influ,function(x) x$influential.species$influ.sp.estimate))
  influ.sp.estimate <- as.data.frame(as.matrix(influ.sp.estimate))
  names(influ.sp.estimate) <- "influ.sp.estimate"
  influ.sp.estimate$tree<-row.names(influ.sp.estimate)
  
  #influ species intercept
  influ.sp.intercept <- (lapply(tree.influ,function(x) x$influential.species$influ.sp.intercept))
  influ.sp.intercept <- as.data.frame(as.matrix(influ.sp.intercept))
  names(influ.sp.intercept) <- "influ.sp.intercept"
  influ.sp.intercept$tree<-row.names(influ.sp.intercept)
  
  #influ.estimates
  influ.estimates <- recombine(tree.influ, slot1 = 5)
  
  #Generates output:
  res <- list(call = match.call(),
              cutoff=cutoff,
              formula = formula,
              full.model.estimates = full.estimates,
              influential.species = list(influ.sp.estimate=influ.sp.estimate,influ.sp.intercept=influ.sp.intercept),
              sensi.estimates = influ.estimates,
              data = full.data)
  
  
  class(res) <- c("sensiTree_Influ","sensiTree_InfluL")
  
  ### Warnings:
  if (length(res$errors) >0){
    warning("Some species deletion presented errors, please check: output$errors")}
  else {
    res$errors <- "No errors found."
  }
  return(res)
}