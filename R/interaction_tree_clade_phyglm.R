#' Influential clade detection and phylogenetic uncertainty - Phylogenetic Logistic Regression
#'
#' Estimate the impact on model estimates of phylogenetic logistic regression after 
#' removing clades from the analysis and evaluating uncertainty in trees topology.
#'
#' @param formula The model formula
#' @param data Data frame containing species traits with row names matching tips
#' in \code{phy}.
#' @param phy A phylogeny (class 'multiPhylo', see ?\code{ape}).
#' @param clade.col The name of a column in the provided data frame with clades 
#' specification (a character vector with clade names).
#' @param n.species Minimum number of species in the clade in order to include
#' this clade in the leave-one-out deletion analyis. Default is \code{5}.
#' @param times.clade Number of simulations for the randomization test.
#' @param times.tree Number of times to repeat the analysis with n different trees picked 
#' randomly in the multiPhylo file.
#' If NULL, \code{times} = 2
#' @param btol Bound on searching space. For details see \code{phyloglm}
#' @param model The phylogenetic model to use (see Details). Default is \code{lambda}.
#' @param track Print a report tracking function progress (default = TRUE)
#' @param ... Further arguments to be passed to \code{phyloglm}
#' @details
#' Currently only logistic regression using the "logistic_MPLE"-method from
#' \code{phyloglm} is implemented.
#' 
#' This function sequentially removes one clade at a time, fits a phylogenetic
#' logistic regression model using \code{\link[phylolm]{phyloglm}} and stores the
#' results. The impact of of a specific clade on model estimates is calculated by the
#' comparison between the full model (with all species) and the model without 
#' the species belonging to a clade. It repeats this operation using n trees, 
#' randomly picked in a multiPhylo file.
#' 
#' Additionally, to account for the influence of the number of species on each 
#' clade (clade sample size), this function also estimates a null distribution of slopes
#' expected for the number of species in a given clade. This is done by fitting
#' models without the same number of species in the given clade. 
#' The number of simulations to be performed is set by 'times'. To test if the 
#' clade influence differs from the null expectation, a randomization test can
#' be performed using 'summary(x)'. 
#'
#'
#' \code{clade_phyglm} detects influential clades based on
#' difference in intercept and/or slope when removing a given clade compared
#' to the full model including all species. This is done for n trees in the multiphylo file.
#' 
#' Currently, this function can only implements simple logistic models (i.e. \eqn{trait~
#' predictor}). In the future we will implement more complex models.
#'
#' Output can be visualised using \code{sensi_plot}.
#'
#' @return The function \code{clade_phyglm} returns a list with the following
#' components:
#' @return \code{formula}: The formula
#' @return \code{full.model.estimates}: Coefficients, aic and the optimised
#' value of the phylogenetic parameter (e.g. \code{lambda}) for the full model
#' without deleted species.
#' @return \code{clade.model.estimates}: A data frame with all simulation
#' estimates. Each row represents a deleted clade. Columns report the calculated
#' regression intercept (\code{intercept}), difference between simulation
#' intercept and full model intercept (\code{DFintercept}), the percentage of change
#' in intercept compared to the full model (\code{intercept.perc}) and intercept
#' p-value (\code{pval.intercept}). All these parameters are also reported for the regression
#' slope (\code{DFslope} etc.). Additionally, model aic value (\code{AIC}) and
#' the optimised value (\code{optpar}) of the phylogenetic parameter 
#' (e.g. \code{kappa} or \code{lambda}, depending on the phylogenetic model used) 
#' are reported.
#' @return \code{data}: Original full dataset.
#' @return \code{errors}: Clades and/or trees where deletion resulted in errors.
#' @author Gustavo Paterno, Caterina Penone & Gijsbert D.A. Werner
#' @seealso \code{\link[phylolm]{phylolm}}, \code{\link[sensiPhy]{samp_phyglm}},
#'  \code{\link{influ_phyglm}}, \code{\link{sensi_plot}}
#' \code{\link{sensi_plot}}
#' @references Ho, L. S. T. and Ane, C. 2014. "A linear-time algorithm for 
#' Gaussian and non-Gaussian trait evolution models". Systematic Biology 63(3):397-408.
#' @examples 
#' \dontrun{
#' # Simulate Data:
#' set.seed(6987)
#' mphy = rmtree(150, N = 30)
#' x = rTrait(n=1,phy=mphy[[1]])
#' X = cbind(rep(1,150),x)
#' y = rbinTrait(n=1,phy=mphy[[1]], beta=c(-1,0.5), alpha=.7 ,X=X)
#' cla <- rep(c("A","B","C","D","E"), each = 30)
#' dat = data.frame(y, x, cla)
#' # Run sensitivity analysis:
#' clade <- interaction_tree_clade_phyglm(y ~ x, phy = mphy, data = dat, times.tree = 3, times.clade = 5, clade.col = "cla")
#'# To check summary results and most influential clades:
#'summary(clade)
#'# Visual diagnostics for clade removal:
#'sensi_plot(clade)
#'# Specify which clade removal to plot:
#'sensi_plot(clade, "B")
#'sensi_plot(clade, "C")
#'sensi_plot(clade, "D") #The clade with the largest effect on slope and intercept
#'}
#' @export

interaction_tree_clade_phyglm <- function(formula, data, phy, clade.col, n.species = 5, 
                                         times.clade = 100, times.tree = 2, btol=50, track = TRUE,...) {
  # Error checking:
  if(!is.data.frame(data)) stop("data must be class 'data.frame'")
  if(missing(clade.col)) stop("clade.col not defined. Please, define the",
                              " column with clade names.")
  if(class(formula)!="formula") stop("formula must be class 'formula'")
  if(class(phy)!="multiPhylo") stop("phy must be class 'multiPhylo'")
  if(length(phy)<times.tree) stop("'times' must be smaller (or equal) than the number of trees in the 'multiPhylo' object")
  
  #Match data and phy
  data_phy <- match_dataphy(formula, data, phy, ...)
  phy <- data_phy$phy
  full.data <- data_phy$data
  if (is.na(match(clade.col, names(full.data)))) {
    stop("Names column '", clade.col, "' not found in data frame'")
  }
  
  
  # If the class of tree is multiphylo pick n=times random trees
  trees<-sample(length(phy),times.tree,replace=F)
  
  
  # Identify CLADES to use and their sample size 
  all.clades <- levels(full.data[ ,clade.col])
  wc <- table(full.data[ ,clade.col]) > n.species
  uc <- table(full.data[ , clade.col])[wc]
  
  if (length(uc) == 0) stop(paste("There is no clade with more than ",
                                  n.species," species. Change 'n.species' to fix this
                                  problem",sep=""))
  
  #List to store information
  tree.clade <- list ()
  
  #Start tree loop here
  errors <- NULL
  counter = 1
  
  for (j in trees){
    
    #Match data order to tip order
    full.data <- full.data[phy[[j]]$tip.label,]
    
    #Select tree
    tree <- phy[[j]]
    
    tree.clade[[counter]] <- clade_phyglm(formula, data=full.data, phy=tree, btol, track = FALSE,
                                   clade.col, n.species, times.clade, verbose = FALSE, ...)
    
    counter = counter +1
  }
  

  names(tree.clade) <- trees
  
  # Merge lists into data.frames between iterations:
  full.estimates  <- recombine(tree.clade, slot1 = 4)
  clade.estimates <- recombine(tree.clade, slot1 = 5)
  null.dist       <- recombine(tree.clade, slot1 = 6)
  
  #Generate output:
  res <- list(call = match.call(),
              model = model,
              formula = formula,
              full.model.estimates = full.estimates,
              clade.model.estimates = clade.estimates,
              null.dist = null.dist, 
              data = full.data,
              errors = errors,
              clade.col = clade.col)
  
  class(res) <- "sensiTree_Clade"
  
  ### Warnings:
  if (length(res$errors) >0){
    warning("Some clades deletion presented errors, please check: output$errors")}
  else {
    res$errors <- "No errors found."
  }
  return(res)
}

