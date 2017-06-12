#' Influential clade detection and phylogenetic uncertainty - Phylogenetic Linear Regression
#'
#' Estimate the impact on model estimates of phylogenetic linear regression after 
#' removing clades from the analysis and evaluating uncertainty in trees topology.. 
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
#' @param model The phylogenetic model to use (see Details). Default is \code{lambda}.
#' @param track Print a report tracking function progress (default = TRUE)
#' @param ... Further arguments to be passed to \code{phylolm}
#' @details
#' This function sequentially removes one clade at a time, fits a phylogenetic
#' linear regression model using \code{\link[phylolm]{phylolm}} and stores the
#' results. The impact of of a specific clade on model estimates is calculated by the
#' comparison between the full model (with all species) and the model without 
#' the species belonging to a clade. It repeats this operation using n trees, 
#' randomly picked in a multiPhylo file.
#' 
#'  Additionally, to account for the influence of the number of species on each 
#'  clade (clade sample size), this function also estimate a null distribution of slopes
#'  expected for the number of species in a given clade. This is done by fitting
#'  models without the same number of species in the given clade. 
#'  The number of simulations to be performed is set by 'times'. To test if the 
#'  clade influence differs from the null expectation, a randomization test can
#'  be performed using 'summary(x)'. 
#'
#' All phylogenetic models from \code{phylolm} can be used, i.e. \code{BM},
#' \code{OUfixedRoot}, \code{OUrandomRoot}, \code{lambda}, \code{kappa},
#' \code{delta}, \code{EB} and \code{trend}. See ?\code{phylolm} for details.
#'
#' \code{clade_phylm} detects influential clades based on
#' difference in intercept and/or slope when removing a given clade compared
#' to the full model including all species. This is done for n trees in the multiphylo file.
#' 
#' Currently, this function can only implement simple linear models (i.e. 
#' \eqn{y = a + bx}). In the future we will implement more complex models.
#'
#' Output can be visualised using \code{sensi_plot}.
#'
#' @return The function \code{clade_phylm} returns a list with the following
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
#' @author Gustavo Paterno (clade) & Caterina Penone (tree)
#' @seealso \code{\link[phylolm]{phylolm}}, \code{\link[sensiPhy]{samp_phyglm}},
#'  \code{\link{influ_phylm}}, \code{\link{sensi_plot}}
#' \code{\link{sensi_plot}}
#' @references Ho, L. S. T. and Ane, C. 2014. "A linear-time algorithm for 
#' Gaussian and non-Gaussian trait evolution models". Systematic Biology 63(3):397-408.
#' @examples 
#' \dontrun{
#'# Load data:
#'data(primates)
#'# run analysis:
#'clade_tree <- clade_phylm(log(sexMaturity) ~ log(adultMass), 
#'phy = primates$phy, data = primates$data, clade.col = "family", times.clade = 30, times.tree = 5)
#'# To check summary results and most influential clades:
#'summary(clade_tree)
#'# Visual diagnostics for clade removal:
#'sensi_plot(clade_tree)
#'# Specify which clade removal to plot:
#'sensi_plot(clade_tree, "Cercopithecidae")
#'sensi_plot(clade_tree, "Cebidae")
#'}
#' @export

interaction_clade_tree_phylm <- function(formula, data, phy, clade.col, n.species = 5, 
                                         times.clade = 100, times.tree = 2, model = "lambda", track = TRUE,...) {
  # Error checking:
  if(!is.data.frame(data)) stop("data must be class 'data.frame'")
  if(missing(clade.col)) stop("clade.col not defined. Please, define the",
                              " column with clade names.")
  if(class(formula)!="formula") stop("formula must be class 'formula'")
  if(class(phy)!="multiPhylo") stop("phy must be class 'multiPhylo'")
  if(length(phy)<times.tree) stop("'times' must be smaller (or equal) than the number of trees in the 'multiPhylo' object")
  
  #Match data and phy
  data_phy <- match_dataphy(formula, data, phy)
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
  clade.tree <- list ()
  
  #Start tree loop here
  errors <- NULL
  pb <- utils::txtProgressBar(min = 0, max = length(uc)*times.clade*times.tree, style = 1)
  
  for (j in trees){
    
    #Match data order to tip order
    full.data <- full.data[phy[[j]]$tip.label,]
    
    #Select tree
    tree <- phy[[j]]
    
    clade.tree[[j]] <- clade_phylm(formula, data=full.data, phy=tree, model = "lambda", track = FALSE,
                         clade.col, n.species = 5, times.clade, ...)
    
  }
  
  on.exit(close(pb))
  
  #Do something with clade.tree
  
  #Generates output:
  res <- list()
  
  class(res) <- "sensiClade_Tree"
  ### Warnings:
  if (length(res$errors) >0){
    warning("Some clades deletion presented errors, please check: output$errors")}
  else {
    res$errors <- "No errors found."
  }
  return(res)
}

