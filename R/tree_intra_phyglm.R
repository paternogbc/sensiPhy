#' Interaction between phylogenetic uncertainty and intraspecific variability - Phylogenetic logistic Regression
#'
#' Performs Phylogenetic logistic regression evaluating
#' intraspecific variability in response and/or predictor variables
#' and uncertainty in trees topology.
#'
#' @param formula The model formula: \code{response~predictor}. 
#' @param data Data frame containing species traits and species names as row names.
#' @param phy A phylogeny (class 'phylo', see ?\code{ape}).
#' @param Vx Name of the column containing the standard deviation or the standard error of the predictor 
#' variable. When information is not available for one taxon, the value can be 0 or \code{NA}
#' @param x.transf Transformation for the predictor variable (e.g. \code{log} or \code{sqrt}). Please use this 
#' argument instead of transforming data in the formula directly (see also details below).
#' @param n.intra Number of times to repeat the analysis generating a random value for response and/or predictor variables.
#' If NULL, \code{n.intra} = 30
#' @param n.tree Number of times to repeat the analysis with n different trees picked 
#' randomly in the multiPhylo file.
#' If NULL, \code{n.tree} = 2
#' @param distrib A character string indicating which distribution to use to generate a random value for the response 
#' and/or predictor variables. Default is normal distribution: "normal" (function \code{\link{rnorm}}).
#' Uniform distribution: "uniform" (\code{\link{runif}})
#' Warning: we recommend to use normal distribution with Vx or Vy = standard deviation of the mean.
#' @param btol Bound on searching space. For details see \code{phyloglm}
#' @param track Print a report tracking function progress (default = TRUE)
#' @param ... Further arguments to be passed to \code{phyloglm}
#' @details
#' This function fits a phylogenetic logistic regression model using \code{\link[phylolm]{phyloglm}} to n trees (\code{n.tree}), 
#' randomly picked in a multiPhylo file. The regression is also repeated \code{n.intra} times.
#' At each iteration the function generates a random value for each row in the dataset using the standard deviation 
#' or errors supplied and assuming a normal or uniform distribution. To calculate means and se for your raw data, 
#' you can use the \code{summarySE} function from the package \code{Rmisc}.
#' 
#' #' All phylogenetic models from \code{phyloglm} can be used, i.e. \code{BM},
#' \code{OUfixedRoot}, \code{OUrandomRoot}, \code{lambda}, \code{kappa},
#' \code{delta}, \code{EB} and \code{trend}. See ?\code{phyloglm} for details.
#'
#' Currently, this function can only implement simple logistic models (i.e. \eqn{trait~
#' predictor}). In the future we will implement more complex models.
#'
#' Output can be visualised using \code{sensi_plot}.
#' 
#' @section Warning:  
#' When Vy or Vx exceed Y or X, respectively, negative (or null) values can be generated, this might cause problems
#' for data transformation (e.g. log-transformation). In these cases, the function will skip the simulation. This problem can
#' be solved by increasing \code{times}, changing the transformation type and/or checking the target species in output$sp.pb.
#'  
#' @return The function \code{tree_intra_phylm} returns a list with the following
#' components:
#' @return \code{formula}: The formula
#' @return \code{data}: Original full dataset
#' @return \code{sensi.estimates}: Coefficients, aic and the optimised value of the phylogenetic 
#' parameter (e.g. \code{lambda}) for each regression using a value in the interval of variation and 
#' a different phylogenetic tree.
#' @return \code{N.obs}: Size of the dataset after matching it with tree tips and removing NA's.
#' @return \code{stats}: Main statistics for model parameters.\code{CI_low} and \code{CI_high} are the lower 
#' and upper limits of the 95% confidence interval.
#' @return \code{all.stats}: Complete statistics for model parameters.
#' Fields coded using \code{all} describe statistics due to both intraspecific variation and phylogenetic uncertainty.
#' Fields coded using \code{intra} describe statistics due to intraspecific variation only.
#' Fields coded using \code{tree} describe statistics due to phylogenetic uncertainty only.
#' \code{sd} is the standard deviation. \code{CI_low} and \code{CI_high} are the lower and upper limits 
#' of the 95% confidence interval.
#' @return \code{sp.pb}: Species that caused problems with data transformation (see details above).
#' 
#' @author Caterina Penone & Pablo Ariel Martinez
#' @seealso \code{\link[phylolm]{phyloglm}}, \code{\link{tree_phyglm}}, \code{\link{intra_phyglm}},
#' \code{\link{tree_intra_phylm}}, \code{\link{sensi_plot}}
#' @references 
#' 
#' Paterno, G. B., Penone, C. Werner, G. D. A. 
#' \href{http://doi.wiley.com/10.1111/2041-210X.12990}{sensiPhy: 
#' An r-package for sensitivity analysis in phylogenetic 
#' comparative methods.} Methods in Ecology and Evolution 
#' 2018, 9(6):1461-1467
#'
#' Martinez, P. a., Zurano, J.P., Amado, T.F., Penone, C., Betancur-R, R., 
#' Bidau, C.J. & Jacobina, U.P. (2015). Chromosomal diversity in tropical reef 
#' fishes is related to body size and depth range. Molecular Phylogenetics and 
#' Evolution, 93, 1-4
#' 
#' Ho, L. S. T. and Ane, C. 2014. "A linear-time algorithm for 
#' Gaussian and non-Gaussian trait evolution models". Systematic Biology 63(3):397-408.
#' @examples 
#'# Simulate data
#'set.seed(6987)
#'mphy = ape::rmtree(150, N = 30)
#'x = phylolm::rTrait(n=1,phy=mphy[[1]])
#'x_sd = rnorm(150,mean = 0.8,sd=0.2)
#'X = cbind(rep(1,150),x)
#'y = rbinTrait(n=1,phy=mphy[[1]], beta=c(-1,0.5), alpha=.7 ,X=X)
#'dat = data.frame(y, x, x_sd)
# Run sensitivity analysis:
#'intra.tree <- tree_intra_phyglm(y ~ x, data = dat, phy = mphy, n.intra = 3, 
#'                                            n.tree = 3, Vx = "x_sd")
#'# summary results:
#'summary(intra.tree)
#'# Visual diagnostics for phylogenetic uncertainty:
#'sensi_plot(intra.tree, uncer.type = "all") #or uncer.type = "tree", uncer.type = "intra"
#' @export


tree_intra_phyglm <- function(formula,
                              data,
                              phy,
                              Vx = NULL,
                              x.transf = NULL,
                              n.intra = 10,
                              n.tree = 2,
                              distrib = "normal",
                              track = TRUE,
                              btol = 50,
                              ...) {
  #Error check
  if (is.null(Vx))
    stop("Vx must be defined")
  if (!inherits(formula, "formula"))
    stop("formula must be class 'formula'")
  if (!inherits(data, "data.frame"))
    stop("data must be class 'data.frame'")
  if (!inherits(phy, "multiPhylo"))
    stop("phy must be class 'multiPhylo'")
  if (formula[[2]] != all.vars(formula)[1] ||
      formula[[3]] != all.vars(formula)[2])
    stop("Please use argument x.transf for data transformation")
  if (length(phy) < n.tree)
    stop("'n.tree' must be smaller (or equal) than the number of trees in the 'multiPhylo' object")
  if (distrib == "normal")
    warning("distrib = normal: make sure that standard deviation is provided for Vx and/or Vy")
  
  
  #Matching tree and phylogeny using utils.R
  datphy <- match_dataphy(formula, data, phy)
  full.data <- datphy[[1]]
  phy <- datphy[[2]]
  
  # If the class of tree is multiphylo pick n=n.tree random trees
  trees <- sample(length(phy), n.tree, replace = F)
  
  #Model calculation
  tree.intra <- list()
  if (track == TRUE)
    pb <- utils::txtProgressBar(min = 0, max = n.tree, style = 3)
  counter = 1
  
  for (j in trees) {
    #Match data order to tip order
    full.data <- full.data[phy[[j]]$tip.label, ]
    
    #Select tree
    tree <- phy[[j]]
    
    #model (remove warnings about standard deviation in intra)
    
    withCallingHandlers(
      tree.intra[[counter]] <-
        intra_phyglm(
          formula = formula,
          data = full.data,
          phy = tree,
          Vx,
          x.transf,
          n.intra = n.intra,
          distrib = distrib,
          btol = btol,
          track = F,
          verbose = F
        ),
      
      warning = function(w) {
        if (grepl("make sure that standard deviation", w$message))
          invokeRestart("muffleWarning")
      }
    )
    
    if (track == TRUE)
      utils::setTxtProgressBar(pb, counter)
    
    counter = counter + 1
    
  }
  
  if (track == TRUE)
    close(pb)
  names(tree.intra) <- trees
  
  mod_results <- recombine(tree.intra, slot1 = 4)
  mod_results$info <- NULL
  names(mod_results)[1] <- "n.tree"
  mod_results$n.tree <- as.numeric(mod_results$n.tree)
  
  #calculate mean and sd for each parameter
  #variation due to intraspecific variability
  mean_by_randomval <-
    stats::aggregate(. ~ n.intra, data = mod_results, mean)
  
  #variation due to tree choice
  mean_by_tree <- stats::aggregate(. ~ n.tree, data = mod_results, mean)
  
  statresults <- data.frame(
    min.all = apply(mod_results, 2, min),
    max.all = apply(mod_results, 2, max),
    mean.all = apply(mod_results, 2, mean),
    sd_all = apply(mod_results, 2, stats::sd),
    
    min.intra = apply(mean_by_randomval, 2, min),
    max.intra = apply(mean_by_randomval, 2, max),
    mean.intra = apply(mean_by_randomval, 2, mean),
    sd_intra = apply(mean_by_randomval, 2, stats::sd),
    
    min.tree = apply(mean_by_tree, 2, min),
    max.tree = apply(mean_by_tree, 2, max),
    mean.tree = apply(mean_by_tree, 2, mean),
    sd_tree = apply(mean_by_tree, 2, stats::sd)
  )[-(1:2),]
  
  statresults$CI_low_all    <-
    statresults$mean.all - stats::qt(0.975, df = n.intra * n.tree - 1) * statresults$sd_all / sqrt(n.intra *
                                                                                                     n.tree)
  statresults$CI_low_intra  <-
    statresults$mean.intra - stats::qt(0.975, df = n.intra - 1) * statresults$sd_intra / sqrt(n.intra)
  statresults$CI_low_tree   <-
    statresults$mean.tree - stats::qt(0.975, df = n.tree - 1) * statresults$sd_intra / sqrt(n.tree)
  
  statresults$CI_high_all    <-
    statresults$mean.all + stats::qt(0.975, df = n.intra * n.tree - 1) * statresults$sd_all / sqrt(n.intra *
                                                                                                     n.tree)
  statresults$CI_high_intra  <-
    statresults$mean.intra + stats::qt(0.975, df = n.intra - 1) * statresults$sd_intra / sqrt(n.intra)
  statresults$CI_high_tree   <-
    statresults$mean.tree + stats::qt(0.975, df = n.tree - 1) * statresults$sd_intra / sqrt(n.tree)
  
  #reoder to later match sensi_plot for the single functions
  statresults <-
    statresults[, c(
      "min.all",
      "max.all",
      "mean.all",
      "sd_all",
      "CI_low_all",
      "CI_high_all",
      "min.intra",
      "max.intra",
      "mean.intra",
      "sd_intra",
      "CI_low_intra",
      "CI_high_intra",
      "min.tree",
      "max.tree",
      "mean.tree",
      "sd_tree",
      "CI_low_tree",
      "CI_high_tree"
    )]
  
  
  #species with transformation problems
  nr <- n.tree * n.intra - nrow(mod_results)
  sp.pb <- unique(unlist(lapply(tree.intra, function(x)
    x$sp.pb)))
  
  if (length(sp.pb) > 0)
    warning (
      paste(
        "in",
        nr,
        "simulations, data transformations generated NAs, please consider using another function
                   for x.transf and check output$sp.pb",
        sep = " "
      )
    )
  
  
  res <- list(
    call = match.call(),
    formula = formula,
    x.transf = x.transf,
    data = full.data,
    sensi.estimates = mod_results,
    N.obs = tree.intra[[1]]$N.obs,
    stats = round(statresults[c(1:6), c(3, 13, 16, 7, 14, 17, 11, 15, 18)], digits =
                    3),
    all.stats = statresults,
    sp.pb = sp.pb
  )
  class(res) <- c("sensiTree_Intra", "sensiTree_IntraL")
  return(res)
}