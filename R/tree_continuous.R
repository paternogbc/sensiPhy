#' Phylogenetic uncertainty - Trait Evolution Continuous Characters
#' 
#' Fits models for trait evolution of continuous characters, 
#' evaluating phylogenetic uncertainty. 
#'
#' @param data Data vector for a single continuous trait, with names matching tips in \code{phy}.
#' @param phy Phylogenies (class 'multiPhylo', see ?\code{ape}).
#' @param n.tree Number of times to repeat the analysis with n different trees picked 
#' randomly in the multiPhylo file. If NULL, \code{n.tree} = 10
#' @param model The evolutionary model (see Details). 
#' @param bounds settings to constrain parameter estimates. See \code{\link[geiger]{fitContinuous}}
#' @param n.cores number of cores to use. If 'NULL', number of cores is detected.
#' @param track Print a report tracking function progress (default = TRUE)
#' @param ... Further arguments to be passed to \code{\link[geiger]{fitContinuous}}
#' @details
#' This function fits different models of continuous character evolution using \code{\link[geiger]{fitContinuous}}
#' to n trees, randomly picked in a multiPhylo file. 
#'
#' Different evolutionary models from \code{fitContinuous} can be used, i.e. \code{BM},\code{OU},
#' \code{EB}, \code{trend}, \code{lambda}, \code{kappa}, \code{delta} and \code{drift}.
#' 
#' See \code{\link[geiger]{fitContinuous}} for more details on character models and tree transformations. 
#' 
#' Output can be visualised using \code{sensi_plot}.
#'
#' @return The function \code{tree_continuous} returns a list with the following
#' components:
#' @return \code{call}: The function call
#' @return \code{data}: The original full data vector
#' @return \code{sensi.estimates}:  (rate of evolution \code{sigsq}, 
#' root state \code{z0} and where applicable \code{optpar}),
#' AICc and the optimised value of the phylogenetic transformation parameter (e.g. \code{lambda}) 
#' for each analysis with a different phylogenetic tree.
#' @return \code{N.tree}: Number of trees \code{n.tree} analysed
#' @return \code{stats}: Main statistics for model parameters, i.e. minimum, maximum, mean, median and sd-values
#' @return \code{optpar}: Evolutionary model used (e.g. \code{lambda}, \code{kappa} etc.)
#' @author Gijsbert Werner & Caterina Penone
#' @seealso \code{\link[geiger]{fitContinuous}}
#' @references 
#' 
#' Paterno, G. B., Penone, C. Werner, G. D. A. 
#' \href{http://doi.wiley.com/10.1111/2041-210X.12990}{sensiPhy: 
#' An r-package for sensitivity analysis in phylogenetic 
#' comparative methods.} Methods in Ecology and Evolution 
#' 2018, 9(6):1461-1467
#'
#' Yang Z. 2006. Computational Molecular Evolution. Oxford University Press: Oxford. 
#' 
#' Harmon Luke J, Jason T Weir, Chad D Brock, Richard E Glor, and Wendell Challenger. 2008.
#' GEIGER: investigating evolutionary radiations. Bioinformatics 24:129-131.
#' 
#' @examples 
#' \dontshow{
#' #Load data:
#' data("primates")
#' #Model trait evolution accounting for phylogenetic uncertainty
#' adultMass<-primates$data$adultMass
#' names(adultMass)<-rownames(primates$data)
#' tree_cont<-tree_continuous(data = adultMass,phy = primates$phy,
#' model = "OU",n.tree=1,n.cores = 2,track = TRUE)
#' }
#' \dontrun{
#' #Load data:
#' data("primates")
#' #Model trait evolution accounting for phylogenetic uncertainty
#' adultMass<-primates$data$adultMass
#' names(adultMass)<-rownames(primates$data)
#' tree_cont<-tree_continuous(data = adultMass,phy = primates$phy,
#' model = "OU",n.tree=30,n.cores = 2,track = TRUE)
#' #Print summary statistics
#' summary(tree_cont)
#' sensi_plot(tree_cont)
#' sensi_plot(tree_cont,graphs="sigsq")
#' sensi_plot(tree_cont,graphs="optpar")
#' #Use a different evolutionary model 
#' tree_cont2<-tree_continuous(data = adultMass,phy = primates$phy,
#' model = "delta",n.tree=30,n.cores = 2,track = TRUE)
#' summary(tree_cont2)
#' sensi_plot(tree_cont2)
#' }
#' @export


tree_continuous <- function(data,
                            phy,
                            n.tree = 10,
                            model,
                            bounds = list(),
                            n.cores = NULL,
                            track = TRUE,
                            ...) {
  #Error check
  if (is.null(model))
    stop("model must be specified, e.g. 'OU' or 'lambda'")
  if (!inherits(data, "numeric") |
      is.null(names(data)))
    stop("data must supplied as a numeric vector with species as names")
  if (!inherits(phy, "multiPhylo"))
    stop("phy must be class 'multiPhylo'")
  if (length(phy) < n.tree)
    stop("'n.tree' must be smaller (or equal) than the number of trees in the 'multiPhylo' object")
  if (model == "white")
    stop("the white-noise (non-phylogenetic) model is not allowed")
  else
    
    #Matching tree and phylogeny using utils.R
    full.data <- data
  phy <- phy
  
  # If the class of tree is multiphylo pick n=n.tree random trees
  trees <- sample(length(phy), n.tree, replace = F)
  
  #Create the results data.frame
  sensi.estimates <-
    data.frame(
      "n.tree" = numeric(),
      "sigsq" = numeric(),
      "optpar" = numeric(),
      "z0" = numeric(),
      "aicc" = numeric()
    )
  
  #Model calculation
  counter = 1
  errors <- NULL
  c.data <- list()
  if (track == TRUE)
    pb <- utils::txtProgressBar(min = 0, max = n.tree, style = 3)
  for (j in trees) {
    #Match data order to tip order
    full.data <- full.data[phy[[j]]$tip.label]
    
    #phylolm model
    mod = try(geiger::fitContinuous(
      phy = phy[[j]],
      dat = full.data,
      model = model,
      bounds = bounds,
      ncores = n.cores,
      ...
    ),
    FALSE)
    
    if (isTRUE(class(mod) == "try-error")) {
      error <- j
      names(error) <- rownames(c.data$full.data)[j]
      errors <- c(errors, error)
      next
    }
    
    else{
      sigsq               <- mod$opt$sigsq
      z0                  <- mod$opt$z0
      aicc                <- mod$opt$aicc
      
      if (model == "BM") {
        optpar <- NA
      }
      if (model == "OU") {
        optpar        <- mod$opt$alpha
      }
      if (model == "EB") {
        optpar               <- mod$opt$a
      }
      if (model == "trend") {
        optpar               <- mod$opt$slope
      }
      if (model == "lambda") {
        optpar               <- mod$opt$lambda
      }
      if (model == "kappa") {
        optpar               <- mod$opt$kappa
      }
      if (model == "delta") {
        optpar               <- mod$opt$delta
      }
      if (model == "drift") {
        optpar              <- mod$opt$drift
      }
      
      if (track == TRUE)
        utils::setTxtProgressBar(pb, counter)
      
      #write in a table
      estim.simu <- data.frame(j, sigsq, optpar, z0, aicc,
                               stringsAsFactors = F)
      sensi.estimates[counter,]  <- estim.simu
      counter = counter + 1
      
    }
  }
  if (track == TRUE)
    on.exit(close(pb))
  #calculate mean and sd for each parameter
  
  statresults <- data.frame(
    min = apply(sensi.estimates, 2, min),
    max = apply(sensi.estimates, 2, max),
    mean = apply(sensi.estimates, 2, mean),
    median = apply(sensi.estimates, 2, median),
    sd = apply(sensi.estimates, 2, sd)
  )[-1, ]
  #Output
  res <- list(
    call = match.call(),
    data = full.data,
    sensi.estimates = sensi.estimates,
    N.tree = n.tree,
    stats = statresults[c(1:4), ],
    optpar = model
  )
  class(res) <- "sensiTree.TraitEvol"
  return(res)
}
