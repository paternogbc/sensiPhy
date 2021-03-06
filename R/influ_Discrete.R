#' Influential Species Detection - Trait Evolution Discrete Characters
#' 
#' Fits models for trait evolution of discrete (binary) characters, 
#' detecting influential species. 
#'
#' @param data Data vector for a single binary trait, with names matching tips in \code{phy}.
#' @param phy A phylogeny (class 'phylo') matching \code{data}.
#' @param model The Mkn model to use (see Details). 
#' @param transform The evolutionary model to transform the tree (see Details). Default is \code{none}.
#' @param cutoff The cut-off parameter for influential species (see Details). 
#' @param bounds settings to constrain parameter estimates. See \code{\link[geiger]{fitDiscrete}}
#' @param n.cores number of cores to use. If 'NULL', number of cores is detected.
#' @param track Print a report tracking function progress (default = TRUE)
#' @param ... Further arguments to be passed to \code{\link[geiger]{fitDiscrete}}
#' @details
#' This function sequentially removes one species at a time,  
#' fits a model of discrete character evolution using \code{\link[geiger]{fitDiscrete}}, 
#' stores the results and calculates the effects on model parameters. Currently, only 
#' binary discrete traits are supported.
#' 
#' \code{influ_discrete} detects influential species based on the standardised
#' difference in q12 or q21 when removing a given species compared
#' to the full model including all species. Species with a standardised difference
#' above the value of \code{cutoff} are identified as influential. 
#' 
#' Different character model from \code{fitDiscrete} can be used, including \code{ER} (equal-rates), 
#' \code{SYM} (symmetric), \code{ARD} (all-rates-different) and \code{meristic} (stepwise fashion). 
#'
#' Different transformations to the phylogenetic tree from \code{fitDiscrete} can be used, i.e. \code{none},
#' \code{EB}, \code{lambda}, \code{kappa} and\code{delta}.
#' 
#' See \code{\link[geiger]{fitDiscrete}} for more details on character models and tree transformations. 
#'
#' @return The function \code{tree_discrete} returns a list with the following
#' components:
#' @return \code{call}: The function call
#' @return \code{cutoff}: The value selected for \code{cutoff}
#' @return \code{data}: The original full data vector
#' @return \code{optpar}: Transformation parameter used (e.g. \code{lambda}, \code{kappa} etc.)
#' @return \code{full.model.estimates}: Parameter estimates (transition rates q12 and q21), 
#' AICc and the optimised value of the phylogenetic transformation parameter (e.g. \code{lambda}) 
#' for the full model.
#' @return \code{influential_species}: List of influential species, based on standardised 
#' difference in estimates for q12 and q21. Species are ordered from most influential to 
#' less influential and only include species with a standardised difference > \code{cutoff}.
#' @return \code{sensi.estimates}: Parameter estimates (transition rates q12 and q21),,(percentual) difference 
#' in parameter estimate compared to the full model (DIFq12, sigsq.q12,sDIFq12, DIFq21, optpar.q21,sDIFq21),  
#' AICc and the optimised value of the phylogenetic transformation parameter (e.g. \code{lambda}) 
#' for each analysis with a species deleted.
#' @author Gijsbert Werner & Gustavo Paterno
#' @seealso \code{\link[geiger]{fitDiscrete}}
#' @references 
#' 
#' Paterno, G. B., Penone, C. Werner, G. D. A. 
#' \href{http://doi.wiley.com/10.1111/2041-210X.12990}{sensiPhy: 
#' An r-package for sensitivity analysis in phylogenetic 
#' comparative methods.} Methods in Ecology and Evolution 
#' 2018, 9(6):1461-1467.  
#'
#' Yang Z. 2006. Computational Molecular Evolution. Oxford University Press: Oxford. 
#' 
#' Harmon Luke J, Jason T Weir, Chad D Brock, Richard E Glor, and Wendell Challenger. 2008.
#' GEIGER: investigating evolutionary radiations. Bioinformatics 24:129-131.
#' 
#' @examples 
#' \dontrun{
#' #Load data:
#' data("primates")
#' #Create a binary trait factor 
#' adultMass_binary<-ifelse(primates$data$adultMass > 7350, "big", "small")
#' adultMass_binary<-as.factor(as.factor(adultMass_binary))
#' names(adultMass_binary)<-rownames(primates$data)
#' #Model trait evolution accounting for influential species
#' influ_binary<-influ_discrete(data = adultMass_binary,phy = primates$phy[[1]],
#' model = "SYM",transform = "none",cutoff = 2,n.cores = 2,track = TRUE)
#' #Print summary statistics
#' summary(influ_binary)
#' sensi_plot(influ_binary) #q12 and q21 are, as expected, exactly the same in symmetrical model. 
#' #Use a different evolutionary model. 
#' influ_binary2<-influ_discrete(data = adultMass_binary,phy = primates$phy[[1]],
#' model = "SYM",transform = "delta",n.cores = 2,track = TRUE)
#' summary(influ_binary2)
#' sensi_plot(influ_binary2)
#' #Or change the cutoff and transformation
#' influ_binary3<-influ_discrete(data = adultMass_binary,phy = primates$phy[[1]],
#' model = "ARD",transform = "none",cutoff = 1.2,n.cores = 2,track = TRUE)
#' summary(influ_binary3)
#' sensi_plot(influ_binary3) 
#' }
#' @export


influ_discrete <- function(data,
                           phy,
                           model,
                           transform = "none",
                           bounds = list(),
                           cutoff = 2,
                           n.cores = NULL,
                           track = TRUE,
                           ...) {
  #Error check
  if (is.null(model))
    stop("model must be specified (e.g. 'ARD' or 'SYM'")
  if (!inherits(data, "factor"))
    stop("data must supplied as a factor with species as names. Consider as.factor()")
  if (length(levels(data)) > 2)
    stop("discrete data can have maximal two levels")
  if (!inherits(phy, "phylo"))
    stop("phy must be class 'phylo'")
  if (transform == "white")
    stop("the white-noise (non-phylogenetic) model is not allowed")
  else
    
    #Matching tree
    full.data <- data
  phy <- phy
  
  #Calculates the full model, extracts model parameters
  N                   <- length(full.data)
  mod.0               <-
    geiger::fitDiscrete(
      phy = phy,
      dat = full.data,
      model = model,
      transform = transform,
      bounds = bounds,
      ncores = n.cores,
      ...
    )
  q12.0               <- mod.0$opt$q12
  q21.0               <- mod.0$opt$q21
  aicc.0              <- mod.0$opt$aicc
  if (transform == "none") {
    optpar.0 <- NA
  }
  if (transform == "EB") {
    optpar.0               <- mod.0$opt$a
  }
  if (transform == "lambda") {
    optpar.0               <- mod.0$opt$lambda
  }
  if (transform == "kappa") {
    optpar.0               <- mod.0$opt$kappa
  }
  if (transform == "delta") {
    optpar.0               <- mod.0$opt$delta
  }
  
  #Creates empty data frame to store model outputs
  sensi.estimates <- data.frame(
    "species" = numeric(),
    "q12" = numeric(),
    "DIFq12" = numeric(),
    "q12.perc" = numeric(),
    "q21" = numeric(),
    "DIFq21" = numeric(),
    "q21.perc" = numeric(),
    "aicc" = numeric(),
    "optpar" = numeric()
  )
  
  #Loops over all species, and removes each one individually
  counter <- 1
  errors <- NULL
  if (track == TRUE)
    pb <- utils::txtProgressBar(min = 0, max = N, style = 3)
  for (i in 1:N) {
    crop.data <- full.data[c(1:N)[-i]]
    crop.phy <-
      ape::drop.tip(phy, setdiff(phy$tip.label, names(crop.data)))
    
    mod = try(geiger::fitDiscrete(
      phy = crop.phy,
      dat = crop.data,
      model = model,
      transform = transform,
      bounds = bounds,
      ncores = n.cores,
      ...
    ),
    TRUE)
    if (isTRUE(class(mod) == "try-error")) {
      error <- i
      names(error) <- rownames(full.data$data)[i]
      errors <- c(errors, error)
      next
    }
    else {
      sp                   <- phy$tip.label[i]
      q12               <- mod$opt$q12
      q21               <- mod$opt$q21
      DIFq12            <- q12 - q12.0
      DIFq21            <- q21 - q21.0
      q12.perc          <-
        round((abs(DIFq12 / q12.0)) * 100,
              digits = 1)
      q21.perc          <-
        round((abs(DIFq21 / q21.0)) * 100,
              digits = 1)
      aicc              <- mod$opt$aicc
      if (transform == "none") {
        optpar <- NA
      }
      if (transform == "EB") {
        optpar               <- mod$opt$a
      }
      if (transform == "lambda") {
        optpar               <- mod$opt$lambda
      }
      if (transform == "kappa") {
        optpar               <- mod$opt$kappa
      }
      if (transform == "delta") {
        optpar               <- mod$opt$delta
      }
      
      if (track == TRUE)
        utils::setTxtProgressBar(pb, i)
      # Stores values for each simulation
      # Store reduced model parameters:
      estim.simu <- data.frame(sp,
                               q12,
                               DIFq12,
                               q12.perc,
                               q21,
                               DIFq21,
                               q21.perc,
                               aicc,
                               optpar,
                               stringsAsFactors = F)
      sensi.estimates[counter,]  <- estim.simu
      counter = counter + 1
    }
  }
  if (track == TRUE)
    on.exit(close(pb))
  #Calculates Standardized DFbeta and DIFq12
  sDIFq12 <- sensi.estimates$DIFq12 /
    stats::sd(sensi.estimates$DIFq12)
  sDIFq21     <- sensi.estimates$DIFq21 /
    stats::sd(sensi.estimates$DIFq21)
  
  sensi.estimates$sDIFq21     <- sDIFq21
  sensi.estimates$sDIFq12     <- sDIFq12
  
  #Creates a list with full model estimates:
  #full model estimates:
  param0 <- list(
    q12 = q12.0,
    q21 = q21.0,
    aicc = aicc.0,
    optpar = optpar.0
  )
  
  #Identifies influencital species (sDF > cutoff) and orders by influence
  reorder.on.q21         <- sensi.estimates[order(abs(sensi.estimates$sDIFq21), decreasing =
                                                    T), c("species", "sDIFq21")]
  influ.sp.q21           <-
    as.character(reorder.on.q21$species[abs(reorder.on.q21$sDIFq21) > cutoff])
  reorder.on.q12     <- sensi.estimates[order(abs(sensi.estimates$sDIFq12), decreasing =
                                                T), c("species", "sDIFq12")]
  influ.sp.q12       <-
    as.character(reorder.on.q12$species[abs(reorder.on.q12$sDIFq12) > cutoff])
  
  #Generates output:
  res <- list(
    call = match.call(),
    cutoff = cutoff,
    data = full.data,
    optpar = transform,
    full.model.estimates = param0,
    influential.species = list(influ.sp.q12 = influ.sp.q12,
                               influ.sp.q21 = influ.sp.q21),
    sensi.estimates = sensi.estimates,
    errors = errors
  )
  class(res) <- "sensiInflu.TraitEvol"
  ### Warnings:
  if (length(res$errors) > 0) {
    warning("Some species deletion presented errors, please check: output$errors")
  }
  else {
    res$errors <- "No errors found."
  }
  
  return(res)
  
}
