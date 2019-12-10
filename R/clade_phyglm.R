#'
#' Estimate the impact on model estimates of phylogenetic logistic regression after 
#' removing clades from the analysis. 
#'
#' @param formula The model formula
#' @param data Data frame containing species traits with row names matching tips
#' in \code{phy}.
#' @param phy A phylogeny (class 'phylo') matching \code{data}.
#' @param btol Bound on searching space. For details see \code{phyloglm}
#' @param track Print a report tracking function progress (default = TRUE)
#' @param clade.col The column in the provided data frame which specifies the
#' clades (a character vector with clade names).
#' @param n.species Minimum number of species in a clade for the clade to be
#' included in the leave-one-out deletion analysis. Default is \code{5}.
#' @param n.sim Number of simulations for the randomization test.
#' @param ... Further arguments to be passed to \code{phyloglm}
#' 
#' @details
#' This function sequentially removes one clade at a time, fits a phylogenetic
#' logistic regression model using \code{\link[phylolm]{phyloglm}} and stores the
#' results. The impact of of a specific clade on model estimates is calculated by a
#' comparison between the full model (with all species) and the model without 
#' the species belonging to a clade.
#' 
#' To account for the influence of the number of species on each 
#' clade (clade sample size), this function also estimates a null distribution
#' expected for the number of species in a given clade. This is done by fitting
#'  models without the same number of species as in the given clade. 
#'  The number of simulations to be performed is set by 'n.sim'. To test if the 
#'  clade influence differs from the null expectation for a clade of that size, 
#'  a randomization test can be performed using 'summary(x)'. 
#'
#' Currently, only logistic regression using the "logistic_MPLE"-method from
#' \code{phyloglm} is implemented.
#'
#' \code{clade_phyglm} detects influential clades based on
#' difference in intercept and/or estimate when removing a given clade compared
#' to the full model including all species.
#' 
#' Additionally, to account for the influence of the number of species on each 
#' clade (clade sample size), this function also estimates a null distribution 
#' expected for the number of species in a given clade. This is done by fitting
#' models without the same number of species in the given clade. 
#'  The number of simulations to be performed is set by 'n.sim'. To test if the 
#'  clade influence differs from the null expectation for a clade of that size, 
#'  a randomization test can be performed using 'summary(x)'. 
#' 
#' Currently, this function can only implement simple logistic models (i.e. \eqn{trait~
#' predictor}). In the future we will implement more complex models.
#'
#' Output can be visualised using \code{sensi_plot}.
#'
#' @return The function \code{clade_phyglm} returns a list with the following
#' components:
#' @return \code{formula}: The formula
#' @return \code{full.model.estimates}: Coefficients, aic and the optimised
#' value of the phylogenetic parameter (e.g. \code{alpha}) for the full model
#' without deleted species.
#' @return \code{sensi.estimates}: A data frame with all simulation
#' estimates. Each row represents a deleted clade. Columns report the calculated
#' regression intercept (\code{intercept}), difference between simulation
#' intercept and full model intercept (\code{DIFintercept}), the percentage of change
#' in intercept compared to the full model (\code{intercept.perc}) and intercept
#' p-value (\code{pval.intercept}). All these parameters are also reported for the regression
#' slope (\code{DIFestimate} etc.). Additionally, model aic value (\code{AIC}) and
#' the optimised value (\code{optpar}) of the phylogenetic parameter are
#' reported.
#' @return \code{null.dist}: A data frame with estimates for the null distributions
#' for all clades analysed.
#' @return \code{data}: Original full dataset.
#' @return \code{errors}: Clades where deletion resulted in errors.
#' @return \code{clade.col}: Which column was used to specify the clades?
#' @author Gustavo Paterno & Gijsbert Werner
#' @seealso \code{\link[phylolm]{phyloglm}}, \code{\link[sensiPhy]{clade_phylm}},
#'  \code{\link{influ_phyglm}}, \code{\link{sensi_plot}}
#' @references 
#' 
#' Paterno, G. B., Penone, C. Werner, G. D. A. 
#' \href{http://doi.wiley.com/10.1111/2041-210X.12990}{sensiPhy: 
#' An r-package for sensitivity analysis in phylogenetic 
#' comparative methods.} Methods in Ecology and Evolution 
#' 2018, 9(6):1461-1467
#' 
#' Ho, L. S. T. and Ane, C. 2014. "A linear-time algorithm for 
#' Gaussian and non-Gaussian trait evolution models". Systematic Biology 63(3):397-408.
#' 
#' @examples 
#' \dontrun{
#'# Simulate Data:
#'set.seed(6987)
#'phy = rtree(150)
#'x = rTrait(n=1,phy=phy)
#'X = cbind(rep(1,150),x)
#'y = rbinTrait(n=1,phy=phy, beta=c(-1,0.5), alpha=.7 ,X=X)
#'cla <- rep(c("A","B","C","D","E"), each = 30)
#'dat = data.frame(y, x, cla)
#'# Run sensitivity analysis:
#'clade <- clade_phyglm(y ~ x, phy = phy, data = dat, n.sim = 30, clade.col = "cla")
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


clade_phyglm <- function(formula,
                         data,
                         phy,
                         btol = 50,
                         track = TRUE,
                         clade.col,
                         n.species = 5,
                         n.sim = 100,
                         ...) {
  # To check summary results and most influential clades:
  
  
  # Error checking:
  if (!inherits(data, "data.frame"))
    stop("data must be class 'data.frame'")
  if (missing(clade.col))
    stop("clade.col not defined. Please, define the",
         " column with clade names.")
  if (!inherits(phy, "phylo"))
    stop("phy must be class 'phylo'")
  
  #Calculates the full model, extracts model parameters
  data_phy <- match_dataphy(formula, data, phy, ...)
  phy <- data_phy$phy
  full.data <- data_phy$data
  if (is.na(match(clade.col, names(full.data)))) {
    stop("Names column '", clade.col, "' not found in data frame'")
  }
  
  # Identify CLADES to use and their sample size
  wc <- table(full.data[, clade.col]) > n.species
  uc <- table(full.data[, clade.col])[wc]
  
  #k <- names(which(table(full.data[,clade.col]) > n.species ))
  if (length(uc) == 0)
    stop(
      paste(
        "There is no clade with more than ",
        n.species,
        " species. Change 'n.species' to fix this
                                  problem",
        sep = ""
      )
    )
  
  # FULL MODEL PARAMETERS:
  N                <- nrow(full.data)
  mod.0            <- phylolm::phyloglm(
    formula,
    data = full.data,
    phy = phy,
    method = "logistic_MPLE",
    btol = btol
  )
  intercept.0      <- mod.0$coefficients[[1]]
  estimate.0       <- mod.0$coefficients[[2]]
  
  if (isTRUE(mod.0$convergence != 0))
    stop("Full model failed to converge,
                                        consider changing btol. See ?phyloglm")
  
  #Create dataframe to store estmates for each clade
  sensi.estimates <-
    data.frame(
      "clade" = I(as.character()),
      "N.species" = numeric(),
      "intercept" = numeric(),
      "DIFintercept" = numeric(),
      "intercept.perc" = numeric(),
      "pval.intercept" = numeric(),
      "estimate" = numeric(),
      "DIFestimate" = numeric(),
      "estimate.perc" = numeric(),
      "pval.estimate" = numeric(),
      "AIC" = numeric(),
      "optpar" = numeric()
    )
  
  # Create dataframe store simulations (null distribution)
  null.dist <- data.frame(
    "clade" = rep(names(uc), each = n.sim),
    "intercept" = numeric(length(uc) * n.sim),
    "estimate" = numeric(length(uc) * n.sim),
    "DIFintercept" = numeric(length(uc) * n.sim),
    "DIFestimate" = numeric(length(uc) * n.sim)
  )
  
  
  ### START LOOP between CLADES:
  # counters:
  aa <- 1
  bb <- 1
  errors <- NULL
  
  if (track == TRUE)
    pb <- utils::txtProgressBar(min = 0,
                                max = length(uc) * n.sim,
                                style = 3)
  for (A in names(uc)) {
    ### Number of species in clade A
    cN  <- as.numeric(uc[names(uc) == A])
    
    ### Fit reduced model (without clade)
    crop.data <- full.data[!full.data[, clade.col] %in% A, ]
    crop.sp <-   which(full.data[, clade.col] %in% A)
    crop.phy <-  ape::drop.tip(phy, phy$tip.label[crop.sp])
    mod = try(phylolm::phyloglm(
      formula,
      data = crop.data,
      phy = crop.phy,
      method = "logistic_MPLE",
      btol = btol
    ),
    TRUE)
    intercept           <- mod$coefficients[[1]]
    estimate            <- mod$coefficients[[2]]
    DIFintercept        <- intercept - intercept.0
    DIFestimate         <- estimate - estimate.0
    intercept.perc      <-
      round((abs(DIFintercept / intercept.0)) * 100,
            digits = 1)
    estimate.perc       <-
      round((abs(DIFestimate / estimate.0)) * 100,
            digits = 1)
    pval.intercept      <-
      phylolm::summary.phylolm(mod)$coefficients[[1, 4]]
    pval.estimate       <-
      phylolm::summary.phylolm(mod)$coefficients[[2, 4]]
    aic.mod             <- mod$aic
    optpar              <- mod$alpha
    
    
    # Store reduced model parameters:
    estim.simu <-
      data.frame(
        A,
        cN,
        intercept,
        DIFintercept,
        intercept.perc,
        pval.intercept,
        estimate,
        DIFestimate,
        estimate.perc,
        pval.estimate,
        aic.mod,
        optpar,
        stringsAsFactors = F
      )
    sensi.estimates[aa,]  <- estim.simu
    
    ### START LOOP FOR NULL DIST:
    # number of species in clade A:
    for (i in 1:n.sim) {
      exclude <- sample(1:N, cN)
      crop.data <- full.data[-exclude, ]
      crop.phy <-  ape::drop.tip(phy, phy$tip.label[exclude])
      mod = try(phylolm::phyloglm(
        formula,
        data = crop.data,
        phy = crop.phy,
        method = "logistic_MPLE",
        btol = btol
      ),
      TRUE)
      intercept         <- mod$coefficients[[1]]
      estimate          <- mod$coefficients[[2]]
      DIFintercept      <- intercept - intercept.0
      DIFestimate       <- estimate - estimate.0
      
      null.dist[bb,]   <- data.frame(clade = as.character(A),
                                     intercept,
                                     estimate,
                                     DIFintercept,
                                     DIFestimate)
      
      if (track == TRUE)
        utils::setTxtProgressBar(pb, bb)
      bb <- bb + 1
    }
    aa <- aa + 1
  }
  if (track == TRUE)
    on.exit(close(pb))
  
  #OUTPUT
  #full model estimates:
  param0 <- list(
    coef = phylolm::summary.phylolm(mod.0)$coefficients,
    aic = phylolm::summary.phylolm(mod.0)$aic,
    optpar = mod.0$optpar
  )
  
  #Generates output:
  res <- list(
    call = match.call(),
    formula = formula,
    full.model.estimates = param0,
    sensi.estimates = sensi.estimates,
    null.dist = null.dist,
    data = full.data,
    errors = errors,
    clade.col = clade.col
  )
  class(res) <- c("sensiClade", "sensiCladeL")
  ### Warnings:
  if (length(res$errors) > 0) {
    warning("Some clades deletion presented errors, please check: output$errors")
  }
  else {
    res$errors <- "No errors found."
  }
  return(res)
}

