#' Influential species detection - Phylogenetic Logistic Regression
#'
#' Performs leave-one-out deletion analysis for phylogenetic logistic regression,
#' and detects influential species.
#'
#' @param formula The model formula
#' @param data Data frame containing species traits with row names matching tips
#' in \code{phy}.
#' @param phy A phylogeny (class 'phylo') matching \code{data}.
#' @param btol Bound on searching space. For details see \code{phyloglm}
#' @param cutoff The cutoff value used to identify for influential species
#' (see Details)
#' @param track Print a report tracking function progress (default = TRUE)
#' @param ... Further arguments to be passed to \code{phyloglm}
#' @details
#' This function sequentially removes one species at a time, fits a phylogenetic
#' logistic regression model using \code{\link[phylolm]{phyloglm}}, stores the
#' results and detects influential species.
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
#' based on standardised difference in intercept and in the slope of the
#' regression. Species are ordered from most influential to less influential and
#' only include species with a standardised difference > \code{cutoff}.
#' @return \code{sensi.estimates}: A data frame with all simulation
#' estimates. Each row represents a deleted clade. Columns report the calculated
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
#' @author Gustavo Paterno & Gijsbert D.A. Werner
#' @seealso \code{\link[phylolm]{phyloglm}}, \code{\link{samp_phyglm}},
#' \code{\link{influ_phylm}}, \code{\link{sensi_plot}}
#' @references 
#' 
#' Paterno, G. B., Penone, C. Werner, G. D. A. 
#' \href{http://doi.wiley.com/10.1111/2041-210X.12990}{sensiPhy: 
#' An r-package for sensitivity analysis in phylogenetic 
#' comparative methods.} Methods in Ecology and Evolution 
#' 2018, 9(6):1461-1467.  
#' 
#' Ho, L. S. T. and Ane, C. 2014. "A linear-time algorithm for 
#' Gaussian and non-Gaussian trait evolution models". Systematic Biology 63(3):397-408.
#' @examples 
#'# Simulate Data:
#'set.seed(6987)
#'phy = rtree(100)
#'x = rTrait(n=1,phy=phy)
#'X = cbind(rep(1,100),x)
#'y = rbinTrait(n=1,phy=phy, beta=c(-1,0.5), alpha=.7 ,X=X)
#'dat = data.frame(y, x)
#'# Run sensitivity analysis:
#'influ <- influ_phyglm(y ~ x, data = dat, phy = phy) 
#'# To check summary results and most influential species:
#'summary(influ)
#'# Visual diagnostics for clade removal:
#'sensi_plot(influ)
#' @export



influ_phyglm <-         function(formula,
                                 data,
                                 phy,
                                 btol = 50,
                                 cutoff = 2,
                                 track = TRUE,
                                 ...) {
  if (!inherits(formula, "formula"))
    stop("formula must be class 'formula'")
  if (!inherits(data, "data.frame"))
    stop("data must be class 'data.frame'")
  if (!inherits(phy, "phylo"))
    stop("phy must be class 'phylo'")
  else
    
    # Check match between data and phy
    data_phy <-
      match_dataphy(formula, data, phy, ...)
  #Calculates the full model, extracts model parameters
  full.data <- data_phy$data
  phy <- data_phy$phy
  #Calculates the full model, extracts model parameters
  N               <- nrow(full.data)
  mod.0           <-
    phylolm::phyloglm(
      formula,
      data = full.data,
      phy = phy,
      method = "logistic_MPLE",
      btol = btol
    )
  intercept.0      <- mod.0$coefficients[[1]]
  estimate.0       <- mod.0$coefficients[[2]]
  pval.intercept.0 <-
    phylolm::summary.phyloglm(mod.0)$coefficients[[1, 4]]
  pval.estimate.0  <-
    phylolm::summary.phyloglm(mod.0)$coefficients[[2, 4]]
  optpar.0         <- mod.0$alpha
  if (isTRUE(mod.0$convergence != 0))
    stop("Full model failed to converge,
                                        consider changing btol. See ?phyloglm")
  else
    
    #Creates empty data frame to store model outputs
    sensi.estimates <-
    data.frame(
      "species" = numeric(),
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
  
  #Loops over all species, and removes each one individually
  counter <- 1
  errors <- NULL
  if (track == TRUE)
    pb <- utils::txtProgressBar(min = 0,
                                max = N,
                                style = 3)
  for (i in 1:N) {
    crop.data <- full.data[c(1:N)[-i],]
    crop.phy <-
      ape::drop.tip(phy, phy$tip.label[i])
    mod = try(phylolm::phyloglm(
      formula,
      data = crop.data,
      phy = crop.phy,
      method = "logistic_MPLE",
      btol = btol
    ),
    TRUE)
    if (isTRUE(class(mod) == "try-error")) {
      error <- i
      names(error) <-
        rownames(full.data$data)[i]
      errors <- c(errors, error)
      next
    }
    else {
      sp                   <- phy$tip.label[i]
      intercept            <-
        mod$coefficients[[1]]
      estimate             <-
        mod$coefficients[[2]]
      DIFintercept         <-
        intercept - intercept.0
      DIFestimate          <-
        estimate - estimate.0
      intercept.perc       <-
        round((abs(DIFintercept / intercept.0)) * 100, digits = 1)
      estimate.perc        <-
        round((abs(DIFestimate / estimate.0)) * 100, digits = 1)
      pval.intercept       <-
        phylolm::summary.phyloglm(mod)$coefficients[[1, 4]]
      pval.estimate        <-
        phylolm::summary.phyloglm(mod)$coefficients[[2, 4]]
      aic.mod              <- mod$aic
      optpar               <- mod$alpha
      
      if (track == TRUE)
        utils::setTxtProgressBar(pb, i)
      
      #Stores values for eacht simulation
      estim.simu <-
        data.frame(
          sp,
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
      sensi.estimates[counter, ]  <-
        estim.simu
      counter = counter + 1
    }
  }
  if (track == TRUE)
    on.exit(close(pb))
  #Calculates standardized DFbeta and DIFintercept
  sDIFintercept <- sensi.estimates$DIFintercept /
    stats::sd(sensi.estimates$DIFintercept)
  sDIFestimate     <- sensi.estimates$DIFestimate /
    stats::sd(sensi.estimates$DIFestimate)
  
  sensi.estimates$sDIFestimate     <- sDIFestimate
  sensi.estimates$sDIFintercept <- sDIFintercept
  
  #Creates a list with full model estimates:
  param0 <-
    list(
      coef = phylolm::summary.phyloglm(mod.0)$coefficients,
      aic = phylolm::summary.phyloglm(mod.0)$aic,
      optpar = phylolm::summary.phyloglm(mod.0)$alpha
    )
  
  #Identifies influencital species (sDF > cutoff) and orders by influence
  reorder.on.estimate         <-
    sensi.estimates[order(abs(sensi.estimates$sDIFestimate), decreasing =
                            T), c("species", "sDIFestimate")]
  influ.sp.estimate           <-
    as.character(reorder.on.estimate$species[abs(reorder.on.estimate$sDIFestimate) >
                                               cutoff])
  reorder.on.intercept     <-
    sensi.estimates[order(abs(sensi.estimates$sDIFintercept), decreasing =
                            T), c("species", "sDIFintercept")]
  influ.sp.intercept       <-
    as.character(reorder.on.intercept$species[abs(reorder.on.intercept$sDIFintercept) >
                                                cutoff])
  
  #Generates output:
  res <- list(
    cutoff = cutoff,
    formula = formula,
    full.model.estimates = param0,
    influential.species = list(
      influ.sp.estimate = influ.sp.estimate,
      influ.sp.intercept = influ.sp.intercept
    ),
    sensi.estimates = sensi.estimates,
    data = full.data,
    errors = errors
  )
  class(res) <- "sensiInflu"
  ### Warnings:
  if (length(res$errors) > 0) {
    warning("Some species deletion presented errors, please check: output$errors")
  }
  else {
    res$errors <- "No errors found."
  }
  
  return(res)
  
}


