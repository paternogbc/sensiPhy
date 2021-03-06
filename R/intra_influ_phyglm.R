#' Interaction between intraspecific variability and influential species - Phylogenetic Logistic Regression
#'
#' Performs leave-one-out deletion analysis for phylogenetic logistic regression,
#' and detects influential species, while taking into account potential
#' interactions with intraspecific variability.
#'
#' @param formula The model formula:  \code{response~predictor}. 
#' @param data Data frame containing species traits with row names matching tips
#' in \code{phy}.
#' @param phy A phylogeny (class 'phylo') matching \code{data}.
#' @param cutoff The cutoff value used to identify for influential species
#' (see Details)
#' @param n.intra Number of datasets resimulated taking into account intraspecific variation (see: \code{"intra_phylgm"})
#' @param Vx Name of the column containing the standard deviation or the standard error of the predictor 
#' variable. When information is not available for one taxon, the value can be 0 or \code{NA}
#' @param x.transf Transformation for the predictor variable (e.g. \code{"log"} or \code{"sqrt"}). Please use this 
#' argument instead of transforming data in the formula directly (see also details below).
#' @param distrib A character string indicating which distribution to use to generate a random value for the response 
#' and/or predictor variables. Default is normal distribution: "normal" (function \code{\link{rnorm}}).
#' Uniform distribution: "uniform" (\code{\link{runif}})
#' Warning: we recommend to use normal distribution with Vx = standard deviation of the mean.
#' @param btol Bound on searching space. For details see \code{phyloglm}
#' @param track Print a report tracking function progress (default = TRUE)
#' @param ... Further arguments to be passed to \code{phylolm}
#' @details
#' This function fits a phylogenetic linear regression model using \code{\link[phylolm]{phylolm}}, and detects
#' influential species by sequentially deleting one at a time. The regression is repeated \code{n.intra} times for 
#' simulated values of the dataset, taking into account intraspecific variation. At each iteration, the function 
#' generates a random value for each row in the dataset using the standard deviation or errors supplied, and 
#' detect the influential species within that iteration. 
#'
#' \code{influ_phylm} detects influential species based on the standardised
#' difference in intercept and/or slope when removing a given species compared
#' to the full model including all species. Species with a standardised difference
#' above the value of \code{cutoff} are identified as influential. The default
#' value for the cutoff is 2 standardised differences change.
#'
#' Currently, this function can only implement simple models (i.e. \eqn{trait~
#' predictor}). In the future we will implement more complex models.
#'
#' Output can be visualised using \code{sensi_plot}.
#' 
#' @section Warning:  
#' When Vx exceeds X negative (or null) values can be generated, this might cause problems
#' for data transformation (e.g. log-transformation). In these cases, the function will skip the simulation. This problem can
#' be solved by increasing \code{n.intra}, changing the transformation type and/or checking the target species in output$sp.pb.
#' 
#' Setting \code{n.intra} at high values can take a long time to execute, since the total number of iterations equals \code{n.intra * nrow(data)}.
#' 
#' @return The function \code{intra_influ_phylm} returns a list with the following
#' components:
#' @return \code{cutoff}: The value selected for \code{cutoff}
#' @return \code{formula}: The formula
#' @return \code{full.model.estimates}: Coefficients, aic and the optimised
#' value of the phylogenetic parameter (e.g. \code{lambda}) for the full model
#' without deleted species.
#' @return \code{influential_species}: List of influential species, both
#' based on standardised difference in intercept and in the slope of the
#' regression. Species are ordered from most influential to less influential and
#' only include species with a standardised difference > \code{cutoff}.
#' @return \code{sensi.estimates}: A data frame with all simulation
#' estimates. Each row represents a deleted clade for an iteration of resimulated
#' data. Columns report the calculated regression intercept (\code{intercept}), 
#' difference between simulation intercept and full model intercept (\code{DIFintercept}), 
#' the standardised difference (\code{sDIFintercept}), the percentage of change in intercept compared
#' to the full model (\code{intercept.perc}) and intercept p-value
#' (\code{pval.intercept}). All these parameters are also reported for the regression
#' slope (\code{DIFestimate} etc.). Additionally, model aic value (\code{AIC}) and
#' the optimised value (\code{optpar}) of the phylogenetic parameter
#' (e.g. \code{kappa} or \code{lambda}, depending on the phylogenetic model used) are
#' reported.
#' @return \code{data}: Original full dataset.
#' @return \code{errors}: Species where deletion resulted in errors. 
#' @author Gustavo Paterno, Caterina Penone & Gijsbert D.A. Werner
#' @seealso \code{\link[phylolm]{phylolm}}, \code{\link{intra_phyglm}},
#' \code{\link{influ_phyglm}},\code{\link{intra_influ_phylm}},\code{\link{sensi_plot}}.
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
#' \dontrun{
#'#Generate data
#'set.seed(6987)
#'phy = rtree(100)
#'x = rTrait(n=1,phy=phy,parameters=list(ancestral.state=2,optimal.value=2,sigma2=1,alpha=1))
#'X = cbind(rep(1,100),x)
#'y = rbinTrait(n=1,phy=phy, beta=c(-1,0.5), alpha=.7 ,X=X)
#'z = rnorm(n = length(x),mean = mean(x),sd = 0.1*mean(x))
#'dat = data.frame(y, x, z)
#'# Run sensitivity analysis:
#'intra_influ <- intra_influ_phyglm(formula = y ~ x, data = dat, phy = phy, 
#'                        Vx = "z", n.intra = 5,track = TRUE,distrib="normal",x.transf=NULL) 
#'# To check summary results and most influential species:
#'summary(intra_influ)
#'# Visual diagnostics for clade removal:
#'sensi_plot(intra_influ)
#'}
#'\dontshow{
#Generate data
#'phy = rtree(100)
#'x = rTrait(n=1,phy=phy,parameters=list(ancestral.state=2,optimal.value=2,sigma2=1,alpha=1))
#'X = cbind(rep(1,100),x)
#'y = rbinTrait(n=1,phy=phy, beta=c(-1,0.5), alpha=.7 ,X=X)
#'z = rnorm(n = length(x),mean = mean(x),sd = 0.1*mean(x))
#'dat = data.frame(y, x, z)
#'# Run sensitivity analysis:
#'intra_influ <- intra_influ_phyglm(formula = y ~ x, data = dat[1:20,], phy = phy, 
#'                                  Vx = "z", n.intra = 2,track = TRUE,
#'                                  distrib="normal",x.transf=NULL) 
#'# To check summary results and most influential species:
#'summary(intra_influ)
#'# Visual diagnostics for clade removal:
#'sensi_plot(intra_influ)
#'}
#' @export


intra_influ_phyglm <- function(formula,
                               data,
                               phy,
                               Vx = NULL,
                               n.intra = 30,
                               x.transf = NULL,
                               distrib = "normal",
                               cutoff = 2,
                               btol = 50,
                               track = TRUE,
                               ...) {
  #Error check
  if (is.null(Vx))
    stop("Vx must be defined")
  if (!inherits(formula, "formula"))
    stop("formula must be class 'formula'")
  if (!inherits(data, "data.frame"))
    stop("data must be class 'data.frame'")
  if (!inherits(phy, "phylo"))
    stop("phy must be class 'phylo'")
  if (formula[[2]] != all.vars(formula)[1] ||
      formula[[3]] != all.vars(formula)[2])
    stop("Please use argument x.transf for data transformation")
  if (distrib == "normal")
    warning ("distrib=normal: make sure that standard deviation is provided for Vx")
  
  #Matching tree and phylogeny using utils.R
  datphy <- match_dataphy(formula, data, phy)
  full.data <- datphy[[1]]
  phy <- datphy[[2]]
  
  resp1 <- all.vars(formula)[1]
  if (length(all.vars(formula)) > 2) {
    resp2 <- all.vars(formula)[2]
  }
  pred <- all.vars(formula)[length(all.vars(formula))]
  
  if (!is.null(Vx) && sum(is.na(full.data[, Vx])) != 0) {
    full.data[is.na(full.data[, Vx]), Vx] <- 0
  }
  
  #Function to pick a random value in the interval
  if (distrib == "normal")
    funr <- function(a, b) {
      stats::rnorm(1, a, b)
    }
  else
    funr <- function(a, b) {
      stats::runif(1, a - b, a + b)
    }
  
  #Start intra loop here
  intra.influ <- list ()
  species.NA <- list()
  errors <- NULL
  if (track == TRUE)
    pb <- utils::txtProgressBar(min = 0, max = n.intra, style = 3)
  counter = 1
  
  for (i in 1:n.intra) {
    ##Set predictor variable
    #Vx is not provided or is not numeric, do not pick random value
    if (!inherits(full.data[, pred], c("numeric", "integer")) ||
        is.null(Vx)) {
      full.data$predV <- full.data[, pred]
    }
    
    #choose a random value in [mean-se,mean+se] if Vx is provided
    if (!is.null(Vx) && is.null(dim(Vx)))
    {
      full.data$predV <-
        apply(full.data[, c(pred, Vx)], 1, function(x)
          funr(x[1], x[2]))
    }
    
    full.data$resp1 <-
      full.data[, resp1] #try to improve this in future
    if (length(all.vars(formula)) > 2) {
      full.data$resp2 <- full.data[, resp2]
    }
    
    #transform Vx if x.transf is provided
    if (!is.null(x.transf))
    {
      suppressWarnings (full.data$predV <- x.transf(full.data$predV))
    }
    
    #skip iteration if there are NA's in the dataset
    species.NA[[i]] <-
      rownames(full.data[with(full.data, is.na(predV)), ])
    if (sum(is.na(full.data[, "predV"])) > 0)
      next
    
    #model
    #Run the model
    if (length(all.vars(formula)) > 2) {
      intra.influ[[i]] <-
        influ_phyglm(
          cbind(resp1, resp2) ~ predV,
          data = full.data,
          phy = phy,
          method = "logistic_MPLE",
          cutoff = cutoff,
          btol = btol,
          track = FALSE,
          verbose = FALSE,
          ...
        )
    } else
      intra.influ[[i]] <-
      influ_phyglm(
        resp1 ~ predV,
        data = full.data,
        phy = phy,
        method = "logistic_MPLE",
        cutoff = cutoff,
        btol = btol,
        track = FALSE,
        verbose = FALSE,
        ...
      )
    
    if (track == TRUE)
      utils::setTxtProgressBar(pb, counter)
    counter = counter + 1
  }
  
  if (track == TRUE)
    close(pb)
  names(intra.influ) <- 1:n.intra
  
  # Merge lists into data.frames between iterations:
  full.estimates  <-
    suppressWarnings(recombine(intra.influ, slot1 = 3, slot2 = 1))
  
  #influ species slope
  influ.sp.estimate <-
    (lapply(intra.influ, function(x)
      x$influential.species$influ.sp.estimate))
  influ.sp.estimate <- as.data.frame(as.matrix(influ.sp.estimate))
  names(influ.sp.estimate) <- "influ.sp.estimate"
  influ.sp.estimate$tree <- row.names(influ.sp.estimate)
  
  #influ species intercept
  influ.sp.intercept <-
    (lapply(intra.influ, function(x)
      x$influential.species$influ.sp.intercept))
  influ.sp.intercept <- as.data.frame(as.matrix(influ.sp.intercept))
  names(influ.sp.intercept) <- "influ.sp.intercept"
  influ.sp.intercept$tree <- row.names(influ.sp.intercept)
  
  #influ.estimates
  influ.estimates <- recombine(intra.influ, slot1 = 5)
  influ.estimates$info <- NULL
  
  #Generates output:
  res <- list(
    call = match.call(),
    cutoff = cutoff,
    formula = formula,
    full.model.estimates = full.estimates,
    influential.species = list(
      influ.sp.estimate = influ.sp.estimate,
      influ.sp.intercept = influ.sp.intercept
    ),
    sensi.estimates = influ.estimates,
    data = full.data
  )
  
  class(res) <- c("sensiIntra_Influ", "sensiIntra_InfluL")
  ### Warnings:
  if (length(res$errors) > 0) {
    warning("Some species deletion presented errors, please check: output$errors")
  }
  else {
    res$errors <- "No errors found."
  }
  
  return(res)
}

    
