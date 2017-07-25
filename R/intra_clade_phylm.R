#' Interaction of intraspecific variability & influential clade detection - Phylogenetic Linear Regression
#'
#' Estimate the impact on model estimates of phylogenetic linear regression after 
#' removing clades from the analysis, while taking into account potential
#' interactions with intraspecific variability.
#'
#' @param formula The model formula
#' @param data Data frame containing species traits with row names matching tips
#' in \code{phy}.
#' @param phy A phylogeny (class 'phylo') matching \code{data}.
#' @param clade.col The name of a column in the provided data frame with clades 
#' specification (a character vector with clade names).
#' @param n.species Minimum number of species in the clade in order to include
#' this clade in the leave-one-out deletion analyis. Default is \code{5}.
#' @param times.clade Number of simulations for the randomization test.
#' @param times.intra Number of datasets resimulated taking into account intraspecific variation (see: \code{"intra_phylm"})
#' @param Vy Name of the column containing the standard deviation or the standard error of the response 
#' variable. When information is not available for one taxon, the value can be 0 or \code{NA}.
#' @param Vx Name of the column containing the standard deviation or the standard error of the predictor 
#' variable. When information is not available for one taxon, the value can be 0 or \code{NA}
#' @param y.transf Transformation for the response variable (e.g. \code{"log"} or \code{"sqrt"}). Please use this 
#' argument instead of transforming data in the formula directly (see also details below).
#' @param x.transf Transformation for the predictor variable (e.g. \code{"log"} or \code{"sqrt"}). Please use this 
#' argument instead of transforming data in the formula directly (see also details below).
#' @param distrib A character string indicating which distribution to use to generate a random value for the response 
#' and/or predictor variables. Default is normal distribution: "normal" (function \code{\link{rnorm}}).
#' Uniform distribution: "uniform" (\code{\link{runif}})
#' Warning: we recommend to use normal distribution with Vx or Vy = standard deviation of the mean.
#' @param model The phylogenetic model to use (see Details). Default is \code{lambda}.
#' @param track Print a report tracking function progress (default = TRUE)
#' @param ... Further arguments to be passed to \code{phylolm}
#' @details
#' This function sequentially removes one clade at a time, fits a phylogenetic
#' linear regression model using \code{\link[phylolm]{phylolm}} and stores the
#' results. The impact of of a specific clade on model estimates is calculated by the
#' comparison between the full model (with all species) and the model without 
#' the species belonging to a clade. This operation is repeated \code{times.intra} times for
#' simulated values of the dataset, taking into account intraspecific variation. At each iteration, the function 
#' generates a random value for each row in the dataset using the standard deviation or errors supplied, and 
#' detect the influential species within that iteration. 
#' 
#' Additionally, to account for the influence of the number of species on each 
#' clade (clade sample size), this function also estimate a null distribution of slopes
#' expected for the number of species in a given clade. This is done by fitting
#' models without the same number of species in the given clade. 
#' The number of simulations to be performed is set by 'times'. To test if the 
#' clade influence differs from the null expectation, a randomization test can
#' be performed using 'summary(x)'. 
#'
#' All phylogenetic models from \code{phylolm} can be used, i.e. \code{BM},
#' \code{OUfixedRoot}, \code{OUrandomRoot}, \code{lambda}, \code{kappa},
#' \code{delta}, \code{EB} and \code{trend}. See ?\code{phylolm} for details.
#'
#' \code{clade_phylm} detects influential clades based on
#' difference in intercept and/or slope when removing a given clade compared
#' to the full model including all species.
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
#' @return \code{errors}: Clades and/or iterations where deletion resulted in errors.
#' @author Gustavo Paterno, Caterina Penone
#' @seealso \code{\link[phylolm]{phylolm}}, \code{\link[sensiPhy]{samp_phyglm}},
#' \code{\link{influ_phylm}}, \code{\link{sensi_plot}}
#' \code{\link{sensi_plot}}, \code{\link{intra_phylm}}
#' @references Ho, L. S. T. and Ane, C. 2014. "A linear-time algorithm for 
#' Gaussian and non-Gaussian trait evolution models". Systematic Biology 63(3):397-408.
#' @examples
#' #load data
#' data(alien)
#' intra_clade <- intra_clade_phylm(gestaLen ~ adultMass, phy = alien$phy[[1]], data = alien$data, 
#' clade.col = "family", times.clade = 30, times.intra = 3, y.transf = log, Vy="SD_gesta")
#' summary(intra_clade)
#' sensi_plot(intra_clade)
#' @export


intra_clade_phylm <- function(formula, data, phy, clade.col, n.species = 5,
                             times.clade = 100, times.intra = 2,
                             Vy = NULL, Vx = NULL, distrib = "normal",
                             y.transf = NULL, x.transf = NULL,
                             model = "lambda", track = TRUE,...) {
  
  # Error checking:
  if(is.null(Vx) & is.null(Vy)) stop("Vx or Vy must be defined")
  if(!is.data.frame(data)) stop("data must be class 'data.frame'")
  if(missing(clade.col)) stop("clade.col not defined. Please, define the column with clade names.")
  if(class(phy)!="phylo") stop("phy must be class 'phylo'")
  if(class(formula)!="formula") stop("formula must be class 'formula'")
  if(formula[[2]]!=all.vars(formula)[1] || formula[[3]]!=all.vars(formula)[2])
    stop("Please use arguments y.transf or x.transf for data transformation")
  if(distrib == "normal") warning ("distrib=normal: make sure that standard deviation 
                                   is provided for Vx and/or Vy")

  #Match data and phy
  data_phy <- match_dataphy(formula, data, phy, ...)
  phy <- data_phy$phy
  full.data <- data_phy$data
  if (is.na(match(clade.col, names(full.data)))) {
    stop("Names column '", clade.col, "' not found in data frame'")
  }
  
  #Prepare data
  resp <- all.vars(formula)[1]
  pred <- all.vars(formula)[2]
  
  if(!is.null(Vy) && sum(is.na(full.data[, Vy])) != 0) {
    full.data[is.na(full.data[, Vy]), Vy] <- 0}
  
  if(!is.null(Vx) && sum(is.na(full.data[, Vx])) != 0) {
    full.data[is.na(full.data[, Vx]), Vx] <- 0}
  
  #Function to pick a random value in the interval
  if (distrib == "normal") funr <- function(a,b) {stats::rnorm(1,a,b)}
  else  funr <- function(a,b) {stats::runif(1, a - b, a + b)}
  
  
  #Identify CLADES to use and their sample size 
  all.clades <- levels(full.data[ ,clade.col])
  wc <- table(full.data[ ,clade.col]) > n.species
  uc <- table(full.data[ , clade.col])[wc]
  
  if (length(uc) == 0) stop(paste("There is no clade with more than ",
                                  n.species," species. Change 'n.species' to fix this
                                  problem",sep=""))
  
  #List to store information
  intra.clade <- list ()
  
  #Start tree loop here
  errors <- NULL
  counter = 1
  
  for (j in 1:times.intra){
    ##Set response and predictor variables
    #Vy is not provided or is not numeric, do not pick random value
    if(!inherits(full.data[,resp], c("numeric","integer")) || is.null(Vy)) 
    {full.data$respV <- stats::model.frame(formula, data = full.data)[,1]}
    
    #choose a random value in [mean-se,mean+se] if Vy is provided
    if (!is.null(Vy))
    {full.data$respV <- apply(full.data[,c(resp,Vy)],1,function(x)funr(x[1],x[2]))}
    
    #Vx is not provided or is not numeric, do not pick random value
    if (!inherits(full.data[,pred], c("numeric","integer")) || is.null(Vx))
    {full.data$predV <- stats::model.frame(formula, data = full.data)[,2]}
    
    #choose a random value in [mean-se,mean+se] if Vx is provided
    if(!is.null(Vx))
    {full.data$predV <- apply(full.data[,c(pred,Vx)],1,function(x)funr(x[1],x[2]))}
    
    #transform Vy and/or Vx if x.transf and/or y.transf are provided
    if(!is.null(y.transf)) 
    {suppressWarnings (full.data$respV <- y.transf(full.data$respV))}
    
    if(!is.null(x.transf)) 
    {suppressWarnings (full.data$predV <- x.transf(full.data$predV))}
    
    intra.clade[[j]] <- clade_phylm(formula, data=full.data, phy=phy, model=model,
                                    clade.col=clade.col, n.species=n.species, times.clade=times.clade, 
                                    track = FALSE, verbose = FALSE,...)
    
    counter = counter + 1
  }
  
  names(intra.clade) <- 1:times.intra
  
  # Merge lists into data.frames between iterations:
  full.estimates  <- suppressWarnings(recombine(intra.clade, slot1 = 4, slot2 = 1))
  clade.estimates <- recombine(intra.clade, slot1 = 5)
  null.dist       <- recombine(intra.clade, slot1 = 6)
  
  #Generates output:
  res <- list(call = match.call(),
              model = model,
              formula = formula,
              full.model.estimates = full.estimates,
              clade.model.estimates = clade.estimates,
              null.dist = null.dist, 
              data = full.data,
              errors = errors,
              clade.col = clade.col)
  
  class(res) <- "sensiIntra_Clade"
  
  ### Warnings:
  if (length(res$errors) >0){
    warning("Some clades deletion presented errors, please check: output$errors")}
  else {
    res$errors <- "No errors found."
  }
  return(res)
}

