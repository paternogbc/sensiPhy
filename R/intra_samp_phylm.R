#' Interaction between intraspecific variability and species sampling - Phylogenetic Linear Regression
#'
#' Performs analyses of sensitivity to species sampling by randomly removing
#' species and detecting the effects on parameter estimates in a phylogenetic
#' linear regression, while taking into account potential
#' interactions with intraspecific variability.
#'
#' @param formula The model formula
#' @param data Data frame containing species traits with row names matching tips
#' in \code{phy}.
#' @param phy A phylogeny (class 'phylo') matching \code{data}.
#' @param n.sim The number of times species are randomly deleted for each
#' \code{break}.
#' @param n.intra Number of datasets resimulated taking into account intraspecific variation (see: \code{"intra_phylm"}) 
#' @param breaks A vector containing the percentages of species to remove.
#' @param model The phylogenetic model to use (see Details). Default is \code{lambda}.
#' #' @param Vy Name of the column containing the standard deviation or the standard error of the response 
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
#' @param track Print a report tracking function progress (default = TRUE)
#' @param ... Further arguments to be passed to \code{phylolm}
#' @details
#'
#' This function randomly removes a given percentage of species (controlled by
#' \code{breaks}) from the full phylogenetic linear regression, fits a phylogenetic
#' linear regression model without these species using \code{\link[phylolm]{phylolm}},
#' repeats this many times (controlled by \code{n.sim}), stores the results and
#' calculates the effects on model parameters. 
#' This operation is repeated \code{n.intra} times for simulated values of the dataset, 
#' taking into account intraspecific variation. At each iteration, the function generates a 
#' random value for each row in the dataset using the standard deviation or errors supplied, and 
#' evaluates the effects of sampling within that iteration.
#'
#' All phylogenetic models from \code{phylolm} can be used, i.e. \code{BM},
#' \code{OUfixedRoot}, \code{OUrandomRoot}, \code{lambda}, \code{kappa},
#' \code{delta}, \code{EB} and \code{trend}. See ?\code{phylolm} for details.
#'
#' Currently, this function can only implement simple linear models (i.e. \eqn{trait~
#' predictor}). In the future we will implement more complex models.
#'
#' Output can be visualised using \code{sensi_plot}.
#' @return The function \code{samp_phylm} returns a list with the following
#' components:
#' @return \code{formula}: The formula
#' @return \code{full.model.estimates}: Coefficients, aic and the optimised
#' value of the phylogenetic parameter (e.g. \code{lambda} or \code{kappa}) for
#' the full model without deleted species.
#' @return \code{sensi.estimates}: A data frame with all simulation
#' estimates. Each row represents a model rerun with a given number of species
#' \code{n.remov} removed, representing \code{n.percent} of the full dataset.
#' Columns report the calculated regression intercept (\code{intercept}),
#' difference between simulation intercept and full model intercept (\code{DIFintercept}),
#' the percentage of change in intercept compared to the full model (\code{intercept.perc})
#' and intercept p-value (\code{pval.intercept}). All these parameters are also reported
#' for the regression slope (\code{DIFestimate} etc.). Additionally, model aic value
#' (\code{AIC}) and the optimised value (\code{optpar}) of the phylogenetic
#' parameter (e.g. \code{kappa} or \code{lambda}, depending on the phylogenetic model
#' used) are reported. Lastly we reported the standardised difference in intercept 
#' (\code{sDIFintercept}) and slope (\code{sDIFestimate}). 
#' @return \code{sign.analysis} For each break (i.e. each percentage of species
#' removed) this reports the percentage of statistically signficant (at p<0.05)
#' intercepts (\code{perc.sign.intercept}) over all repititions as well as the
#' percentage of statisticaly significant (at p<0.05) slopes (\code{perc.sign.estimate}).
#' @return \code{data}: Original full dataset.
#' @note Please be aware that dropping species may reduce power to detect 
#' significant slopes/intercepts and may partially be responsible for a potential 
#' effect of species removal on p-values. Please also consult standardised differences
#' in the (summary) output.
#' @author Gustavo Paterno, Gijsbert D.A. Werner & Caterina Penone
#' @seealso \code{\link[phylolm]{phylolm}},\code{\link{samp_phylm}}, 
#' \code{\link{intra_phylm}}, \code{\link{_intra_samp_phyglm}},
#' \code{\link{sensi_plot}}
#' @references 
#' 
#' Werner, G.D.A., Cornwell, W.K., Sprent, J.I., Kattge, J. & Kiers, E.T. (2014).
#' A single evolutionary innovation drives the deep evolution of symbiotic N2-fixation
#' in angiosperms. Nature Communications, 5, 4087.
#'   
#' Ho, L. S. T. and Ane, C. 2014. "A linear-time algorithm for 
#' Gaussian and non-Gaussian trait evolution models". Systematic Biology 63(3):397-408.
#' 
#' @import ape phylolm
#' 
#' @examples 
#' # Load data:
#' data(alien)
#' # Run analysis:
#' samp <- intra_samp_phylm(gestaLen ~ adultMass, phy = alien$phy[[1]],
#'                          y.transf = log,x.transf = NULL,Vy="SD_gesta",Vx=NULL,
#'                          data = alien$data, n.intra = 5, n.sim=10)
#' summary(samp)
#' head(samp$sensi.estimates)
#' # Visual diagnostics
#' \dontrun{
#' sensi_plot(samp)
#' # You can specify which graph and parameter ("estimate" or "intercept") to print: 
#' sensi_plot(samp, graphs = 1)
#' sensi_plot(samp, graphs = 2)
#' }
#' @export


intra_samp_phylm <- function(formula, data, phy, n.sim=10, n.intra = 3,
                             breaks=seq(.1,.5,.1),model="lambda",
                             Vy = NULL, Vx = NULL, distrib = "normal",
                             y.transf = NULL, x.transf = NULL, track=TRUE,...) { 
  
  #Error check
  if(is.null(Vx) & is.null(Vy)) stop("Vx or Vy must be defined")
  if(class(formula) != "formula") stop("formula must be class 'formula'")
  if(class(data) != "data.frame") stop("data must be class 'data.frame'")
  if(class(phy) != "phylo") stop("phy must be class 'phylo'")
  if(formula[[2]]!=all.vars(formula)[1] || formula[[3]]!=all.vars(formula)[2])
    stop("Please use arguments y.transf or x.transf for data transformation")
  if(distrib == "normal") warning ("distrib=normal: make sure that standard deviation is provided for Vx and/or Vy")
  if(length(breaks) < 2) stop("Please include more than one break, e.g. breaks=c(.3,.5)")
  if((model == "trend") & (sum(is.ultrametric(phy))>1)) 
    stop("Trend is unidentifiable for ultrametric trees., see ?phylolm for details")
  else


  #Matching tree and phylogeny using utils.R
  datphy <- match_dataphy(formula, data, phy,...)
  full.data <- datphy[[1]]
  phy <- datphy[[2]]
  
  resp <- all.vars(formula)[1]
  pred <- all.vars(formula)[2]
  
  if(!is.null(Vy) && sum(is.na(full.data[, Vy])) != 0) {
    full.data[is.na(full.data[, Vy]), Vy] <- 0}
  
  if(!is.null(Vx) && sum(is.na(full.data[, Vx])) != 0) {
    full.data[is.na(full.data[, Vx]), Vx] <- 0}
  
  #Function to pick a random value in the interval
  if (distrib == "normal") funr <- function(a,b) {stats::rnorm(1,a,b)}
  else  funr <- function(a,b) {stats::runif(1, a - b, a + b)}
  
  #List to store information
  intra.samp <- list()
  species.NA <- list()

  #Start intra loop here
  errors <- NULL
  if(track==TRUE) pb <- utils::txtProgressBar(min = 0, max = n.intra, style = 3)
  counter = 1
  
  for (i in 1:n.intra) {
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
    
    #skip iteration if there are NA's in the dataset
    species.NA[[i]]<-rownames(full.data[with(full.data,is.na(predV) | is.na(respV)),])
    if(sum(is.na(full.data[,c("respV","predV")])>0)) next
    
    #Run the model
    intra.samp[[i]] <- samp_phylm(respV~predV, data = full.data, phy=phy, n.sim = n.sim,
                                    model, breaks=breaks, track = FALSE, verbose = FALSE,...)

    if(track==TRUE) utils::setTxtProgressBar(pb, counter)
    counter = counter + 1
  }
  
  if(track==TRUE) close(pb)
  names(intra.samp)<-1:n.intra
  
  # Merge lists into data.frames between iterations:
  full.estimates  <- suppressWarnings(recombine(intra.samp, slot1 = 4, slot2 = 1))
  influ.estimates <- recombine(intra.samp, slot1 = 5)
  influ.estimates$info <- NULL
  perc.sign <- recombine(intra.samp, slot1 = 6)
  perc.sign$info <- NULL
  
  
  #Generates output:
  res <- list(call = match.call(),
              model = model,
              formula = formula,
              full.model.estimates = full.estimates,
              sensi.estimates = influ.estimates,
              sign.analysis = perc.sign,
              data = full.data)
  
  
  class(res) <- "sensiIntra_Samp"
  ### Warnings:
  if (length(res$errors) >0){
    warning("Some species deletion presented errors, please check: output$errors")}
  else {
    res$errors <- "No errors found."
  }
  return(res)
}
  
  