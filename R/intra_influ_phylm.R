#' Interaction of intraspecific variability & influential species - Phylogenetic Linear Regression
#'
#' Performs leave-one-out deletion analyis for phylogenetic linear regression,
#' and detects influential species, while taking into account potential
#' interactions with intraspecific variability.
#'
#' @param formula The model formula:  \code{response~predictor}. 
#' @param data Data frame containing species traits with row names matching tips
#' in \code{phy}.
#' @param phy A phylogeny (class 'phylo') matching \code{data}.
#' @param model The phylogenetic model to use (see Details). Default is \code{lambda}.
#' @param cutoff The cutoff value used to identify for influential species
#' (see Details)
#' @param n.intra Number of datasets resimulated taking into account intraspecific variation (see: \code{"intra_phylm"})
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
#' @param track Print a report tracking function progress (default = TRUE)
#' @param ... Further arguments to be passed to \code{phylolm}
#' @details
#' This function fits a phylogenetic linear regression model using \code{\link[phylolm]{phylolm}}, and detects
#' influential species by sequentially deleting one at a time. The regression is repeated \code{n.intra} times for 
#' simulated values of the dataset, taking into account intraspecific variation. At each iteration, the function 
#' generates a random value for each row in the dataset using the standard deviation or errors supplied, and 
#' detect the influential species within that iteration. 
#'
#' All phylogenetic models from \code{phylolm} can be used, i.e. \code{BM},
#' \code{OUfixedRoot}, \code{OUrandomRoot}, \code{lambda}, \code{kappa},
#' \code{delta}, \code{EB} and \code{trend}. See ?\code{phylolm} for details.
#'
#' \code{influ_phylm} detects influential species based on the standardised
#' difference in intercept and/or slope when removing a given species compared
#' to the full model including all species. Species with a standardised difference
#' above the value of \code{cutoff} are identified as influential. The default
#' value for the cutoff is 2 standardised differences change.
#'
#' Currently, this function can only implement simple linear models (i.e. \eqn{trait~
#' predictor}). In the future we will implement more complex models.
#'
#' Output can be visualised using \code{sensi_plot}.
#' 
#' @section Warning:  
#' When Vy or Vx exceed Y or X, respectively, negative (or null) values can be generated, this might cause problems
#' for data transformation (e.g. log-transformation). In these cases, the function will skip the simulation. 
#' 
#' Setting \code{n.intra} at high values can take a long time to exectue, since the total number of iterations 
#' equals \code{n.intra * nrow(data)}.
#'
#' The function returns a list of \code{sensiInflu}-objects (the output of \code{influ_phylm} and \code{influ_phyglm}). 
#' of length \code{n.intra}. The user can use \code{summary} to evaluate this list. This will give, both for the 
#' regression slope and for the intercept, a table indicating how often across the \code{n.times} simulations a given
#' species was identified as the most influential species, as well as a table listing the mean estimate (slope), DIFestimate, 
#' Percentage change and P-value across all species that occured as most influential species in at least one simulation.
#' 
#' Additionally, users can evaluate each element in the list as a regula \code{sensiInflu}-object. 
#'
#' @author Gustavo Paterno, Caterina Penone & Gijsbert D.A. Werner
#' @seealso \code{\link[phylolm]{phylolm}}, \code{\link{samp_phylm}},
#' \code{\link{influ_phylm}},\code{\link{intra_phylm}},\code{\link{sensi_plot}}.
#' @references Ho, L. S. T. and Ane, C. 2014. "A linear-time algorithm for 
#' Gaussian and non-Gaussian trait evolution models". Systematic Biology 63(3):397-408.
#' @examples 
#' # Load data:
#' data(alien)
#' # run analysis:
#' intra_influ <- intra_influ_phylm(formula = gestaLen ~ adultMass, phy = alien$phy[[1]],
#' data=alien$data,model="lambda",y.transf = log,x.transf = NULL,Vy="SD_gesta",Vx=NULL,
#' n.intra=30,distrib = "normal")
#' summary(intra_influ)
#' sensi_plot(intra_influ)
#'  
#' @export

intra_influ_phylm <- function(formula, data, phy,
                        Vy = NULL, Vx = NULL,
                        y.transf = NULL, x.transf = NULL,
                        n.intra = 10, distrib = "normal",
                        model = "lambda", cutoff=2,
                        track = TRUE, ...){
  #Error check
  if(is.null(Vx) & is.null(Vy)) stop("Vx or Vy must be defined")
  if(class(formula) != "formula") stop("formula must be class 'formula'")
  if(class(data) != "data.frame") stop("data must be class 'data.frame'")
  if(class(phy) != "phylo") stop("phy must be class 'phylo'")
  if(formula[[2]]!=all.vars(formula)[1] || formula[[3]]!=all.vars(formula)[2])
    stop("Please use arguments y.transf or x.transf for data transformation")
  if(distrib == "normal") warning ("distrib=normal: make sure that standard deviation is provided for Vx and/or Vy")
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
  intra.influ <- list ()
  N  <- nrow(full.data)
  
  #Start intra loop here
  species.NA <- list()
  errors <- NULL
  pb <- utils::txtProgressBar(min = 0, max = N*n.intra, style = 3)
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
    intra.influ[[i]] <- influ_phylm(respV~predV, data = full.data, phy=phy, 
                                   model, cutoff, track = FALSE, verbose = FALSE)
    
    if(track==TRUE) utils::setTxtProgressBar(pb, counter)
    counter = counter + N
  }
  
  names(intra.influ)<-1:n.intra
  
  on.exit(close(pb))
  
  # Merge lists into data.frames between iterations:
  full.estimates  <- suppressWarnings(recombine(intra.influ, slot1 = 4, slot2 = 1))
  
  #influ species slope
  influ.sp.estimate <- (lapply(intra.influ,function(x) x$influential.species$influ.sp.estimate))
  influ.sp.estimate <- as.data.frame(as.matrix(influ.sp.estimate))
  names(influ.sp.estimate) <- "influ.sp.estimate"
  influ.sp.estimate$tree<-row.names(influ.sp.estimate)
  
  #influ species intercept
  influ.sp.intercept <- (lapply(intra.influ,function(x) x$influential.species$influ.sp.intercept))
  influ.sp.intercept <- as.data.frame(as.matrix(influ.sp.intercept))
  names(influ.sp.intercept) <- "influ.sp.intercept"
  influ.sp.intercept$tree<-row.names(influ.sp.intercept)
  
  #influ.estimates
  influ.estimates <- recombine(intra.influ, slot1 = 6)
  
  #Generates output:
  res <- list(call = match.call(),
              cutoff=cutoff,
              formula = formula,
              full.model.estimates = full.estimates,
              influential.species = list(influ.sp.estimate=influ.sp.estimate,influ.sp.intercept=influ.sp.intercept),
              sensi.estimates = influ.estimates,
              data = full.data)

  
  class(res) <- "sensiIntra_Influ"
  
  ### Warnings:
  if (length(res$errors) >0){
    warning("Some species deletion presented errors, please check: output$errors")}
  else {
    res$errors <- "No errors found."
  }
  
  return(res)
}
