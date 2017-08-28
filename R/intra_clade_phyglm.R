#' Interaction between intraspecific variability and influential clade detection - Phylogenetic Logistic Regression
#'
#' Estimate the impact on model estimates of phylogenetic logistic regression after 
#' removing clades from the analysis, while taking into account potential
#' interactions with intraspecific variability.
#'
#' @param formula The model formula
#' @param data Data frame containing species traits with row names matching tips
#' in \code{phy}.
#' @param phy A phylogeny (class 'phylo') matching \code{data}.
#' @param clade.col The column in the provided data frame which specifies the
#' clades (a character vector with clade names).
#' @param n.species Minimum number of species in a clade for the clade to be
#' included in the leave-one-out deletion analyis. Default is \code{5}.
#' @param n.sim Number of simulations for the randomization test.
#' @param n.intra Number of datasets resimulated taking into account intraspecific variation (see: \code{"intra_phyglm"})
#' @param Vx Name of the column containing the standard deviation or the standard error of the predictor 
#' variable. When information is not available for one taxon, the value can be 0 or \code{NA}
#' @param x.transf Transformation for the predictor variable (e.g. \code{"log"} or \code{"sqrt"}). Please use this 
#' argument instead of transforming data in the formula directly (see also details below).
#' @param distrib A character string indicating which distribution to use to generate a random value for the response 
#' and/or predictor variables. Default is normal distribution: "normal" (function \code{\link{rnorm}}).
#' Uniform distribution: "uniform" (\code{\link{runif}})
#' Warning: we recommend to use normal distribution with Vx or Vy = standard deviation of the mean.
#' @param btol Bound on searching space. For details see \code{phyloglm}
#' @param track Print a report tracking function progress (default = TRUE)
#' @param ... Further arguments to be passed to \code{phyloglm}
#' @details
#' This function sequentially removes one clade at a time, fits a phylogenetic
#' logistic regression model using \code{\link[phylolm]{phyloglm}} and stores the
#' results. The impact of of a specific clade on model estimates is calculated by the
#' comparison between the full model (with all species) and the model without 
#' the species belonging to a clade. This operation is repeated \code{n.intra} times for
#' simulated values of the dataset, taking into account intraspecific variation. At each iteration, the function 
#' generates a random value for each row in the dataset using the standard deviation or errors supplied, and 
#' detect the influential species within that iteration. 
#' 
#' Additionally, to account for the influence of the number of species on each 
#' clade (clade sample size), this function also estimate a null distribution of slopes
#' expected for the number of species in a given clade. This is done by fitting
#' models without the same number of species in the given clade. 
#'  The number of simulations to be performed is set by 'n.sim'. To test if the 
#'  clade influence differs from the null expectation for a clade of that size, 
#'  a randomization test can be performed using 'summary(x)'. 
#'
#' All phylogenetic models from \code{phyloglm} can be used, i.e. \code{BM},
#' \code{OUfixedRoot}, \code{OUrandomRoot}, \code{lambda}, \code{kappa},
#' \code{delta}, \code{EB} and \code{trend}. See ?\code{phyloglm} for details.
#'
#' \code{clade_phyglm} detects influential clades based on
#' difference in intercept and/or slope when removing a given clade compared
#' to the full model including all species.
#' 
#' Currently, this function can only implement simple logistic models (i.e. 
#' \eqn{y = a + bx}). In the future we will implement more complex models.
#'
#' Output can be visualised using \code{sensi_plot}.
#'
#' @return The function \code{clade_phyglm} returns a list with the following
#' components:
#' @return \code{formula}: The formula
#' @return \code{full.model.estimates}: Coefficients, aic and the optimised
#' value of the phylogenetic parameter (e.g. \code{lambda}) for the full model
#' without deleted species.
#' @return \code{sensi.estimates}: A data frame with all simulation
#' estimates. Each row represents a deleted clade. Columns report the calculated
#' regression intercept (\code{intercept}), difference between simulation
#' intercept and full model intercept (\code{DIFintercept}), the percentage of change
#' in intercept compared to the full model (\code{intercept.perc}) and intercept
#' p-value (\code{pval.intercept}). All these parameters are also reported for the regression
#' slope (\code{DIFestimate} etc.). Additionally, model aic value (\code{AIC}) and
#' the optimised value (\code{optpar}) of the phylogenetic parameter 
#' (e.g. \code{kappa} or \code{lambda}, depending on the phylogenetic model used) 
#' are reported.
#' @return \code{data}: Original full dataset.
#' @return \code{errors}: Clades and/or iterations where deletion resulted in errors.
#' @author Gustavo Paterno, Caterina Penone
#' @seealso \code{\link[phylolm]{phyloglm}}, \code{\link[sensiPhy]{samp_phyglm}},
#' \code{\link{influ_phyglm}}, \code{\link{sensi_plot}}
#' \code{\link{sensi_plot}}, \code{\link{intra_phyglm}}
#' @references Ho, L. S. T. and Ane, C. 2014. "A linear-time algorithm for 
#' Gaussian and non-Gaussian trait evolution models". Systematic Biology 63(3):397-408.
#' @examples
#' set.seed(6987)
#' phy = rtree(100)
#' x = rTrait(n=1,phy=phy,parameters=list(ancestral.state=2,optimal.value=2,sigma2=1,alpha=1))
#' X = cbind(rep(1,100),x)
#' y = rbinTrait(n=1,phy=phy, beta=c(-1,0.5), alpha=.7 ,X=X)
#' z = rnorm(n = length(x),mean = mean(x),sd = 0.1*mean(x))
#' cla <- rep(c("A","B","C","D"), each = 25)
#' dat = data.frame(y, x, z, cla)

#' intra_clade <- intra_clade_phyglm(formula=y ~ x, data = dat, phy = phy,
#'                                   clade.col = "cla", n.sim = 30, n.intra = 3,
#'                                   y.transf = log, Vx = "z", distrib="normal")
#' summary(intra_clade)
#' sensi.plot(intra_clade)

#' @export


intra_clade_phyglm <- function(formula, data, phy, clade.col, n.species = 5,
                              n.sim = 100, n.intra = 2,
                              Vy = NULL, Vx = NULL, distrib = "normal",
                              y.transf = NULL, x.transf = NULL,
                              model = "lambda", track = TRUE,...) {
  
  # Error checking:
  if(is.null(Vx)) stop("Vx must be defined")
  if(!is.data.frame(data)) stop("data must be class 'data.frame'")
  if(missing(clade.col)) stop("clade.col not defined. Please, define the column with clade names.")
  if(class(phy)!="phylo") stop("phy must be class 'phylo'")
  if(class(formula)!="formula") stop("formula must be class 'formula'")
  if(formula[[2]]!=all.vars(formula)[1] || formula[[3]]!=all.vars(formula)[2])
    stop("Please use arguments y.transf or x.transf for data transformation")
  if(distrib=="normal") warning ("distrib=normal: make sure that standard deviation is provided for Vx")
  
  #Match data and phy
  data_phy <- match_dataphy(formula, data, phy,...)
  phy <- data_phy$phy
  full.data <- data_phy$data
  if (is.na(match(clade.col, names(full.data)))) {
    stop("Names column '", clade.col, "' not found in data frame'")
  }
  
  resp1<-all.vars(formula)[1]
  if(length(all.vars(formula))>2){resp2<-all.vars(formula)[2]}
  pred<-all.vars(formula)[length(all.vars(formula))]
  
  if(!is.null(Vx) && sum(is.na(full.data[,Vx]))!=0){
    full.data[is.na(full.data[,Vx]), Vx] <- 0}
  
  #Function to pick a random value in the interval
  if (distrib=="normal") funr <- function(a, b) {stats::rnorm(1,a,b)}
  else  funr <- function(a,b) {stats::runif(1,a-b,a+b)}  
  
  #Identify CLADES to use and their sample size 
  all.clades <- levels(full.data[ ,clade.col])
  wc <- table(full.data[ ,clade.col]) > n.species
  uc <- table(full.data[ , clade.col])[wc]
  
  if (length(uc) == 0) stop(paste("There is no clade with more than ",
                                  n.species," species. Change 'n.species' to fix this
                                  problem",sep=""))
  
  #List to store information
  intra.clade <- list()
  
  #Start clade loop here
  errors <- NULL
  if(track==TRUE) pb <- utils::txtProgressBar(min = 0, max = n.intra, style = 3)
  counter = 1
  
  for (j in 1:n.intra){
    
    ##Set predictor variable
    #Vx is not provided or is not numeric, do not pick random value
    if(!inherits(full.data[,pred], c("numeric","integer")) || is.null(Vx)) {full.data$predV<-full.data[,pred]}
    
    #choose a random value in [mean-se,mean+se] if Vx is provided
    if(!is.null(Vx) && is.null(dim(Vx)))
    {full.data$predV<-apply(full.data[,c(pred,Vx)],1,function(x)funr(x[1],x[2]))}
    
    full.data$resp1<-full.data[,resp1] #try to improve this in future
    if(length(all.vars(formula))>2){full.data$resp2<-full.data[,resp2]}
    
    #transform Vx if x.transf is provided
    if(!is.null(x.transf)) 
    {suppressWarnings (full.data$predV <- x.transf(full.data$predV))}
    
    if(length(all.vars(formula))>2) {
    intra.clade[[j]] <- clade_phyglm(cbind(resp1,resp2)~predV, data=full.data, phy=phy, method="logistic_MPLE",
                                    clade.col=clade.col, n.species=n.species, n.sim=n.sim, 
                                    track = FALSE, verbose = FALSE,...)}
    else 
    intra.clade[[j]] <- clade_phyglm(resp1~predV, data = full.data, phy=phy, method="logistic_MPLE",
                                    clade.col=clade.col, n.species=n.species, n.sim=n.sim, 
                                    track = FALSE, verbose = FALSE,...)
    
    if(track==TRUE) utils::setTxtProgressBar(pb, counter)
    counter = counter + 1
  }
  
  names(intra.clade) <- 1:n.intra
  
  # Merge lists into data.frames between iterations:
  full.estimates  <- suppressWarnings(recombine(intra.clade, slot1 = 3, slot2 = 1))
  clade.estimates <- recombine(intra.clade, slot1 = 4)
  clade.estimates$info <- NULL
  null.dist       <- recombine(intra.clade, slot1 = 5)
  null.dist$info <- NULL
  
  #Generates output:
  res <- list(call = match.call(),
              formula = formula,
              full.model.estimates = full.estimates,
              sensi.estimates = clade.estimates,
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

