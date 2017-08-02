#' Interaction between intraspecific variability and species sampling - Phylogenetic Logistic Regression
#'
#' Performs analyses of sensitivity to species sampling by randomly removing
#' species and detecting the effects on parameter estimates in a phylogenetic
#' logistic regression, while taking into account potential
#' interactions with intraspecific variability.
#'
#' @param formula The model formula
#' @param data Data frame containing species traits with row names matching tips
#' in \code{phy}.
#' @param phy A phylogeny (class 'phylo') matching \code{data}.
#' @param n.sim The number of times species are randomly deleted for each
#' \code{break}.
#' @param n.intra Number of datasets resimulated taking into account intraspecific variation (see: \code{"intra_phyglm"}) 
#' @param breaks A vector containing the percentages of species to remove.
#' @param model The phylogenetic model to use (see Details). Default is \code{lambda}.
#' @param Vx Name of the column containing the standard deviation or the standard error of the predictor 
#' variable. When information is not available for one taxon, the value can be 0 or \code{NA}
#' @param y.transf Transformation for the response variable (e.g. \code{"log"} or \code{"sqrt"}). Please use this 
#' argument instead of transforming data in the formula directly (see also details below).
#' @param x.transf Transformation for the predictor variable (e.g. \code{"log"} or \code{"sqrt"}). Please use this 
#' argument instead of transforming data in the formula directly (see also details below).
#' @param distrib A character string indicating which distribution to use to generate a random value for the response 
#' and/or predictor variables. Default is normal distribution: "normal" (function \code{\link{rnorm}}).
#' Uniform distribution: "uniform" (\code{\link{runif}})
#' Warning: we recommend to use normal distribution with Vx = standard deviation of the mean.
#' @param track Print a report tracking function progress (default = TRUE)
#' @param ... Further arguments to be passed to \code{phylolm}
#' @details
#'
#' This function randomly removes a given percentage of species (controlled by
#' \code{breaks}) from the full phylogenetic logistic regression, fits a phylogenetic
#' logistic regression model without these species using \code{\link[phylolm]{phylolm}},
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
#' Currently, this function can only implement simple logistic models (i.e. \eqn{trait~
#' predictor}). In the future we will implement more complex models.
#'
#' Output can be visualised using \code{sensi_plot}.
#' @return The function \code{samp_phyglm} returns a list with the following
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
#' @seealso \code{\link[phylolm]{phylolm}}, \code{\link{samp_phyglm}},
#' \code{\link{samp_phyglm}},\code{\link{sensi_plot}}
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
#' set.seed(6987)
#' phy = rtree(100)
#' x = rTrait(n=1,phy=phy,parameters=list(ancestral.state=2,optimal.value=2,sigma2=1,alpha=1))
#' X = cbind(rep(1,100),x)
#' y = rbinTrait(n=1,phy=phy, beta=c(-1,0.5), alpha=.7 ,X=X)
#' z = rnorm(n = length(x),mean = mean(x),sd = 0.1*mean(x))
#' dat = data.frame(y, x, z)
#' #Run sensitivity analysis:
#' intra_samp <- intra_samp_phyglm(formula = y ~ x, data = dat, phy = phy, 
#'                                n.sim=10, n.intra = 3,
#'                                breaks=seq(.1,.5,.1),
#'                                Vx = "z", distrib="normal",x.transf=NULL)
#' summary(samp)
#' sensi_plot.sensiSamp(samp)

#' @export


intra_samp_phyglm <- function(formula, data, phy, n.sim=10, n.intra = 3,
                             breaks=seq(.1,.5,.1), 
                             Vx = NULL, distrib = "normal", x.transf = NULL, 
                             btol = 50, track=TRUE,...) { 
  
  #Error check
  if(is.null(Vx)) stop("Vx must be defined")
  if(class(formula) != "formula") stop("formula must be class 'formula'")
  if(class(data) != "data.frame") stop("data must be class 'data.frame'")
  if(class(phy) != "phylo") stop("phy must be class 'phylo'")
  if(formula[[2]]!=all.vars(formula)[1] || formula[[3]]!=all.vars(formula)[2])
    stop("Please use arguments y.transf or x.transf for data transformation")
  if(distrib == "normal") warning ("distrib=normal: make sure that standard deviation is provided for Vx")
  if(length(breaks) < 2) stop("Please include more than one break, e.g. breaks=c(.3,.5)")

  #Matching tree and phylogeny using utils.R
  datphy<-match_dataphy(formula,data,phy, ...)
  full.data<-datphy[[1]]
  phy<-datphy[[2]]
  
  resp1<-all.vars(formula)[1]
  if(length(all.vars(formula))>2){resp2<-all.vars(formula)[2]}
  pred<-all.vars(formula)[length(all.vars(formula))]
  
  if(!is.null(Vx) && sum(is.na(full.data[,Vx]))!=0){
    full.data[is.na(full.data[,Vx]), Vx] <- 0}
  
  #Function to pick a random value in the interval
  if (distrib=="normal") funr <- function(a, b) {stats::rnorm(1,a,b)}
  else  funr <- function(a,b) {stats::runif(1,a-b,a+b)}
  
  
  #List to store information
  intra.samp <- list()
  species.NA <- list()
  
  #Start intra loop here
  errors <- NULL
  counter = 1
  
  for (i in 1:n.intra) {
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
    
    #skip iteration if there are NA's in the dataset
    species.NA[[i]]<-rownames(full.data[with(full.data,is.na(predV)),])
    if(sum(is.na(full.data[,"predV"])>0)) next
    
    #Run the model
    if(length(all.vars(formula))>2) {
      intra.samp[[i]] <- samp_phyglm(cbind(resp1,resp2)~predV, data = full.data, phy=phy, n.sim = n.sim,
                                     breaks=breaks, btol=btol, method="logistic_MPLE",
                                     track = FALSE, verbose = FALSE,...)
    } else 
      intra.samp[[i]] <- samp_phyglm(resp1~predV, data = full.data, phy=phy, n.sim = n.sim,
                                     breaks=breaks, btol=btol, method="logistic_MPLE",
                                     track = FALSE, verbose = FALSE,...)
    
    
    
    counter = counter + 1
  }
  
  names(intra.samp)<-1:n.intra
  
  # Merge lists into data.frames between iterations:
  full.estimates  <- suppressWarnings(recombine(intra.samp, slot1 = 4, slot2 = 1))
  influ.estimates <- recombine(intra.samp, slot1 = 5)
  perc.sign <- recombine(intra.samp, slot1 = 6)
  
  #Generates output:
  res <- list(call = match.call(),
              model = "logistic_MPLE",
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

