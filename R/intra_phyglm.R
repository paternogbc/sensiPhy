#' Intraspecific variability - Phylogenetic Logistic Regression
#'
#' Performs Phylogenetic logistic regression evaluating
#' intraspecific variability in predictor variables.
#'
#' @param formula The model formula: \code{response~predictor}. 
#' @param data Data frame containing species traits with species as row names.
#' @param phy A phylogeny (class 'phylo', see ?\code{ape}).
#' @param Vx Name of the column containing the standard deviation or the standard error of the predictor 
#' variable. When information is not available for one taxon, the value can be 0 or \code{NA}
#' @param x.transf Transformation for the predictor variable (e.g. \code{log} or \code{sqrt}). Please use this 
#' argument instead of transforming data in the formula directly (see also details below).
#' @param times Number of times to repeat the analysis generating a random value for the predictor variable.
#' If NULL, \code{times} = 2
#' @param distrib A character string indicating which distribution to use to generate a random value for the
#' predictor variable. Default is normal distribution: "normal" (function \code{\link{rnorm}}).
#' Uniform distribution: "uniform" (\code{\link{runif}})
#' Warning: we recommend to use normal distribution with Vx = standard deviation of the mean.
#' @param btol Bound on searching space. For details see \code{phyloglm}
#' @param track Print a report tracking function progress (default = TRUE)
#' @param ... Further arguments to be passed to \code{phyloglm}
#' @details
#' This function fits a phylogenetic logistic regression model using \code{\link[phylolm]{phyloglm}}.
#' The regression is repeated \code{times} times. At each iteration the function generates a random value
#' for each row in the dataset using the standard deviation or error supplied and assuming a normal or uniform distribution.
#' To calculate means and se for your raw data, you can use the \code{\link[Rmisc]{summarySE}} function.
#' Warning: when Vx exceeds X, negative (or null) values can be generated, this might cause problems
#' for data transformation (e.g. log-transformation). In these cases, the function will skip the simulation. This problem can
#' be solved by increasing \code{times}, changing the transformation type and/or checking the target species in output$sp.pb.
#'
#' All phylogenetic models from \code{phyloglm} can be used, i.e. \code{BM},
#' \code{OUfixedRoot}, \code{OUrandomRoot}, \code{lambda}, \code{kappa},
#' \code{delta}, \code{EB} and \code{trend}. See ?\code{phyloglm} for details.
#'
#' Currently, this function can only implement simple logistic models (i.e. \eqn{trait~
#' predictor}). In the future we will implement more complex models.
#'
#' Output can be visualised using \code{sensi_plot}.
#'
#' @return The function \code{intra_phyglm} returns a list with the following
#' components:
#' @return \code{formula}: The formula
#' @return \code{data}: Original full dataset
#' @return \code{model_results}: Coefficients, aic and the optimised
#' value of the phylogenetic parameter (e.g. \code{lambda}) for each regression.
#' @return \code{N.obs}: Size of the dataset after matching it with tree tips and removing NA's.
#' @return \code{stats}: Main statistics for model parameters.\code{CI_low} and \code{CI_high} are the lower 
#' and upper limits of the 95% confidence interval.
#' @return \code{all.stats}: Complete statistics for model parameters. \code{sd_intra} is the standard deviation 
#' due to intraspecific variation. \code{CI_low} and \code{CI_high} are the lower and upper limits 
#' of the 95% confidence interval.
#' @return \code{sp.pb}: Species that caused problems with data transformation (see details above).
#' 
#' @author Caterina Penone & Pablo Ariel Martinez
#' @seealso \code{\link[phylolm]{phyloglm}}, \code{\link{sensi_plot}}
#' @references
#' Martinez, P. a., Zurano, J.P., Amado, T.F., Penone, C., Betancur-R, R., 
#' Bidau, C.J. & Jacobina, U.P. (2015). Chromosomal diversity in tropical reef 
#' fishes is related to body size and depth range. Molecular Phylogenetics and 
#' Evolution, 93, 1-4
#' 
#' Ho, L. S. T. and Ane, C. 2014. "A linear-time algorithm for 
#' Gaussian and non-Gaussian trait evolution models". Systematic Biology 63(3):397-408.
#' @export


intra_phyglm <- function(formula, data, phy,
                          Vx=NULL, times = 30,
                          x.transf = NULL,
                          distrib="uniform", btol=50, track=TRUE,...){
  #Error check
  if(is.null(Vx)) stop("Vx must be defined")
  if(class(formula)!="formula") stop("formula must be class 'formula'")
  if(class(data)!="data.frame") stop("data must be class 'data.frame'")
  if(class(phy)!="phylo") stop("phy must be class 'phylo'")
  if(formula[[2]]!=all.vars(formula)[1] || formula[[3]]!=all.vars(formula)[2])
     stop("Please use argument x.transf for data transformation")
  if(distrib=="normal") warning ("distrib=normal: make sure that standard deviation 
                                 is provided for Vx")

  #Matching tree and phylogeny using utils.R
  datphy<-match_dataphy(formula,data,phy)
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
  
  
  #Create the results data.frame
  intra.model.estimates<-data.frame("n.intra"=numeric(),"intercept"=numeric(),"se.intercept"=numeric(),
                                    "pval.intercept"=numeric(),"slope"=numeric(),"se.slope"=numeric(),
                                    "pval.slope"=numeric(),"aic"=numeric(),"optpar"=numeric())
  

  #Model calculation
  counter=1
  errors <- NULL
  species.NA <- list()
  pb <- utils::txtProgressBar(min = 0, max = times, style = 1)
  for (i in 1:times) {
    
    ##Set predictor variable
    #Vx is not provided or is not numeric, do not pick random value
    if(!inherits(full.data[,pred], c("numeric","integer")) || is.null(Vx)) {full.data$predV<-full.data[,pred]}

    #choose a random value in [mean-se,mean+se] if Vx is provided
    if(!is.null(Vx) && is.null(dim(Vx)))
    {full.data$predV<-apply(full.data[,c(pred,Vx)],1,function(x)funr(x[1],x[2]))}
    
    full.data$resp1<-full.data[,resp1] #try to improve this in future
    if(length(all.vars(formula))>2){full.data$resp2<-full.data[,resp2]}
    
    #transform Vy and/or Vx if x.transf and/or y.transf are provided
    if(!is.null(x.transf)) 
    {suppressWarnings (full.data$predV <- x.transf(full.data$predV))}
    
    #skip iteration if there are NA's in the dataset
    species.NA[[i]]<-rownames(full.data[with(full.data,is.na(predV)),])
    if(sum(is.na(full.data[,"predV"])>0)) next
  
    #model
    if(length(all.vars(formula))>2){mod = try(phylolm::phyloglm(cbind(resp1,resp2)~predV, data=full.data, phy=phy,method="logistic_MPLE",btol=btol),FALSE)}
    else
    mod = try(phylolm::phyloglm(resp1~predV, data=full.data, phy=phy,method="logistic_MPLE",btol=btol),FALSE)

    if(isTRUE(class(mod)=="try-error")) {
      error <- i
      errors <- c(errors,error)
      next }
    
    
    else{
      intercept            <- phylolm::summary.phyloglm(mod)$coefficients[[1,1]]
      se.intercept         <- phylolm::summary.phyloglm(mod)$coefficients[[1,2]]
      slope                <- phylolm::summary.phyloglm(mod)$coefficients[[2,1]]
      se.slope             <- phylolm::summary.phyloglm(mod)$coefficients[[2,2]]
      pval.intercept       <- phylolm::summary.phyloglm(mod)$coefficients[[1,4]]
      pval.slope           <- phylolm::summary.phyloglm(mod)$coefficients[[2,4]]
      aic.mod              <- mod$aic
      n                    <- mod$n
      #d                   <- mod$d
      optpar               <- mod$alpha

      
      if(track==TRUE) utils::setTxtProgressBar(pb, i)

      #write in a table
      estim.simu <- data.frame(i, intercept, se.intercept, pval.intercept,
                               slope, se.slope, pval.slope, aic.mod, optpar,
                               stringsAsFactors = F)
      intra.model.estimates[counter, ]  <- estim.simu
      counter=counter+1
      
    }
  }
  on.exit(close(pb))
  
  #calculate mean and sd for each parameter
  #variation due to intraspecific variability
  mean_by_randomval<-stats::aggregate(.~n.intra, data=intra.model.estimates, mean)
  
  statresults<-data.frame(min=apply(intra.model.estimates,2,min),
                          max=apply(intra.model.estimates,2,max),
                          mean=apply(intra.model.estimates,2,mean),
                          sd_intra=apply(mean_by_randomval,2,stats::sd))[-1,]
  
  statresults$CI_low  <- statresults$mean - stats::qt(0.975, df = times-1) * statresults$sd_intra / sqrt(times)
  statresults$CI_high <- statresults$mean + stats::qt(0.975, df = times-1) * statresults$sd_intra / sqrt(times)
  
  #species with transformation problems
  nr <- times - nrow(intra.model.estimates)
  sp.pb <- unique(unlist(species.NA))
  if (length(sp.pb) >0) 
    warning (paste("in", nr,"simulations, data transformations generated NAs, please consider using another function
  for x.transf and check output$sp.pb",sep=" "))
  
  res <- list(call = match.call(),
              formula=formula,
              data=full.data,
              model_results=intra.model.estimates,N.obs=n,
              stats = round(statresults[c(1:6),c(3,5,6)],digits=3),
              all.stats = statresults,sp.pb=sp.pb)
  class(res) <- c("sensiIntra","sensiIntraL")
  return(res)
}

