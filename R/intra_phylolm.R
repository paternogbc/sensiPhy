#' Intraspecific variability - Phylogenetic Linear Regression
#'
#' Performs Phylogenetic linear regression evaluating
#' intraspecific variability.
#'
#' @param formula The model formula: \code{response~predictor}. 
#' If \code{minmax}=TRUE the formula should point to the columns containing the minimum values of the response 
#' and/or predictor variables i.e. \code{min.response~min.predictor}.
#' @param data Data frame containing species traits with species as row names.
#' @param phy A phylogeny (class 'phylo', see ?\code{ape}).
#' @param vari.resp Name of the column containing the standard error, the standard deviation of the response 
#' variable or its maximum value (if \code{minmax}=TRUE).
#' When information is not available for one taxon, the value can be 0 or \code{NA}.
#' @param vari.pred Name of the column containing the standard error, the standard deviation  of the predictor 
#' variable or its maximum value (if \code{minmax}=TRUE).
#' When information is not available for one taxon, the value can be 0 or \code{NA}
#' @param minmax logical; if TRUE, generates a value in the interval [min,max]. If TRUE, please select columns containing minimum
#' values in the formula and columns containig maximum values in \code{vari.resp} and/or \code{vari.pred} (default = FALSE).
#' @param nintra Number of times to repeat the analysis generating a random value for response and/or predictor variables.
#' If NULL, \code{nintra} = 2
#' @param distrib A character string indicating which distribution to use to generate a random value for the response 
#' and/or predictor variables.Default is uniform distribution: "uniform" (\code{\link{runif}})
#' normal distribution is "normal" (function \code{\link{rnorm}}).
#' Warning: normal distribution can be used oly if vari.pred is the standard 
#' deviation of the mean. If minmax=T, only "uniform" distribution is available.
#' @param model The phylogenetic model to use (see Details). Default is \code{lambda}.
#' @param track Print a report tracking function progress (default = TRUE)
#' @param ... Further arguments to be passed to \code{phylolm}
#' @details
#' This function fits a phylogenetic linear regression model using \code{\link[phylolm]{phylolm}}.
#' The regression is repeated \code{nintra} times. At each iteration the functions generates for each row in the dataset
#' a random value in the interval [\code{predictor}-\code{vari.pred},\code{predictor}+\code{vari.pred}] and
#' a random value in the interval [\code{response}-\code{vari.resp},\code{response}+\code{vari.resp}].
#' If \code{minmax}=TRUE the value is randomly generated in the interval [min,max].
#' If you log-transform your predictor and/or response variables, make sure that 
#' vari.pred and/or vari.resp are in a log-scale too.
#'
#' All phylogenetic models from \code{phylolm} can be used, i.e. \code{BM},
#' \code{OUfixedRoot}, \code{OUrandomRoot}, \code{lambda}, \code{kappa},
#' \code{delta}, \code{EB} and \code{trend}. See ?\code{phylolm} for details.
#'
#' Currently, this function can only implement simple linear models (i.e. \eqn{trait~
#' predictor}). In the future we will implement more complex models.
#'
#' Output can be visualised using \code{sensi_plot}.
#'
#' @return The function \code{intra_phylm} returns a list with the following
#' components:
#' @return \code{formula}: The formula
#' @return \code{data}: Original full dataset
#' @return \code{model_results}: Coefficients, aic and the optimised
#' value of the phylogenetic parameter (e.g. \code{lambda}) for each regression with a 
#' different phylogenetic tree.
#' @return \code{N.obs}: Size of the dataset after matching it with tree tips and removing NA's.
#' @return \code{stats}: Statistics for model parameters. \code{sd_tree} is the standard deviation 
#' due to phylogenetic uncertainty.
#' 
#' @author Caterina Penone & Pablo Ariel Martinez
#' @seealso \code{\link{sensi_plot}}
#' @references Here still: reference to phylolm paper + our own?
#' @export


intra_phylm <- function(formula,data,phy,
                          vari.resp=NULL,vari.pred=NULL,minmax=FALSE,
                          nintra=2,distrib="uniform",model="lambda",track=TRUE,...){

  #Error check
  if(class(formula)!="formula") stop("formula must be class 'formula'")
  if(class(data)!="data.frame") stop("data must be class 'data.frame'")
  if(class(phy)!="phylo") stop("phy must be class 'phylo'")
  if (distrib=="normal" && minmax==T)
    stop("Cannot generate normal distribution from min and max values!")
  else
    
  #Matching tree and phylogeny using utils.R
  datphy<-match_dataphy(formula,data,phy)
  full.data<-datphy[[1]]
  phy<-datphy[[2]]
  
  resp<-all.vars(formula)[1]
  pred<-all.vars(formula)[2]
  
  if(!is.null(vari.resp) && sum(is.na(full.data[,vari.resp]))!=0){
    full.data[is.na(full.data[,vari.resp]),] <- 0}
    
    if(!is.null(vari.pred) && sum(is.na(full.data[,vari.pred]))!=0){
      full.data[is.na(full.data[,vari.pred]),] <- 0}



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
  c.data<-list()
  
  for (i in 1:nintra) {

      ##Set response and predictor variables
      #vari.resp is not provided or is not numeric, do not pick random value
      if(!inherits(full.data[,resp], c("numeric","integer")) || is.null(vari.resp)) {full.data$respV<-full.data[,resp]}
      
      #choose a random value in min/max if vari.resp is provided and minmax=T
      if(!is.null(vari.resp) && minmax==T)
      {full.data$respV<-apply(full.data[,c(resp,vari.resp)],1,function(x)stats::runif(1,x[1],x[2]))}
      
      #choose a random value in [mean-se,mean+se] if vari.resp is provided and minmax=F
      if(!is.null(vari.resp) && minmax==F)
      {full.data$respV<-apply(full.data[,c(resp,vari.resp)],1,function(x)funr(x[1],x[2]))}
      
    
      #vari.pred is not provided or is not numeric, do not pick random value
      if(!inherits(full.data[,pred], c("numeric","integer")) || is.null(vari.pred)) {full.data$predV<-full.data[,pred]}
      
      #choose a random value in min/max if vari.pred is provided and minmax=T
      if(!is.null(vari.pred) && minmax==T)
      {full.data$predV<-apply(full.data[,c(pred,vari.pred)],1,function(x)stats::runif(1,x[1],x[2]))}
      
      #choose a random value in [mean-se,mean+se] if vari.pred is provided
      if(!is.null(vari.pred) && is.null(dim(vari.pred)))
      {full.data$predV<-apply(full.data[,c(pred,vari.pred)],1,function(x)funr(x[1],x[2]))}

      #model
      mod = try(phylolm::phylolm(respV~predV, data=full.data, model=model,phy=phy),FALSE)

      if(isTRUE(class(mod)=="try-error")) {
        error <- i
        names(error) <- rownames(c.data$full.data)[i]
        errors <- c(errors,error)
        next }
      
      
      else{
        intercept            <- phylolm::summary.phylolm(mod)$coefficients[[1,1]]
        se.intercept         <- phylolm::summary.phylolm(mod)$coefficients[[1,2]]
        slope                <- phylolm::summary.phylolm(mod)$coefficients[[2,1]]
        se.slope             <- phylolm::summary.phylolm(mod)$coefficients[[2,2]]
        pval.intercept       <- phylolm::summary.phylolm(mod)$coefficients[[1,4]]
        pval.slope           <- phylolm::summary.phylolm(mod)$coefficients[[2,4]]
        aic.mod              <- mod$aic
        n                    <- mod$n
        #d                    <- mod$d
        
        if (model == "BM"){
          optpar <- NA
        }
        if (model != "BM"){
          optpar               <- mod$optpar
        }
        
        if(track==TRUE) print(paste("intra: ",i,sep=""))
        
        #write in a table
        estim.simu <- data.frame(i, intercept, se.intercept, pval.intercept,
                                 slope, se.slope, pval.slope, aic.mod, optpar,
                                 stringsAsFactors = F)
        intra.model.estimates[counter, ]  <- estim.simu
        counter=counter+1
        
      }
     }
  
  #calculate mean and sd for each parameter
  #variation due to intraspecific variability
  mean_by_randomval<-stats::aggregate(.~n.intra, data=intra.model.estimates, mean)
  
  statresults<-data.frame(min=apply(intra.model.estimates,2,min),
                          max=apply(intra.model.estimates,2,max),
                          mean=apply(intra.model.estimates,2,mean),
                          sd_intra=apply(mean_by_randomval,2,stats::sd))[-1,]
  
  
  res <- list(formula=formula,
              datas=full.data,
              model_results=intra.model.estimates,N.obs=n,
              stats=statresults)
  class(res) <- "sensiIntra"
  return(res)
}

