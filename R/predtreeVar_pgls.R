#' PGLS analysis accounting for variability in predictors and trees.
#'
#' \code{predtreeVar_pgls} performs \code{pgls} analysis accounting for variability in predictors and trees.
#' It picks n times a random value within an interval for the predictor variable and 
#' it repeats the analysis m times for m trees in a \code{multiPhylo} file. 
#' The output gives a mean and sd of the model parameters.
#' @aliases predtreeVar_pgls
#' @inheritParams ??inheritParams
#' @param resp Numeric vector containing the response variable
#' @param pred Vector containing the predictor variable
#' @param se.pred Vector containing the standard error of \code{pred}  
#' @param tree A tree or list of tree of class \code{Phylo} or \code{multiPhylo}
#' @param ntree If TRUE or class(tree)=multiPhylo, the number of times to repeat the analysis with n different 
#' trees picked randomly in the multiPhylo file. If NULL, \code{ntree}=1
#' @param npred If TRUE, the number of times to repeat the analysis generating a random value in the interval 
#' [\code{pred}-\code{se.pred},\code{pred}+\code{se.pred}]. If NULL, \code{npred}=1
#' @param method A character string indicating which method to use to generate a random value in the interval 
#' [\code{pred} - \code{se.pred}, \code{pred} + \code{se.pred}]. Default is normal distribution: "normal" 
#' (function \code{\link{rnorm}}), uniform distribution is "uniform" (\code{\link{runif}})
#' @param taxa.col A character vector of taxa names that will be used to match data rows to phylogeny tips.
#' @param lambda A value for the lambda transformation. If NULL, \code{lambda}="ML"
#' @details This functions only works for simple linear regression \eqn{y = bx +
#' a}. Future implementation will deal with more complex models.
#' @return The function \code{predtreeVar_pgls} returns a list with the following
#' components:
#' @return \code{model_results} Model parameters for each iteration
#' @return \code{stats} Mean and sd of each parameter of the model
#' @section Warning: This code is not fully checked. Please be aware.
#' If \code{ntree} and \code{npred} are set to 1, the function computes a classic pgls.
#' @seealso \code{\link{pgls}}, \code{\link{comparative.data}}, \code{\link{runif}}, \code{\link{rnorm}}
#' @examples
#' to do :)
#' @export

predtreeVar_pgls <- function(resp,pred,se.pred=NA,tree,ntree=1,npred=1,method=c("normal","uniform"),taxa.col,lambda="ML"){
  
  # Error check
  method <- match.arg(method)
  
  #Function to pick a random value in the interval
  if (method=="uniform") funr <- function(a, b) {runif(1,a-b,a+b)}
  else funr <- function(a, b) {rnorm(1,a,b)}
  
  # If the class of tree is multiphylo pick n=ntree random trees
  if(inherits(tree, "multiPhylo")){trees<-sample(length(tree),ntree,replace=F)}
  else {trees=1}
  
  #Create the results data.frame
  resultados<-data.frame("n.tree"=numeric(),"n.pred"=numeric(),"estim"=numeric(),
                         "rsq"=numeric(),"pval"=numeric(),"Aicc"=numeric(),"Lambda"=numeric())
  
  #Assemble the dataframe
  if(!inherits(pred, c("numeric","integer")) || is.na(se.pred)) {data<-data.frame(taxa.col,resp,pred)}
  else {data<-data.frame(taxa.col,resp,pred,se.pred)}
  
  #Model calculation
  counter=1
  c.data<-list()
  for (i in 1:length(trees)) {
    for (j in 1:npred){
      
      #choose a random value in [mean-sd,mean+sd]
      if(!inherits(pred, c("numeric","integer")) || is.na(se.pred)) {data$variab<-data$pred}
      else {data$variab<-apply(data[,c("pred","se.pred")],1,function(x)funr(x["pred"],x["se.pred"]))}
      
      
      
      #comparative data creation if tree is class=multiphylo
      if (inherits(tree, "multiPhylo")) {
        c.data[[i]]<-comparative.data(data=data, phy=tree[[i]], names.col="taxa.col", vcv=T, vcv.dim=3) ###
        try(ModeloSimple<- caper::pgls(resp~variab, c.data[[i]], lambda='ML'),TRUE)
      }
      
      else {
        c.data<-comparative.data(data=data, phy=tree, names.col="taxa.col", vcv=T, vcv.dim=3) ###
        try(ModeloSimple<- caper::pgls(resp~variab, c.data, lambda='ML'),TRUE)
      }
      
      if(isTRUE(class(ModeloSimple)=="try-error")) { next }
      else {
        
        #extract model coefficients
        estim<-summary(ModeloSimple)$coef[2,1]
        resq<-summary(ModeloSimple)$r.squared
        pval<-summary(ModeloSimple)$coef[2,4]
        Aicc<-ModeloSimple$aicc
        Lambda<-as.numeric(ModeloSimple$param.CI$lambda$opt)
        
        #write in a table
        resultados[counter,1]<- i
        resultados[counter,2]<- j    
        resultados[counter,3]<- estim
        resultados[counter,4]<- resq
        resultados[counter,5]<- pval
        resultados[counter,6]<- Aicc
        resultados[counter,7]<- Lambda
        counter=counter+1
      }
    }  
  }
  
  #calculate mean and sd for each parameter
  #variation due to tree choice
  mean_by_tree<-aggregate(.~n.tree, data=resultados, mean)
  #variation due to continuous trait
  mean_by_randomval<-aggregate(.~n.pred, data=resultados, mean)
  
  statresults<-data.frame(mean=apply(resultados,2,mean),
                          sd_all=apply(resultados,2,sd),
                          sd_tree=apply(mean_by_tree,2,sd),
                          sd_pred=apply(mean_by_randomval,2,sd))[-(1:2),]
  
  
  output <- list(model_results=resultados,stats=statresults)
  
  return(output)  
}

