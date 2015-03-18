#' Pgls analysis accounting for intraspecific variability and phylogenetic uncertainty
#'
#' \code{variPgls} performs phylogenetic \code{gls} analysis accounting for intraspecific variability and phylogenetic 
#' uncertainty in trees topology. 
#' It picks n times a random value within an interval for the predictor variable (intraspecific variability) and 
#' it repeats the analysis m times for m trees from a \code{multiPhylo} file. The trees are chosen randomly in the \code{multiPhylo} file.
#' The output gives a minimum, maximum, mean and sd of the model parameters.
#' @aliases variPgls
#' @param resp Numeric vector containing the response variable
#' @param pred Vector containing the predictor variable
#' @param se.pred Vector containing the standard error of \code{pred}  
#' @param tree A tree or list of tree of class \code{Phylo} or \code{multiPhylo}
#' @param ntree If TRUE or class(tree)=multiPhylo, the number of times to repeat the analysis with n different 
#' trees picked randomly in the multiPhylo file. If NULL, \code{ntree} = 1
#' @param npred If TRUE, the number of times to repeat the analysis generating a random value in the interval 
#' [\code{pred}-\code{se.pred},\code{pred}+\code{se.pred}]. If NULL, \code{npred} = 1
#' @param method A character string indicating which method to use to generate a random value in the interval 
#' [\code{pred} - \code{se.pred}, \code{pred} + \code{se.pred}]. Default is normal distribution: "normal" 
#' (function \code{\link{rnorm}}), uniform distribution is "uniform" (\code{\link{runif}})
#' @param taxa.col A character vector of taxa names that will be used to match data rows to phylogeny tips.
#' @param lambda A value for the lambda transformation. If NULL, \code{lambda}="ML"
#' Note that the model can be weighted by the sample size of each species, see \code{weights} in \code{\link{gls}}
#' @details This functions only works for simple linear regression \eqn{y = bx +
#' a}. Future implementation will deal with more complex models.
#' @return The function \code{variPgls} returns a list with the following
#' components:
#' @return \code{model_results} Model parameters for each iteration
#' parameters: intercept, beta, beta standard error, beta confidence interval, beta pvalue, AIC, lambda
#' @return \code{stats} Mean and sd of each parameter of the model
#' \code{min}, \code{max} and \code{mean} are the minimum, maximum and mean values for each parameter
#' \code{sd_all} is the total standard deviation (sd), \code{sd_tree} is the sd due to phylogenetic uncertainty,
#' \code{sd_pred} is the sd due to intraspecific variation
#' @section Warning: This code is not fully checked. Please be aware.
#' If \code{ntree} and \code{npred} are set to 1, the function computes a classic phylogenetic gls.
#' @seealso \code{\link{gls}}, \code{\link{corPagel}}, \code{\link{runif}}, \code{\link{rnorm}}
#' @examples
#' library(caper);library(phylolm);library(phytools)
#' 
#' # Example with a single tree
#' N <- 100 # Number of species in the tree
#' tree<-pbtree(n=N)
#' # Simulating response variable with phylogenetic signal
#' Ly <- rTrait(n=1, tree, model=c("lambda"),parameters=list(lambda=.8))
#' # Simulating explanatory variable and its standard error
#' Lx <- data.frame(xmean=Ly + rnorm(N,mean(Ly),1),xse=rep(0.5,100))
#' # Including Species names
#' sp <- tree$tip.label
#' # variPgls analysis
#' variPgls(resp=Ly,pred=Lx$xmean,se.pred=Lx$xse,tree,npred=2,method="normal",taxa.col=sp)
#' 
#' # Example with a set of trees of class multiphylo
#' N <- 100 # Number of species in the trees
#' multitree <- rmtree(N=50,n=100)
#' # Simulating response variable with phylogenetic signal
#' Ly <- rTrait(n=1, multitree[[1]], model=c("lambda"),parameters=list(lambda=.8))
#' # Simulating explanatory variable and its standard error
#' Lx <- data.frame(xmean=Ly + rnorm(N,mean(Ly),1),xse=rep(0.5,100))
#' # Including Species names
#' sp <- multitree[[1]]$tip.label               
#' # variPgls analysis
#' variPgls(resp=Ly,pred=Lx$xmean,se.pred=Lx$xse,tree=multitree,ntree=2,npred=2,method="normal",taxa.col=sp) 
#' @export

variPgls <- function(resp,pred,se.pred=NA,tree,ntree=1,npred=1,method=c("normal","uniform"),taxa.col,lambda="ML"){
  
  # Error check
  method <- match.arg(method)
  
  #Function to pick a random value in the interval
  if (method=="uniform") funr <- function(a, b) {runif(1,a-b,a+b)}
  else funr <- function(a, b) {rnorm(1,a,b)}
  
  # If the class of tree is multiphylo pick n=ntree random trees
  if(inherits(tree, "multiPhylo")){trees<-sample(length(tree),ntree,replace=F)}
  else {trees=1}
  
  #Create the results data.frame
  resultados<-data.frame("n.pred"=numeric(),"n.tree"=numeric(),"intercept"=numeric(),
                         "beta"=numeric(),"beta.se"=numeric(),"beta.ic"=numeric(),
                         "pval"=numeric(),"AIC"=numeric(),"Lambda"=numeric())

  #Assemble the dataframe
  if(!inherits(pred, c("numeric","integer")) || is.na(se.pred)) {data<-data.frame(taxa.col,resp,pred)}
  else {data<-data.frame(taxa.col,resp,pred,se.pred)}
  
  #Model calculation
  counter=1
  c.data<-list()
  for (i in 1:npred) {
    for (j in trees){
      
      #choose a random value in [mean-sd,mean+sd]
      if(!inherits(pred, c("numeric","integer")) || is.na(se.pred)) {data$variab<-data$pred}
      else {data$variab<-apply(data[,c("pred","se.pred")],1,function(x)funr(x["pred"],x["se.pred"]))}
      
      c.data[[i]]<-data.frame(data)
      
      #comparative data creation if tree is class=multiphylo
      if (inherits(tree, "multiPhylo")) {
        cor.0 <- ape::corPagel(1,phy=tree[[i]],fixed=F)
        try(mod.0 <- nlme::gls(resp~variab, data=c.data[[i]],method="ML",correlation=cor.0),TRUE)
      }
      
      else {
        cor.0 <- ape::corPagel(1,phy=tree,fixed=F)
        try(mod.0 <- nlme::gls(resp~variab, data=c.data,method="ML",correlation=cor.0),TRUE)
      }
        
      
      if(isTRUE(class(mod.0)=="try-error")) { next }
      else {
      
      
        #extract model coefficients
        sumMod <- as.data.frame(summary(mod.0)$tTable)
        df.0<-mod.0$dims$N-mod.0$dims$p
        intercept<-sumMod[1,1]
        beta<-sumMod[2,1]
        beta.se<-sumMod[2,2]
        pval<-sumMod[2,4]
        beta.IC <- qt(0.975,df.0)*beta.se

        AIC<-summary(mod.0)$AIC
        Lambda<-attr(mod.0$apVar,"Pars")["corStruct"]
        
        #write in a table
        resultados[counter,1]<- i
        resultados[counter,2]<- j    
        resultados[counter,3]<- intercept
        resultados[counter,4]<- beta
        resultados[counter,5]<- beta.se
        resultados[counter,6]<- beta.IC
        resultados[counter,7]<- pval
        resultados[counter,8]<- AIC
        resultados[counter,9]<- Lambda
      
        counter=counter+1
      }
    }  
  }
  
  #calculate mean and sd for each parameter
  #variation due to tree choice
  mean_by_tree<-aggregate(.~n.tree, data=resultados, mean)
  #variation due to continuous trait
  mean_by_randomval<-aggregate(.~n.pred, data=resultados, mean)
  
  statresults<-data.frame(min=apply(resultados,2,min),
                          max=apply(resultados,2,max),
                          mean=apply(resultados,2,mean),
                          sd_all=apply(resultados,2,sd),
                          sd_tree=apply(mean_by_tree,2,sd),
                          sd_pred=apply(mean_by_randomval,2,sd))[-(1:2),]
  
  
  output <- list(model_results=resultados,paste("Residual degrees of freedom:",df.0),
                 stats=statresults)
  
  return(output)  
}

