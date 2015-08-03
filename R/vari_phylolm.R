#' Phylogenetic linear regression accounting for intraspecific variability and phylogenetic uncertainty
#'
#' \code{vari_phylolm} performs Phylogenetic linear regression analysis accounting 
#' for intraspecific variability and phylogenetic
#' uncertainty in trees topology.
#' It picks n times a random value within an interval for the predictor or response
#'  variable (intraspecific variability) and
#' it repeats the analysis m times for m trees from a \code{multiPhylo} file. 
#' The trees are chosen randomly in the \code{multiPhylo} file.
#' The output gives a minimum, maximum, mean and sd of the model parameters.
#' @aliases vari_phylolm
#' @param resp Numeric vector containing the response variable
#' @param pred Vector containing the predictor variable
#' @param vari.resp Vector containing the standard error, the standard deviation 
#' or the minimum and maximum (cbind(min,max)) of \code{resp}
#' When information is not available for one species, the value can be 0 or \code{NA}
#' @param vari.pred Vector containing the standard error, the standard deviation 
#' or the minimum and maximum (cbind(min,max)) of \code{pred}
#' When information is not available for one species, the value can be 0 or \code{NA}
#' @param tree A tree or list of tree of class \code{Phylo} or \code{multiPhylo}
#' @param ntree If TRUE or class(tree)=multiPhylo, the number of times to repeat
#'  the analysis with n different trees picked randomly in the multiPhylo file. 
#'  If NULL, \code{ntree} = 1
#' @param nintra If TRUE, the number of times to repeat the analysis generating 
#' a random value in the interval [\code{pred}-\code{vari.pred},
#' \code{pred}+\code{vari.pred}] or in the interval [min,max] 
#' (if vari.pred=(cbind(min,max)))
#' and a random value in the interval
#' [\code{resp}-\code{vari.resp},\code{resp}+\code{vari.resp}] or in the interval 
#' [min,max] (if vari.resp=(cbind(min,max))).
#' If NULL, \code{nintra} = 1
#' @param method A character string indicating which method to use to generate a
#'  random value in the interval
#' [\code{pred} - \code{vari.pred}, \code{pred} + \code{vari.pred}]. Default is 
#' uniform distribution: "uniform"
#' (function \code{\link{rnorm}}), normal distribution is "normal" (\code{\link{runif}})
#' Warning: normal distribution can be used oly if vari.pred is the standard 
#' deviation of the mean.
#' @param taxa.nam A character vector of taxa names that will be used to match 
#' data rows to phylogeny tips.
#' @param lambda A value for the lambda transformation. If NULL, \code{lambda}="ML"
#' Note that the model can be weighted by the sample size of each species, see
#'  \code{weights} in \code{\link{gls}}
#' @inheritParams influ_phylolm
#' @details This functions only works for simple linear regression \eqn{y = bx +a}.
#' Future implementation will deal with more complex models.
#' If you log-transform your predictor variable, make sure that your vari.pred 
#' is in a log-scale too
#' to build your input data.
#' @return The function \code{variPgls} returns a list with the following components:
#' @return \code{model_results} Model parameters for each iteration
#' parameters: intercept, intercept standard error and pvalue, beta, beta
#'  standard error and pvalue and confidence interval, AIC, lambda
#' @return \code{stats} Statistics for model parameters.
#' @return Residual degrees of freedom and number of models that converged properly
#' \code{min}, \code{max} and \code{mean} are the minimum, maximum and mean values
#'  for each parameter
#' \code{sd_all} is the total standard deviation (sd), \code{sd_tree} is the sd 
#' due to phylogenetic uncertainty,
#' \code{sd_pred} is the sd due to intraspecific variation
#' @section Warning: This code is not fully checked. Please be aware.
#' If \code{ntree} and \code{nintra} are set to 1, the function computes a 
#' classic phylogenetic gls.
#' @seealso \code{\link{gls}}, \code{\link{corPagel}}, \code{\link{runif}},
#'  \code{\link{rnorm}}
#' @examples
#' \dontrun{
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
#' variPgls(resp=Ly,pred=Lx$xmean,vari.pred=Lx$xse,tree,nintra=2,
#' method="normal",taxa.nam=sp)
#'
#' # Example with a set of trees of class multiphylo
#' N <- 100 # Number of species in the trees
#' multitree <- rmtree(N=50,n=100)
#' # Simulating response variable with phylogenetic signal
#' Ly <- rTrait(n=1, multitree[[1]], model=c("lambda"),
#' parameters=list(lambda=.8))
#' # Simulating explanatory variable and its standard error
#' Lx <- data.frame(xmean=Ly + rnorm(N,mean(Ly),1),xse=rep(0.5,100))
#' # Including Species names
#' sp <- multitree[[1]]$tip.label
#' # variPgls analysis
#' variPgls(resp=Ly,pred=Lx$xmean,vari.pred=Lx$xse,tree=multitree,ntree=2,
#' nintra=2,method="normal",taxa.nam=sp)
#' }
#'  @export

vari_phylolm <- function(resp,pred,vari.resp=NA,vari.pred=NA,taxa.nam,tree,
                         nintra=1,ntree=1,model="lambda",method="normal",cutoff=2,track=TRUE){
  
  #Error check
  if (!inherits(tree, "phylo") & !inherits(tree, "multiPhylo"))
    stop("'", deparse(substitute(tree)), "' not of class 'phylo' or 'multiPhylo'")
  
  #Remove NA's before matching data and tips
  data<-data.frame(taxa.nam,resp,pred)
  
  #Rename names vari.resp and vari.pred columns if 2 columns are provided (min and max case)
  if(exists("vari.resp") && !is.null(dim(vari.resp))){colnames(vari.resp)<-c("resp.min","resp.max")}
  if(exists("vari.pred") && !is.null(dim(vari.pred))){colnames(vari.pred)<-c("pred.min","pred.max")}
  if(exists("vari.resp")) {data<-cbind(data,vari.resp)}
  if(exists("vari.pred")) {data<-cbind(data,vari.pred)}
  else

  if (sum(is.na(resp))!=0 || sum(is.na(pred))!=0)
  {data<-data[!is.na(data$resp),]
   data<-data[!is.na(data$pred),]
   warning("NA's in response or predictor, rows with NA's were removed")}

  #Match data and phylogeny in comparative.data style
  if(inherits(tree, "multiPhylo")){  
    tree1<-tree[[1]]}
  else
    tree1<-tree
  
  tiplabl<-tree1$tip.label
  taxa.nam<-as.character(data$taxa.nam)
  
  in.both <- intersect(taxa.nam, tiplabl)
  
  if (length(in.both) == 0)
    stop("No tips are common to the dataset and phylogeny")
  
  mismatch<-union(setdiff(tiplabl,taxa.nam),setdiff(taxa.nam,tiplabl))
  if (length(mismatch) != 0)   warning("Phylogeny tips do not match the species list,
                                       species were dropped from phylogeny or species list")
  
  #Drop species from tree
  tree<-lapply(tree,ape::drop.tip,tip=mismatch)
  class(tree)<-"multiPhylo"
  
  #transform NA's in SE columns into zeros
  data[is.na(data)] <- 0
  
  #Reorder rows according to tip labels
  if(inherits(tree, "multiPhylo")){  
    tree1<-tree[[1]]}
  else
    tree1<-tree
  
  rownames(data)<-data$taxa.nam
  tip.order <- match(tree1$tip.label, rownames(data))
  if (any(is.na(tip.order)))
    stop("Problem with sorting data frame: mismatch between tip labels and data frame labels")
  data <- data[tip.order, , drop = FALSE]
  rownames(data) <- data$taxa.nam
  
    
  #Function to pick a random value in the interval
  if (method=="normal") funr <- function(a, b) {rnorm(1,a,b)}
  else  funr <- function(a,b) {runif(1,a-b,a+b)}
  
  # If the class of tree is multiphylo pick n=ntree random trees
  if(inherits(tree, "multiPhylo")){trees<-sample(length(tree),ntree,replace=F)}
  else {trees=1}

  #Create the results data.frame
  resultados<-data.frame("n.intra"=numeric(),"n.tree"=numeric(),
                         "intercept"=numeric(),"se.intercept"=numeric(),"pval.intercept"=numeric(),
                         "slope"=numeric(),"se.slope"=numeric(),"pval.slope"=numeric(),
                         "aic"=numeric(),"optpar"=numeric())
  
  #Model calculation
  counter=1
  errors <- NULL
  c.data<-list()
  for (i in 1:nintra) {
    for (j in trees){
        
        ##Set response and predictor variables
        #vari.resp is not provided
        if(!inherits(resp, c("numeric","integer")) || !exists("vari.resp")) {data$respV<-data$resp}
        
        #choose a random value in min/max if vari.resp is provided
        if(exists("vari.resp") && !is.null(dim(vari.resp)))
        {data$respV<-apply(data[,c("resp.min","resp.max")],1,function(x)runif(1,x[1],x[2]))}
        
        #choose a random value in [mean-se,mean+se] if vari.resp is provided
        if(exists("vari.resp") && is.null(dim(vari.resp)))
        {data$respV<-apply(data[,c("resp","vari.resp")],1,function(x)funr(x[1],x[2]))}
          
        #vari.pred is not provided
        if(!inherits(pred, c("numeric","integer")) || !exists("vari.pred")) {data$predV<-data$pred}
        
        #choose a random value in min/max if vari.pred is provided
        if(exists("vari.pred") && !is.null(dim(vari.pred)))
        {data$predV<-apply(data[,c("pred.min","pred.max")],1,function(x)runif(1,x[1],x[2]))}
        
        #choose a random value in [mean-se,mean+se] if vari.pred is provided
        if(exists("vari.pred") && is.null(dim(vari.pred)))
        {data$predV<-apply(data[,c("pred","vari.pred")],1,function(x)funr(x[1],x[2]))}
        else
  
        row.names(data)<-data$taxa.nam
        
        #comparative data creation if tree is class=multiphylo
        if (inherits(tree, "multiPhylo")) {
          mod <- try(phylolm::phylolm(respV~predV, data=data, model=model,phy=tree[[j]]),TRUE)
         }
        
        else {
          mod <- try(phylolm::phylolm(respV~predV, data=data, model=model,phy=tree),TRUE)
        }
        
        
        if(isTRUE(class(mod)=="try-error")) {
          error <- i
          names(error) <- rownames(c.data$data)[i]
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
          d                    <- mod$d
          if (model == "BM"){
            optpar <- NA
          }
          if (model != "BM"){
            optpar               <- mod$optpar
          }
          
          if(track==TRUE) print(paste("intra: ",i,", tree: ",j,sep=""))
        
        #write in a table
        resultados[counter,1] <- i
        resultados[counter,2] <- j
        resultados[counter,3] <- intercept
        resultados[counter,4] <- se.intercept
        resultados[counter,5] <- pval.intercept
        resultados[counter,6] <- slope
        resultados[counter,7] <- se.slope
        resultados[counter,8] <- pval.slope
        resultados[counter,9] <- aic.mod
        resultados[counter,10]<- optpar
        
        
        counter=counter+1
        
    }
  }
  }
  
  #calculate mean and sd for each parameter
  #variation due to tree choice
  mean_by_tree<-aggregate(.~n.tree, data=resultados, mean)
  #variation due to continuous trait
  mean_by_randomval<-aggregate(.~n.intra, data=resultados, mean)
  
  statresults<-data.frame(min=apply(resultados,2,min),
                          max=apply(resultados,2,max),
                          mean=apply(resultados,2,mean),
                          sd_all=apply(resultados,2,sd),
                          sd_tree=apply(mean_by_tree,2,sd),
                          sd_intra=apply(mean_by_randomval,2,sd))[-(1:2),]
  
  
  output <- list(output="variPgls",model_results=resultados,paste("Number of observations:",n),
                 paste("Number of dependent variables:",d),
                 paste("Result based on ",nrow(resultados),"Models that converged properly"),
                 stats=statresults)
  
  return(output)
} 

