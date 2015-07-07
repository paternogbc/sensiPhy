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
#' @param vari.resp Vector containing the standard error, the standard deviation or the minimum and maximum (cbind(min,max)) of \code{resp}
#' When information is not available for one species, the value can be 0 or \code{NA}
#' @param vari.pred Vector containing the standard error, the standard deviation or the minimum and maximum (cbind(min,max)) of \code{pred}
#' When information is not available for one species, the value can be 0 or \code{NA}
#' @param tree A tree or list of tree of class \code{Phylo} or \code{multiPhylo}
#' @param ntree If TRUE or class(tree)=multiPhylo, the number of times to repeat the analysis with n different
#' trees picked randomly in the multiPhylo file. If NULL, \code{ntree} = 1
#' @param nvari If TRUE, the number of times to repeat the analysis generating a random value in the interval
#' [\code{pred}-\code{vari.pred},\code{pred}+\code{vari.pred}] or in the interval [min,max] (if vari.pred=(cbind(min,max)))
#' and a random value in the interval
#' [\code{resp}-\code{vari.resp},\code{resp}+\code{vari.resp}] or in the interval [min,max] (if vari.resp=(cbind(min,max))).
#' If NULL, \code{nvari} = 1
#' @param method A character string indicating which method to use to generate a random value in the interval
#' [\code{pred} - \code{vari.pred}, \code{pred} + \code{vari.pred}]. Default is uniform distribution: "uniform"
#' (function \code{\link{rnorm}}), normal distribution is "normal" (\code{\link{runif}})
#' Warning: normal distribution can be used oly if vari.pred is the standard deviation of the mean.
#' @param taxa.nam A character vector of taxa names that will be used to match data rows to phylogeny tips.
#' @param lambda A value for the lambda transformation. If NULL, \code{lambda}="ML"
#' Note that the model can be weighted by the sample size of each species, see \code{weights} in \code{\link{gls}}
#' @details This functions only works for simple linear regression \eqn{y = bx +a}.
#' Future implementation will deal with more complex models.
#' If you log-transform your predictor variable, make sure that your vari.pred is in a log-scale too or use ##### function
#' to build your input data.
#' @return The function \code{variPgls} returns a list with the following components:
#' @return \code{model_results} Model parameters for each iteration
#' parameters: intercept, intercept standard error and pvalue, beta, beta standard error and pvalue and confidence interval, AIC, lambda
#' @return \code{stats} Statistics for model parameters.
#' @return Residual degrees of freedom and number of models that converged properly
#' \code{min}, \code{max} and \code{mean} are the minimum, maximum and mean values for each parameter
#' \code{sd_all} is the total standard deviation (sd), \code{sd_tree} is the sd due to phylogenetic uncertainty,
#' \code{sd_pred} is the sd due to intraspecific variation
#' @section Warning: This code is not fully checked. Please be aware.
#' If \code{ntree} and \code{nvari} are set to 1, the function computes a classic phylogenetic gls.
#' @seealso \code{\link{gls}}, \code{\link{corPagel}}, \code{\link{runif}}, \code{\link{rnorm}}
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
#' variPgls(resp=Ly,pred=Lx$xmean,vari.pred=Lx$xse,tree,nvari=2,method="normal",taxa.nam=sp)
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
#' variPgls(resp=Ly,pred=Lx$xmean,vari.pred=Lx$xse,tree=multitree,ntree=2,nvari=2,method="normal",taxa.nam=sp)
#' }
#'  @export

variPgls <- function(resp,pred,vari.resp=NA,vari.pred=NA,taxa.nam,tree,ntree=1,nvari=1,method="normal",lambda="ML"){

  #Error check
  if (!inherits(tree, "phylo") & !inherits(tree, "multiPhylo"))
    stop("'", deparse(substitute(tree)), "' not of class 'phylo' or 'multiPhylo'")

  if(inherits(tree, "multiPhylo")){
  tree1<-tree[[1]]}
  else
  tree1<-tree

  #Error check for phylo
    if (!is.rooted(tree1)) {
    if (force.root) {
      tree$root.edge <- 1
    }
    else {
      stop("'", deparse(substitute(tree1)), "' is not rooted or has a basal polytomy.")
    }
  }
  if (any(duplicated(tree1$tip.label)))
    stop("Duplicate tip labels present in phylogeny")
  if (any(duplicated(c(tree1$tip.label, tree1$node.label))))
    stop("Labels duplicated between tips and nodes in phylogeny")


  #Match data and phylogeny in comparative.data style
  tiplabl<-tree1$tip.label
  spnames<-as.data.frame(taxa.nam)

  in.both <- intersect(taxa.nam, tiplabl)
  if (length(in.both) == 0)
    stop("No tips are common to the dataset and phylogeny")

  mismatch<-union(setdiff(tiplabl,taxa.nam),setdiff(taxa.nam,tiplabl))
  if (length(mismatch) != 0)   warning("Phylogeny tips not not match the species list,
  species were dropped from phylogeny or species list")

  #Drop species from tree
  tree<-lapply(tree,ape::drop.tip,tip=mismatch)
  class(tree)<-"multiPhylo"

  #Assemble the dataframe and drop species if needed
  data<-data.frame(taxa.nam,resp,pred)

      #Rename names vari.resp and vari.pred columns if 2 columns are provided (min and max case)
      if(exists("vari.resp") && !is.null(dim(vari.resp))){colnames(vari.resp)<-c("resp.min","resp.max")}
      if(exists("vari.pred") && !is.null(dim(vari.pred))){colnames(vari.pred)<-c("pred.min","pred.max")}

  if(exists("vari.resp")) {data<-cbind(data,vari.resp)}
  if(exists("vari.pred")){data<-cbind(data,vari.pred)}
  else
  data<-data[!taxa.nam %in% as.factor(mismatch), ]


  #NA's check
  if (sum(is.na(resp))!=0 || sum(is.na(pred))!=0)
    {data<-data[!is.na(data$resp),]
     data<-data[!is.na(data$pred),]
     warning("NA's in response or predictor, row with NA's were removed")}
  else

  #transform NA's in SE columns into zeros
  data[is.na(data)] <- 0

  if(inherits(tree, "multiPhylo")){
    tree1<-tree[[1]]}
  else
    tree1<-tree

  rownames(data)<-data$taxa.nam
  tip.order <- match(tree1$tip.label, rownames(data))
  if (any(is.na(tip.order)))
    stop("Problem with sorting data frame: mismatch between tip labels and data frame labels")
  data <- data[tip.order, , drop = FALSE]
  rownames(data) <- tree1$tip.label

  #Function to pick a random value in the interval
  if (method=="normal") funr <- function(a, b) {rnorm(1,a,b)
      warning("With method=normal vari.pred must be the standard deviation of pred")}
  else  funr <- function(a,b) {runif(1,a-b,a+b)}

  # If the class of tree is multiphylo pick n=ntree random trees
  if(inherits(tree, "multiPhylo")){trees<-sample(length(tree),ntree,replace=F)}
  else {trees=1}

  #Create the results data.frame
  resultados<-data.frame("n.pred"=numeric(),"n.tree"=numeric(),
                         "intercept"=numeric(),"interc.se"=numeric(),"interc.pval"=numeric(),
                         "beta"=numeric(),"beta.se"=numeric(),"beta.pval"=numeric(),
                         "beta.ic"=numeric(),"AIC"=numeric(),"Lambda"=numeric())

  #Model calculation
  counter=1
  c.data<-list()
  for (i in 1:nvari) {
    for (j in trees){
      tryCatch({

      ##Set response and predictor variables
        #vari.resp and vari.resp are not provided
        if(!inherits(resp, c("numeric","integer")) || !exists("vari.resp") || !exists("vari.resp")) {data$respV<-data$resp}

        #choose a random value in min/max if vari.resp is provided
        if(exists("vari.resp") && !is.null(dim(vari.resp)))
          {data$respV<-apply(data[,c("resp.min","resp.max")],1,function(x)runif(1,x[1],x[2]))}

        #choose a random value in [mean-se,mean+se] if vari.resp is provided
        if(exists("vari.resp") && !is.null(dim(vari.resp)))
          {data$respV<-apply(data[,c("resp","vari.resp")],1,function(x)funr(x[1],x[2]))}

        else
        #vari.pred and vari.pred are not provided
        if(!inherits(pred, c("numeric","integer")) || !exists("vari.pred") || !exists("vari.pred")) {data$predV<-data$pred}

        #choose a random value in min/max if vari.pred is provided
        if(exists("vari.pred") && !is.null(dim(vari.pred)))
          {data$predV<-apply(data[,c("pred.min","pred.max")],1,function(x)runif(1,x[1],x[2]))}

        #choose a random value in [mean-se,mean+se] if vari.pred is provided
        if(exists("vari.pred") && !is.null(dim(vari.pred)))
          {data$predV<-apply(data[,c("pred","vari.pred")],1,function(x)funr(x[1],x[2]))}

        else

      c.data[[i]]<-data.frame(data)

      #comparative data creation if tree is class=multiphylo
      if (inherits(tree, "multiPhylo")) {
        cor.0 <- ape::corPagel(1,phy=tree[[j]],fixed=F)
        mod.0 <- nlme::gls(respV~predV, data=c.data[[i]],method="ML",correlation=cor.0,na.action=na.omit) #control=nlme::glsControl(singular.ok=TRUE)
      }

      else {
        cor.0 <- ape::corPagel(1,phy=tree,fixed=F)
        mod.0 <- nlme::gls(respV~predV, data=c.data,method="ML",correlation=cor.0,na.action=na.omit)
      }


      #if(isTRUE(class(mod.0)=="try-error")) { next }
      #else {


        #extract model coefficients
        sumMod <- as.data.frame(summary(mod.0)$tTable)
        df.0<-mod.0$dims$N-mod.0$dims$p
        intercept<-sumMod[1,1]
        interc.SE<-sumMod[1,2]
        interc.Pval<-sumMod[1,4]
        beta<-sumMod[2,1]
        beta.se<-sumMod[2,2]
        beta.pval<-sumMod[2,4]
        beta.IC <- qt(0.975,df.0)*beta.se

        AIC<-summary(mod.0)$AIC
        Lambda<-attr(mod.0$apVar,"Pars")["corStruct"]

        #write in a table
        resultados[counter,1] <- i
        resultados[counter,2] <- j
        resultados[counter,3] <- intercept
        resultados[counter,4] <- interc.SE
        resultados[counter,5] <- interc.Pval
        resultados[counter,6] <- beta
        resultados[counter,7] <- beta.se
        resultados[counter,8] <- beta.pval
        resultados[counter,9] <- beta.IC
        resultados[counter,10]<- AIC
        resultados[counter,11]<- Lambda

        counter=counter+1
        print(paste("pred: ",i,", tree: ",j,sep=""))
        }, error=function(e){cat("Warning :",conditionMessage(e), "\n")})
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


  output <- list(output="variPgls",model_results=resultados,paste("Residual degrees of freedom:",df.0),
                 paste("Result based on ",nrow(resultados),"models that converged properly"),
                 stats=statresults)

  return(output)
}

