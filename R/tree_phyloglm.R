#' Phylogenetic uncertainty - Phylogenetic Logistic Regression
#'
#' Performs Phylogenetic logistic regression evaluating
#' uncertainty in trees topology.
#'
#' @param formula The model formula
#' @param data Data frame containing species traits with species as row names.
#' @param phy A phylogeny (class 'multiPhylo', see ?\code{ape}).
#' @param ntree Number of times to repeat the analysis with n different trees picked 
#' randomly in the multiPhylo file.
#' If NULL, \code{ntree} = 2
#' @param btol Bound on searching space. For details see \code{phyloglm}.
#' @param track Print a report tracking function progress (default = TRUE)
#' @param ... Further arguments to be passed to \code{phyloglm}
#' @details
#' This function fits a phylogenetic linear regression model using \code{\link[phylolm]{phyloglm}}
#' to n trees, randomly picked in a multiPhylo file.
#'
#' Currently, this function can only implement simple logistic models (i.e. \eqn{trait~
#' predictor}). In the future we will implement more complex models.
#'
#' Output can be visualised using \code{sensi_plot}.
#'
#' @return The function \code{tree_phyloglm} returns a list with the following
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
#' #which example should we put here?
#' @author Caterina Penone & Pablo Ariel Martinez
#' @seealso \code{\link{sensi_plot}}
#' @references Here still: reference to phylolm paper + our own?
#' @export

tree_phyloglm <- function(formula,data,phy,
                         ntree=2,btol=50,track=TRUE,...){
  
  #Error check
  if(class(formula)!="formula") stop("formula must be class 'formula'")
  if(class(data)!="data.frame") stop("data must be class 'data.frame'")
  if(class(phy)!="multiPhylo") stop("phy must be class 'multiPhylo'")
  else
    
    #Matching tree and phylogeny using utils.R
    datphy<-match_dataphy(formula,data,phy)
  full.data<-datphy[[1]]
  phy<-datphy[[2]]
  
  # If the class of tree is multiphylo pick n=ntree random trees
  trees<-sample(length(phy),ntree,replace=F)
  
  #Create the results data.frame
  tree.model.estimates<-data.frame("n.tree"=numeric(),"intercept"=numeric(),"se.intercept"=numeric(),
                                   "pval.intercept"=numeric(),"slope"=numeric(),"se.slope"=numeric(),
                                   "pval.slope"=numeric(),"aic"=numeric(),"optpar"=numeric())
  
  #Model calculation
  counter=1
  errors <- NULL
  c.data<-list()
  for (j in trees){
    
    #phyloglm model
    mod = try(phylolm::phyloglm(formula, data=full.data,phy=phy[[j]],method="logistic_MPLE",btol=btol),FALSE)

    if(isTRUE(class(mod)=="try-error")) {
      error <- j
      names(error) <- rownames(c.data$full.data)[j]
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
      
      if(track==TRUE) print(paste("tree: ",j,sep=""))
      
      #write in a table
      estim.simu <- data.frame(j, intercept, se.intercept, pval.intercept,
                               slope, se.slope, pval.slope, aic.mod, optpar,
                               stringsAsFactors = F)
      tree.model.estimates[counter, ]  <- estim.simu
      counter=counter+1
      
    }
  }
  
  #calculate mean and sd for each parameter
  #variation due to tree choice
  mean_by_tree<-stats::aggregate(.~n.tree, data=tree.model.estimates, mean)
  
  statresults<-data.frame(min=apply(tree.model.estimates,2,min),
                          max=apply(tree.model.estimates,2,max),
                          mean=apply(tree.model.estimates,2,mean),
                          sd_tree=apply(mean_by_tree,2,stats::sd))[-1,]
  
  
  res <- list(formula=formula,
              datas=full.data,
              model_results=tree.model.estimates,N.obs=n,
              stats=statresults)
  class(res) <- "sensiTree"
  return(res)
}
