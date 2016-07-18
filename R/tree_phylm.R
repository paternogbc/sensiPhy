#' Phylogenetic uncertainty - Phylogenetic Linear Regression
#'
#' Performs Phylogenetic linear regression evaluating
#' uncertainty in trees topology.
#'
#' @param formula The model formula
#' @param data Data frame containing species traits with species as row names.
#' @param phy A phylogeny (class 'multiPhylo', see ?\code{ape}).
#' @param times Number of times to repeat the analysis with n different trees picked 
#' randomly in the multiPhylo file.
#' If NULL, \code{times} = 2
#' @param model The phylogenetic model to use (see Details). Default is \code{lambda}.
#' @param track Print a report tracking function progress (default = TRUE)
#' @param ... Further arguments to be passed to \code{phylolm}
#' @details
#' This function fits a phylogenetic linear regression model using \code{\link[phylolm]{phylolm}}
#' to n trees, randomly picked in a multiPhylo file.
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
#' @return The function \code{tree_phylm} returns a list with the following
#' components:
#' @return \code{formula}: The formula
#' @return \code{data}: Original full dataset
#' @return \code{model_results}: Coefficients, aic and the optimised
#' value of the phylogenetic parameter (e.g. \code{lambda}) for each regression with a 
#' different phylogenetic tree.
#' @return \code{N.obs}: Size of the dataset after matching it with tree tips and removing NA's.
#' @return \code{stats}: Main statistics for model parameters.\code{CI_low} and \code{CI_high} are the lower 
#' and upper limits of the 95% confidence interval.
#' @return \code{all.stats}: Complete statistics for model parameters. \code{sd_intra} is the standard deviation 
#' due to intraspecific variation. \code{CI_low} and \code{CI_high} are the lower and upper limits 
#' of the 95% confidence interval.
#' @author Caterina Penone & Pablo Ariel Martinez
#' @seealso \code{\link[phylolm]{phylolm}}, \code{\link{sensi_plot}}
#' @references 
#' Donoghue, M.J. & Ackerly, D.D. (1996). Phylogenetic Uncertainties and 
#' Sensitivity Analyses in Comparative Biology. Philosophical Transactions:
#'  Biological Sciences, pp. 1241-1249.
#'  
#' Ho, L. S. T. and Ane, C. 2014. "A linear-time algorithm for 
#' Gaussian and non-Gaussian trait evolution models". Systematic Biology 63(3):397-408.
#' @examples 
#'# Load data:
#'data(alien)
#'# This analysis needs a multiphylo file:
#'class(alien$phy)
#'alien$phy
#'# run PGLS accounting for phylogenetic uncertain:
#'tree <- tree_phylm(log(gestaLen) ~ log(adultMass), phy = alien$phy, 
#'data = alien$data, times = 30)
#'# To check summary results:
#'summary(tree)
#'# Visual diagnostics
#'sensi_plot(tree)
#'# You can specify which graph to print: 
#'sensi_plot(tree, graphs = 3)
#' @export

tree_phylm <- function(formula,data,phy,
                         times=2,model="lambda",track=TRUE,...){
  #Error check
  if(class(formula)!="formula") stop("formula must be class 'formula'")
  if(class(data)!="data.frame") stop("data must be class 'data.frame'")
  if(class(phy)!="multiPhylo") stop("phy must be class 'multiPhylo'")
  if(length(phy)<times) stop("'times' must be smaller (or equal) than the number of trees in the 'multiPhylo' object")
  else
  
  #Matching tree and phylogeny using utils.R
  datphy<-match_dataphy(formula,data,phy)
  full.data<-datphy[[1]]
  phy<-datphy[[2]]

  # If the class of tree is multiphylo pick n=times random trees
  trees<-sample(length(phy),times,replace=F)

  #Create the results data.frame
  tree.model.estimates<-data.frame("n.tree"=numeric(),"intercept"=numeric(),"se.intercept"=numeric(),
                         "pval.intercept"=numeric(),"slope"=numeric(),"se.slope"=numeric(),
                         "pval.slope"=numeric(),"aic"=numeric(),"optpar"=numeric())

  #Model calculation
  counter=1
  errors <- NULL
  c.data<-list()
  pb <- utils::txtProgressBar(min = 0, max = times, style = 1)
  for (j in trees){
      
      #phylolm model
      mod = try(phylolm::phylolm(formula, data=full.data,model=model,phy=phy[[j]]),FALSE)

      
      if(isTRUE(class(mod)=="try-error")) {
        error <- j
        names(error) <- rownames(c.data$full.data)[j]
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
        
        if(track==TRUE) utils::setTxtProgressBar(pb, counter)
        
        #write in a table
        estim.simu <- data.frame(j, intercept, se.intercept, pval.intercept,
                                 slope, se.slope, pval.slope, aic.mod, optpar,
                                 stringsAsFactors = F)
        tree.model.estimates[counter, ]  <- estim.simu
        counter=counter+1
        
      }
    }
  on.exit(close(pb))
  #calculate mean and sd for each parameter
  #variation due to tree choice
  mean_by_tree<-stats::aggregate(.~n.tree, data=tree.model.estimates, mean)

  statresults<-data.frame(min=apply(tree.model.estimates,2,min),
                          max=apply(tree.model.estimates,2,max),
                          mean=apply(tree.model.estimates,2,mean),
                          sd_tree=apply(mean_by_tree,2,stats::sd))[-1,]
  
  statresults$CI_low  <- statresults$mean - qt(0.975, df = times-1) * statresults$sd_tree / sqrt(times)
  statresults$CI_high <- statresults$mean + qt(0.975, df = times-1) * statresults$sd_tree / sqrt(times)
  
  res <- list(   call = match.call(),
                 formula=formula,
                 data=full.data,
                 model_results=tree.model.estimates,N.obs=n,
                 stats = round(statresults[c(1:6),c(3,5,6)],digits=3),
                 all.stats = statresults)
  class(res) <- "sensiTree"
  return(res)
}
