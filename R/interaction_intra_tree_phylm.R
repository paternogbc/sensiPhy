#' Interaction of intraspecific variability & Phylogenetic uncertainty - Phylogenetic Linear Regression
#'
#' Performs Phylogenetic linear regression evaluating
#' intraspecific variability in response and/or predictor variables
#' and uncertainty in trees topology.
#'
#' @param formula The model formula: \code{response~predictor}. 
#' @param data Data frame containing species traits and species names as row names.
#' @param phy A phylogeny (class 'phylo', see ?\code{ape}).
#' @param Vy Name of the column containing the standard deviation or the standard error of the response 
#' variable. When information is not available for one taxon, the value can be 0 or \code{NA}.
#' @param Vx Name of the column containing the standard deviation or the standard error of the predictor 
#' variable. When information is not available for one taxon, the value can be 0 or \code{NA}
#' @param y.transf Transformation for the response variable (e.g. \code{log} or \code{sqrt}). Please use this 
#' argument instead of transforming data in the formula directly (see also details below).
#' @param x.transf Transformation for the predictor variable (e.g. \code{log} or \code{sqrt}). Please use this 
#' argument instead of transforming data in the formula directly (see also details below).
#' @param times.intra Number of times to repeat the analysis generating a random value for response and/or predictor variables.
#' If NULL, \code{times.intra} = 30
#' @param times.tree Number of times to repeat the analysis with n different trees picked 
#' randomly in the multiPhylo file.
#' If NULL, \code{times.tree} = 2
#' @param distrib A character string indicating which distribution to use to generate a random value for the response 
#' and/or predictor variables. Default is normal distribution: "normal" (function \code{\link{rnorm}}).
#' Uniform distribution: "uniform" (\code{\link{runif}})
#' Warning: we recommend to use normal distribution with Vx or Vy = standard deviation of the mean.
#' @param model The phylogenetic model to use (see Details). Default is \code{lambda}.
#' @param track Print a report tracking function progress (default = TRUE)
#' @param ... Further arguments to be passed to \code{phylolm}
#' @details
#' This function fits a phylogenetic linear regression model using \code{\link[phylolm]{phylolm}} to n trees (\code{times.tree}), 
#' randomly picked in a multiPhylo file. The regression is also repeated \code{times.intra} times.
#' At each iteration the function generates a random value for each row in the dataset using the standard deviation 
#' or errors supplied and assuming a normal or uniform distribution. To calculate means and se for your raw data, 
#' you can use the \code{summarySE} function from the package \code{Rmisc}.
#' 
#' #' All phylogenetic models from \code{phylolm} can be used, i.e. \code{BM},
#' \code{OUfixedRoot}, \code{OUrandomRoot}, \code{lambda}, \code{kappa},
#' \code{delta}, \code{EB} and \code{trend}. See ?\code{phylolm} for details.
#'
#' Currently, this function can only implement simple linear models (i.e. \eqn{trait~
#' predictor}). In the future we will implement more complex models.
#'
#' Output can be visualised using \code{sensi_plot}.
#' 
#' @section Warning:  
#' When Vy or Vx exceed Y or X, respectively, negative (or null) values can be generated, this might cause problems
#' for data transformation (e.g. log-transformation). In these cases, the function will skip the simulation. This problem can
#' be solved by increasing \code{times}, changing the transformation type and/or checking the target species in output$sp.pb.
#'  
#' @return The function \code{interaction_intra_tree_phylm} returns a list with the following
#' components:
#' @return \code{formula}: The formula
#' @return \code{data}: Original full dataset
#' @return \code{model_results}: Coefficients, aic and the optimised value of the phylogenetic 
#' parameter (e.g. \code{lambda}) for each regression using a value in the interval of variation and 
#' a different phylogenetic tree.
#' @return \code{N.obs}: Size of the dataset after matching it with tree tips and removing NA's.
#' @return \code{stats}: Main statistics for model parameters.\code{CI_low} and \code{CI_high} are the lower 
#' and upper limits of the 95% confidence interval.
#' @return \code{all.stats}: Complete statistics for model parameters.
#' Fields coded using \code{all} describe statistics due to both intraspecific variation and phylogenetic uncertainty.
#' Fields coded using \code{intra} describe statistics due to intraspecific variation only.
#' Fields coded using \code{tree} describe statistics due to phylogenetic uncertainty only.
#' \code{sd} is the standard deviation. \code{CI_low} and \code{CI_high} are the lower and upper limits 
#' of the 95% confidence interval.
#' @return \code{sp.pb}: Species that caused problems with data transformation (see details above).
#' 
#' @author Caterina Penone & Pablo Ariel Martinez
#' @seealso \code{\link[phylolm]{phylolm}}, \code{\link{sensi_plot}}
#' @references
#' Martinez, P. a., Zurano, J.P., Amado, T.F., Penone, C., Betancur-R, R., 
#' Bidau, C.J. & Jacobina, U.P. (2015). Chromosomal diversity in tropical reef 
#' fishes is related to body size and depth range. Molecular Phylogenetics and 
#' Evolution, 93, 1-4
#' 
#' Ho, L. S. T. and Ane, C. 2014. "A linear-time algorithm for 
#' Gaussian and non-Gaussian trait evolution models". Systematic Biology 63(3):397-408.
#' @examples 
#'# Load data:
#'data(alien)
#'# run PGLS accounting for intraspecific variation:
#'intra.tree <- interaction_intra_tree_phylm(gestaLen ~ adultMass, data = alien$data, phy = alien$phy,
#'Vy = "SD_gesta", times.intra = 10, times.tree = 10, y.transf = log, x.transf = log)
#'# To check summary results:
#'summary(intra.tree)
#'# Visual diagnostics
#'sensi_plot(intra.tree, uncer.type = "all") #or uncer.type = "tree", uncer.type = "intra"
#'
#' @export


interaction_intra_tree_phylm <- function(formula, data, phy,
                        Vy = NULL, Vx = NULL,
                        y.transf = NULL, x.transf = NULL,
                        times.intra = 10, times.tree = 2, 
                        distrib = "normal", model = "lambda", 
                        track = TRUE, ...){
  
  #Error check
  if(is.null(Vx) & is.null(Vy)) stop("Vx or Vy must be defined")
  if(class(formula) != "formula") stop("formula must be class 'formula'")
  if(class(data) != "data.frame") stop("data must be class 'data.frame'")
  if(class(phy)!="multiPhylo") stop("phy must be class 'multiPhylo'")
  if(formula[[2]]!=all.vars(formula)[1] || formula[[3]]!=all.vars(formula)[2])
    stop("Please use arguments y.transf or x.transf for data transformation")
  if(length(phy)<times.tree) stop("'times.tree' must be smaller (or equal) than the number of trees in the 'multiPhylo' object")
  if(distrib == "normal") warning ("distrib = normal: make sure that standard deviation is provided for Vx and/or Vy")


  #Matching tree and phylogeny using utils.R
  datphy <- match_dataphy(formula, data, phy)
  full.data <- datphy[[1]]
  phy <- datphy[[2]]
  
  resp <- all.vars(formula)[1]
  pred <- all.vars(formula)[2]
  
  if(!is.null(Vy) && sum(is.na(full.data[, Vy])) != 0) {
    full.data[is.na(full.data[, Vy]), Vy] <- 0}
  
  if(!is.null(Vx) && sum(is.na(full.data[, Vx])) != 0) {
    full.data[is.na(full.data[, Vx]), Vx] <- 0}
  

  #Function to pick a random value in the interval
  if (distrib == "normal") funr <- function(a,b) {stats::rnorm(1,a,b)}
  else  funr <- function(a,b) {stats::runif(1, a - b, a + b)}
  
  # If the class of tree is multiphylo pick n=times.tree random trees
  trees<-sample(length(phy),times.tree,replace=F)
  
  #Create the results data.frame
  tree.intra.model.estimates <- data.frame("n.tree"=numeric(), "n.intra"=numeric(),
                                          "intercept"=numeric(),"se.intercept"=numeric(),
                                          "pval.intercept"=numeric(),
                                          "estimate"=numeric(), "se.estimate"=numeric(),
                                          "pval.estimate"=numeric(),
                                          "aic"=numeric(), "optpar"=numeric())

  #Model calculation
  counter = 1
  errors <- NULL
  species.NA <- list()
  pb <- utils::txtProgressBar(min = 0, max = (times.tree*times.intra), style = 1)
  
  for (j in 1:times.tree) {
      for (i in 1:times.intra) {
    
        ##Set response and predictor variables
    #Vy is not provided or is not numeric, do not pick random value
    if(!inherits(full.data[,resp], c("numeric","integer")) || is.null(Vy)) 
    {full.data$respV <- stats::model.frame(formula, data = full.data)[,1]}
    
    #choose a random value in [mean-se,mean+se] if Vy is provided
    if (!is.null(Vy))
    {full.data$respV <- apply(full.data[,c(resp,Vy)],1,function(x)funr(x[1],x[2]))}
    
    #Vx is not provided or is not numeric, do not pick random value
    if (!inherits(full.data[,pred], c("numeric","integer")) || is.null(Vx))
    {full.data$predV <- stats::model.frame(formula, data = full.data)[,2]}
    
    #choose a random value in [mean-se,mean+se] if Vx is provided
    if(!is.null(Vx))
    {full.data$predV <- apply(full.data[,c(pred,Vx)],1,function(x)funr(x[1],x[2]))}
    
    #transform Vy and/or Vx if x.transf and/or y.transf are provided
    if(!is.null(y.transf)) 
    {suppressWarnings (full.data$respV <- y.transf(full.data$respV))}
    
    if(!is.null(x.transf)) 
    {suppressWarnings (full.data$predV <- x.transf(full.data$predV))}
    
    #skip iteration if there are NA's in the dataset
    species.NA[[i]]<-rownames(full.data[with(full.data,is.na(predV) | is.na(respV)),])
    if(sum(is.na(full.data[,c("respV","predV")])>0)) next
    
    #model
    mod = try(phylolm::phylolm(respV ~ predV, data = full.data, model = model,
                               phy = phy[[j]]), FALSE)
    
    if(isTRUE(class(mod) == "try-error")) {
      error <- i
      errors <- c(errors,error)
      next }
    
    
    else{
      intercept            <- phylolm::summary.phylolm(mod)$coefficients[[1,1]]
      se.intercept         <- phylolm::summary.phylolm(mod)$coefficients[[1,2]]
      estimate             <- phylolm::summary.phylolm(mod)$coefficients[[2,1]]
      se.estimate          <- phylolm::summary.phylolm(mod)$coefficients[[2,2]]
      pval.intercept       <- phylolm::summary.phylolm(mod)$coefficients[[1,4]]
      pval.estimate        <- phylolm::summary.phylolm(mod)$coefficients[[2,4]]
      aic.mod              <- mod$aic
      n                    <- mod$n
      
      if (model == "BM") {
        optpar <- NA
      }
      if (model != "BM") {
        optpar <- mod$optpar
      }
      
      #write in a table
      estim.simu <- data.frame(j, i, intercept, se.intercept, pval.intercept,
                               estimate, se.estimate, pval.estimate, aic.mod, optpar,
                               stringsAsFactors = F)
      tree.intra.model.estimates[counter, ]  <- estim.simu
      counter=counter+1
      
    }
      }
    if(track == TRUE) utils::setTxtProgressBar(pb, i*j)
        }
  on.exit(close(pb))
  
  #calculate mean and sd for each parameter
  #variation due to intraspecific variability
  mean_by_randomval <- stats::aggregate(.~n.intra, data = tree.intra.model.estimates,mean)
  
  #variation due to tree choice
  mean_by_tree<-stats::aggregate(.~n.tree, data=tree.intra.model.estimates, mean)
  
  statresults <- data.frame(min.all = apply(tree.intra.model.estimates, 2, min),
                            max.all = apply(tree.intra.model.estimates, 2, max),
                            mean.all = apply(tree.intra.model.estimates, 2, mean),
                            sd_all = apply(tree.intra.model.estimates, 2, stats::sd),
                            
                            min.intra = apply(mean_by_randomval, 2, min),
                            max.intra = apply(mean_by_randomval, 2, max),
                            mean.intra = apply(mean_by_randomval, 2, mean),
                            sd_intra = apply(mean_by_randomval, 2, stats::sd),
                            
                            min.tree = apply(mean_by_tree, 2, min),
                            max.tree = apply(mean_by_tree, 2, max),
                            mean.tree = apply(mean_by_tree, 2, mean),
                            sd_tree = apply(mean_by_tree, 2, stats::sd))[-(1:2), ]
  
  statresults$CI_low_all    <- statresults$mean.all - stats::qt(0.975, df = times.intra*times.tree-1) * statresults$sd_all / sqrt(times.intra*times.tree)
  statresults$CI_low_intra  <- statresults$mean.intra - stats::qt(0.975, df = times.intra-1) * statresults$sd_intra / sqrt(times.intra)
  statresults$CI_low_tree   <- statresults$mean.tree - stats::qt(0.975, df = times.tree-1) * statresults$sd_intra / sqrt(times.tree)
  
  statresults$CI_high_all    <- statresults$mean.all + stats::qt(0.975, df = times.intra*times.tree-1) * statresults$sd_all / sqrt(times.intra*times.tree)
  statresults$CI_high_intra  <- statresults$mean.intra + stats::qt(0.975, df = times.intra-1) * statresults$sd_intra / sqrt(times.intra)
  statresults$CI_high_tree   <- statresults$mean.tree + stats::qt(0.975, df = times.tree-1) * statresults$sd_intra / sqrt(times.tree)
  
  #reoder to later match sensi_plot for the single functions
  statresults <- statresults[,c("min.all","max.all","mean.all","sd_all","CI_low_all","CI_high_all",
                                "min.intra","max.intra","mean.intra","sd_intra","CI_low_intra","CI_high_intra",
                                "min.tree","max.tree","mean.tree","sd_tree","CI_low_tree","CI_high_tree")]
  

  #species with transformation problems
  nr <- times.tree*times.intra - nrow(tree.intra.model.estimates)
  sp.pb <- unique(unlist(species.NA))
  if (length(sp.pb) >0) 
    warning (paste("in", nr,"simulations, data transformations generated NAs, please consider using another function
                   for x.transf or y.transf and check output$sp.pb",sep=" "))
  
  
  res <- list(call = match.call(),
              formula = formula,
              y.transf = y.transf, 
              x.transf = x.transf,
              data = full.data,
              model_results = tree.intra.model.estimates, N.obs = n,
              stats = round(statresults[c(1:6),c(3,13,16,7,14,17,11,15,18)],digits=3),
              all.stats = statresults,sp.pb=sp.pb)
  class(res) <- "sensiIntra_Tree"
  return(res)
}

