#' Intraspecific variability - Phylogenetic Linear Regression
#'
#' Performs Phylogenetic linear regression evaluating
#' intraspecific variability in response and/or predictor variables.
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
#' @param n.intra Number of times to repeat the analysis generating a random value for response and/or predictor variables.
#' If NULL, \code{n.intra} = 30
#' @param distrib A character string indicating which distribution to use to generate a random value for the response 
#' and/or predictor variables. Default is normal distribution: "normal" (function \code{\link{rnorm}}).
#' Uniform distribution: "uniform" (\code{\link{runif}})
#' Warning: we recommend to use normal distribution with Vx or Vy = standard deviation of the mean.
#' @param model The phylogenetic model to use (see Details). Default is \code{lambda}.
#' @param track Print a report tracking function progress (default = TRUE)
#' @param ... Further arguments to be passed to \code{phylolm}
#' @details
#' This function fits a phylogenetic linear regression model using \code{\link[phylolm]{phylolm}}.
#' The regression is repeated \code{n.intra} times. At each iteration the function generates a random value
#' for each row in the dataset using the standard deviation or errors supplied and assuming a normal or uniform distribution.
#' To calculate means and se for your raw data, you can use the \code{summarySE} function from the 
#' package \code{Rmisc}.
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
#' be solved by increasing \code{n.intra}, changing the transformation type and/or checking the target species in output$sp.pb.
#'  
#' @return The function \code{intra_phylm} returns a list with the following
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
#'intra <- intra_phylm(gestaLen ~ adultMass, y.transf = log, x.transf = log, 
#'phy = alien$phy[[1]], data = alien$data, Vy = "SD_gesta", n.intra = 30)
#'# To check summary results:
#'summary(intra)
#'# Visual diagnostics
#'sensi_plot(intra)
#'
#' @export


intra_phylm <- function(formula, data, phy,
                        Vy = NULL, Vx = NULL,
                        y.transf = NULL, x.transf = NULL,
                        n.intra = 30, distrib = "normal",
                        model = "lambda", track = TRUE, ...){
  #Error check
  if(is.null(Vx) & is.null(Vy)) stop("Vx or Vy must be defined")
  if(class(formula) != "formula") stop("formula must be class 'formula'")
  if(class(data) != "data.frame") stop("data must be class 'data.frame'")
  if(class(phy) != "phylo") stop("phy must be class 'phylo'")
    if(formula[[2]]!=all.vars(formula)[1] || formula[[3]]!=all.vars(formula)[2])
    stop("Please use arguments y.transf or x.transf for data transformation")
  if(distrib == "normal") warning ("distrib=normal: make sure that standard deviation 
                                 is provided for Vx and/or Vy")
  

    
  #Matching tree and phylogeny using utils.R
  datphy <- match_dataphy(formula, data, phy, ...)
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
  
  #Create the results data.frame
  intra.model.estimates <- data.frame("n.intra" = numeric(),"intercept" = numeric(),
                                    "se.intercept" = numeric(), 
                                    "pval.intercept" = numeric(),
                                    "estimate" = numeric(),
                                    "se.estimate" = numeric(),
                                    "pval.estimate" = numeric(),"aic" = numeric(),
                                    "optpar" = numeric())
  #Model calculation
  counter = 1
  errors <- NULL
  species.NA <- list()
  if(track == TRUE) pb <- utils::txtProgressBar(min = 0, max = n.intra, style = 1)
  for (i in 1:n.intra) {
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
                               phy = phy), FALSE)
    
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
      if(track == TRUE) utils::setTxtProgressBar(pb, i)
      #write in a table
      estim.simu <- data.frame(i, intercept, se.intercept, pval.intercept,
                               estimate, se.estimate, pval.estimate, aic.mod, optpar,
                               stringsAsFactors = F)
      intra.model.estimates[counter, ]  <- estim.simu
      counter=counter+1
      
    }
  }
  if(track == TRUE) on.exit(close(pb))
  
  #calculate mean and sd for each parameter
  #variation due to intraspecific variability
  mean_by_randomval <- stats::aggregate(.~n.intra, data = intra.model.estimates,
                                        mean)
  
  statresults <- data.frame(min = apply(intra.model.estimates, 2, min),
                          max = apply(intra.model.estimates, 2, max),
                          mean = apply(intra.model.estimates, 2, mean),
                          sd_intra = apply(mean_by_randomval, 2, stats::sd))[-1, ]
  
  statresults$CI_low  <- statresults$mean - stats::qt(0.975, df = n.intra-1) * statresults$sd_intra / sqrt(n.intra)
  statresults$CI_high <- statresults$mean + stats::qt(0.975, df = n.intra-1) * statresults$sd_intra / sqrt(n.intra)
  
  #species with transformation problems
  nr <- n.intra - nrow(intra.model.estimates)
  sp.pb <- unique(unlist(species.NA))
  if (length(sp.pb) >0) 
  warning (paste("in", nr,"simulations, data transformations generated NAs, please consider using another function
  for x.transf or y.transf and check output$sp.pb",sep=" "))
  
                                       
  res <- list(call = match.call(),
              formula = formula,
              y.transf = y.transf, 
              x.transf = x.transf,
              data = full.data,
              model_results = intra.model.estimates, N.obs = n,
              stats = round(statresults[c(1:6),c(3,5,6)],digits=3),
              all.stats = statresults,sp.pb=sp.pb)
  class(res) <- "sensiIntra"
  return(res)
}

