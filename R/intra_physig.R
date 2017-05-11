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
#' @param times Number of times to repeat the analysis generating a random value for response and/or predictor variables.
#' If NULL, \code{times} = 2
#' @param distrib A character string indicating which distribution to use to generate a random value for the response 
#' and/or predictor variables. Default is normal distribution: "normal" (function \code{\link{rnorm}}).
#' Uniform distribution: "uniform" (\code{\link{runif}})
#' Warning: we recommend to use normal distribution with Vx or Vy = standard deviation of the mean.
#' @param model The phylogenetic model to use (see Details). Default is \code{lambda}.
#' @param track Print a report tracking function progress (default = TRUE)
#' @param ... Further arguments to be passed to \code{phylolm}
#' @details
#' This function fits a phylogenetic linear regression model using \code{\link[phylolm]{phylolm}}.
#' The regression is repeated \code{times} times. At each iteration the function generates a random value
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
#' be solved by increasing \code{times}, changing the transformation type and/or checking the target species in output$sp.pb.
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
#'phy = alien$phy[[1]], data = alien$data, Vy = "SD_gesta", times = 30)
#'# To check summary results:
#'summary(intra)
#'# Visual diagnostics
#'sensi_plot(intra)
#'
#' @export


intra_physig <- function(trait.col, data, phy,
                        V = NULL, times = 30, distrib = "normal",
                        method = "K", track = TRUE, ...){
  #Error check
  if(is.null(V)) stop("V must be defined")
  if(class(data) != "data.frame") stop("data must be class 'data.frame'")
  if(distrib == "normal") message (paste("distrib = normal: make sure that standard deviation", 
                                   "is provided for V", sep = " "))
  if(!inherits(phy,"phylo"))
    stop("tree should be an object of class \"phylo\".")
  if(!inherits(trait.col,"character")) 
    stop("trait.col should be an object of class \"character\".")
  
  # Check match between data and phy 
  datphy <- match_dataphy(get(trait.col) ~ 1, data, phy)
  full.data <- datphy$data
  phy <- datphy$phy
  trait     <- full.data[[trait.col]]
  names(trait)  <- phy$tip.label
  N <- nrow(full.data)
  
  if(!is.null(V) && sum(is.na(full.data[, V])) != 0) {
    full.data[is.na(full.data[, V]), V] <- 0}
  
  #Function to pick a random value in the interval
  if (distrib == "normal") funr <- function(a,b) {stats::rnorm(1,a,b)}
  else  funr <- function(a,b) {stats::runif(1, a - b, a + b)}
  
  #Create the results data.frame
  intra.physig.estimates <- data.frame("n.intra" = numeric(),"estimate" = numeric(),
                                      "pval" = numeric())
  counter = 1
  pb <- utils::txtProgressBar(min = 0, max = times, style = 1)
  for (i in 1:times) {

    #choose a random value in [mean-se,mean+se] if Vx is provided
    predV <- apply(full.data[,c(trait.col,V)],1,function(x)funr(x[1],x[2]))
    
    #model
    mod.s    <- phytools::phylosig(tree = phy, x = predV, method = "K", test = TRUE)
    estimate <- mod.s[[1]] 
    pval     <- mod.s$P
    
    if(track == TRUE) utils::setTxtProgressBar(pb, i)
    #write in a table
    estim.simu <- data.frame(i, estimate, pval)
    intra.physig.estimates[counter, ]  <- estim.simu
    counter = counter + 1
      
    }
  on.exit(close(pb))
  
  #calculate mean and sd for each parameter
  #variation due to intraspecific variability
  mean_by_randomval <- stats::aggregate(.~n.intra, data = intra.physig.estimates,
                                        mean)
  
  statresults <- data.frame(min = apply(intra.physig.estimates, 2, min),
                            max = apply(intra.physig.estimates, 2, max),
                            mean = apply(intra.physig.estimates, 2, mean),
                            sd_intra = apply(intra.physig.estimates, 2, stats::sd))[-1, ]
  
  statresults$CI_low  <- statresults$mean - stats::qt(0.975, df = times-1) * statresults$sd_intra / sqrt(times)
  statresults$CI_high <- statresults$mean + stats::qt(0.975, df = times-1) * statresults$sd_intra / sqrt(times)
  
  stats <- round(statresults[c(1:2),c(3,5,6,1,2)],digits=5)
  
  cl <- match.call()
  res <- list(   call = cl,
                 Trait = trait.col,
                 physig_results = intra.physig.estimates,
                 N.obs = N,
                 stats = stats,
                 data = full.data)
  class(res) <- "intra.physig"
  return(res)
}

