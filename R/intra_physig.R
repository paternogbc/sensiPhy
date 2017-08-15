#' Intraspecific variability - Phylogenetic Linear Regression
#'
#' Performs Phylogenetic siganl estimates evaluating
#' trait intraspecific variability
#'
#' @param trait.col The name of a column in the provided data frame with trait 
#'  to be analyzed  (e.g. "Body_mass").#' @param data Data frame containing species traits and species names as row names.
#' @param phy A phylogeny (class 'phylo', see ?\code{ape}).
#' @param V Name of the column containing the standard deviation or the standard error of the trait 
#' variable. When information is not available for one taxon, the value can be 0 or \code{NA}.
#' @param n.intra Number of times to repeat the analysis generating a random trait value.
#' If NULL, \code{n.intra} = 30
#' @param distrib A character string indicating which distribution to use to generate a random value for the response 
#' and/or predictor variables. Default is normal distribution: "normal" (function \code{\link{rnorm}}).
#' Uniform distribution: "uniform" (\code{\link{runif}})
#' Warning: we recommend to use normal distribution with Vx or Vy = standard deviation of the mean.
#' @param track Print a report tracking function progress (default = TRUE)
#' @param ... Further arguments to be passed to \code{phytools::physig}
#' @details
#' This function estimates phylogenetic signal using \code{\link[phytools]{phylosig}}.
#' The analysis is repeated \code{n.intra} times. At each iteration the function generates a random value
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
#' @return The function \code{intra_phylm} returns a list with the following
#' components:
#' @return \code{Trait}: Column name of the trait analysed
#' @return \code{data}: Original full dataset
#' @return \code{intra.physig.estimates}: Run number, phylogenetic signal estimate 
#' (lambda or K) and the p-value for each run with a different simulated datset.
#' @return \code{N.obs}: Size of the dataset after matching it with tree tips and removing NA's.
#' @return \code{stats}: Main statistics for signal estimate\code{CI_low} and \code{CI_high} are the lower 
#' and upper limits of the 95% confidence interval.
#' @author Caterina Penone & Pablo Ariel Martinez & Gustavo Paterno
#' @seealso \code{\link[phylolm]{phylolm}}, \code{\link{sensi_plot}}
#' @references
#' Martinez, P. a., Zurano, J.P., Amado, T.F., Penone, C., Betancur-R, R., 
#' Bidau, C.J. & Jacobina, U.P. (2015). Chromosomal diversity in tropical reef 
#' fishes is related to body size and depth range. Molecular Phylogenetics and 
#' Evolution, 93, 1-4
#' 
#' Blomberg, S. P., T. Garland Jr., A. R. Ives (2003) 
#' Testing for phylogenetic signal in comparative data: 
#' Behavioral traits are more labile. Evolution, 57, 717-745.
#' 
#' Pagel, M. (1999) Inferring the historical patterns of biological evolution. 
#' Nature, 401, 877-884.
#' 
#' Kamilar, J. M., & Cooper, N. (2013). Phylogenetic signal in primate behaviour,
#'  ecology and life history. Philosophical Transactions of the Royal Society 
#'  B: Biological Sciences, 368: 20120341.
#'  
#' @examples 
#'data(alien)
#'# Run sensitivity analysis:
#'intra <- intra_physig(trait.col = "gestaLen", V = "SD_gesta" , data = alien.data, phy = alien.phy[[1]])
#'summary(intra)
#'sensi_plot(intra)
#'sensi_plot(intra, graphs = 1)
#'sensi_plot(intra, graphs = 2)
#' @export


intra_physig <- function(trait.col, data, phy,
                        V = NULL, n.intra = 100, distrib = "normal",
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
  if(track == TRUE) pb <- utils::txtProgressBar(min = 0, max = n.intra, style = 1)
  for (i in 1:n.intra) {

    #choose a random value in [mean-se,mean+se] if Vx is provided
    predV <- apply(full.data[,c(trait.col,V)],1,function(x)funr(x[1],x[2]))
    
    #model
    mod.s    <- phytools::phylosig(tree = phy, x = predV, method = method, test = TRUE, ...)
    estimate <- mod.s[[1]] 
    pval     <- mod.s$P
    
    if(track == TRUE) utils::setTxtProgressBar(pb, i)
    #write in a table
    estim.simu <- data.frame(i, estimate, pval)
    intra.physig.estimates[counter, ]  <- estim.simu
    counter = counter + 1
      
    }
  on.exit(close(pb))
  
  statresults <- data.frame(min = apply(intra.physig.estimates, 2, min),
                            max = apply(intra.physig.estimates, 2, max),
                            mean = apply(intra.physig.estimates, 2, mean),
                            sd_intra = apply(intra.physig.estimates, 2, stats::sd))[-1, ]
  
  statresults$CI_low  <- statresults$mean - stats::qt(0.975, df = n.intra-1) * statresults$sd_intra / sqrt(n.intra)
  statresults$CI_high <- statresults$mean + stats::qt(0.975, df = n.intra-1) * statresults$sd_intra / sqrt(n.intra)
  
  stats <- round(statresults[c(1:2),c(3,5,6,1,2)],digits=5)
  
  cl <- match.call()
  res <- list(   call = cl,
                 Trait = trait.col,
                 intra.physig.estimates = intra.physig.estimates,
                 N.obs = N,
                 stats = stats,
                 data = full.data)
  class(res) <- "intra.physig"
  return(res)
}

