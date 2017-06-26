#' Interaction of intraspecific variability & influential species - Phylogenetic Logistic Regression
#'
#' Performs leave-one-out deletion analyis for phylogenetic logistic regression,
#' and detects influential species, while taking into account potential
#' interactions with intraspecific variability.
#'
#' @param formula The model formula:  \code{response~predictor}. 
#' @param data Data frame containing species traits with row names matching tips
#' in \code{phy}.
#' @param phy A phylogeny (class 'phylo') matching \code{data}.
#' @param cutoff The cutoff value used to identify for influential species
#' (see Details)
#' @param n.intra Number of datasets resimulated taking into account intraspecific variation (see: \code{"intra_phylgm"})
#' @param Vx Name of the column containing the standard deviation or the standard error of the predictor 
#' variable. When information is not available for one taxon, the value can be 0 or \code{NA}
#' @param x.transf Transformation for the predictor variable (e.g. \code{"log"} or \code{"sqrt"}). Please use this 
#' argument instead of transforming data in the formula directly (see also details below).
#' @param distrib A character string indicating which distribution to use to generate a random value for the response 
#' and/or predictor variables. Default is normal distribution: "normal" (function \code{\link{rnorm}}).
#' Uniform distribution: "uniform" (\code{\link{runif}})
#' Warning: we recommend to use normal distribution with Vx or Vy = standard deviation of the mean.
#' @param btol Bound on searching space. For details see \code{phyloglm}
#' @param track Print a report tracking function progress (default = TRUE)
#' @param ... Further arguments to be passed to \code{phylolm}
#' @details
#' This function fits a phylogenetic linear regression model using \code{\link[phylolm]{phylolm}}, and detects
#' influential species by sequentially deleting one at a time. The regression is repeated \code{n.intra} times for 
#' simulated values of the dataset, taking into account intraspecific variation. At each iteration, the function 
#' generates a random value for each row in the dataset using the standard deviation or errors supplied, and 
#' detect the influential species within that iteration. 
#'
#' \code{influ_phylm} detects influential species based on the standardised
#' difference in intercept and/or slope when removing a given species compared
#' to the full model including all species. Species with a standardised difference
#' above the value of \code{cutoff} are identified as influential. The default
#' value for the cutoff is 2 standardised differences change.
#'
#' Currently, this function can only implement simple models (i.e. \eqn{trait~
#' predictor}). In the future we will implement more complex models.
#'
#' Output can be visualised using \code{sensi_plot}.
#' 
#' @section Warning:  
#' When Vx exceeds X negative (or null) values can be generated, this might cause problems
#' for data transformation (e.g. log-transformation). In these cases, the function will skip the simulation. This problem can
#' be solved by increasing \code{n.intra}, changing the transformation type and/or checking the target species in output$sp.pb.
#' 
#' Setting \code{n.intra} at high values can take a long time to exectue, since the total number of iterations equals \code{n.intra * nrow(data)}.
#' 
#' The function returns a list of \code{sensiInflu}-objects (the output of \code{influ_phylm} and \code{influ_phyglm}). 
#' of length \code{n.intra}. The user can use \code{summary} to evaluate this list. This will give, both for the 
#' regression slope and for the intercept, a table indicating how often across the \code{n.times} simulations a given
#' species was identified as the most influential species, as well as a table listing the mean slope, DFslope, 
#' Percentage change and P-value across all species that occured as most influential species in at least one simulation.
#' 
#' Additionally, users can evaluate each element in the list as a regula \code{sensiInflu}-object. 
#' 
#' @author Gustavo Paterno, Caterina Penone & Gijsbert D.A. Werner
#' @seealso \code{\link[phylolm]{phylolm}}, \code{\link{samp_phylm}},
#' \code{\link{influ_phylm}},\code{\link{intra_phylm}},\code{\link{sensi_plot}}.
#' @references Ho, L. S. T. and Ane, C. 2014. "A linear-time algorithm for 
#' Gaussian and non-Gaussian trait evolution models". Systematic Biology 63(3):397-408.
#' @examples 
#' # Load data:
#' data(alien)
#' # run analysis:
#' intra_influ <- interaction_intra_influ_phylm(formula = gestaLen ~ adultMass, phy = alien$phy[[1]],
#' data=alien$data,model="lambda",y.transf = "log",x.transf = NULL,Vy="SD_gesta",Vx=NULL,
#' times=3,distrib = "normal")
#' # To check summary results:
#'summary(intra_influ)
#'# Most influential speciesL
#'intra_influ$influential.species
#'# Visual diagnostics
#'sensi_plot(intra_influ)
#'# You can specify which graph and parameter ("slope" or "intercept") to print: 
#'sensi_plot(intra_influ, param = "slope", graphs = 2)
#'
#'set.seed(6987)
#'phy = rtree(100)
#'x = rTrait(n=1,phy=phy,parameters=list(ancestral.state=2,optimal.value=2,sigma2=1,alpha=1))
#'X = cbind(rep(1,100),x)
#'y = rbinTrait(n=1,phy=phy, beta=c(-1,0.5), alpha=.7 ,X=X)
#'z = rnorm(n = length(x),mean = mean(x),sd = 0.1*mean(x))
#'dat = data.frame(y, x, z)
#'# Run sensitivity analysis:
#'influ_test <- interaction_intra_influ_phyglm(formula = y ~ x, data = dat, phy = phy, Vx = "z", 
#'                                             times = 3,track = TRUE,distrib="normal",x.transf=NULL) 
#'# To check summary results and most influential species:
#'summary(influ_test)
#'# Visual diagnostics for clade removal:
#'sensi_plot(influ_test)
#' @export