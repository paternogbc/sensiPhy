#' Sensitivity Analysis Species Sampling  - Phylogenetic Logistic Regression
#'
#' Performs analyses of sensitivity to species sampling by randomly removing
#' species and detecting the effects on parameter estimates in phylogenetic
#' logistic regression.
#'
#' @param formula The model formula
#' @param data Data frame containing species traits with row names matching tips
#' in \code{phy}.
#' @param phy A phylogeny (class 'phylo') matching \code{data}.
#' @param n.sim The number of times species are randomly deleted for each
#' \code{break}.
#' @param breaks A vector containing the percentages of species to remove.
#' @param btol Bound on searching space. For details see \code{phyloglm}
#' @param track Print a report tracking function progress (default = TRUE)
#' @param ... Further arguments to be passed to \code{phyloglm}
#' @details
#'
#' This function randomly removes a given percentage of species (controlled by
#' \code{breaks}) from the full phylogenetic logistic regression, fits a phylogenetic
#' logistic regression model without these species using \code{\link[phylolm]{phyloglm}},
#' repeats this many times (controlled by \code{n.sim}), stores the results and
#' calculates the effects on model parameters.
#'
#' Only logistic regression using the "logistic_MPLE"-method from
#' \code{phyloglm} is implemented.
#'
#' Currently, this function can only implement simple logistic models (i.e. \eqn{trait~
#' predictor}). In the future we will implement more complex models.
#'
#' Output can be visualised using \code{sensi_plot}.
#' @return The function \code{samp_phylm} returns a list with the following
#' components:
#' @return \code{formula}: The formula
#' @return \code{full.model.estimates}: Coefficients, aic and the optimised
#' value of the phylogenetic parameter (e.g. \code{lambda} or \code{kappa}) for
#' the full model without deleted species.
#' @return \code{sensi.estimates}: A data frame with all simulation
#' estimates. Each row represents a model rerun with a given number of species
#' \code{n.remov} removed, representing \code{n.percent} of the full dataset.
#' Columns report the calculated regression intercept (\code{intercept}),
#' difference between simulation intercept and full model intercept (\code{DIFintercept}),
#' the percentage of change in intercept compared to the full model (\code{intercept.perc})
#' and intercept p-value (\code{pval.intercept}). All these parameters are also reported
#' for the regression slope (\code{DIFestimate} etc.). Additionally, model aic value
#' (\code{AIC}) and the optimised value (\code{optpar}) of the phylogenetic
#' parameter (e.g. \code{kappa} or \code{lambda}, depending on the phylogenetic model
#' used) are reported. Lastly we reported the standardised difference in intercept 
#' (\code{sDIFintercept}) and slope (\code{sDIFestimate}). 
#' @return \code{sign.analysis} For each break (i.e. each percentage of species
#' removed) this reports the percentage of statistically signficant (at p<0.05)
#' intercepts (\code{perc.sign.intercept}) over all repititions as well as the
#' percentage of statisticaly significant (at p<0.05) slopes (\code{perc.sign.estimate}).
#' @return \code{data}: Original full dataset.
#' @note Please be aware that dropping species may reduce power to detect 
#' significant slopes/intercepts and may partially be responsible for a potential 
#' effect of species removal on p-values. Please also consult standardised differences
#' in the (summary) output. 
#' @author Gustavo Paterno & Gijsbert D.A. Werner
#' @seealso \code{\link[phylolm]{phyloglm}}, \code{\link{samp_phylm}},
#' \code{\link{influ_phyglm}}, \code{\link{sensi_plot}}
#' @references 
#' Werner, G.D.A., Cornwell, W.K., Sprent, J.I., Kattge, J. & Kiers, E.T. (2014).
#'  A single evolutionary innovation drives the deep evolution of symbiotic N2-fixation
#'   in angiosperms. Nature Communications, 5, 4087.
#'   
#' #' Ho, L. S. T. and Ane, C. 2014. "A linear-time algorithm for 
#' Gaussian and non-Gaussian trait evolution models". Systematic Biology 63(3):397-408.
#' @examples
#'# Simulate Data:
#'set.seed(6987)
#'phy = rtree(100)
#'x = rTrait(n=1,phy=phy)
#'X = cbind(rep(1,100),x)
#'y = rbinTrait(n=1,phy=phy, beta=c(-1,0.5), alpha=.7 ,X=X)
#'dat = data.frame(y, x)
#'# Run sensitivity analysis:
#'samp <- samp_phyglm(y ~ x, data = dat, phy = phy, n.sim = 10) 
#'# To check summary results and most influential species:
#'summary(samp)
#'\dontrun{
#'# Visual diagnostics for clade removal:
#'sensi_plot(samp)
#'}
#' @export

samp_phyglm <- function(formula, data, phy, n.sim = 30,
                         breaks=seq(.1, .5, .1), btol = 50, track = TRUE, ...)
{
    
if(class(formula) != "formula") stop("formula must be class 'formula'")
if(class(data) != "data.frame") stop("data must be class 'data.frame'")
if(class(phy) != "phylo") stop("phy must be class 'phylo'")
if(length(breaks) < 2) stop("Please include more than one break,
                          e.g. breaks=c(.3,.5)")
else

# Check match between data and phy 
data_phy <- match_dataphy(formula, data, phy, ...)

full.data <- data_phy$data
phy <- data_phy$phy
N <- nrow(full.data)

mod.0 <- phylolm::phyloglm(formula, data = full.data,
                           phy = phy, method = "logistic_MPLE", btol = btol)

if(isTRUE(mod.0$convergence != 0)) stop("Full model failed to converge,
                                              consider changing btol. See ?phyloglm")
intercept.0             <- mod.0$coefficients[[1]]
estimate.0                 <- mod.0$coefficients[[2]]
optpar.0                <- mod.0$alpha
pval.intercept.0        <- phylolm::summary.phyloglm(mod.0)$coefficients[[1,4]]
pval.estimate.0            <- phylolm::summary.phyloglm(mod.0)$coefficients[[2,4]]
aic.0                   <- mod.0$aic

#Creates empty data frame to store model outputs
sensi.estimates<-
        data.frame("n.remov" = numeric(), "n.percent" = numeric(),
                   "intercept" = numeric(), "DIFintercept" = numeric(),
                   "intercept.perc" = numeric(), "pval.intercept" = numeric(),
                   "estimate" = numeric(), "DIFestimate" = numeric(),
                   "estimate.perc" = numeric(), "pval.estimate" = numeric(),
                   "AIC" = numeric(),"optpar" = numeric())

#Loops over breaks, remove percentage of species determined by 'breaks
#and repeat determined by 'n.sim'.
counter=1
limit <- sort(round((breaks) * nrow(full.data), digits = 0))
NL <- length(breaks) * n.sim
if(track==TRUE) pb <- utils::txtProgressBar(min = 0, max = NL, style = 1)
for (i in limit){
    for (j in 1:n.sim){
            exclude <- sample(1:N, i)
            crop.data <- full.data[-exclude, ]
            crop.phy <-  ape::drop.tip(phy,phy$tip.label[exclude])
            mod <- try(phylolm::phyloglm(formula, data = crop.data,
                            phy = crop.phy, method = "logistic_MPLE",
                            btol = btol), TRUE)
            
            if(isTRUE(class(mod)=="try-error")) { next }
            else { 
                intercept             <- mod$coefficients[[1]]
                estimate                 <- mod$coefficients[[2]]
                pval.intercept        <-
                phylolm::summary.phyloglm(mod)$coefficients[[1,4]]
                pval.estimate <- phylolm::summary.phyloglm(mod)$coefficients[[2,4]]
                aic <- mod$aic
                DIFintercept <- intercept - intercept.0
                DIFestimate <- estimate - estimate.0
                intercept.perc  <- round((abs(
                    DIFintercept / intercept.0)) * 100, digits = 1)
                estimate.perc <- round((abs(
                    DIFestimate / estimate.0)) * 100, digits = 1)
                aic <- mod$aic
                optpar <- mod$alpha
                n.remov <- i
                n.percent <- round((n.remov/N)*100,digits=0)
                rep <- j
                
                if(track==TRUE) utils::setTxtProgressBar(pb, counter)
                    # Stores values for each simulation
                    estim.simu <- data.frame(n.remov, n.percent, intercept, 
                                             DIFintercept, intercept.perc,
                                             pval.intercept, estimate,
                                             DIFestimate, estimate.perc,
                                             pval.estimate, aic, optpar,
                                             stringsAsFactors = F)
                    sensi.estimates[counter, ]  <- estim.simu
                    counter=counter+1
            }
    }
}
if(track==TRUE) on.exit(close(pb))

#Calculates Standardized DFbeta and DIFintercept
sDIFintercept <- sensi.estimates$DIFintercept/
  stats::sd(sensi.estimates$DIFintercept)
sDIFestimate     <- sensi.estimates$DIFestimate/
  stats::sd(sensi.estimates$DIFestimate)

sensi.estimates$sDIFintercept <- sDIFintercept
sensi.estimates$sDIFestimate     <- sDIFestimate

#Calculates percentages of signficant intercepts & slopes within breaks.
res                     <- sensi.estimates
n.sim                   <- table(res$n.remov)
breaks                  <- unique(res$n.percent)
sign.intercept          <- res$pval.intercept > .05
sign.estimate              <- res$pval.estimate > .05
res$sign.intercept      <- sign.intercept
res$sign.estimate          <- sign.estimate
perc.sign.intercept <- 1-(with(res,tapply(sign.intercept,n.remov,sum))) / n.sim
perc.sign.estimate     <- 1-(with(res,tapply(sign.estimate,n.remov,sum))) / n.sim
mean.sDIFestimate       <- with(res,tapply(sDIFestimate,n.remov,mean))
mean.sDIFintercept   <- with(res,tapply(sDIFintercept,n.remov,mean))
mean.perc.intercept <- with(res,tapply(intercept.perc,n.remov,mean))
mean.perc.estimate     <- with(res,tapply(estimate.perc,n.remov,mean))
perc.sign.tab       <- data.frame(percent_sp_removed=breaks,
                                  perc.sign.intercept = as.numeric(perc.sign.intercept),
                                  mean.perc.intercept = as.numeric(mean.perc.intercept),
                                  mean.sDIFintercept = as.numeric(mean.sDIFintercept),
                                  perc.sign.estimate = as.numeric(perc.sign.estimate),
                                  mean.perc.estimate = as.numeric(mean.perc.estimate),
                                  mean.sDIFestimate = as.numeric(mean.sDIFestimate))

#Creates a list with full model estimates:
param0<-list(coef = phylolm::summary.phyloglm(mod.0)$coefficients,
             aic = phylolm::summary.phyloglm(mod.0)$aic,
             optpar = mod.0$alpha)

#Generates output:
res <- list(call = match.call(),
            model = "logistic_MPLE",
            formula = formula,
            full.model.estimates = param0,
            sensi.estimates = sensi.estimates,
            sign.analysis = perc.sign.tab,
            data = full.data)
class(res) <- "sensiSamp"
return(res)
        
}
