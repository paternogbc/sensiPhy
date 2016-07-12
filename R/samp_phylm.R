#' Sensitivity Analysis Species Sampling  - Phylogenetic Linear Regression
#'
#' Performs analyses of sensitivity to species sampling by randomly removing
#' species and detecting the effects on parameter estimates in a phylogenetic
#' linear regression.
#'
#' @param formula The model formula
#' @param data Data frame containing species traits with row names matching tips
#' in \code{phy}.
#' @param phy A phylogeny (class 'phylo') matching \code{data}.
#' @param times The number of times species are randomly deleted for each
#' \code{break}.
#' @param breaks A vector containing the percentages of species to remove.
#' @param model The phylogenetic model to use (see Details). Default is \code{lambda}.
#' @param track Print a report tracking function progress (default = TRUE)
#' @param ... Further arguments to be passed to \code{phylolm}
#' @details
#'
#' This function randomly removes a given percentage of species (controlled by
#' \code{breaks}) from the full phylogenetic linear regression, fits a phylogenetic
#' linear regression model without these species using \code{\link[phylolm]{phylolm}},
#' repeats this many times (controlled by \code{times}), stores the results and
#' calculates the effects on model parameters.
#'
#' All phylogenetic models from \code{phylolm} can be used, i.e. \code{BM},
#' \code{OUfixedRoot}, \code{OUrandomRoot}, \code{lambda}, \code{kappa},
#' \code{delta}, \code{EB} and \code{trend}. See ?\code{phylolm} for details.
#'
#' Currently, this function can only implement simple linear models (i.e. \eqn{trait~
#' predictor}). In the future we will implement more complex models.
#'
#' Output can be visualised using \code{sensi_plot}.
#' @return The function \code{samp_phylm} returns a list with the following
#' components:
#' @return \code{formula}: The formula
#' @return \code{full.model.estimates}: Coefficients, aic and the optimised
#' value of the phylogenetic parameter (e.g. \code{lambda} or \code{kappa}) for
#' the full model without deleted species.
#' @return \code{samp.model.estimates}: A data frame with all simulation
#' estimates. Each row represents a model rerun with a given number of species
#' \code{n.remov} removed, representing \code{n.percent} of the full dataset.
#' Columns report the calculated regression intercept (\code{intercept}),
#' difference between simulation intercept and full model intercept (\code{DFintercept}),
#' the percentage of change in intercept compared to the full model (\code{intercept.perc})
#' and intercept p-value (\code{pval.intercept}). All these parameters are also reported
#' for the regression slope (\code{DFslope} etc.). Additionally, model aic value
#' (\code{AIC}) and the optimised value (\code{optpar}) of the phylogenetic
#' parameter (e.g. \code{kappa} or \code{lambda}, depending on the phylogenetic model
#' used) are reported.
#' @return \code{sign.analysis} For each break (i.e. each percentage of species
#' removed) this reports the percentage of statistically signficant (at p<0.05)
#' intercepts (\code{perc.sign.intercept}) over all repititions as well as the
#' percentage of statisticaly significant (at p<0.05) slopes (\code{perc.sign.slope}).
#' @return \code{data}: Original full dataset.
#' @author Gustavo Paterno & Gijsbert D.A. Werner
#' @seealso \code{\link[phylolm]{phylolm}}, \code{\link{samp_phyglm}},
#' \code{\link{influ_phylm}},\code{\link{sensi_plot}}
#' @references 
#' #' Ho, L. S. T. and Ane, C. 2014. "A linear-time algorithm for 
#' Gaussian and non-Gaussian trait evolution models". Systematic Biology 63(3):397-408.
#' @import ape phylolm
#' 
#' @examples 
#' \dontrun{
#' # Load data:
#' data(alien)
#' # Run analysis:
#' samp <- samp_phylm(log(gestaLen) ~ log(adultMass), phy = alien$phy[[1]], 
#' data = alien$data)
#' summary(samp)
#' head(samp$samp.model.estimates)
#' # Visual diagnostics
#' sensi_plot(samp)
#' # You can specify which graph and parameter ("slope" or "intercept") to print: 
#' sensi_plot(samp, graphs = 1, param = "slope")
#' sensi_plot(samp, graphs = 2, param = "intercept")
#' }
#' @export

samp_phylm <- function(formula,data,phy,times=30,
                         breaks=seq(.1,.5,.1),model="lambda",track=TRUE,...){
# Basic error checking:
if(class(formula) != "formula") 
    stop("formula must be class 'formula'")
if(class(data) != "data.frame") 
    stop("data must be class 'data.frame'")
if(class(phy) != "phylo") 
    stop("phy must be class 'phylo'")
if(length(breaks) < 2) 
    stop("Please include more than one break, e.g. breaks=c(.3,.5)")
if ( (model == "trend") & (ape::is.ultrametric(phy)))
    stop("Trend is unidentifiable for ultrametric trees.,
         see ?phylolm for details")
else

# Check match between data and phy 
data_phy <- match_dataphy(formula, data, phy)
#Calculates the full model, extracts model parameters
full.data <- data_phy$data
phy <- data_phy$phy
N       <- nrow(full.data)
mod.0   <- phylolm::phylolm(formula, data = full.data, model = model, phy = phy)
intercept.0      <- mod.0$coefficients[[1]]
slope.0          <- mod.0$coefficients[[2]]
pval.intercept.0 <- phylolm::summary.phylolm(mod.0)$coefficients[[1,4]]
pval.slope.0     <- phylolm::summary.phylolm(mod.0)$coefficients[[2,4]]
optpar.0         <- mod.0$optpar
aic.0            <- mod.0$aic

#Creates empty data frame to store model outputs
samp.model.estimates <-
    data.frame("n.remov" = numeric(), "n.percent"= numeric(),
               "intercept"= numeric(),"DFintercept"= numeric(),
               "intercept.perc"= numeric(),"pval.intercept"=numeric(),
               "slope"= numeric(),"DFslope"= numeric(),
               "slope.perc"= numeric(),"pval.slope"= numeric(),
               "AIC"= numeric(),"optpar" = numeric())

#Loops over breaks, remove percentage of species determined by 'breaks
#and repeat determined by 'times'.
counter <- 1
limit <- sort(round( (breaks) * nrow(full.data),digits=0))
NL <- length(breaks) * times
pb <- utils::txtProgressBar(min = 0, max = NL, style = 1)
for (i in limit){
    for (j in 1:times){
        exclude <- sample(1:N,i)
        crop.data <- full.data[-exclude,]
        crop.phy <-  ape::drop.tip(phy,phy$tip.label[exclude])
        mod <- try(phylolm::phylolm(formula, data = crop.data,
                                    model = model,phy = crop.phy),TRUE)
        if(isTRUE(class(mod) == "try-error")) {next}
        else {  
        intercept      <- mod$coefficients[[1]]
        slope          <- mod$coefficients[[2]]
        optpar         <- mod$optpar
        pval.intercept <- phylolm::summary.phylolm(mod)$coefficients[[1,4]]
        pval.slope     <- phylolm::summary.phylolm(mod)$coefficients[[2,4]]
        aic            <- mod$aic
        DFintercept    <- intercept - intercept.0
        DFslope        <- slope - slope.0
        intercept.perc <- round( (abs(DFintercept / intercept.0)) * 100, digits = 1)
        slope.perc     <- round( (abs(DFslope / slope.0)) * 100, digits = 1)
        aic            <- mod$aic
        
        if (model == "BM" | model == "trend"){
            optpar <- NA
        }
        if (model != "BM" & model != "trend" ){
            optpar               <- mod$optpar
        }
        
        n.remov <- i
        n.percent <- round( (n.remov / N) * 100,digits = 0)
        #rep <- j
        
        if(track == TRUE) (
          utils::setTxtProgressBar(pb, counter))
        # Stores values for each simulation
        estim.simu <- data.frame(n.remov, n.percent, intercept, 
                                 DFintercept, intercept.perc,
                                 pval.intercept, slope,
                                 DFslope, slope.perc,
                                 pval.slope, aic, optpar,
                                 stringsAsFactors = F)
        samp.model.estimates[counter, ]  <- estim.simu
        counter <- counter + 1
        }
    }
}
close(pb)
#Calculates percentages of signficant intercepts & slopes within breaks.
res                 <- samp.model.estimates
times               <- table(res$n.remov)
breaks              <- unique(res$n.percent)
sign.intercept      <- res$pval.intercept > .05
sign.slope          <- res$pval.slope > .05
res$sign.intercept  <- sign.intercept
res$sign.slope      <- sign.slope
perc.sign.intercept <- 1-(with(res,tapply(sign.intercept,n.remov,sum))) / times
perc.sign.slope     <- 1-(with(res,tapply(sign.slope,n.remov,sum))) / times
perc.sign.tab       <- data.frame(percent_sp_removed=breaks,
                        perc.sign.intercept = as.numeric(perc.sign.intercept),
                            perc.sign.slope = as.numeric(perc.sign.slope))

#Creates a list with full model estimates:
param0 <- list(coef=phylolm::summary.phylolm(mod.0)$coefficients,
    aic = phylolm::summary.phylolm(mod.0)$aic,
    optpar = optpar.0)

#Generates output:
res <- list(model = model,
            formula = formula,
            full.model.estimates = param0,
            samp.model.estimates = samp.model.estimates,
            sign.analysis = perc.sign.tab,
            data = full.data)
class(res) <- "sensiSamp"
return(res)

}
