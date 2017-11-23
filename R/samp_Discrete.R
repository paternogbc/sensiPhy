
#' @export

samp_Discrete <- function(data,phy,n.sim=30,
                          breaks=seq(.1,.5,.1),
                          model="ARD",transform="none",
                          bounds=list(),track=TRUE,...){
  
  #Error check
  if(class(data)!="factor") stop("data must supplied as a factor with species as names. Consider as.factor()")
  if(length(levels(data))>2) stop("discrete data can have maximal two levels")
  if(class(phy)!="phylo") stop("phy must be class 'phylo'")
  if(transform=="white") stop("the white-noise (non-phylogenetic) model is not allowed")
  if(length(breaks) < 2) 
    stop("Please include more than one break, e.g. breaks=c(.3,.5)")
  else

# Check match between data and phy 
data_phy <- match_dataphy(formula, data, phy, ...)
#Calculates the full model, extracts model parameters
full.data         <- data_phy$data
phy               <- data_phy$phy
N                 <- nrow(full.data)
mod.0             <- phylolm::phylolm(formula, data = full.data, model = model, phy = phy)
intercept.0       <- mod.0$coefficients[[1]]
estimate.0        <- mod.0$coefficients[[2]]
pval.intercept.0  <- phylolm::summary.phylolm(mod.0)$coefficients[[1,4]]
pval.estimate.0   <- phylolm::summary.phylolm(mod.0)$coefficients[[2,4]]
optpar.0          <- mod.0$optpar
aic.0             <- mod.0$aic

#Creates empty data frame to store model outputs
sensi.estimates <-
    data.frame("n.remov" = numeric(), "n.percent"= numeric(),
               "intercept"= numeric(),"DIFintercept"= numeric(),
               "intercept.perc"= numeric(),"pval.intercept"=numeric(),
               "estimate"= numeric(),"DIFestimate"= numeric(),
               "estimate.perc"= numeric(),"pval.estimate"= numeric(),
               "AIC"= numeric(),"optpar" = numeric())

#Loops over breaks, remove percentage of species determined by 'breaks
#and repeat determined by 'n.sim'.
counter <- 1
limit <- sort(round( (breaks) * nrow(full.data),digits=0))
NL <- length(breaks) * n.sim
if(track==TRUE) pb <- utils::txtProgressBar(min = 0, max = NL, style = 3)
for (i in limit){
    for (j in 1:n.sim){
        exclude <- sample(1:N,i)
        crop.data <- full.data[-exclude,]
        crop.phy <-  ape::drop.tip(phy,phy$tip.label[exclude])
        mod <- try(phylolm::phylolm(formula, data = crop.data,
                                    model = model,phy = crop.phy),TRUE)
        if(isTRUE(class(mod) == "try-error")) {next}
        else {  
        intercept       <- mod$coefficients[[1]]
        estimate        <- mod$coefficients[[2]]
        optpar          <- mod$optpar
        pval.intercept  <- phylolm::summary.phylolm(mod)$coefficients[[1,4]]
        pval.estimate   <- phylolm::summary.phylolm(mod)$coefficients[[2,4]]
        aic             <- mod$aic
        DIFintercept    <- intercept - intercept.0
        DIFestimate     <- estimate - estimate.0
        intercept.perc  <- round( (abs(DIFintercept / intercept.0)) * 100, digits = 1)
        estimate.perc   <- round( (abs(DIFestimate / estimate.0)) * 100, digits = 1)
        aic             <- mod$aic
        
        if (model == "BM" | model == "trend"){
            optpar <- NA
        }
        if (model != "BM" & model != "trend" ){
            optpar               <- mod$optpar
        }
        
        n.remov <- i
        n.percent <- round( (n.remov / N) * 100,digits = 0)
        #rep <- j
        
        if(track == TRUE) (utils::setTxtProgressBar(pb, counter))
        # Stores values for each simulation
        estim.simu <- data.frame(n.remov, n.percent, intercept, 
                                 DIFintercept, intercept.perc,
                                 pval.intercept, estimate,
                                 DIFestimate, estimate.perc,
                                 pval.estimate, aic, optpar,
                                 stringsAsFactors = F)
        sensi.estimates[counter, ]  <- estim.simu
        counter <- counter + 1
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
res                 <- sensi.estimates
n.sim               <- table(res$n.remov)
breaks              <- unique(res$n.percent)
sign.intercept      <- res$pval.intercept > .05
sign.estimate       <- res$pval.estimate > .05
res$sign.intercept  <- sign.intercept
res$sign.estimate   <- sign.estimate
perc.sign.intercept <- 1-(with(res,tapply(sign.intercept,n.remov,sum))) / n.sim
perc.sign.estimate  <- 1-(with(res,tapply(sign.estimate,n.remov,sum))) / n.sim
mean.sDIFestimate   <- with(res,tapply(sDIFestimate,n.remov,mean))
mean.sDIFintercept  <- with(res,tapply(sDIFintercept,n.remov,mean))
mean.perc.intercept <- with(res,tapply(intercept.perc,n.remov,mean))
mean.perc.estimate  <- with(res,tapply(estimate.perc,n.remov,mean))
perc.sign.tab       <- data.frame(percent_sp_removed=breaks,
                                  perc.sign.intercept = as.numeric(perc.sign.intercept),
                                  mean.perc.intercept = as.numeric(mean.perc.intercept),
                                  mean.sDIFintercept = as.numeric(mean.sDIFintercept),
                                  perc.sign.estimate = as.numeric(perc.sign.estimate),
                                  mean.perc.estimate = as.numeric(mean.perc.estimate),
                                  mean.sDIFestimate = as.numeric(mean.sDIFestimate))

#Creates a list with full model estimates:
param0 <- list(coef=phylolm::summary.phylolm(mod.0)$coefficients,
    aic = phylolm::summary.phylolm(mod.0)$aic,
    optpar = optpar.0)

#Generates output:
res <- list(call = match.call(),
            model = model,
            formula = formula,
            full.model.estimates = param0,
            sensi.estimates = sensi.estimates,
            sign.analysis = perc.sign.tab,
            data = full.data)
class(res) <- "sensiSamp"
return(res)

}
