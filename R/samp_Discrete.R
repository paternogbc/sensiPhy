
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

    #Matching tree and phylogeny using utils.R
    full.data<-data
    phy<-phy
    
    #Calculates the full model, extracts model parameters
    N                   <- length(full.data)
    mod.0               <- geiger::fitDiscrete(phy = phy,dat = full.data,
                                               model = model,transform = transform,
                                               bounds = bounds,ncores = NULL,...)
    q12.0               <- mod.0$opt$q12
    q21.0               <- mod.0$opt$q12
    aicc.0              <- mod.0$opt$aicc
    if (transform == "none"){
      optpar.0 <- NA
    }
    if (transform == "EB"){
      optpar.0               <- mod.0$opt$a
    }
    if (transform == "lambda"){
      optpar.0               <- mod.0$opt$lambda
    }
    if (transform == "kappa"){
      optpar.0               <- mod.0$opt$kappa
    }
    if (transform == "delta"){
      optpar.0               <- mod.0$opt$delta
    }
    
    
    #Creates empty data frame to store model outputs
    sensi.estimates<-data.frame("n.remov" = numeric(), "n.percent"= numeric(),
                                "q12"=numeric(),"DIFq12"= numeric(),"q12.perc"= numeric(),
                                "q21"=numeric(),"DIFq21"= numeric(),"q21.perc"= numeric(),
                                "aicc"=numeric(),"optpar"=numeric()) 

        
        #Loops over breaks, remove percentage of species determined by 'breaks
        #and repeat determined by 'n.sim'.
        counter <- 1
        limit <- sort(round( (breaks) * nrow(full.data),digits=0))
        NL <- length(breaks) * n.sim
        if(track==TRUE) pb <- utils::txtProgressBar(min = 0, max = NL, style = 3)
        for (i in limit){
            for (j in 1:n.sim){
              #Prep simulation data
                exclude <- sample(1:N,i)
                crop.data <- full.data[-exclude]
                crop.phy <-  ape::drop.tip(phy,setdiff(phy$tip.label,rownames(crop.data)))
              #Run the model
                mod = try(geiger::fitDiscrete(phy = crop.phy,dat = crop.data,
                                              model = model,transform = transform,
                                              bounds = bounds,ncores = NULL,...),TRUE)
                if(isTRUE(class(mod) == "try-error")) {next}
                else {  
                  q12               <- mod$opt$q12
                  q21               <- mod$opt$q12
                  DIFq12            <- q12 - q12.0
                  DIFq21            <- q21 - q21.0
                  q12.perc          <- round((abs(DIFq12 / q12.0)) * 100,
                                             digits = 1)
                  q21.perc          <- round((abs(DIFq21 / q21.0)) * 100,
                                             digits = 1)
                  aicc              <- mod$opt$aicc
                  if (transform == "none"){
                    optpar <- NA
                  }
                  if (transform == "EB"){
                    optpar               <- mod$opt$a
                  }
                  if (transform == "lambda"){
                    optpar               <- mod$opt$lambda
                  }
                  if (transform == "kappa"){
                    optpar               <- mod$opt$kappa
                  }
                  if (transform == "delta"){
                    optpar               <- mod$opt$delta
                  }
                  
                n.remov <- i
                n.percent <- round( (n.remov / N) * 100,digits = 0)
                #rep <- j
                
                if(track == TRUE) (utils::setTxtProgressBar(pb, counter))
                # Stores values for each simulation
                # Stores values for each simulation
                # Store reduced model parameters: 
                estim.simu <- data.frame(n.remov, n.percent, q12, DIFq12,q12.perc,
                                         q21, DIFq21,q21.perc,
                                         aicc, optpar,
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
