#' @examples 
#' #Load data:
#' data("primates")
#' #Create a binary trait factor 
#' adultMass_binary<-ifelse(primates$data$adultMass > 7350, "big", "small")
#' adultMass_binary<-as.factor(as.factor(adultMass_binary))
#' names(adultMass_binary)<-rownames(primates$data)
#' #Model trait evolution accounting for phylogenetic uncertainty
#' samp_binary<-samp_Discrete(data = adultMass_binary,phy = primates$phy[[1]],
#' n.sim=25,breaks=seq(.1,.3,.1),model = "SYM",transform = "none",track = T)
#' #Print summary statistics for the transitions rates, aic-values and (if applicable) optimisation parameter
#' summary(samp_binary)
#' #Use a different evolutionary model or transformation, 
#' e.g. all-rates-different, with an Early Burst (EB) model of trait evolution
#' samp_binary_ARD_EB<-samp_Discrete(data = adultMass_binary,phy = primates$phy[[1]],
#' n.sim=25,breaks=seq(.1,.3,.1),model = "ARD",transform = "EB",track = T)
#' summary(samp_binary_ARD_EB)
#' @export

samp_Discrete <- function(data,phy,n.sim=30,
                          breaks=seq(.1,.5,.1),
                          model,transform="none",
                          bounds=list(),track=TRUE,...){
  
  #Error check
  if(is.null(model)) stop("model must be specified (e.g. 'ARD' or 'SYM'")
  if(class(data)!="factor") stop("data must supplied as a factor with species as names. Consider as.factor()")
  if(length(levels(data))>2) stop("discrete data can have maximal two levels")
  if(class(phy)!="phylo") stop("phy must be class 'phylo'")
  if(transform=="white") stop("the white-noise (non-phylogenetic) model is not allowed")
  if(length(breaks) < 2) 
    stop("Please include more than one break, e.g. breaks=c(.3,.5)")
  else

    #Matching tree and phylogeny
    full.data<-data
    phy<-phy
    
    #Calculates the full model, extracts model parameters
    N                   <- length(full.data)
    mod.0               <- geiger::fitDiscrete(phy = phy,dat = full.data,
                                               model = model,transform = transform,
                                               bounds = bounds,ncores = NULL,...)
    q12.0               <- mod.0$opt$q12
    q21.0               <- mod.0$opt$q21
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
        limit <- sort(round( (breaks) * length(full.data),digits=0))
        NL <- length(breaks) * n.sim
        if(track==TRUE) pb <- utils::txtProgressBar(min = 0, max = NL, style = 3)
        for (i in limit){
            for (j in 1:n.sim){
              #Prep simulation data
                exclude <- sample(1:N,i)
                crop.data <- full.data[-exclude]
                crop.phy <-  ape::drop.tip(phy,setdiff(phy$tip.label,names(crop.data)))
              #Run the model
                mod = try(geiger::fitDiscrete(phy = crop.phy,dat = crop.data,
                                              model = model,transform = transform,
                                              bounds = bounds,ncores = NULL,...),TRUE)
                if(isTRUE(class(mod) == "try-error")) {next}
                else {  
                  q12               <- mod$opt$q12
                  q21               <- mod$opt$q21
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
                
        #Calculates Standardized DFbeta and DIFq12
        sDIFq12 <- sensi.estimates$DIFq12/
          stats::sd(sensi.estimates$DIFq12)
        sDIFq21     <- sensi.estimates$DIFq21/
          stats::sd(sensi.estimates$DIFq21)
        
        sensi.estimates$sDIFq21     <- sDIFq21
        sensi.estimates$sDIFq12     <- sDIFq12
        
        #Calculates stats
        res                 <- sensi.estimates
        n.sim               <- table(res$n.remov)
        breaks              <- unique(res$n.percent)
        mean.sDIFq12   <- with(res,tapply(sDIFq12,n.remov,mean))
        mean.sDIFq21  <- with(res,tapply(sDIFq21,n.remov,mean))
        mean.perc.q21 <- with(res,tapply(q21.perc,n.remov,mean))
        mean.perc.q12  <- with(res,tapply(q12.perc,n.remov,mean))
        median.sDIFq12   <- with(res,tapply(sDIFq12,n.remov,median))
        median.sDIFq21  <- with(res,tapply(sDIFq21,n.remov,median))
        breaks.summary.tab       <- data.frame(percent_sp_removed=breaks,
                                          mean.perc.q12 = as.numeric(mean.perc.q12),
                                          mean.sDIFq12 = as.numeric(mean.sDIFq12),
                                          median.sDIFq12 = as.numeric(median.sDIFq12),
                                          mean.perc.q21 = as.numeric(mean.perc.q21),
                                          mean.sDIFq21 = as.numeric(mean.sDIFq21),
                                          median.sDIFq21 = as.numeric(median.sDIFq21))
        
        #Creates a list with full model estimates:
        param0 <- list(q12=q12.0,q21=q21.0,
                       aicc=aicc.0,
                       optpar=optpar.0)
        
        #Generates output:
        res <- list(   call = match.call(),
                       data = full.data,
                       optpar = transform,
                       full.model.estimates = param0,
                       breaks.summary.tab = breaks.summary.tab,
                       sensi.estimates=sensi.estimates)
        
        class(res) <- "sensiSamp.TraitEvol"
        return(res)
        
        }
