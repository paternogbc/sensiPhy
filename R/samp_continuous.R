#' Species Sampling uncertainty - Trait Evolution Continuous Characters
#' 
#' Fits models for trait evolution of continuous characters, 
#' evaluating sampling uncertainty. 
#'
#' @param data Data vector for a single binary trait, with names matching tips in \code{phy}.
#' @param phy A phylogeny (class 'phylo') matching \code{data}.
#' @param n.sim The number of times species are randomly deleted for each \code{break}.
#' @param breaks A vector containing the percentages of species to remove.
#' @param model The evolutionary model (see Details). 
#' @param bounds settings to contstrain parameter estimates. See \code{\link[geiger]{fitContinuous}}
#' @param n.cores number of cores to use. If 'NULL', number of cores is detected.
#' @param track Print a report tracking function progress (default = TRUE)
#' @param ... Further arguments to be passed to \code{\link[geiger]{fitContinuous}}
#' @details
#' This function randomly removes a given percentage of species (controlled by \code{breaks}),  
#' fits different models of continuous character evolution using \code{\link[geiger]{fitContinuous}}, 
#' repeats this this many times (controlled by \code{n.sim}), stores the results and calculates 
#' the effects on model parameters.
#' 
#' Different evolutionary models from \code{fitContinuous} can be used, i.e. \code{BM},\code{OU},
#' \code{EB}, \code{trend}, \code{lambda}, \code{kappa}, \code{delta} and \code{drift}.
#' 
#' See \code{\link[geiger]{fitContinuous}} for more details on character models and tree transformations. 
#' 
#' Output can be visualised using \code{sensi_plot}. [Not yet!]
#'
#' @return The function \code{tree_continuous} returns a list with the following
#' components:
#' @return \code{call}: The function call
#' @return \code{data}: The original full data vector
#' @return \code{optpar}: Transformation parameter used (e.g. \code{lambda}, \code{kappa} etc.)
#' \code{full.model.estimates}: Parameter estimates (rate of evolution \code{sigsq}, 
#' root state \code{z0} and where applicable \code{optpar}), 
#' AICc and the optimised value of the phylogenetic transformation parameter (e.g. \code{lambda}) 
#' for the full model without deleted species.
#' @return \code{break.summary.tab}: Summary per \code{break} of the mean and median effects 
#' of species removal on percentage and absolute change parameter estimates. 
#' @return \code{sensi.estimates}: Parameter estimates, 
#' AICc and the optimised value of the phylogenetic transformation parameter (e.g. \code{lambda}) 
#' for each analysis with a different phylogenetic tree.
#' @return \code{N.tree}: Number of trees \code{n.tree} analysed
#' @return \code{stats}: Main statistics for model parameters, i.e. minimum, maximum, mean, median and sd-values
#' @return \code{optpar}: Transformation parameter used (e.g. \code{lambda}, \code{kappa} etc.)
#' @author Gijsbert Werner & Gustavo Paterno
#' @seealso \code{\link[geiger]{fitContinuous}}
#' @references 
#' Yang Z. 2006. Computational Molecular Evolution. Oxford University Press: Oxford. 
#' 
#' Harmon Luke J, Jason T Weir, Chad D Brock, Richard E Glor, and Wendell Challenger. 2008.
#' GEIGER: investigating evolutionary radiations. Bioinformatics 24:129-131.
#' 
#' Werner, G.D.A., Cornwell, W.K., Sprent, J.I., Kattge, J. & Kiers, E.T. (2014). 
#' A single evolutionary innovation drives the deep evolution of symbiotic N2-fixation in angiosperms. Nature Communications, 5, 4087.
#' 
#' @examples 
#' \dontshow{
#' #Load data:
#' data("primates")
#' #Model trait evolution accounting for phylogenetic uncertainty
#' adultMass<-primates$data$adultMass
#' names(adultMass)<-rownames(primates$data)
#' samp_cont<-samp_continuous(data = adultMass,phy = primates$phy[[1]],
#' model = "OU",n.sim=1,breaks=c(.1,.2),n.cores = 2, track = TRUE)
#' }
#' \dontrun{
#' #Load data:
#' data("primates")
#' #Model trait evolution accounting for phylogenetic uncertainty
#' adultMass<-primates$data$adultMass
#' names(adultMass)<-rownames(primates$data)
#' samp_cont<-samp_continuous(data = adultMass,phy = primates$phy[[1]],
#' model = "OU",n.sim=25,breaks=seq(.05,.2,.05),n.cores = 2, track = TRUE)
#' #Print summary statistics
#' summary(samp_cont)
#' sensi_plot(samp_cont)
#' sensi_plot(samp_cont, graphs = 1)
#' #Use a different evolutionary model 
#' samp_cont2<-samp_continuous(data = adultMass,phy = primates$phy[[1]],
#' model = "kappa",n.sim=25,breaks=seq(.05,.2,.05),n.cores = 2,track = TRUE)
#' summary(samp_cont2)
#' sensi_plot(samp_cont2)
#' sensi_plot(samp_cont2, graphs = 2)
#' samp_cont3<-samp_continuous(data = adultMass,phy = primates$phy[[1]],
#' model = "BM",n.sim=25,breaks=seq(.05,.2,.05),n.cores = 2,track = TRUE)
#' summary(samp_cont3)
#' }
#' @export

samp_continuous <- function(data,phy,n.sim=30,
                          breaks=seq(.1,.5,.1),
                          model,n.cores = NULL,
                          bounds=list(),track=TRUE,...){
  
  #Error check
  if(is.null(model)) stop("model must be specified, e.g. 'OU' or 'lambda'")
  if(class(data)!="numeric" | is.null(names(data))) stop("data must supplied as a numeric vector with species as names")
  if(class(phy)!="phylo") stop("phy must be class 'phylo'")
  if(model=="white") stop("the white-noise (non-phylogenetic) model is not allowed")
  if ( (model == "drift") & (ape::is.ultrametric(phy))) stop("A drift model is unidentifiable for ultrametric trees., see ?fitContinuous for details")
  if(length(breaks) < 2) 
    stop("Please include more than one break, e.g. breaks=c(.3,.5)")
  else

    #Matching tree and phylogeny
    full.data<-data
    phy<-phy
    
    #Calculates the full model, extracts model parameters
    N                   <- length(full.data)
    mod.0               <- geiger::fitContinuous(phy = phy,dat = full.data,
                                               model = model,
                                               bounds = bounds,ncores = n.cores,...)
    sigsq.0               <- mod.0$opt$sigsq
    z0.0                  <- mod.0$opt$z0
    aicc.0              <- mod.0$opt$aicc
    if (model == "BM"){
      optpar.0 <- NA
    }
    if (model == "OU"){
      optpar.0        <- mod.0$opt$alpha
    }
    if (model == "EB"){
      optpar.0               <- mod.0$opt$a
    }
    if (model == "trend"){
      optpar.0               <- mod.0$opt$slope
    }
    if (model == "lambda"){
      optpar.0               <- mod.0$opt$lambda
    }
    if (model == "kappa"){
      optpar.0               <- mod.0$opt$kappa
    }
    if (model == "delta"){
      optpar.0               <- mod.0$opt$delta
    }
    if (model == "drift"){
      optpar.0               <- mod.0$opt$drift
    }
    
    #Creates empty data frame to store model outputs
    sensi.estimates<-data.frame("n.remov" = numeric(), "n.percent"= numeric(),
                                "sigsq"=numeric(),"DIFsigsq"= numeric(),"sigsq.perc"= numeric(),
                                "optpar"=numeric(),"DIFoptpar"=numeric(),"optpar.perc"=numeric(),
                                "z0"=numeric(),
                                "aicc"=numeric())
        
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
                mod = try(geiger::fitContinuous(phy = crop.phy,dat = crop.data,
                                              model = model,
                                              bounds = bounds,ncores = n.cores,...),TRUE)
                if(isTRUE(class(mod) == "try-error")) {next}
                else {
                  sigsq               <- mod$opt$sigsq
                  z0                  <- mod$opt$z0
                  aicc              <- mod$opt$aicc
                  DIFsigsq            <- sigsq - sigsq.0
                  sigsq.perc          <- round((abs(DIFsigsq / sigsq.0)) * 100,
                                               digits = 1)
                  if (model == "BM"){
                    optpar <- NA
                  }
                  if (model == "OU"){
                    optpar        <- mod$opt$alpha
                  }
                  if (model == "EB"){
                    optpar               <- mod$opt$a
                  }
                  if (model == "trend"){
                    optpar               <- mod$opt$slope
                  }
                  if (model == "lambda"){
                    optpar               <- mod$opt$lambda
                  }
                  if (model == "kappa"){
                    optpar               <- mod$opt$kappa
                  }
                  if (model == "delta"){
                    optpar               <- mod$opt$delta
                  }
                  if (model == "drift"){
                    optpar              <- mod$opt$drift
                  }
                  DIFoptpar            <- optpar - optpar.0
                  optpar.perc        <- round((abs(DIFoptpar / optpar.0)) * 100,
                                              digits = 1)
                  
                n.remov <- i
                n.percent <- round( (n.remov / N) * 100,digits = 0)
                #rep <- j
                
                if(track == TRUE) (utils::setTxtProgressBar(pb, counter))
                # Store reduced model parameters: 
                estim.simu <- data.frame(n.remov, n.percent,
                                         sigsq, DIFsigsq,sigsq.perc,
                                         optpar,DIFoptpar,optpar.perc,
                                         z0,
                                         aicc,
                                         stringsAsFactors = F)
                sensi.estimates[counter, ]  <- estim.simu
                counter <- counter + 1
                }
            }
        }
        if(track==TRUE) on.exit(close(pb))
                
        #Calculates Standardized DFs
        sDIFsigsq <- sensi.estimates$DIFsigsq/
          stats::sd(sensi.estimates$DIFsigsq)
        sDIFoptpar     <- sensi.estimates$DIFoptpar/
          stats::sd(sensi.estimates$DIFoptpar)
        
        sensi.estimates$sDIFsigsq     <- sDIFsigsq
        sensi.estimates$sDIFoptpar    <- sDIFoptpar
        
        #Calculates stats
        res                 <- sensi.estimates
        n.sim               <- table(res$n.remov)
        breaks              <- unique(res$n.percent)
        mean.sDIFsigsq   <- with(res,tapply(sDIFsigsq,n.remov,mean))
        mean.sDIFoptpar  <- with(res,tapply(sDIFoptpar,n.remov,mean))
        mean.perc.optpar <- with(res,tapply(optpar.perc,n.remov,mean))
        mean.perc.sigsq  <- with(res,tapply(sigsq.perc,n.remov,mean))
        median.sDIFsigsq   <- with(res,tapply(sDIFsigsq,n.remov,median))
        median.sDIFoptpar  <- with(res,tapply(sDIFoptpar,n.remov,median))
        breaks.summary.tab       <- data.frame(percent_sp_removed=breaks,
                                          mean.perc.sigsq = as.numeric(mean.perc.sigsq),
                                          mean.sDIFsigsq = as.numeric(mean.sDIFsigsq),
                                          median.sDIFsigsq = as.numeric(median.sDIFsigsq),
                                          mean.perc.optpar = as.numeric(mean.perc.optpar),
                                          mean.sDIFoptpar = as.numeric(mean.sDIFoptpar),
                                          median.sDIFoptpar = as.numeric(median.sDIFoptpar))
        
        #Creates a list with full model estimates:
        param0 <- list(sigsq=sigsq.0,optpar=optpar.0,
                       z0=z0.0,
                       aicc=aicc.0)
        
        #Generates output:
        res <- list(   call = match.call(),
                       data = full.data,
                       optpar = model,
                       full.model.estimates = param0,
                       breaks.summary.tab = breaks.summary.tab,
                       sensi.estimates=sensi.estimates)
        
        class(res) <- "sensiSamp.TraitEvol"
        return(res)
        
        }
