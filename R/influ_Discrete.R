#' @examples 
#' #Load data:
#' data("primates")
#' #Create a binary trait factor 
#' adultMass_binary<-ifelse(primates$data$adultMass > 7350, "big", "small")
#' adultMass_binary<-as.factor(as.factor(adultMass_binary))
#' names(adultMass_binary)<-rownames(primates$data)
#' #Model trait evolution accounting for phylogenetic uncertainty
#' influ_binary<-influ_Discrete(data = adultMass_binary,phy = primates$phy[[1]],
#' model = "ARD",transform = "none",cutoff = 2,track = T)
#' #Print summary statistics for the transitions rates, aic-values and (if applicable) optimisation parameter
#' summary(influ_binary)
#' #Use a different evolutionary model or transformation, 
#' e.g. symmetrical rates, with an Early Burst (EB) model of trait evolution
#' influ_binary_SYM_EB<-influ_Discrete(data = adultMass_binary,phy = primates$phy[[1]],
#' model = "SYM",transform = "EB",n.tree = 30,track = T)
#' summary(influ_binary_SYM_EB)
#' #Or change the cutoff
#' influ_binary<-influ_Discrete(data = adultMass_binary,phy = primates$phy[[1]],
#' model = "ARD",transform = "none",cutoff = 1.2,track = T)

#' @export

influ_Discrete <- function(data,phy,model,
                           transform = "none",bounds = list(),
                           cutoff=2,track=TRUE,...){
  
  #Error check
  if(is.null(model)) stop("model must be specified (e.g. 'ARD' or 'SYM'")
  if(class(data)!="factor") stop("data must supplied as a factor with species as names. Consider as.factor()")
  if(length(levels(data))>2) stop("discrete data can have maximal two levels")
  if(class(phy)!="phylo") stop("phy must be class 'phylo'")
  if(transform=="white") stop("the white-noise (non-phylogenetic) model is not allowed")
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
  sensi.estimates<-data.frame("species" =numeric(),
                              "q12"=numeric(),"DIFq12"= numeric(),"q12.perc"= numeric(),
                              "q21"=numeric(),"DIFq21"= numeric(),"q21.perc"= numeric(),
                              "aicc"=numeric(),"optpar"=numeric()) 
  
  #Loops over all species, and removes each one individually
  counter <- 1
  errors <- NULL
  if(track==TRUE) pb <- utils::txtProgressBar(min = 0, max = N, style = 3)
  for (i in 1:N){
    
    crop.data <- full.data[c(1:N)[-i]]
    crop.phy <-  ape::drop.tip(phy,setdiff(phy$tip.label,names(crop.data)))
    
    mod = try(geiger::fitDiscrete(phy = crop.phy,dat = crop.data,
                                  model = model,transform = transform,
                                  bounds = bounds,ncores = NULL,...),TRUE)
    if(isTRUE(class(mod)=="try-error")) {
      error <- i
      names(error) <- rownames(full.data$data)[i]
      errors <- c(errors,error)
      next }
    else {  sp                   <- phy$tip.label[i]
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
    
    if(track==TRUE) utils::setTxtProgressBar(pb, i)
    # Stores values for each simulation
    # Store reduced model parameters: 
    estim.simu <- data.frame(sp, q12, DIFq12,q12.perc,
                             q21, DIFq21,q21.perc,
                             aicc, optpar,
                             stringsAsFactors = F)
    sensi.estimates[counter, ]  <- estim.simu
    counter=counter+1
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
  
  #Creates a list with full model estimates:
  #full model estimates:
  param0 <- list(q12=q12.0,q21=q21.0,
                 aicc=aicc.0,
                 optpar=optpar.0)
  
  #Identifies influencital species (sDF > cutoff) and orders by influence
  reorder.on.q21         <-sensi.estimates[order(abs(
    sensi.estimates$sDIFq21),decreasing=T),c("species","sDIFq21")]
  influ.sp.q21           <-as.character(reorder.on.q21$species[abs(
    reorder.on.q21$sDIFq21)>cutoff])
  reorder.on.q12     <-sensi.estimates[order(abs(
    sensi.estimates$sDIFq12),decreasing=T),c("species","sDIFq12")]
  influ.sp.q12       <-as.character(reorder.on.q12$species[abs(
    reorder.on.q12$sDIFq12)>cutoff])
  
  #Generates output:
  res <- list(   call = match.call(),
                 cutoff=cutoff,
                 data = full.data,
                 optpar = transform,
                 full.model.estimates = param0,
                 influential.species= list(influ.sp.q12=influ.sp.q12,
                                           influ.sp.q21=influ.sp.q21),
                 sensi.estimates=sensi.estimates,
                 errors = errors)
  class(res) <- "sensiInflu.TrailEvol"
  ### Warnings:
  if (length(res$errors) >0){
    warning("Some species deletion presented errors, please check: output$errors")}
  else {
    res$errors <- "No errors found."
  }
  
  return(res)
  
}

