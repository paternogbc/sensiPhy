#' @examples 
#' #Load data:
#' data("primates")
#' #Create a binary trait factor 
#' primates$data$adultMass_binary<-ifelse(primates$data$adultMass > 7350, "big", "small")
#' primate_phy_pruned<-drop.tip(phy=primates$phy[[1]],
#' tip=setdiff(primates$phy$tip.label,rownames(primates$data)))
#' clade_test<-clade_Discrete(data=primates$data,phy = primate_phy_pruned,
#' trait.col = "adultMass_binary",clade.col="family",nsim=5)
#' summary(clade_test)
#' @export

clade_Discrete <- function(data, phy, model,transform = "none",
                           trait.col,clade.col,n.species = 5, n.sim = 20,
                           bounds = list(), track=TRUE, ...) {
  # Error checking:
  if(is.null(model)) stop("model must be specified (e.g. 'ARD' or 'SYM'")
  if(!is.data.frame(data)) stop("data must be class 'data.frame'")
  if(missing(clade.col)) stop("clade.col not defined. Please, define the",
                              " column with clade names.")
  if(class(phy)!="phylo") stop("phy must be class 'phylo'")
  if(transform=="white") stop("the white-noise (non-phylogenetic) model is not allowed")
  if(length(which(!phy$tip.label %in% rownames(data)))>0) stop("not all tips are present in data, prune tree")
  if(length(which(!rownames(data) %in% phy$tip.label))>0) stop("not all data species are present in tree, remove superfluous data points")
  else
    
    #Calculates the full model, extracts model parameters
    full.data<-data
    phy <- phy
  if (is.na(match(clade.col, names(full.data)))) {
    stop("Names column '", clade.col, "' not found in data frame'")
  }
  
  # Identify CLADES to use and their sample size 
  all.clades <- levels(full.data[ ,clade.col])
  wc <- table(full.data[ ,clade.col]) > n.species
  uc <- table(full.data[ , clade.col])[wc]
  
  #k <- names(which(table(full.data[,clade.col]) > n.species ))
  if (length(uc) == 0) stop(paste("There is no clade with more than ",
                                  n.species," species. Change 'n.species' to fix this
                                  problem",sep=""))
  
  # FULL MODEL PARAMETERS:
  trait_vec_full<-full.data[,trait.col]
  trait_vec_full<-as.factor(trait_vec_full)
  if(length(levels(trait_vec_full))>2) stop("discrete data can have maximal two levels")
  names(trait_vec_full)<-rownames(full.data)
  
  N                   <- nrow(full.data)
  mod.0               <- geiger::fitDiscrete(phy = phy,dat = trait_vec_full,
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
  
  #Create dataframe to store estmates for each clade
  sensi.estimates<-data.frame("clade" =I(as.character()),"N.species" = numeric(),
                              "q12"=numeric(),"DIFq12"= numeric(),"q12.perc"= numeric(),
                              "q21"=numeric(),"DIFq21"= numeric(),"q21.perc"= numeric(),
                              "aicc"=numeric(),"optpar"=numeric()) 
  
  # Create dataframe store simulations (null distribution)
  null.dist <- data.frame("clade" = rep(names(uc), each = n.sim),
                          "q12"= numeric(length(uc)*n.sim),
                          "DIFq12"= numeric(length(uc)*n.sim),
                          "q21" = numeric(length(uc)*n.sim),
                          "DIFq21" = numeric(length(uc)*n.sim))
  
  
  ### START LOOP between CLADES:
  # counters:
  aa <- 1; bb <- 1
  errors <- NULL
  
  if(track==TRUE) pb <- utils::txtProgressBar(min = 0, max = length(uc)*n.sim,
                                              style = 3)
  for (A in names(uc)){
    
    ### Number of species in clade A
    cN  <- as.numeric(uc[names(uc) == A])
    
    ### Fit reduced model (without clade)
    crop.data <- full.data[!full.data[ ,clade.col] %in% A,]
    crop.phy <-  ape::drop.tip(phy,setdiff(phy$tip.label,rownames(crop.data)))
    crop.trait_vec<-crop.data[,trait.col]
    crop.trait_vec<-as.factor(crop.trait_vec)
    names(crop.trait_vec)<-rownames(crop.data)
    mod = try(geiger::fitDiscrete(phy = crop.phy,dat = crop.trait_vec,
                                  model = model,transform = transform,
                                  bounds = bounds,ncores = NULL,...),TRUE)
    q12               <- mod$opt$q12
    q21               <- mod$opt$q12
    DIFq12            <- q12 - q12.0
    DIFq21            <- q21 - q21.0
    q12.perc      <- round((abs(DIFq12 / q12.0)) * 100,
                           digits = 1)
    q21.perc       <- round((abs(DIFq21 / q21.0)) * 100,
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
    
    # Store reduced model parameters: 
    estim.simu <- data.frame(A, cN, q12, DIFq12,q12.perc,
                             q21, DIFq21,q21.perc,
                             aicc, optpar,
                             stringsAsFactors = F)
    sensi.estimates[aa, ]  <- estim.simu
    
    ### START LOOP FOR NULL DIST:
    # number of species in clade A:
    for (i in 1:n.sim) {
      exclude <- sample(1:N, cN)
      crop.data <- full.data[-exclude,]
      crop.phy <-  ape::drop.tip(phy,setdiff(phy$tip.label,rownames(crop.data)))
      crop.trait_vec<-crop.data[,trait.col]
      crop.trait_vec<-as.factor(crop.trait_vec)
      names(crop.trait_vec)<-rownames(crop.data)
      mod = try(geiger::fitDiscrete(phy = crop.phy,dat = crop.trait_vec,
                                    model = model,transform = transform,
                                    bounds = bounds,ncores = NULL,...),TRUE)
      
      if(isTRUE(class(mod)=="try-error")) {
        error <- i
        names(error) <- rownames(full.data$data)[i]
        errors <- c(errors,error)
        next }
      else 
        
        q12               <- mod$opt$q12
      q21               <- mod$opt$q12
      aicc              <- mod$opt$aicc
      DIFq12            <- q12 - q12.0
      DIFq21            <- q21 - q21.0
      
      null.dist[bb, ] <- data.frame(clade = as.character(A), 
                                    q12,
                                    DIFq12,
                                    q21,
                                    DIFq21)
      
      if(track==TRUE) utils::setTxtProgressBar(pb, bb)
      bb <- bb + 1
    }
    aa <- aa + 1
  }
  if(track==TRUE) on.exit(close(pb))
  
  #OUTPUT
  #full model estimates:
  param0 <- list(q12=q12.0,q21=q21.0,
                 aicc=aicc.0,
                 optpar=optpar.0)
  
  #Generates output:
  res <- list(   call = match.call(),
                 data = full.data,
                 full.model.estimates = param0,
                 sensi.estimates=sensi.estimates,
                 null.dist = null.dist,
                 errors = errors,
                 optpar = transform,
                 clade.col = clade.col)
  class(res) <- "sensiClade.TraitEvol"
  ### Warnings:
  if (length(res$errors) >0){
    warning("Some clades deletion presented errors, please check: output$errors")}
  else {
    res$errors <- "No errors found."
  }
  return(res)
}

