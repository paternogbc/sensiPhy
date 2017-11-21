
#' @export

clade_Discrete <- function(data, phy, model = "ARD",transform = "none",
                              track = TRUE,clade.col, n.species = 5, n.sim = 100, ...) {
  # Error checking:
  if(class(data)!="factor") stop("data must supplied as a factor with species as names. Consider as.factor(")
  if(length(levels(data))>2) stop("discrete data can have maximal two levels")
  if(missing(clade.col)) stop("clade.col not defined. Please, define the",
                              " column with clade names.")
  if(class(phy)!="phylo") stop("phy must be class 'phylo'")
  if(transform=="white") stop("the white-noise (non-phylogenetic) model is not allowed")
  else
  
  #Calculates the full model, extracts model parameters
  data_phy <- match_dataphy(formula, data, phy, ...)
  phy <- data_phy$phy
  full.data <- data_phy$data
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
  N               <- nrow(full.data)
  mod.0           <- phylolm::phylolm(formula, data=full.data,
                                      model=model,phy=phy)
  intercept.0      <- mod.0$coefficients[[1]]
  estimate.0          <- mod.0$coefficients[[2]]
  pval.intercept.0 <- phylolm::summary.phylolm(mod.0)$coefficients[[1,4]]
  pval.estimate.0     <- phylolm::summary.phylolm(mod.0)$coefficients[[2,4]]
  optpar.0 <- mod.0$optpar
  
  #Create dataframe to store estmates for each clade
  sensi.estimates <-
    data.frame("clade" =I(as.character()), 
               "N.species" = numeric(),"intercept"=numeric(),
               "DIFintercept"=numeric(),"intercept.perc"=numeric(),
               "pval.intercept"=numeric(),"estimate"=numeric(),
               "DIFestimate"=numeric(),"estimate.perc"=numeric(),
               "pval.estimate"=numeric(),"AIC"=numeric(),
               "optpar" = numeric())
  
  # Create dataframe store simulations (null distribution)
  null.dist <- data.frame("clade" = rep(names(uc), each = n.sim),
                          "intercept"= numeric(length(uc)*n.sim),
                          "estimate" = numeric(length(uc)*n.sim),
                          "DIFintercept"=numeric(length(uc)*n.sim),
                          "DIFestimate"=numeric(length(uc)*n.sim))
  
  
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
    crop.sp <-   which(full.data[ ,clade.col] %in% A)
    crop.phy <-  ape::drop.tip(phy,phy$tip.label[crop.sp])
    mod=try(phylolm::phylolm(formula, data=crop.data,model=model,
                             phy=crop.phy),TRUE)
    intercept           <- mod$coefficients[[1]]
    estimate            <- mod$coefficients[[2]]
    DIFintercept        <- intercept - intercept.0
    DIFestimate         <- estimate - estimate.0
    intercept.perc      <- round((abs(DIFintercept / intercept.0)) * 100,
                                  digits = 1)
    estimate.perc       <- round((abs(DIFestimate / estimate.0)) * 100,
                                  digits = 1)
    pval.intercept      <- phylolm::summary.phylolm(mod)$coefficients[[1,4]]
    pval.estimate       <- phylolm::summary.phylolm(mod)$coefficients[[2,4]]
    aic.mod             <- mod$aic
    if (model == "BM" | model == "trend"){
      optpar <- NA
    }
    if (model != "BM" & model != "trend" ){
      optpar               <- mod$optpar
    }
    
    # Store reduced model parameters: 
    estim.simu <- data.frame(A, cN, intercept, DIFintercept, intercept.perc,
                             pval.intercept, estimate, DIFestimate, estimate.perc,
                             pval.estimate, aic.mod, optpar,
                             stringsAsFactors = F)
    sensi.estimates[aa, ]  <- estim.simu
    
    ### START LOOP FOR NULL DIST:
    # number of species in clade A:
    for (i in 1:n.sim) {
      exclude <- sample(1:N, cN)
      crop.data <- full.data[-exclude,]
      crop.phy <-  ape::drop.tip(phy,phy$tip.label[exclude])
      mod <- try(phylolm::phylolm(formula, data = crop.data,
                                  model = model,phy = crop.phy),TRUE)
      intercept       <- mod$coefficients[[1]]
      estimate        <- mod$coefficients[[2]]
      DIFintercept    <- intercept - intercept.0
      DIFestimate     <- estimate - estimate.0
      
      null.dist[bb, ] <- data.frame(clade = as.character(A), 
                                     intercept,
                                     estimate,
                                     DIFintercept,
                                     DIFestimate)
      
      if(track==TRUE) utils::setTxtProgressBar(pb, bb)
      bb <- bb + 1
    }
    aa <- aa + 1
  }
  if(track==TRUE) on.exit(close(pb))
  
  #OUTPUT
  #full model estimates:
  param0 <- list(coef=phylolm::summary.phylolm(mod.0)$coefficients,
                 aic=phylolm::summary.phylolm(mod.0)$aic,
                 optpar=mod.0$optpar)
  
  #Generates output:
  res <- list(call = match.call(),
              model = model,
              formula = formula,
              full.model.estimates = param0,
              sensi.estimates = sensi.estimates,
              null.dist = null.dist, 
              data = full.data,
              errors = errors,
              clade.col = clade.col)
  class(res) <- "sensiClade"
  ### Warnings:
  if (length(res$errors) >0){
    warning("Some clades deletion presented errors, please check: output$errors")}
  else {
    res$errors <- "No errors found."
  }
  return(res)
}

