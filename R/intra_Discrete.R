#' @export

tree_Discrete <- function(data,phy,n.tree=10,model = "ARD",
                          transform = "none",bounds = list(),
                          track=TRUE,...){
  #Error check
  if(class(data)!="factor") stop("data must supplied as a factor with species as names. Consider as.factor()")
  if(length(levels(data))>2) stop("discrete data can have maximal two levels")
  if(class(phy)!="multiPhylo") stop("phy must be class 'multiPhylo'")
  if(length(phy)<n.tree) stop("'n.tree' must be smaller (or equal) than the number of trees in the 'multiPhylo' object")
  if(transform=="white") stop("the white-noise (non-phylogenetic) model is not allowed")
  else

intra_phylm <- function(data, phy,
                        model = "ARD", transform = "none",
                        Vx = NULL,
                        y.transf = NULL, x.transf = NULL,
                        n.intra = 30, distrib = "normal",
                        bounds = list(),track = TRUE, ...){
  #Error check
  if(class(data)!="factor") stop("data must supplied as a factor with species as names. Consider as.factor()")
  if(length(levels(data))>2) stop("discrete data can have maximal two levels")
  if(class(phy)!="multiPhylo") stop("phy must be class 'multiPhylo'")
  if(transform=="white") stop("the white-noise (non-phylogenetic) model is not allowed")
  if(is.null(Vx)) stop("Vx must be defined")
  if(class(phy) != "phylo") stop("phy must be class 'phylo'")
  if(distrib == "normal") warning ("distrib=normal: make sure that standard deviation 
                                 is provided for Vx")
  else
  

    
  #Matching tree and phylogeny using utils.R
  datphy <- match_dataphy(formula, data, phy, ...)
  full.data <- datphy[[1]]
  phy <- datphy[[2]]
  
  resp <- all.vars(formula)[1]
  pred <- all.vars(formula)[2]
  
  if(!is.null(Vy) && sum(is.na(full.data[, Vy])) != 0) {
    full.data[is.na(full.data[, Vy]), Vy] <- 0}
  
  if(!is.null(Vx) && sum(is.na(full.data[, Vx])) != 0) {
    full.data[is.na(full.data[, Vx]), Vx] <- 0}
  
  
  
  #Function to pick a random value in the interval
  if (distrib == "normal") funr <- function(a,b) {stats::rnorm(1,a,b)}
  else  funr <- function(a,b) {stats::runif(1, a - b, a + b)}
  
  #Create the results data.frame
  sensi.estimates <- data.frame("n.intra" = numeric(),"intercept" = numeric(),
                                    "se.intercept" = numeric(), 
                                    "pval.intercept" = numeric(),
                                    "estimate" = numeric(),
                                    "se.estimate" = numeric(),
                                    "pval.estimate" = numeric(),"aic" = numeric(),
                                    "optpar" = numeric())
  #Model calculation
  counter = 1
  errors <- NULL
  species.NA <- list()
  if(track == TRUE) pb <- utils::txtProgressBar(min = 0, max = n.intra, style = 3)
  for (i in 1:n.intra) {
    ##Set response and predictor variables
    #Vy is not provided or is not numeric, do not pick random value
    if(!inherits(full.data[,resp], c("numeric","integer")) || is.null(Vy)) 
        {full.data$respV <- stats::model.frame(formula, data = full.data)[,1]}
    
    #choose a random value in [mean-se,mean+se] if Vy is provided
    if (!is.null(Vy))
    {full.data$respV <- apply(full.data[,c(resp,Vy)],1,function(x)funr(x[1],x[2]))}
    
    #Vx is not provided or is not numeric, do not pick random value
    if (!inherits(full.data[,pred], c("numeric","integer")) || is.null(Vx))
        {full.data$predV <- stats::model.frame(formula, data = full.data)[,2]}
    
    #choose a random value in [mean-se,mean+se] if Vx is provided
    if(!is.null(Vx))
    {full.data$predV <- apply(full.data[,c(pred,Vx)],1,function(x)funr(x[1],x[2]))}
    
    #transform Vy and/or Vx if x.transf and/or y.transf are provided
    if(!is.null(y.transf)) 
    {suppressWarnings (full.data$respV <- y.transf(full.data$respV))}
    
    if(!is.null(x.transf)) 
    {suppressWarnings (full.data$predV <- x.transf(full.data$predV))}
    
    #skip iteration if there are NA's in the dataset
    species.NA[[i]]<-rownames(full.data[with(full.data,is.na(predV) | is.na(respV)),])
    if(sum(is.na(full.data[,c("respV","predV")])>0)) next
    
    #model
    mod = try(phylolm::phylolm(respV ~ predV, data = full.data, model = model,
                               phy = phy), FALSE)
    
    if(isTRUE(class(mod) == "try-error")) {
      error <- i
      errors <- c(errors,error)
      next }
    
    
    else{
      intercept            <- phylolm::summary.phylolm(mod)$coefficients[[1,1]]
      se.intercept         <- phylolm::summary.phylolm(mod)$coefficients[[1,2]]
      estimate             <- phylolm::summary.phylolm(mod)$coefficients[[2,1]]
      se.estimate          <- phylolm::summary.phylolm(mod)$coefficients[[2,2]]
      pval.intercept       <- phylolm::summary.phylolm(mod)$coefficients[[1,4]]
      pval.estimate        <- phylolm::summary.phylolm(mod)$coefficients[[2,4]]
      aic.mod              <- mod$aic
      n                    <- mod$n

      if (model == "BM") {
        optpar <- NA
      }
      if (model != "BM") {
        optpar <- mod$optpar
      }
      if(track == TRUE) utils::setTxtProgressBar(pb, i)
      #write in a table
      estim.simu <- data.frame(i, intercept, se.intercept, pval.intercept,
                               estimate, se.estimate, pval.estimate, aic.mod, optpar,
                               stringsAsFactors = F)
      sensi.estimates[counter, ]  <- estim.simu
      counter=counter+1
      
    }
  }
  if(track == TRUE) on.exit(close(pb))
  
  #calculate mean and sd for each parameter
  #variation due to intraspecific variability
  mean_by_randomval <- stats::aggregate(.~n.intra, data = sensi.estimates,
                                        mean)
  
  statresults <- data.frame(min = apply(sensi.estimates, 2, min),
                          max = apply(sensi.estimates, 2, max),
                          mean = apply(sensi.estimates, 2, mean),
                          sd_intra = apply(mean_by_randomval, 2, stats::sd))[-1, ]
  
  statresults$CI_low  <- statresults$mean - stats::qt(0.975, df = n.intra-1) * statresults$sd_intra / sqrt(n.intra)
  statresults$CI_high <- statresults$mean + stats::qt(0.975, df = n.intra-1) * statresults$sd_intra / sqrt(n.intra)
  
  #species with transformation problems
  nr <- n.intra - nrow(sensi.estimates)
  sp.pb <- unique(unlist(species.NA))
  if (length(sp.pb) >0) 
  warning (paste("in", nr,"simulations, data transformations generated NAs, please consider using another function
  for x.transf or y.transf and check output$sp.pb",sep=" "))
  
                                       
  res <- list(call = match.call(),
              formula = formula,
              y.transf = y.transf, 
              x.transf = x.transf,
              data = full.data,
              sensi.estimates = sensi.estimates, N.obs = n,
              stats = round(statresults[c(1:6),c(3,5,6)],digits=3),
              all.stats = statresults,sp.pb=sp.pb)
  class(res) <- "sensiIntra"
  return(res)
}

