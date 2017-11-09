#' Phylogenetic uncertainty - Trait evolution (Binary traits)
#'

#' @examples 
#' #Load data:
#' data("primates")
#' #Create a binary trait factor 
#' adultMass_binary<-ifelse(primates.data$adultMass > 7350, "big", "small")
#' adultMass_binary<-as.factor(as.factor(adultMass_binary))
#' names(adultMass_binary)<-rownames(primates.data)
#' #Model trait evolution accounting for phylogenetic uncertainty
#' tree_binary<-tree_Discrete(data = adultMass_binary,phy = primates.phy,model = "ER",transform = "none",n.tree = 30,track = T)
#'  #Use a different model or transformation
#'  tree_binary_ARD<-tree_Discrete(data = adultMass_binary,phy = primates.phy,model = "ARD",transform = "none",n.tree = 10,track = T)
#'  tree_binary_ARD_lambda<-tree_Discrete(data = adultMass_binary,phy = primates.phy,model = "ARD",transform = "lambda",n.tree = 10,track = T)
#'#Visual and other diagnostics to do still 
#' @export

tree_Discrete <- function(data,phy,model = "ARD",transform = "none",bounds = list(),
                         n.tree=2,track=TRUE,...){
  #Error check
  if(class(data)!="factor") stop("data must supplied as a factor with species as names. Consider as.factor(")
  if(length(levels(data))>2) stop("discrete data can have maximal two levels")
  if(class(phy)!="multiPhylo") stop("phy must be class 'multiPhylo'")
  if(length(phy)<n.tree) stop("'n.tree' must be smaller (or equal) than the number of trees in the 'multiPhylo' object")
  if(transform=="white") stop("the white-noise (non-phylogenetic) model is not allowed")
  else

  
  #Matching tree and phylogeny using utils.R
  #datphy<-match_dataphy(formula,data,phy, ...) #We could update this to work for what is needed for fitDiscrete
  full.data<-data
  phy<-phy

  # If the class of tree is multiphylo pick n=n.tree random trees
  trees<-sample(length(phy),n.tree,replace=F)

  #Create the results data.frame
  sensi.estimates<-data.frame("n.tree"=numeric(),"q12"=numeric(),"q21"=numeric(),
                              "aicc"=numeric(),"optpar"=numeric()) 
  #Model calculation
  counter=1
  errors <- NULL
  c.data<-list()
  if(track==TRUE) pb <- utils::txtProgressBar(min = 0, max = n.tree, style = 3)
  for (j in trees){
      
      #Match data order to tip order
      full.data <- full.data[phy[[j]]$tip.label]
    
      #phylolm model
      mod = try(geiger::fitDiscrete(phy = phy[[j]],dat = full.data,model = model,transform = transform,
                                    bounds = bounds,ncores = NULL),FALSE)

      if(isTRUE(class(mod)=="try-error")) {
        error <- j
        names(error) <- rownames(c.data$full.data)[j]
        errors <- c(errors,error)
        next }
      
      else{
        q12           <- mod$opt$q12
        q21           <- mod$opt$q21
        aicc          <- mod$opt$aicc

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
        
        if(track==TRUE) utils::setTxtProgressBar(pb, counter)
        
        #write in a table
        estim.simu <- data.frame(j, q12, q21, aicc, optpar,
                                 stringsAsFactors = F)
        sensi.estimates[counter, ]  <- estim.simu
        counter=counter+1
        
      }
    }
  if(track==TRUE) on.exit(close(pb))
  #calculate mean and sd for each parameter

  statresults<-data.frame(min=apply(sensi.estimates,2,min),
                          max=apply(sensi.estimates,2,max),
                          mean=apply(sensi.estimates,2,mean),
                          median=apply(sensi.estimates,2,median),
                          sd=apply(sensi.estimates,2,sd))[-1,]
  
  res <- list(   call = match.call(),
                 data=full.data,
                 sensi.estimates=sensi.estimates,N.obs=n,
                 stats = round(statresults[c(1:4),],digits=6))
  class(res) <- "sensiTree.TraitEvol"
  return(res)
}
