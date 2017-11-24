#' Phylogenetic uncertainty - Trait Evolution Discrete Characters
#' 
#' Fits models for trait evolution of discrete (binary) characters, 
#' evaluating phylogenetic uncertainty. 
#'
#' @param data Data vector for a single binary trait, with names matching tips in \code{phy}.
#' @param phy Phylogenies (class 'multiPhylo', see ?\code{ape}).
#' @param n.tree Number of times to repeat the analysis with n different trees picked 
#' randomly in the multiPhylo file. If NULL, \code{n.tree} = 10
#' @param model The Mkn model to use (see Details). 
#' @param transform The evolutionary model to transform the tree (see Details). Default is \code{none}.
#' @param bounds settings to contstrain parameter estimates. See \code{\link[geiger]{fitDiscrete}}
#' @param track Print a report tracking function progress (default = TRUE)
#' @param ... Further arguments to be passed to \code{\link[geiger]{fitDiscrete}}
#' @details
#' This function fits different models of discrete character evolution using \code{\link[geiger]{fitDiscrete}}
#' to n trees, randomly picked in a multiPhylo file. Currenlty, only binary discrete traits are supported
#'
#' Different character model from \code{fitDiscrete} can be used, including \code{ER} (equal-rates), 
#' \code{SYM} (symmetric), \code{ARD} (all-rates-different) and \code{meristic} (stepwise fashion). 
#'
#' All transformations to the phylogenetic tree from \code{fitDiscrete} can be used, i.e. \code{none},
#' \code{EB}, \code{lambda}, \code{kappa} and\code{delta}.
#' 
#' See \code{\link[geiger]{fitDiscrete}} for more details on character models and tree transformations. 
#' 
#' Output can be visualised using \code{sensi_plot}.
#'
#' @return The function \code{tree_discrete} returns a list with the following
#' components:
#' @return \code{call}: The function call
#' @return \code{data}: The original full data vector
#' @return \code{sensi.estimates}: Parameter estimates (transition rates q12 and q21), 
#' AICc and the optimised value of the phylogenetic transformation parameter (e.g. \code{lambda}) 
#' for each analysis with a different phylogenetic tree.
#' @return \code{N.tree}: Number of trees \code{n.tree} analysed
#' @return \code{stats}: Main statistics for model parameters, i.e. minimum, maximum, mean, median and sd-values
#' @return \code{optpar}: Transformation parameter used (e.g. \code{lambda}, \code{kappa} etc.)
#' @author Gijsbert Werner & Caterina Penone
#' @seealso \code{\link[geiger]{fitDiscrete}}
#' @references 
#' Yang Z. 2006. Computational Molecular Evolution. Oxford University Press: Oxford. 
#' 
#' Harmon Luke J, Jason T Weir, Chad D Brock, Richard E Glor, and Wendell Challenger. 2008.
#' GEIGER: investigating evolutionary radiations. Bioinformatics 24:129-131.
#' 
#' @examples 
#' #Load data:
#' data("primates")
#' #Model trait evolution accounting for phylogenetic uncertainty
#' adultMass<-primates$data$adultMass
#' names(adultMass)<-rownames(primates$data)
#' tree_cont<-tree_continuous(data = adultMass,phy = primates$phy,
#' model = "OU",n.tree=30,track = TRUE)
#' #Print summary statistics for the transitions rates, aic-values and (if applicable) optimisation parameter
#' summary(tree_cont)
#' #Use a different evolutionary model 
#' tree_cont2<-tree_continuous(data = adultMass,phy = primates$phy,
#' model = "lambda",n.tree=30,track = TRUE)
#' @export

tree_continuous <- function(data,phy,n.tree=10,model,
                          bounds = list(),
                         track=TRUE,...){
  #Error check
  if(is.null(model)) stop("model must be specified, e.g. 'OU' or 'lambda'")
  if(class(data)!="numeric" | is.null(names(data))) stop("data must supplied as a numeric vector with species as names")
  if(class(phy)!="multiPhylo") stop("phy must be class 'multiPhylo'")
  if ((model == "drift") & (ape::is.ultrametric(phy))) stop("A drift model is unidentifiable for ultrametric trees., see ?fitContinuous for details")
  if(length(phy)<n.tree) stop("'n.tree' must be smaller (or equal) than the number of trees in the 'multiPhylo' object")
  if(model=="white") stop("the white-noise (non-phylogenetic) model is not allowed")
  else

  #Matching tree and phylogeny using utils.R
  full.data<-data
  phy<-phy

  # If the class of tree is multiphylo pick n=n.tree random trees
  trees<-sample(length(phy),n.tree,replace=F)

  #Create the results data.frame
  sensi.estimates<-data.frame("n.tree"=numeric(),"sigsq"=numeric(),"optpar"=numeric(),
                              "z0"=numeric(),"aicc"=numeric()) 
  
  #Model calculation
  counter=1
  errors <- NULL
  c.data<-list()
  if(track==TRUE) pb <- utils::txtProgressBar(min = 0, max = n.tree, style = 3)
  for (j in trees){
      
      #Match data order to tip order
      full.data <- full.data[phy[[j]]$tip.label]
    
      #phylolm model
      mod = try(geiger::fitContinuous(phy = phy[[j]],dat = full.data,model = model,
                                    bounds = bounds,ncores = NULL,...),FALSE)

      if(isTRUE(class(mod)=="try-error")) {
        error <- j
        names(error) <- rownames(c.data$full.data)[j]
        errors <- c(errors,error)
        next }
      
      else{
        sigsq               <- mod$opt$sigsq
        z0                  <- mod$opt$z0
        aicc                <- mod$opt$aicc
        
        if (model == "none"){
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
        
        if(track==TRUE) utils::setTxtProgressBar(pb, counter)
        
        #write in a table
        estim.simu <- data.frame(j, sigsq, optpar, z0, aicc,
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
  #Output
  res <- list(   call = match.call(),
                 data=full.data,
                 sensi.estimates=sensi.estimates,N.tree=n.tree,
                 stats = statresults[c(1:4),],
                 optpar = transform)
  class(res) <- "sensiTree.TraitEvol"
  return(res)
}
