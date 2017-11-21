#' Phylogenetic uncertainty - Trait Evolution Discrete Characters
#' 
#' Fits models for trait evolution of discrete (binary) characters, 
#' evaluating phylogenetic uncertainty. 
#'
#' @param data Data vector for a single binary trait, with names matching tips in \code{phy}.
#' @param phy Phylogenies (class 'multiPhylo', see ?\code{ape}).
#' @param n.tree Number of times to repeat the analysis with n different trees picked 
#' randomly in the multiPhylo file. If NULL, \code{n.tree} = 5
#' @param model The Mkn model to use (see Details). Default is \code{ARD}.
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
#' @return The function \code{tree_Discrete} returns a list with the following
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
#' data("alien")
#' #Create a binary trait factor 
#' adultMass_binary<-ifelse(alien$data$adultMass > 50000, "big", "small")
#' adultMass_binary<-as.factor(as.factor(adultMass_binary))
#' names(adultMass_binary)<-rownames(alien$data)
#' #Model trait evolution accounting for phylogenetic uncertainty
#' tree_binary<-tree_Discrete(data = adultMass_binary,phy = alien$phy,model = "ARD",transform = "none",n.tree = 30,track = T)
#' #Print summary statistics for the transitions rates, aic-values and (if applicable) optimisation parameter
#' summary(tree_binary)
#' sensi_plot(tree_binary)
#' sensi_plot(tree_binary,graphs="q12")
#' sensi_plot(tree_binary,graphs="q21")
#' #Use a different evolutionary model or transformation.
#' tree_binary_lambda<-tree_Discrete(data = adultMass_binary,phy = alien$phy,model = "ARD",transform = "lambda",n.tree = 30,track = T)
#' summary(tree_binary_lambda) #Using Pagel's Lambda
#' sensi_plot(tree_binary_lamda)  
#' #Symmetrical rates, with an Early Burst (EB) model of trait evolution
#' tree_binary_SYM_EB<-tree_Discrete(data = adultMass_binary,phy = alien$phy,model = "SYM",transform = "EB",n.tree = 30,track = T)
#' summary(tree_binary_SYM_EB)
#' sensi_plot(tree_binary_lamda) 
#' sensi_plot(tree_binary_lamda,graphs="optpar") 
#' @export

tree_Discrete <- function(data,phy,n.tree=5,model = "ARD",
                          transform = "none",bounds = list(),
                         track=TRUE,...){
  #Error check
  if(class(data)!="factor") stop("data must supplied as a factor with species as names. Consider as.factor(")
  if(length(levels(data))>2) stop("discrete data can have maximal two levels")
  if(class(phy)!="multiPhylo") stop("phy must be class 'multiPhylo'")
  if(length(phy)<n.tree) stop("'n.tree' must be smaller (or equal) than the number of trees in the 'multiPhylo' object")
  if(transform=="white") stop("the white-noise (non-phylogenetic) model is not allowed")
  else

  
  #Matching tree and phylogeny using utils.R
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
                                    bounds = bounds,ncores = NULL,...),FALSE)

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
  #Output
  res <- list(   call = match.call(),
                 data=full.data,
                 sensi.estimates=sensi.estimates,N.tree=n.tree,
                 stats = statresults[c(1:4),],
                 optpar = transform)
  class(res) <- "sensiTree.TraitEvol"
  return(res)
}
