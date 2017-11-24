#' Influential Species Detection - Trait Evolution Continous Characters
#' Change help still 
#' Fits models for trait evolution of discrete (binary) characters, 
#' detecting influential species. 
#'
#' @param data Data vector for a single binary trait, with names matching tips in \code{phy}.
#' @param phy A phylogeny (class 'phylo') matching \code{data}.
#' @param model The Mkn model to use (see Details). 
#' @param transform The evolutionary model to transform the tree (see Details). Default is \code{none}.
#' @param bounds settings to contstrain parameter estimates. See \code{\link[geiger]{fitDiscrete}}
#' @param track Print a report tracking function progress (default = TRUE)
#' @param ... Further arguments to be passed to \code{\link[geiger]{fitDiscrete}}
#' @details
#' This function sequentially removes one species at a time,  
#' fits different models of discrete character evolution using \code{\link[geiger]{fitDiscrete}}, 
#' repeats this this many times (controlled by \code{n.sim}), stores the results and calculates 
#' the effects on model parameters Currently, only binary discrete traits are supported.
#' 
#' #' \code{influ_discrete} detects influential species based on the standardised
#' difference in q12 or q21 when removing a given species compared
#' to the full model including all species. Species with a standardised difference
#' above the value of \code{cutoff} are identified as influential. 
#' 
#' Different character model from \code{fitDiscrete} can be used, including \code{ER} (equal-rates), 
#' \code{SYM} (symmetric), \code{ARD} (all-rates-different) and \code{meristic} (stepwise fashion). 
#'
#' All transformations to the phylogenetic tree from \code{fitDiscrete} can be used, i.e. \code{none},
#' \code{EB}, \code{lambda}, \code{kappa} and\code{delta}.
#' 
#' See \code{\link[geiger]{fitDiscrete}} for more details on character models and tree transformations. 
#' 
#' Output can be visualised using \code{sensi_plot}. [Not yet!]
#'
#' @return The function \code{tree_discrete} returns a list with the following
#' components:
#' @return \code{call}: The function call
#' @return \code{cutoff}: The value selected for \code{cutoff}
#' @return \code{data}: The original full data vector
#' @return \code{optpar}: Transformation parameter used (e.g. \code{lambda}, \code{kappa} etc.)
#' @return \code{full.model.estimates}: Parameter estimates (transition rates q12 and q21), 
#' AICc and the optimised value of the phylogenetic transformation parameter (e.g. \code{lambda}) 
#' for the full model without deleted species.
#' @return \code{influential_species}: List of influential species, based on standardised 
#' difference in estimates for q12 and q21. Species are ordered from most influential to 
#' less influential and only include species with a standardised difference > \code{cutoff}.
#' @return \code{sensi.estimates}: Parameter estimates (transition rates q12 and q21), 
#' AICc and the optimised value of the phylogenetic transformation parameter (e.g. \code{lambda}) 
#' for each analysis with a different phylogenetic tree.
#' @author Gijsbert Werner & Gustavo Paterno
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
#' #Create a binary trait factor 
#' adultMass_binary<-ifelse(primates$data$adultMass > 7350, "big", "small")
#' adultMass_binary<-as.factor(as.factor(adultMass_binary))
#' names(adultMass_binary)<-rownames(primates$data)
#' #Model trait evolution accounting for phylogenetic uncertainty
#' influ_binary<-influ_discrete(data = adultMass_binary,phy = primates$phy[[1]],
#' model = "ARD",transform = "none",cutoff = 2,track = T)
#' #Print summary statistics for the transitions rates, aic-values and (if applicable) optimisation parameter
#' summary(influ_binary)
#' #Use a different evolutionary model or transformation, 
#' e.g. symmetrical rates, with an Early Burst (EB) model of trait evolution
#' influ_binary_SYM_EB<-influ_discrete(data = adultMass_binary,phy = primates$phy[[1]],
#' model = "SYM",transform = "EB",n.tree = 30,track = T)
#' summary(influ_binary_SYM_EB)
#' #Or change the cutoff
#' influ_binary<-influ_discrete(data = adultMass_binary,phy = primates$phy[[1]],
#' model = "ARD",transform = "none",cutoff = 1.2,track = T)

#' @export

influ_continuous <- function(data,phy,model,
                           bounds = list(),
                           cutoff=2,track=TRUE,...){
  
            #Error check
            if(is.null(model)) stop("model must be specified (e.g. 'ARD' or 'SYM'")
            if(class(data)!="factor") stop("data must supplied as a factor with species as names. Consider as.factor()")
            if(length(levels(data))>2) stop("discrete data can have maximal two levels")
            if(class(phy)!="phylo") stop("phy must be class 'phylo'")
            if(transform=="white") stop("the white-noise (non-phylogenetic) model is not allowed")
            if ( (model == "drift") & (ape::is.ultrametric(phy))) stop("A drift model is unidentifiable for ultrametric trees., see ?fitContinuous for details")
            else
              
            #Matching tree
            full.data<-data
            phy<-phy
            
            #Calculates the full model, extracts model parameters
            N                   <- length(full.data)
            mod.0               <- geiger::fitDiscrete(phy = phy,dat = full.data,
                                                       model = model,transform = transform,
                                                       bounds = bounds,ncores = NULL,...)
            sigsq.0               <- mod.0$opt$sigsq
            z0.0                  <- mod.0$opt$z0
            aicc.0              <- mod.0$opt$aicc
            if (model == "BM"){
              optpar.0 <- NA
            }
            if (model == "OU"){
              optpar.0        <- mod.0$opt$alpha
            }
            if (transform == "EB"){
              optpar.0               <- mod.0$opt$a
            }
            if (transform == "trend"){
              optpar.0               <- mod.0$opt$slope
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
            if (transform == "drift"){
              optpar.0               <- mod.0$opt$drift
            }
            
            #Creates empty data frame to store model outputs
            sensi.estimates<-data.frame("species" =numeric(),
                                        "sigsq"=numeric(),"DIFsigsq"= numeric(),"sigsq.perc"= numeric(),
                                        "optpar"=numeric(),"DIFoptpar"=numeric(),"optpar.perc"=numeric(),
                                        "z0"=numeric(),
                                        "aicc"=numeric()) 
            
            #Loops over all species, and removes each one individually
            counter <- 1
            errors <- NULL
            if(track==TRUE) pb <- utils::txtProgressBar(min = 0, max = N, style = 3)
            for (i in 1:N){
              
              crop.data <- full.data[c(1:N)[-i]]
              crop.phy <-  ape::drop.tip(phy,setdiff(phy$tip.label,names(crop.data)))
              
              mod = try(geiger::fitDiscrete(phy = crop.phy,dat = crop.data,
                                            model = model,
                                            bounds = bounds,ncores = NULL,...),TRUE)
              if(isTRUE(class(mod)=="try-error")) {
                error <- i
                names(error) <- rownames(full.data$data)[i]
                errors <- c(errors,error)
                next }
              else {  sp                   <- phy$tip.label[i]
              sigsq               <- mod$opt$sigsq
              z0                  <- mod$opt$z0
              aicc              <- mod$opt$aicc
              DIFsigsq            <- sigsq - sigsq.0
              sigsq.perc          <- round((abs(sigsq / sigsq.0)) * 100,
                                         digits = 1)
              aicc              <- mod$opt$aicc
              if (transform == "none"){
                optpar <- NA
              }
              if (model == "OU"){
                optpar        <- mod$opt$alpha
              }
              if (transform == "EB"){
                optpar               <- mod$opt$a
              }
              if (transform == "trend"){
                optpar               <- mod$opt$slope
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
              if (transform == "drift"){
                optpar              <- mod$opt$drift
              }
              
              DIFoptpar            <- optpar - optpar.0
              optpar.perc        <- round((abs(optpar / optpar.0)) * 100,
                                           digits = 1)
              
              if(track==TRUE) utils::setTxtProgressBar(pb, i)
              # Stores values for each simulation
              # Store reduced model parameters: 
              estim.simu <- data.frame(sp, 
                                       sigsq, DIFsigsq,sigsq.perc,
                                       optpar,DIFoptpar,optpar.perc,
                                       z0,
                                       aicc,
                                       stringsAsFactors = F)
              sensi.estimates[counter, ]  <- estim.simu
              counter=counter+1
              }
            }
            if(track==TRUE) on.exit(close(pb))
            #Calculates Standardized DFbeta and DIFq12
            sDIFsigsq <- sensi.estimates$DIFsigsq/
              stats::sd(sensi.estimates$DIFsigsq)
            sDIFoptpar     <- sensi.estimates$DIFoptpar/
              stats::sd(sensi.estimates$DIFz0)
            
            sensi.estimates$sDIFsigsq     <- sDIFsigsq
            sensi.estimates$sDIFoptpar     <- sDIFoptpar
            
            #Creates a list with full model estimates:
            #full model estimates:
            param0 <- list(sigsq=sigsq.0,
                           optpar=optpar.0,
                           z0=z0.0,
                           aicc=aicc.0)
            
            #Identifies influencital species (sDF > cutoff) and orders by influence
            reorder.on.sigsq         <-sensi.estimates[order(abs(
              sensi.estimates$sDIFsigsq),decreasing=T),c("species","sDIFsigsq")]
            influ.sp.sigsq           <-as.character(reorder.on.sigsq$species[abs(
              reorder.on.sigsq$sDIFsigsq)>cutoff])
            reorder.on.optpar     <-sensi.estimates[order(abs(
              sensi.estimates$sDIFoptpar),decreasing=T),c("species","sDIFoptpar")]
            influ.sp.optpar       <-as.character(reorder.on.optpar$species[abs(
              reorder.on.optpar$sDIFoptpar)>cutoff])
            
            #Generates output:
            res <- list(   call = match.call(),
                           cutoff=cutoff,
                           data = full.data,
                           optpar = transform,
                           full.model.estimates = param0,
                           influential.species= list(influ.sp.sigsq=influ.sp.sigsq,
                                                     influ.sp.optpar=influ.sp.optpar),
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
          
