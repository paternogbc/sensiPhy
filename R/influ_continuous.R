#' Influential Species Detection - Trait Evolution Continous Characters
#'
#' Fits models for trait evolution of continuous characters, 
#' detecting influential species. 
#'
#' @param data Data vector for a single continuous trait, with names matching tips in \code{phy}.
#' @param phy A phylogeny (class 'phylo') matching \code{data}.
#' @param model The evolutionary model (see Details). 
#' @param cutoff The cut-off parameter for influential species (see Details). 
#' @param bounds settings to contstrain parameter estimates. See \code{\link[geiger]{fitContinuous}}
#' @param track Print a report tracking function progress (default = TRUE)
#' @param ... Further arguments to be passed to \code{\link[geiger]{fitContinuous}}
#' @details
#' This function sequentially removes one species at a time,  
#' fits different models of continuous character evolution using \code{\link[geiger]{fitContinuous}}, 
#' repeats this this many times (controlled by \code{n.sim}), stores the results and calculates 
#' the effects on model parameters.
#' 
#' #' \code{influ_continuous} detects influential species based on the standardised
#' difference in the rate parameter \code{sigsq} and the optimisation parameter \code{optpar} 
#' (e.g. lamda, kappa, alpha, depending on which \code{model} is set), when removing 
#' a given species compared to the full model including all species. 
#' Species with a standardised difference above the value of 
#' \code{cutoff} are identified as influential. 
#'
#' All evolutionary models from \code{fitContinuous} can be used, i.e. \code{BM},\code{OU},
#' \code{EB}, \code{trend}, \code{lambda}, \code{kappa}, \code{delta} and \code{drift}.
#' 
#' See \code{\link[geiger]{fitContinuous}} for more details on evolutionary models. 
#' 
#' Output can be visualised using \code{sensi_plot}. [Not yet!]
#'
#' @return The function \code{tree_discrete} returns a list with the following
#' components:
#' @return \code{call}: The function call
#' @return \code{cutoff}: The value selected for \code{cutoff}
#' @return \code{data}: The original full data vector
#' @return \code{optpar}: Transformation parameter used (e.g. \code{lambda}, \code{kappa} etc.)
#' @return \code{full.model.estimates}: Parameter estimates (rate of evolution \code{sigsq}, 
#' root state \code{z0} and where applicable \code{optpar}), 
#' AICc and the optimised value of the phylogenetic transformation parameter (e.g. \code{lambda}) 
#' for the full model without deleted species.
#' @return \code{influential_species}: List of influential species, based on standardised 
#' difference in estimates for sigsq and optpar. Species are ordered from most influential to 
#' less influential and only include species with a standardised difference > \code{cutoff}.
#' @return \code{sensi.estimates}: Parameter estimates, AICc and the optimised value of 
#' the phylogenetic transformation parameter (e.g. \code{lambda}) for each analysis 
#' with a different phylogenetic tree.
#' @author Gijsbert Werner & Gustavo Paterno
#' @seealso \code{\link[geiger]{fitContinuous}}
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
#' influ_cont<-influ_continuous(data = adultMass,phy = primates$phy[[1]],
#' model = "OU",cutoff = 2,track = T)
#' #Print summary statistics for the transitions rates, aic-values and (if applicable) optimisation parameter
#' summary(influ_cont)
#' #Use a different evolutionary model or cutoff 
#' influ_cont2<-influ_continuous(data = adultMass,phy = primates$phy[[1]],
#' model = "lambda",cutoff = 1.2,track = T)
#' summary(influ_cont2)
#' influ_cont3<-influ_continuous(data = adultMass,phy = primates$phy[[1]],
#' model = "BM",cutoff = 2,track = T)
#' summary(influ_cont3)
#' @export

influ_continuous <- function(data,phy,model,
                           bounds = list(),
                           cutoff=2,track=TRUE,...){
  
            #Error check
            if(is.null(model)) stop("model must be specified, e.g. 'OU' or 'lambda'")
            if(class(data)!="numeric" | is.null(names(data))) stop("data must supplied as a numeric vector with species as names")
            if(class(phy)!="phylo") stop("phy must be class 'phylo'")
            if(model=="white") stop("the white-noise (non-phylogenetic) model is not allowed")
            if ( (model == "drift") & (ape::is.ultrametric(phy))) stop("A drift model is unidentifiable for ultrametric trees., see ?fitContinuous for details")
            else
              
            #Matching tree
            full.data<-data
            phy<-phy
            
            #Calculates the full model, extracts model parameters
            N                   <- length(full.data)
            mod.0               <- geiger::fitContinuous(phy = phy,dat = full.data,
                                                       model = model,
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
              
              mod = try(geiger::fitContinuous(phy = crop.phy,dat = crop.data,
                                            model = model,
                                            bounds = bounds,ncores = NULL,...),TRUE)
              if(isTRUE(class(mod)=="try-error")) {
                error <- i
                names(error) <- rownames(full.data$data)[i]
                errors <- c(errors,error)
                next }
              else { 
              sp                   <- phy$tip.label[i]
              sigsq               <- mod$opt$sigsq
              z0                  <- mod$opt$z0
              aicc              <- mod$opt$aicc
              DIFsigsq            <- sigsq - sigsq.0
              sigsq.perc          <- round((abs(DIFsigsq / sigsq.0)) * 100,
                                         digits = 1)
              aicc              <- mod$opt$aicc
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
            if(stats::sd(sensi.estimates$DIFoptpar)==0){
              sDIFoptpar<-NA
            } else{
            sDIFoptpar     <- sensi.estimates$DIFoptpar/
              stats::sd(sensi.estimates$DIFoptpar)}
            
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

            if(model=="BM"){
              influ.sp.optpar<-"No optpar calculated for BM-model. Influential species not calculated"
            }  
            if(stats::sd(sensi.estimates$DIFoptpar)==0){
              influ.sp.optpar<-"No variation in optpar. Influential species not calculated"
            } 
            else{
            reorder.on.optpar     <-sensi.estimates[order(abs(
              sensi.estimates$sDIFoptpar),decreasing=T),c("species","sDIFoptpar")]
            influ.sp.optpar       <-as.character(reorder.on.optpar$species[abs(
              reorder.on.optpar$sDIFoptpar)>cutoff])
            }
            
            #Generates output:
            res <- list(   call = match.call(),
                           cutoff=cutoff,
                           data = full.data,
                           optpar = model,
                           full.model.estimates = param0,
                           influential.species= list(influ.sp.sigsq=influ.sp.sigsq,
                                                     influ.sp.optpar=influ.sp.optpar),
                           sensi.estimates=sensi.estimates,
                           errors = errors)
            class(res) <- "sensiInflu.TraitEvol"
            ### Warnings:
            if (length(res$errors) >0){
              warning("Some species deletion presented errors, please check: output$errors")}
            else {
              res$errors <- "No errors found."
            }
            
            return(res)
            
          }
          