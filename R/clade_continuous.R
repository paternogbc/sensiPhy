#' Influential Clade Detection - Trait Evolution Continuous Characters
#' 
#' Fits models for trait evolution of continuous characters, 
#' detecting influential clades
#' 
#' @param data Data frame containing species traits with row names matching tips
#' in \code{phy}.
#' @param phy A phylogeny (class 'phylo') matching \code{data}.
#' @param model The evolutionary model (see Details). 
#' @param trait.col The column in the provided data frame which specifies the
#' trait to analyse (which should be a factor with two level)
#' @param clade.col The column in the provided data frame which specifies the
#' clades (a character vector with clade names).
#' @param n.species Minimum number of species in a clade for the clade to be
#' included in the leave-one-out deletion analysis. Default is \code{5}.
#' @param n.sim Number of simulations for the randomization test.
#' @param bounds settings to contstrain parameter estimates. See \code{\link[geiger]{fitContinuous}}
#' @param track Print a report tracking function progress (default = TRUE)
#' @param ... Further arguments to be passed to \code{\link[geiger]{fitContinuous}}
#' @details
#' #' This function sequentially removes one clade at a time,
#' fits different models of continuous character evolution using \code{\link[geiger]{fitContinuous}}, 
#' repeats this this many times (controlled by \code{n.sim}), stores the results and calculates 
#' the effects on model parameters Currently, only binary continuous traits are supported. 
#' 
#' Additionally, to account for the influence of the number of species on each 
#' clade (clade sample size), this function also estimates a null distribution
#' expected for the number of species in a given clade. This is done by fitting
#'  models without the same number of species as in the given clade.The number of 
#'  simulations to be performed is set by 'n.sim'. To test if the 
#'  clade influence differs from the null expectation for a clade of that size, 
#'  a randomization test can be performed using 'summary(x)'. 
#'
#' Different evolutionary models from \code{fitContinuous} can be used, i.e. \code{BM},\code{OU},
#' \code{EB}, \code{trend}, \code{lambda}, \code{kappa}, \code{delta} and \code{drift}.
#' 
#' See \code{\link[geiger]{fitContinuous}} for more details on evolutionary models. 
#' 
#' Output can be visualised using \code{sensi_plot}.
#'
#' @return The function \code{tree_continuous} returns a list with the following
#' components:
#' @return \code{call}: The function call
#' @return \code{data}: The original full data frame. 
#' @return \code{full.model.estimates}: Parameter estimates (rate of evolution \code{sigsq}, 
#' root state \code{z0} and where applicable \code{optpar}), 
#' AICc and the optimised value of the phylogenetic transformation parameter (e.g. \code{lambda}) 
#' for the full model without deleted clades.
#' @return \code{sensi.estimates}: Parameter estimates, 
#' AICc and the optimised value of the phylogenetic transformation parameter (e.g. \code{lambda}) 
#' for each repeat with a clade removed.
#' @return \code{null.dist}: A data frame with estimates for the null distributions
#' for all clades analysed.
#' @return \code{errors}: Clades where deletion resulted in errors.
#' @return \code{clade.col}: Which column was used to specify the clades?
#' @return \code{optpar}: Transformation parameter used (e.g. \code{lambda}, \code{kappa} etc.)
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
#' clade_cont<-clade_continuous(data=primates$data,phy = primates$phy[[1]],model="OU",
#' trait.col = "adultMass",clade.col="family",n.sim=10,n.species=10,track=TRUE)
#' #Print summary statistics for the transitions rates, aic-values and (if applicable) optimisation parameter
#' summary(clade_cont)
#' #Change the evolutionary model, tree transformation or minimum numher of species per clade
#' clade_cont2<-clade_continuous(data=primates$data,phy = primates$phy[[1]],model="BM",
#' trait.col = "adultMass",clade.col="family",n.sim=10,n.species=5,track=TRUE)
#' summary(clade_cont2)
#' 
#' @export

clade_continuous <- function(data, phy, model,
                           trait.col,clade.col,n.species = 5, n.sim = 20,
                           bounds = list(), track=TRUE, ...) {
          # Error checking:
          if(is.null(model)) stop("model must be specified, e.g. 'OU' or 'lambda'")
          if(!is.data.frame(data)) stop("data must be class 'data.frame'")
          if(missing(clade.col)) stop("clade.col not defined. Please, define the",
                                      " column with clade names.")
          if(class(phy)!="phylo") stop("phy must be class 'phylo'")
          if(model=="white") stop("the white-noise (non-phylogenetic) model is not allowed")
          if ( model == "drift" & ape::is.ultrametric(phy)) stop("A drift model is unidentifiable for ultrametric trees., see ?fitContinuous for details")
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
          names(trait_vec_full)<-rownames(full.data)
          
          N                   <- nrow(full.data)
          mod.0               <- geiger::fitContinuous(phy = phy,dat = trait_vec_full,
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
          
          #Create dataframe to store estmates for each clade
          sensi.estimates<-data.frame("clade" =I(as.character()),"N.species" = numeric(),
                                      "sigsq"=numeric(),"DIFsigsq"= numeric(),"sigsq.perc"= numeric(),
                                      "optpar"=numeric(),"DIFoptpar"=numeric(),"optpar.perc"=numeric(),
                                      "z0"=numeric(),
                                      "aicc"=numeric()) 
          
          # Create dataframe store simulations (null distribution)
          null.dist <- data.frame("clade" = rep(names(uc), each = n.sim),
                                  "sigsq"= numeric(length(uc)*n.sim),
                                  "DIFsigsq"= numeric(length(uc)*n.sim),
                                  "optpar" = numeric(length(uc)*n.sim),
                                  "DIFoptpar" = numeric(length(uc)*n.sim))

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
            names(crop.trait_vec)<-rownames(crop.data)
            mod = try(geiger::fitContinuous(phy = crop.phy,dat = crop.trait_vec,
                                          model = model,
                                          bounds = bounds,ncores = NULL,...),TRUE)
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
            
            # Store reduced model parameters: 
            estim.simu <- data.frame(A, cN,
                                     sigsq, DIFsigsq,sigsq.perc,
                                     optpar,DIFoptpar,optpar.perc,
                                     z0,
                                     aicc,
                                     stringsAsFactors = F)
            sensi.estimates[aa, ]  <- estim.simu
            
            ### START LOOP FOR NULL DIST:
            # number of species in clade A:
            for (i in 1:n.sim) {
              exclude <- sample(1:N, cN)
              crop.data <- full.data[-exclude,]
              crop.phy <-  ape::drop.tip(phy,setdiff(phy$tip.label,rownames(crop.data)))
              crop.trait_vec<-crop.data[,trait.col]
              names(crop.trait_vec)<-rownames(crop.data)
              mod = try(geiger::fitContinuous(phy = crop.phy,dat = crop.trait_vec,
                                            model = model,
                                            bounds = bounds,ncores = NULL,...),TRUE)
              
              if(isTRUE(class(mod)=="try-error")) {
                error <- i
                names(error) <- rownames(full.data$data)[i]
                errors <- c(errors,error)
                next }
              else 
              sigsq               <- mod$opt$sigsq
              aicc              <- mod$opt$aicc
              DIFsigsq            <- sigsq - sigsq.0
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
              
              null.dist[bb, ] <- data.frame(clade = as.character(A), 
                                            sigsq,
                                            DIFsigsq,
                                            optpar,
                                            DIFoptpar)
              
              if(track==TRUE) utils::setTxtProgressBar(pb, bb)
              bb <- bb + 1
            }
            aa <- aa + 1
          }
          if(track==TRUE) on.exit(close(pb))
          
          #OUTPUT
          #full model estimates:
          param0 <- list(sigsq=sigsq.0,optpar=optpar.0,
                         optpar=optpar.0,
                         aicc=aicc.0)
          
          #Generates output:
          res <- list(   call = match.call(),
                         data = full.data,
                         full.model.estimates = param0,
                         sensi.estimates=sensi.estimates,
                         null.dist = null.dist,
                         errors = errors,
                         optpar = model,
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
        
