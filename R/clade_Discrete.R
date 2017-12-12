#' Influential Clade Detection - Trait Evolution Discrete Characters
#' 
#' Fits models for trait evolution of discrete (binary) characters, 
#' detecting influential clades
#'
#' @param data Data frame containing species traits with row names matching tips
#' in \code{phy}.
#' @param phy A phylogeny (class 'phylo') matching \code{data}.
#' @param model The Mkn model to use (see Details). 
#' @param transform The evolutionary model to transform the tree (see Details). Default is \code{none}.
#' @param trait.col The column in the provided data frame which specifies the
#' trait to analyse (which should be a factor with two level)
#' @param clade.col The column in the provided data frame which specifies the
#' clades (a character vector with clade names).
#' @param n.species Minimum number of species in a clade for the clade to be
#' included in the leave-one-out deletion analysis. Default is \code{5}.
#' @param n.sim Number of simulations for the randomization test.
#' @param bounds settings to constrain parameter estimates. See \code{\link[geiger]{fitDiscrete}}
#' @param n.cores number of cores to use. If 'NULL', number of cores is detected.
#' @param track Print a report tracking function progress (default = TRUE)
#' @param ... Further arguments to be passed to \code{\link[geiger]{fitDiscrete}}
#' @details
#' This function sequentially removes one clade at a time,
#' fits different models of discrete character evolution using \code{\link[geiger]{fitDiscrete}}, 
#' repeats this this many times (controlled by \code{n.sim}), stores the results and calculates 
#' the effects on model parameters. Currently, only binary discrete traits are supported. 
#' 
#' Additionally, to account for the influence of the number of species on each 
#' clade (clade sample size), this function also estimates a null distribution
#' expected for the number of species in a given clade. This is done by fitting
#'  models without the same number of species as in the given clade.The number of 
#'  simulations to be performed is set by 'n.sim'. To test if the 
#'  clade influence differs from the null expectation for a clade of that size, 
#'  a randomization test can be performed using 'summary(x)'. 
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
#' @return \code{data}: The original full data frame. 
#' @return \code{full.model.estimates}: Parameter estimates (transition rates q12 and q21), 
#' AICc and the optimised value of the phylogenetic transformation parameter (e.g. \code{lambda}) 
#' for the full model without deleted clades.
#' @return \code{sensi.estimates}: Parameter estimates (transition rates q12 and q21), 
#' AICc and the optimised value of the phylogenetic transformation parameter (e.g. \code{lambda}) 
#' for each repeat with a clade removed.
#' @return \code{null.dist}: A data frame with estimates for the null distributions
#' for all clades analysed.
#' @return \code{errors}: Clades where deletion resulted in errors.
#' @return \code{clade.col}: Which column was used to specify the clades?
#' @return \code{optpar}: Transformation parameter used (e.g. \code{lambda}, \code{kappa} etc.)
#' @author Gijsbert Werner & Gustavo Paterno
#' @seealso \code{\link[geiger]{fitDiscrete}}
#' @references 
#' Yang Z. 2006. Computational Molecular Evolution. Oxford University Press: Oxford. 
#' 
#' Harmon Luke J, Jason T Weir, Chad D Brock, Richard E Glor, and Wendell Challenger. 2008.
#' GEIGER: investigating evolutionary radiations. Bioinformatics 24:129-131.
#' 
#' @examples 
#' \dontrun{
#' #Load data:
#' data("primates")
#' #Create a binary trait factor 
#' primates$data$adultMass_binary<-ifelse(primates$data$adultMass > 7350, "big", "small")
#' clade_disc<-clade_discrete(data=primates$data,phy = primates$phy[[1]],model="SYM",
#' trait.col = "adultMass_binary",clade.col="family",nsim=1,n.species=10,n.cores = 2)
#' summary(clade_disc)
#' sensi_plot(clade_disc)
#' sensi_plot(clade_disc, clade = "Cebidae", graph = "q12")
#' #Change the evolutionary model, tree transformation or minimum numher of species per clade
#' clade_disc_2<-clade_discrete(data=primates$data,phy = primates$phy[[1]],
#' model="ARD",transform="kappa",
#' trait.col = "adultMass_binary",clade.col="family",nsim=30,
#' n.species=8,n.cores = 2)
#' summary(clade_disc_2)
#' sensi_plot(clade_disc_2)
#' sensi_plot(clade_disc_2, graph = "q12")
#' sensi_plot(clade_disc_2, graph = "q21")
#' }
#' @export

clade_discrete <- function(data, phy, model,transform = "none",
                           trait.col,clade.col,n.species = 5, n.sim = 20,
                           bounds = list(), n.cores = NULL,track=TRUE, ...) {
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
                                                     bounds = bounds,ncores = n.cores,...)
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
                                          bounds = bounds,ncores = n.cores,...),TRUE)
            q12               <- mod$opt$q12
            q21               <- mod$opt$q21
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
                                            bounds = bounds,ncores = n.cores,...),TRUE)
              
              if(isTRUE(class(mod)=="try-error")) {
                error <- i
                names(error) <- rownames(full.data$data)[i]
                errors <- c(errors,error)
                next }
              else 
                
                q12               <- mod$opt$q12
              q21               <- mod$opt$q21
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
        
