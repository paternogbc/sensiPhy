#' Interaction of intraspecific variability & influential species - Phylogenetic Linear Regression
#'
#' Performs leave-one-out deletion analyis for phylogenetic linear regression,
#' and detects influential species, while taking into account potential
#' interactions with intraspecific variability.
#'
#' @param formula The model formula:  \code{response~predictor}. 
#' @param data Data frame containing species traits with row names matching tips
#' in \code{phy}.
#' @param phy A phylogeny (class 'phylo') matching \code{data}.
#' @param model The phylogenetic model to use (see Details). Default is \code{lambda}.
#' @param cutoff The cutoff value used to identify for influential species
#' (see Details)
#' @param Vy Name of the column containing the standard deviation or the standard error of the response 
#' variable. When information is not available for one taxon, the value can be 0 or \code{NA}.
#' @param Vx Name of the column containing the standard deviation or the standard error of the predictor 
#' variable. When information is not available for one taxon, the value can be 0 or \code{NA}
#' @param y.transf Transformation for the response variable (e.g. \code{"log"} or \code{"sqrt"}). Please use this 
#' argument instead of transforming data in the formula directly (see also details below).
#' @param x.transf Transformation for the predictor variable (e.g. \code{"log"} or \code{"sqrt"}). Please use this 
#' argument instead of transforming data in the formula directly (see also details below).
#' @param distrib A character string indicating which distribution to use to generate a random value for the response 
#' and/or predictor variables. Default is normal distribution: "normal" (function \code{\link{rnorm}}).
#' Uniform distribution: "uniform" (\code{\link{runif}})
#' Warning: we recommend to use normal distribution with Vx or Vy = standard deviation of the mean.
#' @param track Print a report tracking function progress (default = TRUE)
#' @param ... Further arguments to be passed to \code{phylolm}
#' @details
#' This function fits a phylogenetic linear regression model using \code{\link[phylolm]{phylolm}}, and detects
#' influential species by sequentially deleting one at a time. The regression is repeated \code{times} times for 
#' simulated values of the dataset, taking into account intraspecific variation. At each iteration, the function 
#' generates a random value for each row in the dataset using the standard deviation or errors supplied, and 
#' detect the influential species within that iteration. 
#'
#' All phylogenetic models from \code{phylolm} can be used, i.e. \code{BM},
#' \code{OUfixedRoot}, \code{OUrandomRoot}, \code{lambda}, \code{kappa},
#' \code{delta}, \code{EB} and \code{trend}. See ?\code{phylolm} for details.
#'
#' \code{influ_phylm} detects influential species based on the standardised
#' difference in intercept and/or slope when removing a given species compared
#' to the full model including all species. Species with a standardised difference
#' above the value of \code{cutoff} are identified as influential. The default
#' value for the cutoff is 2 standardised differences change.
#'
#' Currently, this function can only implement simple linear models (i.e. \eqn{trait~
#' predictor}). In the future we will implement more complex models.
#'
#' Output can be visualised using \code{sensi_plot}.
#' 
#' @section Warning:  
#' When Vy or Vx exceed Y or X, respectively, negative (or null) values can be generated, this might cause problems
#' for data transformation (e.g. log-transformation). In these cases, the function will skip the simulation. This problem can
#' be solved by increasing \code{times}, changing the transformation type and/or checking the target species in output$sp.pb.
#' 
#' Setting \code{times} at high values can take a long time to exectue, since the total number of iterations equals \code{times * nrow(data)}.
#'
#' @return The function \code{interaction_intra_influ_phylm} returns a list with the following
#' components:
#' @return \code{cutoff}: The value selected for \code{cutoff}
#' @return \code{formula}: The model formula
#' @return \code{full.model.estimates}: Coefficients, aic and the optimised
#' value of the phylogenetic parameter (e.g. \code{lambda}) for the full model
#' without deleted species.
#' @return \code{influential_species}: List of influential species, both
#' based on standardised difference in interecept and in the slope of the
#' regression. Species are ordered from the species that was influential in most
#' of the iterations to the ones least commonly included influential,
#' only including species with a standardised difference > \code{cutoff}.
#' @return \code{influ.model.estimates}: A data frame with all simulation
#' estimates. Each row represents a deleted clade. #' Columns report the calculated
#' regression intercept (\code{intercept}), difference between simulation
#' intercept and full model intercept (\code{DFintercept}), the standardised
#' difference (\code{sDFintercept}), the percentage of change in intercept compared
#' to the full model (\code{intercept.perc}) and intercept p-value
#' (\code{pval.intercept}). All these parameters are also reported for the regression
#' slope (\code{DFslope} etc.). Additionally, model aic value (\code{AIC}) and
#' the optimised value (\code{optpar}) of the phylogenetic parameter
#' (e.g. \code{kappa} or \code{lambda}, depending on the phylogenetic model used) are
#' reported.
#' @return \code{data}: Original full dataset.
#' @return \code{errors}: Species where deletion resulted in errors.
#' @author Gustavo Paterno, Caterina Penone & Gijsbert D.A. Werner
#' @seealso \code{\link[phylolm]{phylolm}}, \code{\link{samp_phylm}},
#' \code{\link{influ_phylm}},\code{\link{sensi_plot}}
#' @references Ho, L. S. T. and Ane, C. 2014. "A linear-time algorithm for 
#' Gaussian and non-Gaussian trait evolution models". Systematic Biology 63(3):397-408.
#' @examples 
#' # Load data:
#' data(alien)
#' # run analysis:
#' influ <- interaction_intra_influ_phylm(formula = gestaLen ~ adultMass, phy = alien$phy[[1]],
#' data=alien$data,model="lambda",y.transf = "log",NULL,Vy="SD_gesta",Vx=NULL,times=10,
#' distrib = "normal")
#' # To check summary results:
#'summary(influ)
#'# Most influential speciesL
#'influ$influential.species
#'# Visual diagnostics
#'sensi_plot(influ)
#'# You can specify which graph and parameter ("slope" or "intercept") to print: 
#'sensi_plot(influ, param = "slope", graphs = 2)
#' @export

interaction_intra_influ_phylm <- function(formula,data,phy,model="lambda",cutoff=2,Vy = NULL, Vx = NULL,
                        y.transf = NULL, x.transf = NULL,
                        times = 10, distrib = "normal",
                        track=TRUE,...){
  
  #Error check
  if(is.null(Vx) & is.null(Vy)) stop("Vx or Vy must be defined")
  if(class(formula) != "formula") stop("formula must be class 'formula'")
  if(class(data) != "data.frame") stop("data must be class 'data.frame'")
  if(class(phy) != "phylo") stop("phy must be class 'phylo'")
  if(formula[[2]]!=all.vars(formula)[1] || formula[[3]]!=all.vars(formula)[2])
    stop("Please use arguments y.transf or x.transf for data transformation")
  if(distrib == "normal") warning ("distrib=normal: make sure that standard deviation 
                                   is provided for Vx and/or Vy")
  if(class(formula)!="formula") stop("formula must be class 'formula'")
  if(class(data)!="data.frame") stop("data must be class 'data.frame'")
  if(class(phy)!="phylo") stop("phy must be class 'phylo'")
  else
    
    
    #Matching tree and phylogeny using utils.R
    datphy <- match_dataphy(formula, data, phy)
    full.data <- datphy[[1]]
    phy <- datphy[[2]]
    
  #To loop over iterations of the dataset based on intraspecific variation. 
  resp <- all.vars(formula)[1]
  pred <- all.vars(formula)[2]
  
  if(!is.null(Vy) && sum(is.na(full.data[, Vy])) != 0) {
    full.data[is.na(full.data[, Vy]), Vy] <- 0}
  
  if(!is.null(Vx) && sum(is.na(full.data[, Vx])) != 0) {
    full.data[is.na(full.data[, Vx]), Vx] <- 0}
  
  #Function to pick a random value in the interval
  if (distrib == "normal") funr <- function(a,b) {stats::rnorm(1,a,b)}
  else  funr <- function(a,b) {stats::runif(1, a - b, a + b)}
  
  #Caculate the null model, i.e. no species deleted and no data uncertainty considered. 
  #transform if x.transf and/or y.transf are provided
  if(is.null(y.transf) & is.null(x.transf))
    {formula.0<-formula} ####This is what happens when there is no transformations. Solve for the other two cases too. 
  
  if(!is.null(y.transf) & is.null(x.transf)) #When there is only a transformation for the response variable
  {formula.0<-as.formula(paste0(y.transf,"(",resp,") ~ ",pred))}
  
  if(is.null(y.transf) & !is.null(x.transf)) #When there is only a transformation for the predictor variable
  {formula.0<-as.formula(paste0(resp," ~ ",x.transf,"(",pred,")"))}
  
  if(!is.null(y.transf) & !is.null(x.transf)) #When they both have a transformation.
  {formula.0<-as.formula(paste0(y.transf,"(",resp,") ~ ",x.transf,"(",pred,")"))}
  
  N               <- nrow(full.data)
  mod.0           <- phylolm::phylolm(formula.0, data=full.data,
                                      model=model,phy=phy)
  intercept.0      <- mod.0$coefficients[[1]]
  slope.0          <- mod.0$coefficients[[2]]
  pval.intercept.0 <- phylolm::summary.phylolm(mod.0)$coefficients[[1,4]]
  pval.slope.0     <- phylolm::summary.phylolm(mod.0)$coefficients[[2,4]]
  optpar.0 <- mod.0$optpar
  total_iteration <- N * times #I.e. how often are we going resimulate the dataset, times the # of species to drop. 
  
  
  #Creates empty data frame to store model outputs
  influ.model.estimates<-
    data.frame("species" =numeric(), "intercept"=numeric(),
               "DFintercept"=numeric(),"intercept.perc"=numeric(),
               "pval.intercept"=numeric(),"slope"=numeric(),
               "DFslope"=numeric(),"slope.perc"=numeric(),
               "pval.slope"=numeric(),"AIC"=numeric(),
               "optpar" = numeric())
  

  #Loops over all species, and removes each one individually
  counter <- 1
  errors <- NULL
  species.NA <- list()
  pb <- utils::txtProgressBar(min = 0, max = total_iteration, style = 1)
  
 #Create a nested for-loop. 
    for (i in 1:times) { #First create the new datset, and then drop all the species on that as previously. 
     
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
      {suppressWarnings (full.data$respV <- 
                           do.call(y.transf,list(x=full.data$respV)))}
      
      if(!is.null(x.transf)) 
      {suppressWarnings (full.data$predV <- 
                           do.call(x.transf,list(x=full.data$predV)))}
      
      #skip iteration if there are NA's in the dataset
      species.NA[[i]]<-rownames(full.data[with(full.data,is.na(predV) | is.na(respV)),])
      if(sum(is.na(full.data[,c("respV","predV")])>0)) next
      
      #Here, calculate the null-model for this particular resimulated dataset, 
      #i.e. no species deleted, but within this resimlated dataset / data unceratinty. 
      mod.0.resim           <- phylolm::phylolm(respV ~ predV, data=full.data,
                                          model=model,phy=phy)
      intercept.0.resim      <- mod.0.resim$coefficients[[1]]
      slope.0.resim          <- mod.0.resim$coefficients[[2]]
      pval.intercept.0.resim <- phylolm::summary.phylolm(mod.0.resim)$coefficients[[1,4]]
      pval.slope.0.resim     <- phylolm::summary.phylolm(mod.0.resim)$coefficients[[2,4]]
      optpar.0.resim <- mod.0.resim$optpar
      #Question, do we want to store these values too in the ultimate data frame for the user? 
      
          #Here, go into the species-drop loop:
          for (k in 1:N){
          crop.data <- full.data[c(1:N)[-k],]
          crop.phy <-  ape::drop.tip(phy,phy$tip.label[k])
          mod=try(phylolm::phylolm(respV ~ predV, data=crop.data,model=model,
                                   phy=crop.phy),
                  TRUE)
          if(isTRUE(class(mod)=="try-error")) {
            error <- k
            names(error) <- rownames(full.data$data)[k]
            errors <- c(errors,error)
            next }
          else {  sp                   <- phy$tip.label[k]
          intercept            <- mod$coefficients[[1]]
          slope                <- mod$coefficients[[2]]
          DFintercept          <- intercept - intercept.0.resim
          DFslope              <- slope - slope.0.resim
          intercept.perc       <- round((abs(DFintercept/intercept.0.resim))*100,digits=1)
          slope.perc           <- round((abs(DFslope/slope.0.resim))*100,digits=1)
          pval.intercept       <- phylolm::summary.phylolm(mod)$coefficients[[1,4]]
          pval.slope           <- phylolm::summary.phylolm(mod)$coefficients[[2,4]]
          aic.mod              <- mod$aic
          if (model == "BM" | model == "trend"){
            optpar <- NA
          }
          if (model != "BM" & model != "trend" ){
            optpar               <- mod$optpar
          }
          
          # Stores values for each simulation
          estim.simu <- data.frame(sp, intercept, DFintercept, intercept.perc,
                                   pval.intercept, slope, DFslope, slope.perc,
                                   pval.slope, aic.mod, optpar,
                                   stringsAsFactors = F)
          influ.model.estimates[counter, ]  <- estim.simu
          
          counter=counter+1
          
          if(track==TRUE)
            utils::setTxtProgressBar(pb, counter)
          
          }
      }
}

   
on.exit(close(pb))
  
  
  #Calculates Standardized DFbeta and DFintercept
  sDFintercept <- influ.model.estimates$DFintercept/
    stats::sd(influ.model.estimates$DFintercept)
  sDFslope     <- influ.model.estimates$DFslope/
    stats::sd(influ.model.estimates$DFslope)
  
  influ.model.estimates$sDFslope     <- sDFslope
  influ.model.estimates$sDFintercept <- sDFintercept
  
  
  
  #Creates a list with full model estimates:
  param0 <- list(coef=phylolm::summary.phylolm(mod.0)$coefficients,
                 aic=phylolm::summary.phylolm(mod.0)$aic,
                 optpar=mod.0$optpar)
  
  #Identifies influencital species (sDF > cutoff) and orders by influence
  reorder.on.slope         <-influ.model.estimates[order(abs(
    influ.model.estimates$sDFslope),decreasing=T),c("species","sDFslope")]
  influ.sp.slope           <-as.character(reorder.on.slope$species[abs(
    reorder.on.slope$sDFslope)>cutoff])
  reorder.on.intercept     <-influ.model.estimates[order(abs(
    influ.model.estimates$sDFintercept),decreasing=T),c("species","sDFintercept")]
  influ.sp.intercept       <-as.character(reorder.on.intercept$species[abs(
    reorder.on.intercept$sDFintercept)>cutoff])
  
  #Generates output:
  res <- list(call = match.call(),
              cutoff=cutoff,
              formula=formula,
              full.model.estimates=param0,
              influential.species= list(influ.sp.slope=influ.sp.slope,
                                        influ.sp.intercept=influ.sp.intercept),
              influ.model.estimates=influ.model.estimates,
              data=full.data,errors=errors)
  class(res) <- "sensiInflu_Intra"
  ### Warnings:
  if (length(res$errors) >0){
    warning("Some species deletion presented errors, please check: output$errors")}
  else {
    res$errors <- "No errors found."
  }
  
  return(res)
  
}

