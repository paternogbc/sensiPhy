#' Sensitivity Analysis Species Sampling  - Phylogenetic signal
#'
#' Performs analyses of sensitivity to species sampling by randomly removing
#' species and detecting the effects on phylogenetic 
#' signal estimates
#'
#' @param trait.col The name of a column in the provided data frame with trait 
#'  to be analyzed  (e.g. "Body_mass").
#' @param data Data frame containing species traits with row names matching tips
#' in \code{phy}.
#' @param phy A phylogeny (class 'phylo') matching \code{data}.
#' @param n.sim The number of times to repeat species random removal for each
#' \code{break} interval.
#' @param breaks A vector containing the percentages of species to remove.
#' @param method Method to compute signal: can be "K" or "lambda".
#' @param track Print a report tracking function progress (default = TRUE)
#' @param ... Further arguments to be passed to \code{\link[phytools]{phylosig}}
#' 
#' @details
#'
#' This function randomly removes a given percentage of species (controlled by
#' \code{breaks}) from the full data, estimates phylogenetic
#' signal for a given trait (K or lambda) without these species using 
#' \code{\link[phytools]{phylosig}}, then
#' repeats the analysis many times (controlled by \code{n.sim}), stores the results and
#' calculates the effect of random species removal on phylogenetic siganl estimates.
#'
#' Output can be visualised using \code{sensi_plot}.
#' @return The function \code{samp_phylosig} returns a list with the following
#' components:
#' @return \code{Trait}: Column name of the trait analysed
#' @return \code{full.model.estimates}: Phylogenetic signal (K or lambda) and 
#' p-value  using the full dataset (without deleted species). See 
#' \code{\link[phytools]{phylosig}} for details.
#' @return \code{sensi.estimates}: A data frame with all simulation
#' estimates. Each row represents a rerun with a given number of species
#' \code{n.remov} removed, representing \code{n.percent} of the full dataset.
#' Columns report the calculated signal estimate (\code{estimate}),
#' difference between reduced data signal estimate and full data signal (\code{DF}),
#' the percentage of change in signal compared to the full data estimate (\code{perc})
#' and signal p-value for the reduced data estimate(\code{pval}). 
#' @return \code{sign.analysis} For each break (i.e. each percentage of species
#' removed) this reports the percentage of statistically signficant (at p<0.05)
#' phylogenetic signal over all repititions with reduced data sets.
#' @return \code{data}: Original full dataset used in the analysis.
#' #' @note Please be aware that dropping species may reduce power to detect 
#' significant signal and may partially be responsible for a potential 
#' effect of species removal on p-values. Please also consult standardised differences
#' in the (summary) output.
#' @author Gustavo Paterno & Gijsbert D.A. Werner
#' @seealso
#' \code{\link[phytools]{phylosig}}, \code{\link{influ_phylosig}},\code{\link{sensi_plot}}
#' @references 
#' 
#' Werner, G.D.A., Cornwell, W.K., Sprent, J.I., Kattge, J. & Kiers, E.T. (2014).
#'  A single evolutionary innovation drives the deep evolution of symbiotic N2-fixation
#'   in angiosperms. Nature Communications, 5, 4087.
#'   
#' Blomberg, S. P., T. Garland Jr., A. R. Ives (2003) 
#' Testing for phylogenetic signal in comparative data: 
#' Behavioral traits are more labile. Evolution, 57, 717-745.
#' 
#' Pagel, M. (1999) Inferring the historical patterns of biological evolution. 
#' Nature, 401, 877-884.
#' 
#' Kamilar, J. M., & Cooper, N. (2013). Phylogenetic signal in primate behaviour,
#'  ecology and life history. Philosophical Transactions of the Royal Society 
#'  B: Biological Sciences, 368: 20120341.
#'  
#' @importFrom phytools phylosig
#' 
#' @examples 
#'data(alien)
#'# Logtransform data
#'alien.data$logMass <- log(alien.data$adultMass) 
#'# Run sensitivity analysis:
#'samp <- samp_physig(trait.col = "logMass", data = alien.data, n.sim = 100, phy = alien.phy[[1]])
#'summary(samp)
#'sensi_plot(samp)
#'sensi_plot(samp, graphs = 1)
#'sensi_plot(samp, graphs = 2)
#' @export
samp_physig <- function(trait.col, data , phy, n.sim = 30,
                        breaks=seq(.1,.5,.1), method = "K", track = TRUE, ...){
  ### data prep------
  # Basic error checking:
  if(class(phy) != "phylo") 
    stop("phy must be class 'phylo'")
  if(length(breaks) < 2) 
    stop("Please include more than one break, e.g. breaks=c(.3,.5)")
  
  # Check match between data and phy 
  data_phy <- match_dataphy(get(trait.col) ~ 1, data, phy)
  full.data <- data_phy$data
  phy <- data_phy$phy
  trait <- full.data[, trait.col]
  names(trait) <- phy$tip.label
  N             <- nrow(full.data)
  
  # FULL MODEL PARAMETERS:
  mod.0     <- phytools::phylosig(x = trait, tree = phy, method = method, test = TRUE, ...)
  
  #Creates empty data frame to store model outputs
  samp.physig.estimates <-
    data.frame("n.remov" = numeric(), "n.percent"= numeric(),
               "estimate"= numeric(),"DF"= numeric(),
               "perc"= numeric(),"pval"=numeric())
  
  #Loops over breaks, remove percentage of species determined by 'breaks
  #and repeat determined by 'n.sim'.
  counter <- 1
  limit <- sort(round( (breaks) * nrow(full.data),digits=0))
  NL <- length(breaks) * n.sim
  if(track==TRUE) pb <- utils::txtProgressBar(min = 0, max = NL, style = 1)
  
  
  ##### Loop----
  for (i in limit){
    for (j in 1:n.sim){
      
      exclude <- sample(1:N,i)
      crop.data <- full.data[-exclude,]
      crop.phy <-  ape::drop.tip(phy,phy$tip.label[exclude])
      crop.trait <- crop.data[, trait.col]
      names(crop.trait) <- crop.phy$tip.label
      
      ### Reduced data estimate
      mod.s     <- phytools::phylosig(x = crop.trait, tree = crop.phy, method = method, 
                                      test = TRUE, ...)
      
      ### Metrics:
      estimate  <- mod.s[[1]]
      DF        <- mod.s[[1]] - mod.0[[1]]
      perc      <- round((abs(DF/mod.0[[1]]))*100,digits=1)
      pval    <- mod.s$P
     
      
      n.remov <- i
      n.percent <- round( (n.remov / N) * 100,digits = 0)
      #rep <- j
      
      if(track == TRUE) (
        utils::setTxtProgressBar(pb, counter))
      # Stores values for each simulation
      estim.simu <- data.frame(n.remov, n.percent, estimate, 
                               DF, perc,
                               pval, stringsAsFactors = F)
      samp.physig.estimates[counter, ]  <- estim.simu
      counter <- counter + 1
    }
  }
  
  close(pb)
  
  #Calculates Standardized DFestimate
  sDF <- samp.physig.estimates$DF/
    stats::sd(samp.physig.estimates$DF)
  
  samp.physig.estimates$sDF <- sDF
  samp.physig.estimates <- samp.physig.estimates[, c(1,2,3,4,7,5,6)]
  #colnames(samp.physig.estimates)[3] <- method
  
  res                 <- samp.physig.estimates
  n.sim               <- table(res$n.remov)
  breaks              <- unique(res$n.percent)
  
  ### Significance Table
    sign                <- res$pval > .05
    res$sign            <- sign
    perc.sign           <- 1-(with(res,tapply(sign, n.remov, sum))) / n.sim
    mean.sDF            <- with(res,tapply(sDF,n.remov,mean))
    mean.perc           <- with(res,tapply(perc,n.remov,mean))
    perc.sign.tab       <- data.frame(percent_sp_removed=breaks,
                                      perc.sign = as.numeric(perc.sign),
                                      mean.perc = as.numeric(mean.perc),
                                      mean.sDF = as.numeric(mean.sDF))
  
  #Creates a list with full model estimates:
  param0 <- list(estimate = mod.0[[1]],
                 Pval = mod.0$P)
  names(param0)[1] = method
  
  #Generates output:
  res <- list(call = match.call(),
              Trait = trait.col,
              full.data.estimates = param0,
              sensi.estimates = samp.physig.estimates,
              sign.analysis = perc.sign.tab,
              data = full.data,
              phy = phy)
  class(res) <- "samp.physig"
  return(res)
}  

