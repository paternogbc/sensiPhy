#' Influential species detection - Phylogenetic signal
#'
#' Performs leave-one-out deletion analysis for phylogenetic signal estimates,
#' and detects influential species for K or lambda.
#' @param trait.col The name of a column in the provided data frame with trait 
#'  to be analyzed  (e.g. "Body_mass").
#' @param data Data frame containing species traits with row names matching tips
#' in \code{phy}.
#' @param phy A phylogeny (class 'phylo') matching \code{data}.
#' @param method Method to compute signal: can be "K" or "lambda".
#' @param cutoff The cutoff value used to identify for influential species
#' (see Details)
#' @param track Print a report tracking function progress (default = TRUE)
#' @param ... Further arguments to be passed to \code{\link[phytools]{phylosig}}
#' 
#' @details
#' This function sequentially removes one species at a time, ans estimates phylogenetic
#' signal (K or lambda) using \code{\link[phytools]{phylosig}}, stores the
#' results and detects the most influential species.
#'
#' \code{influ_physig} detects influential species based on the standardised
#' difference in signal estimate (K or lambda) when removing a given species compared
#' to the full data estimate (with all species). Species with a standardised difference
#' above the value of \code{cutoff} are identified as influential. The default
#' value for the cutoff is 2 standardised differences in signal estimate.
#'
#' Output can be visualised using \code{sensi_plot}.
#'
#' @return The function \code{influ_physig} returns a list with the following
#' components:
#' @return \code{cutoff}: The value selected for \code{cutoff}
#' @return \code{trait.col}: Column name of the trait analysed
#' @return \code{full.data.estimates}: Phylogenetic signal estimate (K or lambda)
#' and the P value (for the full data).
#' @return \code{influential_species}: List of influential species,
#' based on standardised difference in K or lambda. 
#' Species are ordered from most influential to less influential and
#' only include species with a standardised difference > \code{cutoff}.
#' @return \code{influ.physig.estimates}: A data frame with all simulation
#' estimates. Each row represents a deleted species 
#' Columns report the calculated signal estimate (\code{k}) or (\code{lambda}),
#'  difference between signal estimation of the reduced and full data 
#'  (\code{DF}), the percentage of change in signal compared
#' to the full data signal (\code{perc}) and p-value for the phylogenetic signal
#' test (\code{pval})
#' @return \code{data}: Original full dataset.
#' @author Gustavo Paterno
#' @seealso \code{\link[phytools]{phylosig}},
#' \code{\link{influ_phylm}},\code{\link{sensi_plot}}
#' @references 
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
#' @examples 
#' \dontshow{
#'# Load data:
#'data(alien)
#'# Logtransform data
#'alien.data$logMass <- log(alien.data$adultMass) 
#'# Run sensitivity analysis:
#'influ <- influ_physig("logMass", data = alien.data[1:20,],
#' phy = alien.phy[[1]])
#'# To check summary results:
#'summary(influ)
#' }
#' \dontrun{
#'# Load data:
#'data(alien)
#'# Logtransform data
#'alien.data$logMass <- log(alien.data$adultMass) 
#'# Run sensitivity analysis:
#'influ <- influ_physig("logMass", data = alien.data, phy = alien.phy[[1]])
#'# To check summary results:
#'summary(influ)
#'# Most influential speciesL
#'influ$influential.species
#'# Visual diagnostics
#'sensi_plot(influ)
#'# You can specify which graph to print: 
#'sensi_plot(influ, graphs = 1)
#'sensi_plot(influ, graphs = 2)
#'}
#' @export
influ_physig <- function(trait.col, data, phy, method = "K", cutoff = 2, track = TRUE, ...){
  
  ### Basic checking
  if(!inherits(phy,"phylo"))
    stop("tree should be an object of class \"phylo\".")
  if(!inherits(trait.col,"character")) 
    stop("trait.col should be an object of class \"character\".")
  
  #Matching tree and phylogeny using utils.R
  datphy <- match_dataphy(get(trait.col) ~ 1, data, phy)
  full.data <- datphy$data
  phy <- datphy$phy
  trait     <- full.data[[trait.col]]
  names(trait)  <- phy$tip.label
  N <- nrow(full.data)
  
  ### Fit full data model
  mod.0     <- phytools::phylosig(x = trait, tree = phy, method = method, test = T, ...)
  
  #Creates empty data frame to store model outputs
  influ.physig.estimates<-
    data.frame("species" =numeric(), "estimate" = numeric(),
               "DF"=numeric(),"perc"=numeric(),
               "pval"=numeric())
  
  ### Start leave-one-out analysis (between all species)
  counter <- 1
  if(track==TRUE) pb <- utils::txtProgressBar(min = 0, max = N, style = 3)
  
  for (i in 1:N){
    # Crop data: remove one species from data and phy:
    crop.data <- trait[-i]
    crop.phy  <-  ape::drop.tip(phy,phy$tip.label[i])
    # Fit reduced data model
    mod.s     <- phytools::phylosig(tree = crop.phy, crop.data, method = method, test = T, ...)
    
    estimate  <- mod.s[[1]]
    sp        <- phy$tip.label[i]
    DF        <- mod.s[[1]] - mod.0[[1]]
    perc      <- round((abs(DF/mod.0[[1]]))*100,digits=1)
    pval      <- mod.s$P
    
    # Stores values for each simulation
    estim.simu <- data.frame(sp, estimate, DF, perc, pval,
                             stringsAsFactors = F)
    
    influ.physig.estimates[counter, ]  <- estim.simu
    counter=counter+1
    if(track==TRUE) (utils::setTxtProgressBar(pb, i))
    
  }
  
  if(track==TRUE) on.exit(close(pb))
  
  
  #Calculates Standardized DFbeta and DFintercept
  sDF <- influ.physig.estimates$DF/
    stats::sd(influ.physig.estimates$DF)
  
  influ.physig.estimates$sDF <- sDF
  
  #Creates a list with full dataset estimates:
  param0 <- list(estimate = mod.0[[1]],
                 pval = mod.0$P)
  
  #Identifies influencital species (sDF > cutoff) and orders by influence
  reorder  <- influ.physig.estimates[order(abs(
    influ.physig.estimates$sDF),decreasing=T),c("species","sDF")]
  influ.sp <- as.character(reorder$species[abs(reorder$sDF)>cutoff]) 
  
  #Generates output:
  res <- list(
    call = match.call(),
    cutoff=cutoff,
    trait.col=trait.col,
    full.data.estimates = param0,
    influential.species = influ.sp,
    influ.physig.estimates = influ.physig.estimates,
    data = full.data,
    phy = phy)
  class(res) <- "influ.physig"
  return(res)
}
