#' Influential clade detection - Phylogenetic signal 
#'
#' Estimate the influence of clade removal on phylogenetic signal estimates
#'
#' @param trait.col The name of a column in the provided data frame with trait 
#'  to be analyzed  (e.g. "Body_mass").
#' @param data Data frame containing species traits with row names matching tips
#' in \code{phy}.
#' @param phy A phylogeny (class 'phylo') matching \code{data}.
#' @param method Method to compute signal: can be "K" or "lambda".
#' @param clade.col The column in the provided data frame which specifies the
#' clades (a character vector with clade names).
#' @param n.species Minimum number of species in a clade for the clade to be
#' included in the leave-one-out deletion analysis. Default is \code{5}.
#' @param n.sim Number of simulations for the randomization test.
#' @param track Print a report tracking function progress (default = TRUE)
#' @param ... Further arguments to be passed to \code{\link[phytools]{phylosig}}
#' 
#' @details
#' This function sequentially removes one clade at a time, estimates phylogenetic
#' signal (K or lambda) using \code{\link[phytools]{phylosig}} and stores the
#' results. The impact of a specific clade on signal estimates is calculated by the
#' comparison between the full data (with all species) and reduced data estimates
#' (without the species belonging to a clade).
#' 
#' To account for the influence of the number of species on each 
#' clade (clade sample size), this function also estimate a null distribution of signal estimates
#' expected by the removal of the same number of species as in a given clade. This is done by estimating
#' phylogenetic signal without the same number of species in the given clade. 
#'  The number of simulations to be performed is set by 'n.sim'. To test if the 
#'  clade influence differs from the null expectation for a clade of that size, 
#'  a randomization test can be performed using 'summary(x)'. 
#'
#' Output can be visualised using \code{sensi_plot}.
#'
#' @return The function \code{clade_physig} returns a list with the following
#' components:
#' @return \code{trait.col}: Column name of the trait analysed
#' @return \code{full.data.estimates}: Phylogenetic signal estimate (K or lambda)
#' and the P value (for the full data).
#' @return \code{sensi.estimates}: A data frame with all simulation
#' estimates. Each row represents a deleted clade. Columns report the calculated
#' phylogenetic signal (K or lambda) (\code{estimate}), difference between simulation
#' signal and full data signal (\code{DF}), the percentage of change
#' in signal compared to the full data estimate (\code{perc}) and 
#' p-value of the phylogenetic signal with the reduced data (\code{pval}).
#' @return \code{null.dist}: A data frame with estimates for the null distribution
#' of phylogenetic signal for all clades analysed.
#' @return \code{data}: Original full dataset.
#' @author Gustavo Paterno
#' 
#' @seealso \code{\link[phytools]{phylosig}}, 
#' \code{\link{clade_phylm}},\code{\link{sensi_plot}}
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
#' @examples 
#'data(alien)
#'# Logtransform data
#'alien.data$logMass <- log(alien.data$adultMass) 
#'# Run sensitivity analysis:
#'clade <- clade_physig(trait.col = "logMass", data = alien.data, n.sim = 20,
#'                  phy = alien.phy[[1]], clade.col = "family", method = "K")
#'summary(clade)
#'sensi_plot(clade, "Bovidae")
#'sensi_plot(clade, "Sciuridae")
#' @export

clade_physig <- function(trait.col, data, phy, clade.col, n.species = 5, n.sim = 100, method = "K",  track = TRUE, ...) {
  # Error checking:
  if(missing(clade.col)) stop("clade.col not defined. Please, define the",
                              " column with clade names.")
  if(class(phy)!="phylo") stop("phy must be class 'phylo'")
  
  #Matching tree and phylogeny using utils.R
  datphy <- match_dataphy(get(trait.col) ~ 1, data, phy)
  full.data <- datphy$data
  phy <- datphy$phy
  trait     <- full.data[[trait.col]]
  names(trait)  <- phy$tip.label
  N <- nrow(full.data)
  
  if (is.na(match(clade.col, names(full.data)))) {
    stop("Names column '", clade.col, "' not found in data frame'")
  }
  
  # Identify CLADES to use and their sample size 
  all.clades <- levels(full.data[ ,clade.col])
  wc <- table(full.data[ ,clade.col]) > n.species
  uc <- table(full.data[ , clade.col])[wc]
  
  if (length(uc) == 0) stop(paste("There is no clade with more than ",
                                  n.species," species. Change 'n.species' to fix this
                                  problem",sep=""))
  
  ### Fit full data model
  mod.0     <- phytools::phylosig(x = trait, tree = phy, method = method, test = TRUE)
  e.0 <- mod.0[[1]]  
  p.0 <- mod.0$P
  
  #Create dataframe to store estmates for each clade
  sensi.clade <-
    data.frame("clade" =I(as.character()), 
               "N.species" = numeric(),"estimate"=numeric(),
               "DF"=numeric(),"perc"=numeric(),
               "pval"=numeric())
  
  # Create dataframe store simulations (null distribution)
  null.dist <- data.frame("clade" = rep(names(uc), each = n.sim),
                          "estimate"= numeric(length(uc)*n.sim),
                          "DF"=numeric(length(uc)*n.sim))
  
  ### START LOOP between CLADES:
  # counters:
  aa <- 1; bb <- 1

  if(track==TRUE) pb <- utils::txtProgressBar(min = 0, max = length(uc)*n.sim,
                              style = 3)
  for (A in names(uc)){
    ### Number of species in clade A
    cN  <- as.numeric(uc[names(uc) == A])
    
    ### Fit reduced model (without clade)
    crop.data <- full.data[!full.data[ ,clade.col] %in% A,]
    crop.sp <-   which(full.data[ ,clade.col] %in% A)
    crop.phy <-  ape::drop.tip(phy,phy$tip.label[crop.sp])
    crop.trait <- crop.data[, trait.col]
    names(crop.trait) <- crop.phy$tip.label
    
    mod.s = phytools::phylosig(x = crop.trait, crop.phy, method = method, test = TRUE)
                               
    ### Raw differance
    DF              <- mod.s[[1]] - e.0
    ### Percentage of differance
    perc           <- round((abs(DF/e.0))*100,digits=1)
    ### Pvalues
    pval           <- mod.s$P
    
    # Store reduced model parameters: 
    estim.simu <- data.frame(A, cN, mod.s[[1]], DF, perc,
                             pval,
                             stringsAsFactors = F)
    sensi.clade[aa, ]  <- estim.simu
    
    ### START LOOP FOR NULL DIST:
    # number of species in clade A:
    for (i in 1:n.sim) {
      exclude <- sample(1:N, cN)
      crop.data <- full.data[-exclude,]
      crop.phy <-  ape::drop.tip(phy,phy$tip.label[exclude])
      crop.trait <- crop.data[, trait.col]
      names(crop.trait) <- crop.phy$tip.label
      
      mod.s = phytools::phylosig(x = crop.trait, crop.phy, method = method, test = FALSE)
      
      ### Raw differance
      DF              <- mod.s[[1]] - e.0
      ### Percentage of differance
      perc           <- round((abs(DF/e.0))*100,digits=1)
      
      null.dist[bb, ]  <- data.frame(clade = as.character(A), 
                                     estimate = mod.s[[1]], DF)
                                     
      
      if(track==TRUE) (utils::setTxtProgressBar(pb, bb))
      bb <- bb + 1
    }
    aa <- aa + 1
  }
  if(track==TRUE) on.exit(close(pb))
  
  #OUTPUT
  #full model estimates:
  param0 <- data.frame(estimate = e.0,
                 Pval = p.0)
  
  #Generates output:
  res <- list(call = match.call(),
              trait = trait.col,
              method = method,
              full.data.estimates = param0,
              sensi.estimates = sensi.clade,
              null.dist = null.dist,
              data = full.data)
  class(res) <- "clade.physig"
  return(res)
}
