#' Alien Mammals dataset: Example dataset for the package sensiPhy
#'
#' A comparative dataset containing traits for 94 alien mammal species 
#' (alien.data) and a multiphylo object with 101 phylogenies matching the 
#' data (alien.phy). Tip labels are the binomial species names and match 
#' with data rownames. Data was taken from (Gonzalez-Surez et al. 2015) and phylogenies from
#' (Fritz et al 2009). 
#' @usage data(alien)
#' @format A data frame with 94 rows and 7 variables:
#' \itemize{
#'   \item family: Taxonomic family
#'   \item Mass: Mean adult body mass (g)
#'   \item gesta: Mean gestation length (days)
#'   \item range: Mean home range (km)
#'   \item SE_mass: Standard deviation (intraspecific) for mean adult body mass (g)
#'   \item SE_gesta: Standard deviation (intraspecific) for mean gestation length (days)
#'   \item SE_range: Standard deviation (intraspecific) for mean home range (km)
#'   }
#' @format A multiphylo containing 101 trees for 94 mammal species.
#' @references Alien mammal data: Gonzalez-Suarez, Manuela, Sven Bacher, and Jonathan M. Jeschke. 
#' "Intraspecific trait variation is correlated with establishment success of alien mammals." 
#' The American Naturalist 185.6 (2015): 737-746  DOI: 10.1086/681105
#' 
#' Downloaded from: Gonzalez-Surez M, Bacher S, Jeschke J (2015) Data from: Intraspecific trait 
#' variation is correlated with establishment success of alien mammals.
#' Dryad Digital Repository. http://dx.doi.org/10.5061/dryad.sp963
#' 
#' Phylogeny: Kuhn, Tyler S., Arne O. Mooers, and Gavin H. Thomas. "A simple polytomy resolver for 
#' dated phylogenies." Methods in Ecology and Evolution 2.5 (2011): 427-436.
#' 
#' Fritz, Susanne A., Olaf RP Bininda-Emonds, and Andy Purvis. "Geographical variation in predictors 
#' of mammalian extinction risk: big is bad, but only in the tropics." Ecology letters 12.6 (2009): 538-549.
#' 
"alien"


