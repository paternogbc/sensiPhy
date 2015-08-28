#' Alien Mammals dataset: Example dataset for the package sensiPhy
#'
#' A comarative dataset containing traits for 94 alien mammal species 
#' (alien.data) and a multiphylo object with 101 phylogenies matching the 
#' data (alien.phy). Tip labels are the binomial species names and match 
#' with data rownames. Data was taken from (REFERENCE) and phylogenies from
#' (REFERENCE). 
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
#' @references INCLUDE DATA and PHY citation here!
"alien"


