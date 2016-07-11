#' Match data and phylogeny based on model formula
#'
#' Combines phylogeny and data to ensure that tips in phylogeny match data and that observations
#' with missing values are removed. This function uses variables provided in the 
#' `formula` argument to:
#' \itemize{
#'  \item{Remove NA`s:   } {Check if there is any row with NA in the variables included in 
#'  the formula. All rows containing NA will be removed from the data}
#'  \item{Match data and phy:   } {Check if tips from phylogeny match rownames in
#'  data. Tips not present in data and phy will be removed from the phylogeny and
#'  data}
#'  \item{Return matched data and phy:   } {The returned data
#'  has no NA in the variables included in `formula` and only rows that match phylogeny
#'  tips. Returned phy has only tips that match data}
#'  }
#'  Used internally in \code{\link{samp_phylm}},\code{\link{samp_phyglm}},\code{\link{clade_phylm}},
#'  \code{\link{clade_phyglm}},\code{\link{intra_phylm}},\code{\link{intra_phyglm}},\code{\link{tree_phylm}}
#'  and \code{\link{clade_tree_phyglm}}. Users can also directly use this function to combine a phylogeny and a dataset. 
#'
#' @param formula The model formula
#' @param data Data frame containing species traits with row names matching tips
#' in \code{phy}.
#' @param phy A phylogeny (class 'phylo' or 'multiphylo')
#' @return The function \code{match_dataphy} returns a list with the following
#' components:
#' @return \code{data}: Cropped dataset matching phylogeny
#' @return \code{phy}: Cropped phylogeny matching data
#' @return \code{dropped}: Species dropped from phylogeny and removed from data.
#' 
#' @details This function uses all variables provided in the `formula` to match
#' data and phylogeny. To avoid cropping the full dataset, `match_dataphy` searches
#' for NA values only on variables provided by formula. Missing values on 
#' other variables, not included in `formula`, will not be removed from data. 
#' 
#' This ensures consistance between data and phylogeny only for the variables 
#' that are being used in the model (set by `formula`).
#' 
#' If phy is a 'multiphylo' object, all phylogenies will be cropped
#' to match data. The returned phyogeny will be a 'multiphylo' object.
#' @note If tips are removed from the phylogeny and data or if rows containing
#' missing values are removed from data, a message will be printed with the 
#' details. Further, the final number of species that match data and phy will
#' always be reported by a message.
#' 
#' @author Caterina Penone & Gustavo Paterno
#' @references This function is largely inspired by the function \code{comparative.data} in caper package
#' David Orme, Rob Freckleton, Gavin Thomas, Thomas Petzoldt, Susanne Fritz, Nick Isaac and Will Pearse
#' (2013). caper: Comparative Analyses of Phylogenetics and Evolution in R. R package version 0.5.2.
#' http://CRAN.R-project.org/package=caper
#' @examples 
#'# Load data:
#'data(alien)
#'head(alien$data)
#'# Match data and phy based on model formula:
#'comp.data <- match_dataphy(gestaLen ~ adultMass, data = alien$data, alien$phy)
#'comp.data2 <- match_dataphy(homeRange ~ adultMass, data = alien$data, alien$phy)
#'# Check data:
#'comp.data$data
#'# Check phy:
#'comp.data$phy
#'# See species dropped from phy or data:
#'comp.data$dropped
#' @export
match_dataphy <- function(formula, data, phy){
    
    # original data set:
    data.0 <- data
    # Cropping data frame by formula variables:
    mf <- stats::model.frame(formula = formula, data = data.0, na.action = stats::na.exclude)
    if (nrow(data.0) > nrow(mf)) warning("NA's in response or predictor,", 
                                         " rows with NA's were removed")
    
    #Match data and phylogeny in comparative.data style
    if(inherits(phy, "multiPhylo")){  
        phy1 <- phy[[1]]}
    else
        phy1<-phy
    
    tiplabl <- phy1$tip.label
    taxa.nam <- as.character(rownames(mf))
    
    in.both <- intersect(taxa.nam, tiplabl)
    
    if (length(in.both) == 0)
        stop("No tips are common to the dataset and phylogeny, 
        Please check if row names of your dataset contain species/tip names")
    
    mismatch <- union(setdiff(tiplabl,taxa.nam),setdiff(taxa.nam,tiplabl))
    if (length(mismatch) != 0)   warning("Some phylogeny tips do not match species in data,",
                                         " species were dropped from phylogeny",
                                         " or data")
    
    #Drop species from tree
    if(inherits(phy, "multiPhylo")){ 
        phy <- lapply(phy, ape::drop.tip,tip = mismatch)
        class(phy)<-"multiPhylo"
        tip.order <- match(phy[[1]]$tip.label, rownames(data))
    }
    if(inherits(phy, "phylo")){ 
        phy <- ape::drop.tip(phy,tip = mismatch)
        class(phy)<-"phylo"
        tip.order <- match(phy$tip.label, rownames(data))
    }
    
    if (any(is.na(tip.order)))
        stop("Problem with sorting data frame: mismatch between tip labels and data frame labels")
    data <- data[tip.order, , drop = FALSE]
    data.out <- data.0[rownames(data),]
    
    message(paste("Used dataset has ",nrow(data.out)," species that match data and phylogeny"))
    return(list(data = data.out, phy = phy, dropped = mismatch))
}
