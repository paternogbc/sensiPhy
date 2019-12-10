#' Phylogenetic signal in missing data
#'
#' Calculates D statistic (Fritz & Purvis 2010), a measure of phylogenetic
#' signal, for missing data. Missingness is recoded into a binary variable 
#' (1=missing, 0=non missing). This function is an adaptation of
#' \code{\link[caper]{phylo.d}} for missing data.
#'
#' @param data Data frame containing species traits with species as row names.
#' @param phy A phylogeny (class 'phylo', see ?\code{ape}).
#' @inheritParams caper::phylo.d
#' \code{permut} Number of permutations to be used in the randomisation test.
#' \code{binvar} The name of the variable in \code{data} holding the variable of interest with
#' missing data.
#' @param ... Further arguments to be passed to \code{phylo.d}.
#' @details
#' This function builds on \code{\link[caper]{phylo.d}} to calculate a phylogenetic signal
#' in missing data. The variable of interest, usually a trait, is recoded into a binary variable
#' (1=missing data, 0=non missing data). Then the \code{\link[caper]{phylo.d}} function tests the estimated 
#' D value for significant departure from both random association and the clumping expected under a Brownian 
#' evolution threshold model (Fritz & Purvis, 2010).
#' 
#' Output can be visualised using \code{print()} and \code{plot()}
#'
#' @return The function \code{miss.phylo.d} returns an object of class "phylo.d" with the following
#' components, for complete list of arguments see \code{\link[caper]{phylo.d}} :
#' @return \code{DEstimate}: The estimated D value
#' @return \code{Pval1}: A p value, giving the result of testing whether D is significantly different from one
#' @return \code{Pval0}: A p value, giving the result of testing whether D is significantly different from zero
#' The function also prints the percentage of missing data per variable in the dataset.

#' @author Caterina Penone & Gustavo Paterno
#' @references 
#' 
#' Paterno, G. B., Penone, C. Werner, G. D. A. 
#' \href{http://doi.wiley.com/10.1111/2041-210X.12990}{sensiPhy: 
#' An r-package for sensitivity analysis in phylogenetic 
#' comparative methods.} Methods in Ecology and Evolution 
#' 2018, 9(6):1461-1467
#'
#' Fritz, S. A. and Purvis, A. (2010). Selectivity in mammalian extinction risk and threat types: a new measure of 
#' phylogenetic signal strength in binary traits. Conservation Biology, 24(4):1042-1051.
#' 
#' David Orme, Rob Freckleton, Gavin Thomas, Thomas Petzoldt, Susanne Fritz, Nick Isaac and Will Pearse (2013).
#' caper: Comparative Analyses of Phylogenetics and Evolution in R. R package version 0.5.2.
#' https://CRAN.R-project.org/package=caper
#' @examples 
#'# Load caper:
#'library(caper)
#'# Load data
#'data(primates)
#'data<-alien$data
#'phy=alien$phy[[1]]
#'
#'# Test phylogenetic signal for missing data:
#'sexNAsig <- miss.phylo.d(data,phy,binvar=homeRange)
#'print(sexNAsig)
#'plot(sexNAsig)
#'
#'massNAsig <- miss.phylo.d(data,phy,binvar=adultMass)
#'print(massNAsig)
#'plot(massNAsig)
#' @export

miss.phylo.d<-function(data, phy,...){
  
  sp.nam <- NULL
  names.col <- NULL
  
  #error check
  if (class(data) != "data.frame") stop("data must be class 'data.frame'")
  if (class(phy) != "phylo") stop("phy must be class 'phylo'")

  #calculate % of NAs per trait
  tot.sp <- nrow(data)
  nNA <- colSums(is.na(data))
  percNA <- round(nNA/tot.sp*100,digits=2)
  print("Percentage of missing data in traits:")
  print(percNA)
  
  #remove factor columns (categorical traits)
  factCols <- sapply(data,is.factor)
  data <- data[,!factCols]
  
  #recode traits with missing data into binary
  data[!is.na(data)] <- 0
  data[is.na(data)] <- 1

  #match with phylogeny
  if(is.null("names.col")){
    compdat<-caper::comparative.data(phy,data,names.col=names.col)}
  
  if(!is.null("names.col")){
    data$sp.nam<-row.names(data)
    compdat<-caper::comparative.data(phy,data,names.col=sp.nam)}

  
  #calculate d statistic using caper::phylo.d
  d.stat<-caper::phylo.d(compdat,...)
}
