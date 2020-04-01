# sensiPhy __v0.8.5__ (2020-03-31)

## Improvements:

1. This new version further adjusts sensiPhy code for R 4.0.0 release. All functions that call `data.frame()` were adjusted to work with stringsAsFactors = FALSE or 
stringsAsFactors = TRUE. 

# sensiPhy __v0.8.4__ (10 Dec 2019)

## BUG FIXES

1. This new version adjusts sensiPhy code for R 4.0.0 release.
In all sensiPhy functions the use of class(.) == was replaced by inherits(., *).

## NOTES

1. Updated sensiPhy citation reference.

## sensiPhy 0.8.3
#### Bug fix

This new version fix a small bug with functions for phylogenetic uncertainty (tree_xxx).
Minor typos were also fixed.

## sensiPhy 0.8.2
#### Bug fix

This new version fix all issues requested by CRAN.
In all sensiPhy functions the use of if() with conditions of length greater than one was corrected to avoid Errors.

## sensiPhy 0.8.1

#### Additions:

The vigentee now include two new sections:  
1. Using sensiPhy to analyse results from other packages  
2. How long does it take?  

Also available at the online tutorial: https://github.com/paternogbc/sensiPhy/wiki

#### Bug fix
* Corrected issue #182: intra_physig problem with se argument.

## sensiPhy 0.8.0

#### Major additions
`sensiPhy` now performs sensitivity analysis for a new class of methods which allows users to perform sensitivity analyses of both continuous and discrete (binary) macro-evolutionary models of trait evolution (e.g. Mkn models for binary traits, OU, BM, lambda etc. for continuous traits).

`sensiPhy` nor performs sensitivity analysis of phylogenetic uncertainty for simple metrics of diversification and speciation rates (Magallon and Sanderson (2000) method) or speciation rate using bd.km (Kendall-Moran method)

#### New functions (trait evolution)
##### Influential species:
* `influ_continuous()`: Performs sensitivity analysis of influential species for 
 models of trait evolution (continuous characters)
 
* `influ_discrete()`: Performs sensitivity analysis of influential species for 
models of trait evolution (binary discrete characters)

##### Influential clades:
*  `clade_continuous()`: Performs sensitivity analysis of influential clades for models of trait evolution (continuous characters)

*  `clade_discrete()`: Performs sensitivity analysis of influential clades for  for 
models of trait evolution (binary discrete characters)

##### Sampling size
*  `samp_continuous()`: Performs sensitivity analysis of species sampling for models of trait evolution (continuous characters)

*  `samp_discrete()`: Performs sensitivity analysis of species sampling for 
models of trait evolution (binary discrete characters)

##### Phylogenetic uncertainty
*  `tree_continuous()`: Performs sensitivity analysis of phylogenetic uncertainty for models of trait evolution (continuous characters)

*  `tree_discrete()`: Performs sensitivity analysis of phylogenetic uncertainty for 
models of trait evolution (binary discrete characters)

#### New functions (diversification rates)
##### Phylogenetic uncertainty
*  `tree_bd()`: Performs estimates of diversification rate evaluating uncertainty in trees topology.

#### New functions (Diagnostic plots and stats)
* `summary()` and `sensi_plot` methods were implemented (for all new functions) to provide a quick and intuitive  overview of results from sensitivite analysis. 

#### Bug fix
* Corrected progressbar bug when track=FALSE in physig functions

## sensiPhy 0.7.0
#### Core changes
`sensiPhy` now imports the package `phytools` 

#### Major additions
* `sensiPhy` now performs sensitivity analysis by interacting two types of uncertainty at the same time (tree and intra against influ, clade and samp methods)
* `sensiPhy` now performs sensitivity analysis for phylogenetic signal

#### New functions
##### Phylogenetic signal
*  `influ_physig()`: Performs sensitivity analysis of influential species for phylogenetic signal estimate (k or lambda)
*  `clade_physig()`: Performs sensitivity analysis of influential clades for phylogenetic signal estimate (k or lambda)
*  `samp_physig()`: Performs sensitivity analysis of species sampling for phylogenetic signal estimate (k or lambda)
*  `tree_physig()`: Performs sensitivity analysis of phylogenetic signal estimate (k or lambda) accounting for phylogenetic uncertainty
*  `intra_physig()`: Performs sensitivity analysis of phylogenetic signal estimate (k or lambda) accounting for intra-specific variation and measurement errors
##### Interactions for phylolm models
*  `tree_intra_phylm()`: Performs sensitivity analysis of interaction between phylogenetic uncertainty and  intraspecific variability for phylolm models (linear regression)
*  `tree_intra_phyglm()`: Performs sensitivity analysis of interaction between phylogenetic uncertainty  and  intraspecific variability for phylolm models (logistic regression)
*  `tree_clade_phylm()`: Performs sensitivity analysis of interaction between phylogenetic uncertainty and sensitivity to species sampling for phylolm models (linear regression)
*  `tree_clade_phyglm()`: Performs sensitivity analysis of interaction between phylogenetic uncertainty  and sensitivity to species sampling for phylolm models (logistic regression)
*  `tree_influ_phylm()`: Performs sensitivity analysis of interaction between phylogenetic uncertainty and influential species detection for phylolm models (linear regression)
*  `tree_influ_phyglm()`: Performs sensitivity analysis of interaction between phylogenetic uncertainty and influential species detection for phylolm models (logistic regression)
*  `tree_samp_phylm()`: Performs sensitivity analysis of interaction between phylogenetic uncertainty and sensitivity to species sampling for phylolm models (linear regression)
*  `tree_samp_phyglm()`: Performs sensitivity analysis of interaction between phylogenetic uncertainty  and sensitivity to species sampling for phylolm models (logistic regression)
*  `intra_clade_phylm()`: Performs sensitivity analysis of interaction between intraspecific variability and influential clades for phylolm models (linear regression)
*  `intra_clade_phyglm()`: Performs sensitivity analysis of interaction between intraspecific variability and influential clades for phylolm models (logistic regression)
*  `intra_influ_phylm()`: Performs sensitivity analysis of interaction between intraspecific variability and influential species detection for phylolm models (linear regression)
*  `intra_influ_phyglm()`: Performs sensitivity analysis of interaction between intraspecific variability and influential species detection for phylolm models (logistic regression)
*  `intra_samp_phylm()`: Performs sensitivity analysis of interaction between intraspecific variability and species sampling for phylolm models (linear regression)
*  `intra_samp_phyglm()`: Performs sensitivity analysis of interaction between intraspecific variability and species sampling for phylolm models (logistic regression)

#### Improvements
* `match_data_phy()` now accepts datasets with no information on species names as row names. If the number of species corresponds to the number of tips a warning informs the user that the function assumes that the dataset and the phylogeny are in the same order.

#### Naming standardization between functions:
* For all `sensiPhy` function the following changes were made:
1. slope -> estimate
2. DF -> DIF
3. model estimates -> sensi.estimates

#### Bug fix
* `Tree `methods: Data order was matching order of the first tree of the multiphylo file only. This bug is now fixed. Data and order matching is now done at each iteration.


## sensiPhy 0.6.0

#### New Functions

* `miss.phylo.d()` - Calculates phylogenetic signal for missing data (D statistic; Fritz & Purvis 2010).
Missingness is recoded into a binary variable. 

#### Improvements:

* The package now includes a __Vignette__ with a quick introduction to all sensiPhy functions.

* `clade_phylm()` and `clade_phyglm()` now account for clade sample size bias.
This is done by estimating a null distribution of intercepts and slopes considering only
the number of species in the clade.

* `summary()` methods for `clade_phylm()` & `clade_phyglm()` now includes a randomization test
to account for the number of species in clades (tests if change in model parameters (without the focal clade) is within the null distribution - one-tailed test).

* `sensi_plot()` for clade analysis now include a histogram with the simulated DFslopes (null distribution).

* `sensi_plot()` for influential species analysis (`influ_phylm` / `influ_phyglm`) now prints the names
of the most influential species on the regression plot.

* `sensi_plot()` now uses font size = 12 for better visualization.

* Packages datasets ("primates", "alien") now loads data and phylogeny in independent objects to
faciliate usage in examples. 

## sensiPhy 0.5.0

First submission to CRAN.
