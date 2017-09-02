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
*  `samp_physig()`: Performs sensitivity analysis of influential species for phylogenetic signal estimate (k or lambda)
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
