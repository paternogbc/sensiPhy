## sensiPhy 0.7.0
* `match_data_phy()` now accepts datasets with no information on species names as row names. If the number of species corresponds to the number of tips a warning informs the user that the function assumes that the dataset and the phylogeny are in the same order.


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

Fixed all issues before submission to CRAN.
