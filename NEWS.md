## sensiPhy 0.5.0.9000

* `clade_phylm()` and `clade_phyglm()` now account for clade sample size bias.
This is done by estimating a null distribution of intercepts and slopes considering only
the number of species in the clade.

* `summary()` methods `clade_phylm()` & `clade_phyglm()` now includes a permutation test
to account for the number of species in clades (tests model parameters (without the focal clade) is within the null distribution - one-tailed permutation test).

* `sensi_plot()` for clade analysis now include a histogram with the simulated DFslopes (null distribution).

## sensiPhy 0.5.0

Fixed all issues before submission to CRAN.
