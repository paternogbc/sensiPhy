## ----setup, include=FALSE, message=FALSE---------------------------------
knitr::opts_chunk$set(fig.height = 6.5, fig.width = 10)

## ----message=T-----------------------------------------------------------
set.seed(12345)
library(sensiPhy)

### Loading data:
data(alien)
data(primates) # se ?alien & ?primates for details about the data.

## ----samp_analysis, echo=T, cache=T, warning=FALSE-----------------------
# run analysis:
samp <- samp_phylm(log(gestaLen) ~ log(adultMass), phy = alien$phy[[1]], 
                   data = alien$data, times = 10, track = F)

# You can change the number of repetitions and break intervals:
samp2 <- samp_phylm(log(gestaLen) ~ log(adultMass), phy = alien$phy[[1]], track = F,
                    data = alien$data, times = 10, breaks = c(0.1, 0.2, 0.3, 0.4))
# You can change the phylogenetic model:
samp <- samp_phylm(log(gestaLen) ~ log(adultMass), phy = alien$phy[[1]], 
                   data = alien$data, model = "kappa", track = F)

# Check results:
knitr::kable(summary(samp))
# Visual diagnostics
sensi_plot(samp2)
# You can specify which graph and parameter ("slope" or "intercept") to print: 
sensi_plot(samp2, graphs = 1)
sensi_plot(samp2, param = "intercept")

## ----influ_analysis, echo=T, cache=T, warning=FALSE----------------------
# run analysis:
influ <- influ_phylm(log(gestaLen) ~ log(adultMass), phy = alien$phy[[1]], 
                     data = alien$data, track = F)
# To check summary results:
summary(influ)
# Most influential speciesL
influ$influential.species
# Visual diagnostics
sensi_plot(influ)

## ----clade_analysis, echo=T, cache=T, warning=FALSE, fig.height=6, fig.width=10----
# run analysis:
clade <- clade_phylm(log(sexMaturity) ~ log(adultMass), phy = primates$phy[[1]],
                     data = primates$data, clade.col = "family", track = F)
# To check summary results and most influential clades:
knitr::kable(summary(clade))
# Visual diagnostics for clade removal:
sensi_plot(clade, "Cercopithecidae")
sensi_plot(clade, "Cebidae")

## ----tree_analysis, echo=T, cache=T, warning=FALSE-----------------------
# This analysis needs a multiphylo file:
class(alien$phy)
alien$phy
# run PGLS accounting for phylogenetic uncertain:
tree <- tree_phylm(log(gestaLen) ~ log(adultMass), phy = alien$phy, 
                   data = alien$data, times = 100, track = F)
# To check summary results:
knitr::kable(summary(tree))
# Visual diagnostics
sensi_plot(tree)

## ----intra_analysis, echo=T, cache=T, warning=FALSE----------------------
# run PGLS accounting for intraspecific variation:
intra <- intra_phylm(gestaLen ~ adultMass, phy = alien$phy[[1]], track = F, 
                     data = alien$data, Vy = "SD_gesta", Vx = "SD_mass",
                     times = 100, x.transf = log, y.transf = log)
# To check summary results:
knitr::kable(summary(intra))
# Visual diagnostics
sensi_plot(intra)

