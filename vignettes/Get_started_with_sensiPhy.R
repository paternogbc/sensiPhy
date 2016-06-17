## ----loading_package, message=FALSE--------------------------------------
library(sensiPhy)
library(knitr)
data("alien")

## ----data, echo=FALSE----------------------------------------------------
kable(head(alien$data))

## ---- crop.data, message=FALSE, warning=FALSE----------------------------
dat <- match_dataphy(log10(homeRange) ~ log10(adultMass), data = alien$data, 
                 phy = alien$phy)

## ------------------------------------------------------------------------
mod.0 <- phylolm(log10(homeRange) ~ log10(adultMass), data = dat$data, 
                 phy = dat$phy[[1]], model = "lambda")
mod.0$coefficients

## ---- message=F, warning=F-----------------------------------------------
samp <- samp_phylm(log10(homeRange) ~ log10(adultMass), data = dat$data, 
                 phy = dat$phy[[1]], model = "lambda", track = F)

## ---- fig.height=5.5, fig.width=7.5, warning=F, eval=F-------------------
#  sensi_plot(samp)

## ---- message=F, warning=F-----------------------------------------------
influ <- influ_phylm(log10(homeRange) ~ log10(adultMass), data = dat$data, 
                 phy = dat$phy[[1]], track = F)

## ---- fig.height=5.5, fig.width=7.5, warning=F, eval=F-------------------
#  sensi_plot(influ)

## ---- message=F, warning=F-----------------------------------------------
clade <- clade_phylm(log10(homeRange) ~ log10(adultMass), data = dat$data, 
                 phy = dat$phy[[1]], clade.col = "family", n.species = 3)
summary(clade)

## ---- fig.height=4, fig.width=7.5, warning=F, eval=F---------------------
#  sensi_plot(clade, clade = "Cervidae")

## ---- message=F, warning=F-----------------------------------------------
trees <- tree_phylm(log(homeRange) ~ log(adultMass), data = dat$data, 
                 phy = dat$phy, times = 30, track = F)
summary(trees)

## ---- fig.height=5.5, fig.width=7.5, warning=F, eval=F-------------------
#  sensi_plot(trees)

