### First install sensiC from Github:
devtools::install_github("paternogbc/sensiC")

## Required packages:
library(caper);library(ggplot2);library(gridExtra)
library(sensiC)

### Loading data:
data(shorebird)

## Organizing comparative data for pgls:
bird.comp <- comparative.data(shorebird.tree, shorebird.data, Species, vcv=TRUE, vcv.dim=3)

### Original Linear regression (PGLS):
mod0 <- pgls(Egg.Mass ~ M.Mass, data=bird.comp,"ML")
summary(mod0)

## Sensitive analysis with sensiC package:

### Example: Estimating sample size bias with `samp_pgls`

samp <- samp_pgls(log(Egg.Mass) ~ log(M.Mass),data=bird.comp$data,phy=bird.comp$phy)

### You can specify number of simulation and break intervals:
samp2 <- samp_pgls(log(Egg.Mass) ~ log(M.Mass),data=bird.comp$data,phy=bird.comp$phy,
                 times= 200, breaks=c(0.1,.2,.3,.4,.5,.6,.7))


### Example: Estimating influential points and parameter bias with `influ_pgls`
influ <- influ_pgls(log(Egg.Mass) ~ log(M.Mass),data=bird.comp$data,phy=bird.comp$phy)
### Estimated parameters:
head(influ$results)
### Most influential species:
influ[[4]]
### Check for species with erros erros:
influ$errors


## Visualizing Results:
sensi_plot(samp)
sensi_plot(samp2)
sensi_plot(influ)


