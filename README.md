[![Travis-CI Build Status](https://travis-ci.org/paternogbc/sensiPhy.svg?branch=master)](https://travis-ci.org/paternogbc/sensiPhy)

sensiPhy
========

R package to perform sensitivity analysis for comparative methods

####Installing sensiPhy from Github:

```{r}
# First install the package `devtools` (you should skip this if you have `devtools`)
install.packages("devtools")

# Install `sensiPhy` from github: 
devtools::install_github("paternogbc/sensiPhy")

# Loading required packages:
library(phylolm);library(ggplot2);library(gridExtra)
library(sensiPhy)
```

Simulating data:
```{r}
set.seed(2468)
    tree <- rtree(100)
    pred<- rTraitCont(tree,root.value=0,sigma=1,model="BM")
    cont_trait1 <- pred + rTraitCont(tree,model="BM",sigma=2.5)
    bin_trait1<-rbinTrait(n=1,tree,beta=c(-1,0.5),alpha=3,
                          X=cbind(rep(1,length(tree$tip.label)),pred))
    dat<-data.frame(pred,cont_trait1,bin_trait1)
```

Original Linear regression (PGLS):
```{r}
mod0 <- phylolm(cont_trait1 ~ pred, data=dat,phy=tree)
summary(mod0)
```

#### Sensitivity analysis with sensiPhy package:

##### Example: Estimating sampling effort bias with `samp_phylolm`

```{r}
# Run sensitive analysis:
samp <- samp_phylolm(cont_trait1 ~ pred,data=dat,phy=tree)

# To check the results:
head(samp$samp.model.estimates)

# You can also specify the number of simulation and break intervals:
samp2 <- samp_phylolm(cont_trait1 ~ pred,data=dat,phy=tree,
                 times= 50, breaks=c(0.1,.2,.3,.4,.5,.6,.7,.8))
```

##### Example: Estimating influential species and parameter bias with `influ_phylolm`

```{r}
# Run influential analysis:
influ <- influ_phylolm(cont_trait1 ~ pred,data=dat,phy=tree)
# Check the results
influ[[6]]
# Check the most influential species:
influ[[5]]
# Check for species that presented errors during simulations:
influ$errors
```
### Visualizing Results with `sensi_plot`
```{r}
sensi_plot(samp)
sensi_plot(samp2)
sensi_plot(influ)

```

### Output `samp_phylolm`:
![Output samp_phylolm](http://i.imgur.com/YyKMEbX.jpg)

### Output `influ_phylolm`:
![Output samp_gls](http://i.imgur.com/gF6GuEH.jpg)
