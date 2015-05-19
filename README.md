sensiC
======

R package to perform model diagnostic for comparative methods

####Installing sensiC from Github:

```{r}
# First install package `devtools` (you should skip this with you have `devtools`)
install.packages("devtools")

devtools::install_github("paternogbc/sensiC")

# Required packages:
library(caper);library(ggplot2);library(gridExtra)
library(sensiC)
```

Loading and organizing data:
```{r}
data(shorebird)

# Organizing comparative data for pgls:
bird.comp <- comparative.data(shorebird.tree, shorebird.data, Species, 
        vcv=TRUE, vcv.dim=3)
```

Original Linear regression (PGLS):
```{r}
mod0 <- pgls(log(Egg.Mass) ~ log(M.Mass), data=bird.comp,"ML")
summary(mod0)
```

#### Model diagnostics with sensiC package:

##### Example: Estimating sampling effort bias with `samp_pgls`

```{r}
# First match the order of species in data and phy:
ord <- match(shorebird.tree$tip.label,shorebird.data$Species)
shorebird.data <- shorebird.data[ord,]

# Run sensitive analysis:
samp <- samp_pgls(log(Egg.Mass) ~ log(M.Mass),data=shorebird.data,phy=shorebird.tree)

# To check the results:
samp$results

# You can also specify the number of simulation and break intervals:
samp2 <- samp_pgls(log(Egg.Mass) ~ log(M.Mass),data=bird.comp$data,phy=bird.comp$phy,
                 times= 50, breaks=c(0.1,.2,.3,.4,.5,.6,.7,.8))
```

##### Example: Estimating influential points and parameter bias with `influ_pgls`

```{r}
# Run influential analysis:
influ <- influ_pgls(log(Egg.Mass) ~ log(M.Mass),data=shorebird.data,phy=shorebird.tree)
# Check the results
influ$results
# Check the most influential species:
influ[[4]]
# Check for species that presented errors during simulations:
influ$errors
```
### Visualizing Results with `sensi_plot`
```{r}
sensi_plot(samp)
sensi_plot(samp2)
sensi_plot(influ)

```

### Output `samp_gls`:
![Output samp_gls](http://i.imgur.com/zp5JXIJ.jpg)

### Output `influ_gls`:
![Output samp_gls](http://i.imgur.com/gF6GuEH.jpg)
