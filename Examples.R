### Examples:

### First install sensiC from Github:
devtools::install_github("paternogbc/sensiC")

### Required packages:
library(caper);library(ggplot2);library(gridExtra)
library(sensiC)

### Loading data:
data(shorebird)

### Organizing comparative data for pgls:
comp.data <- comparative.data(shorebird.tree, shorebird.data, Species, vcv=TRUE, vcv.dim=3)

### Linear regression (PGLS):
mod0 <- pgls(log(Egg.Mass) ~ log(M.Mass), data=comp.data,"ML")
summary(mod0)

### Example: samp_pgls
samp1 <- samp_pgls(log(Egg.Mass) ~ log(M.Mass),data=comp.data)
### You can specify the number of replicates and break intervals:
samp2 <- samp_pgls(log(Egg.Mass) ~ log(M.Mass),data=comp.data,times=20,breaks=c(.1,.3,.5))

### Example: samp_gls,
# First we need to match tip.labels with rownames in data:
sp.ord <- match(shorebird.tree$tip.label, rownames(shorebird.data))
shorebird.data <- shorebird.data[sp.ord,]
samp3 <- samp_gls(log(Egg.Mass) ~ log(M.Mass),data=shorebird.data,phy=shorebird.tree)


### Example: influ_pgls
influ1 <- influ_pgls(log(Egg.Mass) ~ log(M.Mass),data=comp.data)
### Estimated parameters:
head(influ$results)
### Most influential species:
influ[[5]]
### Check for species with erros erros:
influ$errors

### Example influ_gls:
### First we need to match tip.labels with rownames in data:
sp.ord <- match(shorebird.tree$tip.label, rownames(shorebird.data))
shorebird.data <- shorebird.data[sp.ord,]

### Now we can run the function influ_gls:
influ2 <- influ_gls(log(Egg.Mass) ~ log(M.Mass),data=shorebird.data,phy=shorebird.tree)

### Visualizing Results:
sensi_plot(samp1)
sensi_plot(samp2)
sensi_plot(samp3)
sensi_plot(influ1)
sensi_plot(influ2)
