### Examples:

### First install sensiC from Github:
devtools::install_github("paternogbc/sensiC")

### Required packages:
library(phylolm);library(phytools);library(sensiC)

set.seed(111)
N <- 50 # Number of species
### Simulating tree
tree<-pbtree(n=N)
### Simulating response variable with phylogenetic signal
Ly <- rTrait(n=1, tree, model=c("lambda"),parameters=list(lambda=.8))
### Simulating explanatory variable
Lx <- Ly + rnorm(N,mean(Ly),1)
### Including Species names:
sp <- tree$tip.label
regre <- data.frame(sp,Ly,Lx)

### Organizing comparative data for pgls:
comp.data <- comparative.data(data=regre,phy=tree,vcv=T,vcv.dim=3,names.col="sp")

### Linear regression (PGLS):
mod0 <- pgls(Ly ~Lx, data=comp.data,"ML")
summary(mod0)

### Example: samp_pgls
samp1 <- samp_pgls(Ly ~ Lx,data=comp.data)
### You can specify the number of replicates and break intervals:
samp2 <- samp_pgls(Ly ~ Lx,data=comp.data,times=99,breaks=c(.1,.3,.5))

### Example: influ_pgls
influ <- influ_pgls(Ly ~ Lx,data=comp.data)
### Estimated parameters:
head(influ$results)
### Most influential species:
influ[[5]]
### Check for species with erros erros:
influ$errors

### Visualizing Results:
sensi_plot(samp1,method="sampling")
sensi_plot(samp2,method="sampling")
sensi_plot(influ,method="influence")


