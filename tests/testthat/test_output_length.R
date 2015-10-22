library(ape);
library(phylolm)
library(sensiPhy)
context("output length: samp_xxx & influ_xxx")

###: test output length samp_phylm:
test_that("output length is equal to times*breaks", {
    set.seed(2468)
    tree <- rtree(50)
    pred<- rTraitCont(tree,root.value=0,sigma=1,model="BM")
    cont_trait1 <- pred + rTraitCont(tree,model="BM",sigma=0.1)
    bin_trait1<-rbinTrait(n=1,tree,beta=c(-1,0.5),alpha=0.1,
                          X=cbind(rep(1,length(tree$tip.label)),pred))
    dat<-data.frame(pred,cont_trait1,bin_trait1)
    mod.samp <- samp_phylm(cont_trait1 ~ pred,phy = tree,data = dat,model = "lambda",
                             track = F,times = 5,
                             breaks = c(.1,.2,.3))
    expect_equal(nrow(mod.samp$samp.model.estimates),15)
})

###: test output length influ_phylm:
test_that("output length is equal to times*breaks", {
    set.seed(2468)
    tree <- rtree(50)
    pred<- rTraitCont(tree,root.value=0,sigma=1,model="BM")
    cont_trait1 <- pred + rTraitCont(tree,model="BM",sigma=0.1)
    bin_trait1<-rbinTrait(n=1,tree,beta=c(-1,0.5),alpha=0.1,
                         X=cbind(rep(1,length(tree$tip.label)),pred))
    dat<-data.frame(pred,cont_trait1,bin_trait1)
    mod.influ <- influ_phylm(cont_trait1 ~ pred,
                        phy = tree,data = dat,model = "lambda",track = F)
    expect_equal(nrow(mod.influ$influ.model.estimates),50)
})

###: test output length samp_phyglm:
test_that("mod.0 samp_phyglm is equal to phyloglm", {
    set.seed(2468)
    tree <- rtree(50)
    pred<- rTraitCont(tree,root.value=0,sigma=1,model="BM")
    cont_trait1 <- pred + rTraitCont(tree,model="BM",sigma=0.1)
    bin_trait1<-rbinTrait(n=1,tree,beta=c(-1,0.5),alpha=0.1,
                          X=cbind(rep(1,length(tree$tip.label)),pred))
    dat<-data.frame(pred,cont_trait1,bin_trait1)
    mod.samp <- samp_phyglm(bin_trait1 ~ pred,phy = tree,data = dat,
                             track = F,times = 5,
                             breaks = c(.1,.2,.3))
    expect_equal(nrow(mod.samp$samp.model.estimates),15)
})

###: test output length influ_phyglm:
test_that("mod.0 samp_phyglm is equal to phyloglm", {
    set.seed(2468)
    tree <- rtree(50)
    pred<- rTraitCont(tree,root.value=0,sigma=1,model="BM")
    cont_trait1 <- pred + rTraitCont(tree,model="BM",sigma=0.1)
    bin_trait1<-rbinTrait(n=1,tree,beta=c(-1,0.5),alpha=0.1,
                          X=cbind(rep(1,length(tree$tip.label)),pred))
    dat<-data.frame(pred,cont_trait1,bin_trait1)
    mod.influ <- influ_phyglm(bin_trait1 ~ pred,phy = tree,data = dat,
                                track = F)
    expect_equal(nrow(mod.influ$influ.model.estimates),50)
})

