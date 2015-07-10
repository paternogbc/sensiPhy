library(ape);
library(phylolm)
library(sensiPhy)
context("model equality: samp_xxx & influ_xxx")

###: test samp_phylolm:
test_that("mod.0 samp_phylolm is equal to phylolm", {
    set.seed(2468)
    tree <- rtree(50)
    pred<- rTraitCont(tree,root.value=0,sigma=1,model="BM")
    cont_trait1 <- pred + rTraitCont(tree,model="BM",sigma=0.1)
    bin_trait1<-rbinTrait(n=1,tree,beta=c(-1,0.5),alpha=0.1,
                          X=cbind(rep(1,length(tree$tip.label)),pred))
    dat<-data.frame(pred,cont_trait1,bin_trait1)
    original <- phylolm(cont_trait1 ~ pred,phy = tree,data = dat,model = "lambda")
    mod.samp <- samp_phylolm(cont_trait1 ~ pred,phy = tree,data = dat,model = "lambda",
                             track = F,times = 1)
    expect_equal(summary.phylolm(original)$coefficients,
                 mod.samp$full.model.estimates$coef)
})

###: test samp_phyloglm:
test_that("mod.0 samp_phyloglm is equal to phyloglm", {
    set.seed(2468)
    tree <- rtree(50)
    pred<- rTraitCont(tree,root.value=0,sigma=1,model="BM")
    cont_trait1 <- pred + rTraitCont(tree,model="BM",sigma=0.1)
    bin_trait1<-rbinTrait(n=1,tree,beta=c(-1,0.5),alpha=0.1,
                          X=cbind(rep(1,length(tree$tip.label)),pred))
    dat<-data.frame(pred,cont_trait1,bin_trait1)
    original <- phyloglm(bin_trait1~pred,data = dat,phy = tree)
    mod.samp <- samp_phyloglm(bin_trait1~pred,data = dat,phy = tree,track = F,times = 1)
    expect_equal(summary.phyloglm(original)$coefficients,
                 mod.samp$full.model.estimates$coef)
})

###: test influ_phylolm:
test_that("mod.0 influ_phylolm is equal to phylolm", {
    set.seed(2468)
    tree <- rtree(50)
    pred<- rTraitCont(tree,root.value=0,sigma=1,model="BM")
    cont_trait1 <- pred + rTraitCont(tree,model="BM",sigma=0.1)
    bin_trait1<-rbinTrait(n=1,tree,beta=c(-1,0.5),alpha=0.1,
                          X=cbind(rep(1,length(tree$tip.label)),pred))
    dat<-data.frame(pred,cont_trait1,bin_trait1)
    original <- phylolm(cont_trait1 ~ pred,phy = tree,data = dat,model = "lambda")
    mod.samp <- influ_phylolm(cont_trait1 ~ pred,phy = tree,data = dat,model = "lambda",
                             track = F,times = 1)
    expect_equal(summary.phylolm(original)$coefficients,
                 mod.samp$full.model.estimates$coef)
})

###: test influ_phyloglm:
test_that("mod.0 influ_phyloglm is equal to phyloglm", {
    set.seed(2468)
    tree <- rtree(50)
    pred<- rTraitCont(tree,root.value=0,sigma=1,model="BM")
    cont_trait1 <- pred + rTraitCont(tree,model="BM",sigma=0.1)
    bin_trait1<-rbinTrait(n=1,tree,beta=c(-1,0.5),alpha=0.1,
                          X=cbind(rep(1,length(tree$tip.label)),pred))
    dat<-data.frame(pred,cont_trait1,bin_trait1)
    original <- phyloglm(bin_trait1~pred,data = dat,phy = tree)
    mod.samp <- influ_phyloglm(bin_trait1~pred,data = dat,phy = tree,track = F)
    expect_equal(summary.phyloglm(original)$coefficients,
                 mod.samp$full.model.estimates$coef)
})


