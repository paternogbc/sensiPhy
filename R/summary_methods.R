### METHODS for: PGLS--------------------------------------------

### Summary method for class: sensiClade:--------------------------------------

#' @export
summary.sensiClade <- function(object, ...){
    ### Permutation test:
    ce <- object$sensi.estimates
    nd <- object$null.dist
    c <- levels(nd$clade)
    
    
    stats.slo <- data.frame("clade removed" = c, 
                            "N.species" = ce$N.species,
                            "estimate" = numeric(length(c)),
                            "DIFestimate" = numeric(length(c)),
                            "change" = numeric(length(c)),
                            "Pval" = numeric(length(c)),
                            "m.null.estimate" = numeric(length(c)),
                            "Pval.randomization" = numeric(length(c)))
    stats.int <- data.frame("clade removed" = c, 
                            "N.species" = ce$N.species,
                            "intercept" = numeric(length(c)),
                            "DIFintercept" = numeric(length(c)),
                            "change" = numeric(length(c)),
                            "Pval" = numeric(length(c)),
                            "m.null.intercept" = numeric(length(c)),
                            "Pval.randomization" = numeric(length(c)))
    aa <- 1
    for(j in c) {
      
      nes <- nd[nd$clade == j, ] # null estimates
      ces <- ce[ce$clade == j, ] # reduced model estimates
      times <- nrow(nes)
      
      ### Permutation test SLOPE:
      if (ces$DIFestimate > 0){
        p.slo <- sum(nes$estimate >= ces$estimate)/times
      }
      if (ces$DIFestimate < 0){
        p.slo <- sum(nes$estimate <= ces$estimate)/times
      }
      
      stats.slo[aa, -c(1:2)] <- data.frame(
            estimate = ces$estimate,
            DIFestimate = ces$DIFestimate,
            ces$estimate.perc,
            Pval = ces$pval.estimate,
            m.null.estimate = mean((nes$estimate)),
            Pval.randomization = p.slo)
      names(stats.slo)[5] <- "Change (%)"      
      
      ### Permutation test intercept:
      if (ces$DIFintercept > 0){
        p.int <- sum(nes$intercept >= ces$intercept)/times
      }
      if (ces$DIFintercept < 0){
        p.int <- sum(nes$intercept <= ces$intercept)/times
      }
      
      stats.int[aa, -c(1:2)] <- data.frame(
          intercept = ces$intercept,
          DIFintercept = ces$DIFintercept,
          ces$intercept.perc,
          Pval = ces$pval.estimate,
          m.null.intercept = mean((nes$intercept)),
          Pval.randomization = p.int)
      
      names(stats.int)[5] <- "Change (%)"
      
      aa <- aa+1
    }
    
    
    ### Sort by % of change:
    ord.slo <- order(object$sensi.estimates$estimate.perc, decreasing = TRUE)
    
    res <- list(stats.slo[ord.slo, ], stats.int[ord.slo, ])
    names(res) <- c("estimate", "Intercept")
    res
}

### Summary method for class: sensiClade.TraitEvol:--------------------------------------

#' @export
summary.sensiClade.TraitEvol <- function(object, ...){
  if(as.character(object$call[[1]])=="clade_discrete"){ #Check what type of TraitEvolution is evaluated
  ### Permutation test:
  ce <- object$sensi.estimates
  nd <- object$null.dist
  c <- levels(nd$clade)
  
  
  stats.q12 <- data.frame("clade removed" = c, 
                          "N.species" = ce$N.species,
                          "q12" = numeric(length(c)),
                          "DIFq12" = numeric(length(c)),
                          "change" = numeric(length(c)),
                          "m.null.q12" = numeric(length(c)),
                          "Pval.randomization" = numeric(length(c)))
  stats.q21 <- data.frame("clade removed" = c, 
                          "N.species" = ce$N.species,
                          "q21" = numeric(length(c)),
                          "DIFq21" = numeric(length(c)),
                          "change" = numeric(length(c)),
                          "m.null.q21" = numeric(length(c)),
                          "Pval.randomization" = numeric(length(c)))
  aa <- 1
  for(j in c) {
    
    nes <- nd[nd$clade == j, ] # null estimates
    ces <- ce[ce$clade == j, ] # reduced model estimates
    times <- nrow(nes)
    
    ### Permutation test q12:
    if (ces$DIFq12 > 0){
      p.q12 <- sum(nes$q12 >= ces$q12)/times
    }
    if (ces$DIFq12 < 0){
      p.q12 <- sum(nes$q12 <= ces$q12)/times
    }
    
    stats.q12[aa, -c(1:2)] <- data.frame(
      q12 = ces$q12,
      DIFq12 = ces$DIFq12,
      ces$q12.perc,
      m.null.q12 = mean((nes$q12)),
      Pval.randomization = p.q12)
    names(stats.q12)[5] <- "Change (%)"      
    
    ### Permutation test q21:
    if (ces$DIFq21 > 0){
      p.q21 <- sum(nes$q21 >= ces$q21)/times
    }
    if (ces$DIFq21 < 0){
      p.q21 <- sum(nes$q21 <= ces$q21)/times
    }
    
    stats.q21[aa, -c(1:2)] <- data.frame(
      q21 = ces$q21,
      DIFq21 = ces$DIFq21,
      ces$q21.perc,
      m.null.q21 = mean((nes$q21)),
      Pval.randomization = p.q21)
    
    names(stats.q21)[5] <- "Change (%)"
    
    aa <- aa+1
  }
  
  
  ### Sort by % of change:
  ord.q12 <- order(object$sensi.estimates$q12.perc, decreasing = TRUE)
  
  res <- list(stats.q12[ord.q12, ], stats.q21[ord.q12, ])
  names(res) <- c("q12", "q21")
  return(res)
  }
  #Check again what type of TraitEvolution is evaluated. If not discrete, it should be continuous. 
    if(as.character(object$call[[1]])=="clade_continuous"){
    ### Permutation test:
    ce <- object$sensi.estimates
    nd <- object$null.dist
    c <- levels(nd$clade)
    
    
    stats.sigsq <- data.frame("clade removed" = c, 
                            "N.species" = ce$N.species,
                            "sigsq" = numeric(length(c)),
                            "DIFsigsq" = numeric(length(c)),
                            "change" = numeric(length(c)),
                            "m.null.sigsq" = numeric(length(c)),
                            "Pval.randomization" = numeric(length(c)))
    stats.optpar <- data.frame("clade removed" = c, 
                            "N.species" = ce$N.species,
                            "optpar" = numeric(length(c)),
                            "DIFoptpar" = numeric(length(c)),
                            "change" = numeric(length(c)),
                            "m.null.optpar" = numeric(length(c)),
                            "Pval.randomization" = numeric(length(c)))
    aa <- 1
    for(j in c) {
      
      nes <- nd[nd$clade == j, ] # null estimates
      ces <- ce[ce$clade == j, ] # reduced model estimates
      times <- nrow(nes)
      
      ### Permutation test sigsq:
      if (ces$DIFsigsq > 0){
        p.sigsq <- sum(nes$sigsq >= ces$sigsq)/times
      }
      if (ces$DIFsigsq < 0){
        p.sigsq <- sum(nes$sigsq <= ces$sigsq)/times
      }
      
      stats.sigsq[aa, -c(1:2)] <- data.frame(
        sigsq = ces$sigsq,
        DIFsigsq = ces$DIFsigsq,
        ces$sigsq.perc,
        m.null.sigsq = mean((nes$sigsq)),
        Pval.randomization = p.sigsq)
      names(stats.sigsq)[5] <- "Change (%)"      
      
      if(object$optpar!="BM"){
      ### Permutation test optpar:
      if (ces$DIFoptpar > 0){
        p.optpar <- sum(nes$optpar >= ces$optpar)/times
      }
      if (ces$DIFoptpar < 0){
        p.optpar <- sum(nes$optpar <= ces$optpar)/times
      }
      
      stats.optpar[aa, -c(1:2)] <- data.frame(
        optpar = ces$optpar,
        DIFoptpar = ces$DIFoptpar,
        ces$optpar.perc,
        m.null.optpar = mean((nes$optpar)),
        Pval.randomization = p.optpar)
      
      names(stats.optpar)[5] <- "Change (%)"
      }
      
      aa <- aa+1
    }
    
    
    ### Sort by % of change:
    ord.sigsq <- order(object$sensi.estimates$sigsq.perc, decreasing = TRUE)
    
    if(object$optpar!="BM") {
    res <- list(stats.sigsq[ord.sigsq, ], stats.optpar[ord.sigsq, ])
    names(res) <- c("sigsq", "optpar")
    return(res)
    }
    
    if(object$optpar=="BM") {
      res <- list(stats.sigsq[ord.sigsq, ])
      names(res) <- c("sigsq")
      return(res)
    }
  }
}

### Summary method for class: sensiIntra_Clade:--------------------------------------

#' @export
summary.sensiIntra_Clade <- function(object, ...){
  ### Permutation test:
  ce <- object$sensi.estimates
  nd <- object$null.dist
  c <- levels(nd$clade)
  it <- unique(nd$iteration)
  
  
  stats.slo <- data.frame("clade removed" = rep(c,each=length(it)), 
                          "N.species" = rep(ce$N.species[1:length(c)],each=length(it)),
                          "estimate" = numeric(length(c)),
                          "DIFestimate" = numeric(length(c)),
                          "change" = numeric(length(c)),
                          "Pval" = numeric(length(c)),
                          "m.null.estimate" = numeric(length(c)),
                          "Pval.randomization" = numeric(length(c)))
  stats.int <- data.frame("clade removed" = rep(c,each=length(it)), 
                          "N.species" = rep(ce$N.species[1:length(c)],each=length(it)),
                          "intercept" = numeric(length(c)),
                          "DIFintercept" = numeric(length(c)),
                          "change" = numeric(length(c)),
                          "Pval" = numeric(length(c)),
                          "m.null.intercept" = numeric(length(c)),
                          "Pval.randomization" = numeric(length(c)))
  aa <- 1
  
  for(j in c) {
    for(i in 1:length(it)){
      
      nes <- nd[nd$clade == j & nd$iteration == it[i], ] # null estimates
      ces <- ce[ce$clade == j & ce$iteration == it[i], ] # reduced model estimates
      times <- nrow(nes)
      
      ### Permutation test SLOPE:
      if (ces$DIFestimate > 0){
        p.slo <- sum(nes$estimate >= ces$estimate)/times
      }
      if (ces$DIFestimate < 0){
        p.slo <- sum(nes$estimate <= ces$estimate)/times
      }
      
      stats.slo[aa, -c(1:2)] <- data.frame(
        estimate = ces$estimate,
        DIFestimate = ces$DIFestimate,
        ces$estimate.perc,
        Pval = ces$pval.estimate,
        m.null.estimate = mean((nes$estimate)),
        Pval.randomization = p.slo)
      
      names(stats.slo)[5] <- "Change (%)" 
      
      ### Permutation test intercept:
      if (ces$DIFintercept > 0){
        p.int <- sum(nes$intercept >= ces$intercept)/times
      }
      if (ces$DIFintercept < 0){
        p.int <- sum(nes$intercept <= ces$intercept)/times
      }
      
      stats.int[aa, -c(1:2)] <- data.frame(
        intercept = ces$intercept,
        DIFintercept = ces$DIFintercept,
        ces$intercept.perc,
        Pval = ces$pval.estimate,
        m.null.intercept = mean((nes$intercept)),
        Pval.randomization = p.int)
      
      names(stats.int)[5] <- "Change (%)"
      
      aa <- aa+1
    }
  }
  
  #calculate means and number of P.vals<0.005 for each clade
  perc.ran.slo <- ((stats::aggregate(stats.slo$Pval.randomization<=0.05,by=list(stats.slo$clade.removed),FUN=sum))$x)/length(it)*100
  stats.slo <- stats::aggregate(.~clade.removed, data=stats.slo, mean)
  stats.slo$perc.ran.slo<-round(perc.ran.slo,1)
  
  names(stats.slo)[9] <-"Significant (%)"
  
  perc.ran.int <- ((stats::aggregate(stats.int$Pval.randomization<=0.05,by=list(stats.int$clade.removed),FUN=sum))$x)/length(it)*100
  stats.int <- stats::aggregate(.~clade.removed, data=stats.int, mean)
  stats.int$perc.ran.int<-round(perc.ran.int,1)
  
  names(stats.int)[9] <-"Significant (%)"
  
  ### Sort by % of change:
  ord.slo <- order(stats.slo$`Change (%)`, decreasing = TRUE)
  
  res <- list(stats.slo[ord.slo, ], stats.int[ord.slo, ])
  names(res) <- c("Estimate", "Intercept")
  res
}



### Summary method for class: sensiTree_Clade:--------------------------------------

#' @export
summary.sensiTree_Clade <- function(object, ...){
  ### Permutation test:
  ce <- object$sensi.estimates
  nd <- object$null.dist
  c <- levels(nd$clade)
  it <- unique(nd$iteration)
  
  
  stats.slo <- data.frame("clade removed" = rep(c,each=length(it)), 
                          "N.species" = rep(ce$N.species[1:length(c)],each=length(it)),
                          "estimate" = numeric(length(c)),
                          "DIFestimate" = numeric(length(c)),
                          "change" = numeric(length(c)),
                          "Pval" = numeric(length(c)),
                          "m.null.estimate" = numeric(length(c)),
                          "Pval.randomization" = numeric(length(c)))
  stats.int <- data.frame("clade removed" = rep(c,each=length(it)), 
                          "N.species" = rep(ce$N.species[1:length(c)],each=length(it)),
                          "intercept" = numeric(length(c)),
                          "DIFintercept" = numeric(length(c)),
                          "change" = numeric(length(c)),
                          "Pval" = numeric(length(c)),
                          "m.null.intercept" = numeric(length(c)),
                          "Pval.randomization" = numeric(length(c)))
  aa <- 1
  
  for(j in c) {
    for(i in 1:length(it)){
      
      nes <- nd[nd$clade == j & nd$iteration == it[i], ] # null estimates
      ces <- ce[ce$clade == j & ce$iteration == it[i], ] # reduced model estimates
      times <- nrow(nes)
      
      ### Permutation test SLOPE:
      if (ces$DIFestimate > 0){
        p.slo <- sum(nes$estimate >= ces$estimate)/times
      }
      if (ces$DIFestimate < 0){
        p.slo <- sum(nes$estimate <= ces$estimate)/times
      }
      
      stats.slo[aa, -c(1:2)] <- data.frame(
        estimate = ces$estimate,
        DIFestimate = ces$DIFestimate,
        ces$estimate.perc,
        Pval = ces$pval.estimate,
        m.null.estimate = mean((nes$estimate)),
        Pval.randomization = p.slo)
      
      names(stats.slo)[5] <- "Change (%)" 
      
      ### Permutation test intercept:
      if (ces$DIFintercept > 0){
        p.int <- sum(nes$intercept >= ces$intercept)/times
      }
      if (ces$DIFintercept < 0){
        p.int <- sum(nes$intercept <= ces$intercept)/times
      }
      
      stats.int[aa, -c(1:2)] <- data.frame(
        intercept = ces$intercept,
        DIFintercept = ces$DIFintercept,
        ces$intercept.perc,
        Pval = ces$pval.estimate,
        m.null.intercept = mean((nes$intercept)),
        Pval.randomization = p.int)
      
      names(stats.int)[5] <- "Change (%)"
      
      aa <- aa+1
    }
  }
  
  #calculate means and number of P.vals<0.005 for each clade
  perc.ran.slo <- ((stats::aggregate(stats.slo$Pval.randomization<=0.05,by=list(stats.slo$clade.removed),FUN=sum))$x)/length(it)*100
  stats.slo <- stats::aggregate(.~clade.removed, data=stats.slo, mean)
  stats.slo$perc.ran.slo<-round(perc.ran.slo,1)
  
  names(stats.slo)[9] <-"Significant (%)"
  
  perc.ran.int <- ((stats::aggregate(stats.int$Pval.randomization<=0.05,by=list(stats.int$clade.removed),FUN=sum))$x)/length(it)*100
  stats.int <- stats::aggregate(.~clade.removed, data=stats.int, mean)
  stats.int$perc.ran.int<-round(perc.ran.int,1)
  
  names(stats.int)[9] <-"Significant (%)"
  
  ### Sort by % of change:
  ord.slo <- order(stats.slo$`Change (%)`, decreasing = TRUE)
  
  res <- list(stats.slo[ord.slo, ], stats.int[ord.slo, ])
  names(res) <- c("Estimate", "Intercept")
  res
}

### Summary method for class: sensiInflu:--------------------------------------

#' @export
summary.sensiInflu <- function(object, ...){
    sp.estimate <- object$influential.species$influ.sp.estimate
    rows.estimate <- match(sp.estimate, object$sensi.estimates$species)
    estimate <- object$sensi.estimates[rows.estimate, c("species","estimate","DIFestimate","estimate.perc","pval.estimate")]
    ord.estimate <- order(estimate$estimate.perc, 
                       decreasing = TRUE)
    estimate <- estimate[ord.estimate, ]
    rownames(estimate) <- NULL
    colnames(estimate) <- c("Species removed", "Estimate", "DIFestimate", "Change(%)", "Pval")
    
    sp.inter <-object$influential.species$influ.sp.intercept
    rows.inter <- match(sp.inter, object$sensi.estimates$species)
    inter <- object$sensi.estimates[rows.inter, c("species","intercept","DIFintercept","intercept.perc","pval.intercept")]
    ord.inter <- order(inter$intercept.perc, 
                       decreasing = TRUE)
    inter <- inter[ord.inter, ]
    rownames(inter) <- NULL
    colnames(inter) <- c("Species removed", "Intercept", "DIFintercept", "Change(%)", "Pval")
    
    res <- list("Influential species for the Estimate" = sp.estimate, "Estimate" = estimate,
                "Influential species for the Intercept" = sp.inter, "Intercept" = inter)
    return(res)
    
}

### Summary method for class: sensiInflu.TraitEvol:--------------------------------------

#' @export
summary.sensiInflu.TraitEvol <- function(object, ...){
  if(as.character(object$call[[1]])=="influ_discrete"){
  sp.q12 <- object$influential.species$influ.sp.q12
  rows.q12 <- match(sp.q12, object$sensi.estimates$species)
  q12 <- object$sensi.estimates[rows.q12, c("species","q12","DIFq12","q12.perc")]
  ord.q12 <- order(q12$q12.perc, 
                        decreasing = TRUE)
  q12 <- q12[ord.q12, ]
  rownames(q12) <- NULL
  colnames(q12) <- c("Species removed", "q12", "DIFq12", "Change(%)")
  
  sp.q21 <-object$influential.species$influ.sp.q21
  rows.q21 <- match(sp.q21, object$sensi.estimates$species)
  q21 <- object$sensi.estimates[rows.q21, c("species","q21","DIFq21","q21.perc")]
  ord.q21 <- order(q21$q21.perc, 
                     decreasing = TRUE)
  q21 <- q21[ord.q21, ]
  rownames(q21) <- NULL
  colnames(q21) <- c("Species removed", "q21", "DIFq21", "Change(%)")
  
  res <- list("Influential species for q12" = sp.q12, "q12" = q12,
              "Influential species for q21" = sp.q21, "q21" = q21)
  return(res)
  }
  if(as.character(object$call[[1]])=="influ_continuous"){
    sp.sigsq <- object$influential.species$influ.sp.sigsq
    rows.sigsq <- match(sp.sigsq, object$sensi.estimates$species)
    sigsq <- object$sensi.estimates[rows.sigsq, c("species","sigsq","DIFsigsq","sigsq.perc")]
    ord.sigsq <- order(sigsq$sigsq.perc, 
                     decreasing = TRUE)
    sigsq <- sigsq[ord.sigsq, ]
    rownames(sigsq) <- NULL
    colnames(sigsq) <- c("Species removed", "sigsq", "DIFsigsq", "Change(%)")
    
    sp.optpar <-object$influential.species$influ.sp.optpar
    rows.optpar <- match(sp.optpar, object$sensi.estimates$species)
    optpar <- object$sensi.estimates[rows.optpar, c("species","optpar","DIFoptpar","optpar.perc")]
    ord.optpar <- order(optpar$optpar.perc, 
                     decreasing = TRUE)
    optpar <- optpar[ord.optpar, ]
    rownames(optpar) <- NULL
    colnames(optpar) <- c("Species removed", "optpar", "DIFoptpar", "Change(%)")
    
    res <- list("Influential species for sigsq" = sp.sigsq, "sigsq" = sigsq,
                "Influential species for optpar" = sp.optpar, "optpar" = optpar)
    return(res)
  }
}

### Summary method for class: sensiIntra_Influ:--------------------------------------

#' @export
summary.sensiIntra_Influ <- function(object, ...){
  n.tree <- length(unique(object$sensi.estimates$iteration))
  
  ## Estimate:
  sp.estimate <- unlist(as.list(object$influential.species$influ.sp.estimate$influ.sp.estimate))
  sp.estimate.tab <- table(sp.estimate)
  sp.estimate <- sp.estimate.tab[order(sp.estimate.tab,decreasing=T)] 
  sp.estimate.tab <- data.frame("Species removed" = names(sp.estimate), 
                                "Significant" = (as.numeric(sp.estimate)/n.tree)*100)
  colnames(sp.estimate.tab) <- c("Species removed", "(%) of iterations")
  if(nrow(sp.estimate.tab) > 20) sp.estimate.tab <- sp.estimate.tab[1:20, ]
  
  sensi.estimates<-object$sensi.estimates
  #rows.estimate <- match(names(sp.estimate), sensi.estimates$species)
  rows.estimate <- sensi.estimates$species %in% sp.estimate.tab$`Species removed`
  estimate <- sensi.estimates[rows.estimate, c("species","estimate","DIFestimate","estimate.perc","pval.estimate")]
  estimate <- stats::aggregate(estimate[,2:5],list(estimate$species),mean)
  names(estimate)[1]<-"species"
  ord.estimate <- order(estimate$estimate.perc,decreasing = TRUE)
  estimate <- estimate[ord.estimate, ]
  rownames(estimate) <- NULL
  colnames(estimate) <- c("Species removed", "Estimate", "DIFestimate", "Change(%)", "Pval")
  
  
  ### Intercept:
  sp.inter <-unlist(as.list(object$influential.species$influ.sp.intercept$influ.sp.intercept))
  sp.inter.tab <- table(sp.inter)
  sp.inter <- sp.inter.tab[order(sp.inter.tab,decreasing=T)] #Consider giving the counts, rather than just order> 
  sp.inter.tab <- data.frame(names(sp.inter), (as.numeric(sp.inter)/n.tree)*100)
  colnames(sp.inter.tab) <- c("Species removed", "(%) of iterations")
  if(nrow(sp.inter.tab) > 20) sp.inter.tab <- sp.inter.tab[1:20, ]
  
  
  rows.inter <- sensi.estimates$species %in% sp.inter.tab$`Species removed`
  inter <- sensi.estimates[rows.inter, c("species","intercept","DIFintercept","intercept.perc","pval.intercept")]
  inter <- stats::aggregate(inter[,2:5],list(inter$species),mean)
  names(inter)[1]<-"species"
  ord.inter <- order(inter$intercept.perc,decreasing = TRUE)
  inter <- inter[ord.inter, ]
  rownames(inter) <- NULL
  colnames(inter) <- c("Species removed", "Intercept", "DIFintercept", "Change(%)", "Pval")
  
  res <- list("Most Common Influential species for the Estimate" = sp.estimate.tab, "Average Estimates" = estimate,
              "Most Common Influential species for the Intercept" = sp.inter.tab, "Average Intercepts" = inter)
  return(res)
  
}

### Summary method for class: sensiTree_Influ:--------------------------------------

#' @export
summary.sensiTree_Influ <- function(object, ...){
  summary.sensiIntra_Influ(object, ...)
}

### Summary method for class: sensiSamp:----------------------------------------

#' @export
summary.sensiSamp <- function(object, ...){
    simu <- nrow(object$sensi.estimates)
    sig <- object$sign.analysis
    sig$perc.sign.intercept <- sig$perc.sign.intercept * 100
    sig$perc.sign.estimate <- sig$perc.sign.estimate * 100
    names(sig) <- c("% Species Removed", 
                    "% Significant Intercepts",
                    "Mean Intercept Change (%)",
                    "Mean sDIFintercept",
                    "% Significant Estimates",
                    "Mean Estimate Change (%)",
                    "Mean sDIFestimate")
    
    message(paste(simu, "simulations saved," ,
                  "see output$sensi.estimates to acess all simulations"))
    return(sig)
}

### Summary method for class: sensiSamp.TraitEvol:----------------------------------------

#' @export
summary.sensiSamp.TraitEvol <- function(object, ...){
  simu <- nrow(object$sensi.estimates)
  sig <- object$breaks.summary.tab
  
  if(as.character(object$call[[1]])=="samp_discrete"){ 
  names(sig) <- c("% Species Removed", 
                  "Mean q12 Change (%)",
                  "Mean sDIFq12",
                  "Median sDIFq12",
                  "Mean q21 Change (%)",
                  "Mean sDIFq21",
                  "Median sDIFq21")
  }
  
  if(as.character(object$call[[1]])=="samp_continuous"){ 
    names(sig) <- c("% Species Removed", 
                    "Mean sigsq Change (%)",
                    "Mean sDIFsigsq",
                    "Median sDIFsigsq",
                    "Mean optpar Change (%)",
                    "Mean sDIFoptpar",
                    "Median sDIFoptpar")
  }
  
  message(paste(simu, "simulations saved," ,
                "see output$sensi.estimates to acess all simulations"))
  return(sig)
}




### Summary method for class: sensiTree_Samp:----------------------------------------

#' @export
summary.sensiTree_Samp <- function(object, ...){
  simu <- nrow(object$sensi.estimates)
  sig <- object$sign.analysis
  sig$perc.sign.intercept <- sig$perc.sign.intercept * 100
  sig$perc.sign.estimate <- sig$perc.sign.estimate * 100
  sig <- stats::aggregate(.~percent_sp_removed, data=sig, mean)
  sig$iteration <- NULL
  
  names(sig) <- c("% Species Removed", 
                  "% Significant Intercepts",
                  "Mean Intercept Change (%)",
                  "Mean sDIFintercept",
                  "% Significant Estimates",
                  "Mean Estimate Change (%)",
                  "Mean sDIFestimate")
  
  message(paste(simu, "simulations saved," ,
                "see output$sensi.estimates to acess all simulations"))
  return(sig)
}

### Summary method for class: sensiIntra_Samp:----------------------------------------

#' @export
summary.sensiIntra_Samp <- function(object, ...){
  simu <- nrow(object$sensi.estimates)
  sig <- object$sign.analysis
  sig$perc.sign.intercept <- sig$perc.sign.intercept * 100
  sig$perc.sign.estimate <- sig$perc.sign.estimate * 100
  sig <- stats::aggregate(.~percent_sp_removed, data=sig, mean)
  sig$iteration <- NULL
  
  names(sig) <- c("% Species Removed", 
                  "% Significant Intercepts",
                  "Mean Intercept Change (%)",
                  "Mean sDIFintercept",
                  "% Significant Estimates",
                  "Mean Estimate Change (%)",
                  "Mean sDIFestimate")
  
  message(paste(simu, "simulations saved," ,
                "see output$sensi.estimates to acess all simulations"))
  return(sig)
}


### Summary method for class: sensiIntra:--------------------------------------

#' @export
summary.sensiIntra <- function(object, ...){
    res <- object$stats
    return(res)
}

### Summary method for class: sensiTree:----------------------------------------

#' @export
summary.sensiTree <- function(object, ...){
    res <- object$stats
    return(res)
}

### Summary method for class: sensiTree_Intra:--------------------------------------

#' @export
summary.sensiTree_Intra <- function(object, ...){
  res <- object$stats
  return(res)
}


### METHODS for: Phylogenetic signal--------------------------------------------
# Summary method for class influ.physig-----------------------------------------
#' @export
summary.influ.physig <- function(object, ...){
  
  method <- object$call$method
  if(is.null(object$call$method)) method <- "K"
  sp <- object$influential.species
  rows <- match(sp, object$influ.physig.estimates$species)
  estim <- object$influ.physig.estimates[rows, -6]
  ord <- order(estim$perc, decreasing = TRUE)
  sta <- estim[ord, ] # Output oredered summay.
  rownames(sta) <- NULL
  colnames(sta) <- c("Species removed", method, "DF", "Change(%)", "Pval")
  
  res.0 <- data.frame(Trait = object$trait, N.species = nrow(object$data), 
                      estimate = object$full.data.estimates[[1]],
                      Pval = object$full.data.estimates[[2]])
  colnames(res.0)[3] <- method
  ### Output list:
  res <- list(res.0, sp, sta)
  names(res) <- c(paste("Full data estimate", sep = " "),
                  paste("Influential species for", method, sep = " "), 
                  paste("Summary for", method, sep = " "))
  return(res)
}

### Summary method for class: clade.physig -------------------------------------

#' @export
summary.clade.physig <- function(object, ...){
  ce <- object$sensi.estimates
  nd <- object$null.dist
  c <- levels(nd$clade)
  
  method <- object$method
  
  stats <- data.frame("clade removed" = c, 
                      "N.species" = ce$N.species,
                      "estimate" = numeric(length(c)),
                      "DIF.est" = numeric(length(c)),
                      "change" = numeric(length(c)),
                      "Pval" = numeric(length(c)),
                      "m.null.estimate" = numeric(length(c)),
                      "Pval.randomization" = numeric(length(c)))
  aa <- 1
  for(j in c) {
    nes <- nd[nd$clade == j, ] # null estimates
    ces <- ce[ce$clade == j, ] # reduced model estimates
    times <- nrow(nes)
    
    ### Permutation test K:
    if (ces$DF > 0){
      p <- sum(nes$estimate  >= ces$estimate)/times
    }
    if (ces$DF < 0){
      p <- sum(nes$estimate <= ces$estimate)/times
    }
    
    stats[aa, -c(1:2)] <- data.frame(
      estimate = ces$estimate,
      DF = ces$DF,
      ces$perc,
      Pval = ces$pval,
      m.null = mean((nes$estimate)),
      Pval.randomization = p)
    names(stats)[5] <- "Change (%)"      
    
    aa <- aa + 1
  }
  
  
  ### Sort by % of change:
  ord <- order(object$sensi.estimates$perc, decreasing = TRUE)
  res.0 <- data.frame(Trait = object$trait, N.species = nrow(object$data), 
                      estimate = object$full.data.estimates$estimate,
                      Pval = object$full.data.estimates$Pval)
  res.1 <- stats[ord, ]
  res <- list(res.0, res.1)
  names(res)[[1]] <- paste("Full data:", method, "phylogenetic signal estimate", sep = " ")
  names(res)[[2]] <- "Summary by clade removal"
  return(res)
}

### Summary method for class: samp.physig:--------------------------------------
#' @export
summary.samp.physig <- function(object, ...){
  method <- object$call$method
  if (is.null(object$call$method)) method <- "K"
  
  simu <- nrow(object$samp.physig.estimates)
  res <- object$sign.analysis
  res$perc.sign <- res$perc.sign * 100
  names(res) <- c("Species Removed (%)", 
                  paste("Significant", method, "(%)"),
                  "Mean Change (%)",
                  "Mean sDFestimate")
  
  message(paste(simu, "simulations saved," ,
                "see output$samp.physig.estimates to acess all simulations"))
  return(res)
}

### Summary method for class: tree.physig:--------------------------------------

#' @export
summary.tree.physig <- function(object, ...){
  res <- list(object$call,
              object$stats)
  names(res) <- c("Call", "Summary")
  return(res)
}

### Summary method for class: intra.physig:--------------------------------------

#' @export
summary.intra.physig <- function(object, ...){
  res <- list(object$call,
              object$stats)
  names(res) <- c("Call", "Summary")
  return(res)
}

### METHODS for: diversification rate--------------------------------------------
### Summary method for class: sensiTree.TraitEvol:--------------------------------------
#' @export
summary.sensiTree.TraitEvol <- function(object, ...){
  res <- list(round(object$stats,4))
  names(res) <- c("Summary")
  return(res)
}

### METHODS for: diversification rate--------------------------------------------
### tree.bd
#' @export

summary.tree.bd <- function(object, ...){
  res <- list(object$call,
              object$stats)
  names(res) <- c("Call", "Summary")
}  


