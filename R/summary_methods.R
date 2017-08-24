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
  #rows.inter <- match(names(sp.inter), sensi.estimates$species)
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
  n.tree <- length(unique(object$sensi.estimates$iteration))
  
  ## Estimate:
  sp.estimate <- unlist(as.list(object$influential.species$influ.sp.estimate$influ.sp.estimate))
  sp.estimate.tab <- table(sp.estimate)
  sp.estimate <- sp.estimate.tab[order(sp.estimate.tab,decreasing=T)] 
  sp.estimate.tab <- data.frame("Species removed" = names(sp.estimate), 
                                "Significant" = (as.numeric(sp.estimate)/n.tree)*100)
  colnames(sp.estimate.tab) <- c("Species removed", "(%) of iterations")
  
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
  #rows.inter <- match(names(sp.inter), sensi.estimates$species)
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
summary.sensiIntra_Tree <- function(object, ...){
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
summary.samp.physig <- function(obejct, ...){
  method <- obejct$call$method
  if (is.null(obejct$call$method)) method <- "K"
  
  simu <- nrow(obejct$samp.physig.estimates)
  res <- obejct$sign.analysis
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
