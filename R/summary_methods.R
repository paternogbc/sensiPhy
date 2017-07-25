### Summary method for class: sensiClade:--------------------------------------

#' @export
summary.sensiClade <- function(object, ...){
    ### Permutation test:
    ce <- object$clade.model.estimates
    nd <- object$null.dist
    c <- levels(nd$clade)
    
    
    stats.slo <- data.frame("clade removed" = c, 
                            "N.species" = ce$N.species,
                            "slope" = numeric(length(c)),
                            "DFslope" = numeric(length(c)),
                            "change" = numeric(length(c)),
                            "Pval" = numeric(length(c)),
                            "m.null.slope" = numeric(length(c)),
                            "Pval.randomization" = numeric(length(c)))
    stats.int <- data.frame("clade removed" = c, 
                            "N.species" = ce$N.species,
                            "intercept" = numeric(length(c)),
                            "DFintercept" = numeric(length(c)),
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
      if (ces$DFslope > 0){
        p.slo <- sum(nes$slope >= ces$slope)/times
      }
      if (ces$DFslope < 0){
        p.slo <- sum(nes$slope <= ces$slope)/times
      }
      
      stats.slo[aa, -c(1:2)] <- data.frame(
            slope = ces$slope,
            DFslope = ces$DFslope,
            ces$slope.perc,
            Pval = ces$pval.slope,
            m.null.slope = mean((nes$slope)),
            Pval.randomization = p.slo)
      names(stats.slo)[5] <- "Change (%)"      
      
      ### Permutation test intercept:
      if (ces$DFintercept > 0){
        p.int <- sum(nes$intercept >= ces$intercept)/times
      }
      if (ces$DFintercept < 0){
        p.int <- sum(nes$intercept <= ces$intercept)/times
      }
      
      stats.int[aa, -c(1:2)] <- data.frame(
          intercept = ces$intercept,
          DFintercept = ces$DFintercept,
          ces$intercept.perc,
          Pval = ces$pval.slope,
          m.null.intercept = mean((nes$intercept)),
          Pval.randomization = p.int)
      
      names(stats.int)[5] <- "Change (%)"
      
      aa <- aa+1
    }
    
    
    ### Sort by % of change:
    ord.slo <- order(object$clade.model.estimates$slope.perc, decreasing = TRUE)
    
    res <- list(stats.slo[ord.slo, ], stats.int[ord.slo, ])
    names(res) <- c("Slope", "Intercept")
    res
    }

### Summary method for class: sensiInflu:--------------------------------------

#' @export
summary.sensiInflu <- function(object, ...){
    sp.slope <- object$influential.species$influ.sp.slope
    rows.slope <- match(sp.slope, object$influ.model.estimates$species)
    slope <- object$influ.model.estimates[rows.slope, c("species","slope","DFslope","slope.perc","pval.slope")]
    ord.slope <- order(slope$slope.perc, 
                       decreasing = TRUE)
    slope <- slope[ord.slope, ]
    rownames(slope) <- NULL
    colnames(slope) <- c("Species removed", "Slope", "DFslope", "Change(%)", "Pval")
    
    sp.inter <-object$influential.species$influ.sp.intercept
    rows.inter <- match(sp.inter, object$influ.model.estimates$species)
    inter <- object$influ.model.estimates[rows.inter, c("species","intercept","DFintercept","intercept.perc","pval.intercept")]
    ord.inter <- order(inter$intercept.perc, 
                       decreasing = TRUE)
    inter <- inter[ord.inter, ]
    rownames(inter) <- NULL
    colnames(inter) <- c("Species removed", "Intercept", "DFintercept", "Change(%)", "Pval")
    
    res <- list("Influential species for the Slope" = sp.slope, "Slope Estimates" = slope,
                "Influential species for the Intercept" = sp.inter, "Intercept Estimates" = inter)
    return(res)
    
}

### Summary method for class: sensiIntra_Influ and sensiTree_Influ:--------------------------------------

#' @export
summary.sensiINTER_Influ <- function(object, ...){
  sp.slope <- unlist(as.list(object$influential.species$influ.sp.slope$influ.sp.slope))
  sp.slope.tab <- table(sp.slope)
  sp.slope <- sp.slope.tab[order(sp.slope.tab,decreasing=T)] 
  influ.model.estimates<-object$influ.model.estimates
  rows.slope <- match(names(sp.slope), influ.model.estimates$species)
  slope <- influ.model.estimates[rows.slope, c("species","slope","DFslope","slope.perc","pval.slope")]
  slope <- aggregate(slope[,2:5],list(slope$species),mean)
  names(slope)[1]<-"species"
  ord.slope <- order(slope$slope.perc,decreasing = TRUE)
  slope <- slope[ord.slope, ]
  rownames(slope) <- NULL
  colnames(slope) <- c("Species removed", "Slope", "DFslope", "Change(%)", "Pval")
  
  sp.inter <-unlist(as.list(object$influential.species$influ.sp.intercept$influ.sp.intercept))
  sp.inter.tab <- table(sp.inter)
  sp.inter <- sp.inter.tab[order(sp.inter.tab,decreasing=T)] #Consider giving the counts, rather than just order> 
  rows.inter <- match(names(sp.inter), influ.model.estimates$species)
  inter <- influ.model.estimates[rows.inter, c("species","intercept","DFintercept","intercept.perc","pval.intercept")]
  inter <- aggregate(inter[,2:5],list(inter$species),mean)
  names(inter)[1]<-"species"
  ord.inter <- order(inter$intercept.perc,decreasing = TRUE)
  inter <- inter[ord.inter, ]
  rownames(inter) <- NULL
  colnames(inter) <- c("Species removed", "Intercept", "DFintercept", "Change(%)", "Pval")
  
  res <- list("Most Common Influential species for the Slope" = sp.slope, "Mean Slope Estimates" = slope,
              "Most Common Influential species for the Intercept" = sp.inter, "Mean Intercept Estimates" = inter)
  return(res)
  
}

### Summary method for class: sensiSamp:----------------------------------------

#' @export
summary.sensiSamp <- function(object, ...){
    simu <- nrow(object$samp.model.estimates)
    sig <- object$sign.analysis
    sig$perc.sign.intercept <- sig$perc.sign.intercept * 100
    sig$perc.sign.slope <- sig$perc.sign.slope * 100
    if (length(intersect(class(object),c("sensiTree_Samp", "sensiTree_SampL","sensiIntra_Samp","sensiIntra_SampL"))) != 0) {
      names(sig) <- c("iteration",
                      "% Species Removed", 
                      "% Significant Intercepts",
                      "Mean Intercept Change (%)",
                      "Mean sDFintercept",
                      "% Significant Slopes",
                      "Mean Slope Change (%)",
                      "Mean sDFslope")}
    
    else
    names(sig) <- c("% Species Removed", 
                    "% Significant Intercepts",
                    "Mean Intercept Change (%)",
                    "Mean sDFintercept",
                    "% Significant Slopes",
                    "Mean Slope Change (%)",
                    "Mean sDFslope")
    message(paste(simu, "simulations saved," ,
                  "see output$samp.model.estimates to acess all simulations"))
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

### Summary method for class: sensiIntra_Tree:--------------------------------------

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
