### Summary method for class: sensiClade:--------------------------------------

#' @export
summary.sensiClade <- function(object, ...){
    ord.inter <- order(object$clade.model.estimates$intercept.perc, 
                       decreasing = TRUE)
    inter <- object$clade.model.estimates[ord.inter, c(1,2,3,4,5, 6)]
    colnames(inter) <- c("Clade removed", "N.species", "Intercept", "DFintercept",
                         "change (%)", "Pval")
    ord.slope <- order(object$clade.model.estimates$slope.perc, 
                       decreasing = TRUE)
    slope <- object$clade.model.estimates[ord.slope, c(1,2,7,8,9,10)]
    colnames(slope) <- c("Clade removed", "N.species", "Slope", "DFslope", 
                         "change (%)", "Pval")
   
    ### Permutation test:
    ce <- object$clade.model.estimates
    nd <- object$null.dist
    c <- levels(nd$clade)
    
    stats.slo <- data.frame(slope,
                            "m.null.slope" = numeric(length(c)),
                            "p.randomization" = numeric(length(c)))
    colnames(stats.slo)[5] <- "change (%)"
    stats.int <- data.frame(inter,
                            "m.null.intercept" = numeric(length(c)),
                            "p.randomization" = numeric(length(c)))
    colnames(stats.int)[5] <- "change (%)"
    
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

    stats.slo[aa, c(7,8)] <- c(
                              null.DFslope = mean((nes$slope)),
                              p.randomization = p.slo)
    names(stats.slo)[5] <- "Change (%)"      

    ### Permutation test intercept:
    if (ces$DFintercept > 0){
      p.int <- sum(nes$intercept >= ces$intercept)/times
    }
    if (ces$DFintercept < 0){
      p.int <- sum(nes$intercept <= ces$intercept)/times
    }
    
    stats.int[aa, c(7,8)] <- c(null.DFintercept = mean((nes$intercept)),
                               p.randomization = p.int)
    colnames(stats.int)[5] <- "Change (%)"
    
    aa <- aa+1
    }
    
    res <- list(stats.slo, stats.int)
    names(res) <- c("Slope", "Intercept")
    res
    }

### Summary method for class: sensiInflu:--------------------------------------

#' @export
summary.sensiInflu <- function(object, ...){
    sp.slope <- object$influential.species$influ.sp.slope
    rows.slope <- match(sp.slope, object$influ.model.estimates$species)
    slope <- object$influ.model.estimates[rows.slope, c(1,6,7,8,9)]
    ord.slope <- order(slope$slope.perc, 
                       decreasing = TRUE)
    slope <- slope[ord.slope, ]
    rownames(slope) <- NULL
    colnames(slope) <- c("Species removed", "Slope", "DFslope", "Change(%)", "Pval")
    
    sp.inter <-object$influential.species$influ.sp.intercept
    rows.inter <- match(sp.inter, object$influ.model.estimates$species)
    inter <- object$influ.model.estimates[rows.inter, c(1,2,3,4,5)]
    ord.inter <- order(inter$intercept.perc, 
                       decreasing = TRUE)
    inter <- inter[ord.inter, ]
    rownames(inter) <- NULL
    colnames(inter) <- c("Species removed", "Intercept", "DFintercept", "Change(%)", "Pval")
    
    res <- list("Influential species for the Slope" = sp.slope, "Slope Estimates" = slope,
                "Influential species for the Intercept" = sp.inter, "Intercept Estimates" = inter)
    return(res)
    
}

### Summary method for class: sensiSamp:----------------------------------------

#' @export
summary.sensiSamp <- function(object, ...){
    simu <- nrow(object$samp.model.estimates)
    sig <- object$sign.analysis
    sig$perc.sign.intercept <- sig$perc.sign.intercept * 100
    sig$perc.sign.slope <- sig$perc.sign.slope * 100
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
