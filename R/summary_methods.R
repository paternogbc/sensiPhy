### Summary method for class: sensiClade:--------------------------------------

#' @export
summary.sensiClade <- function(object, ...){
    ord.inter <- order(object$clade.model.estimates$intercept.perc, 
                       decreasing = TRUE)
    inter <- object$clade.model.estimates[ord.inter, c(1,2,3,4,5)]
    colnames(inter) <- c("Clade removed", "Intercept", "DFintercept",
                         "change (%)", "Pval")
    ord.slope <- order(object$clade.model.estimates$slope.perc, 
                       decreasing = TRUE)
    slope <- object$clade.model.estimates[ord.slope, c(1,6,7,8,9)]
    colnames(slope) <- c("Clade removed", "Slope", "DFslope", 
                         "change (%)", "Pval")
    res <- list(slope, inter)
    names(res) <- c("Summary for Slope", "Summary for Intercept")
    return(res)
}


### Summary method for class: sensiInflu:--------------------------------------

#' @export
summary.sensiInflu <- function(object, ...){
    sp.slope <- object$influential.species$influ.sp.slope
    rows.slope <- match(sp.slope, object$influ.model.estimates$species)
    slope <- object$influ.model.estimates[rows.slope, c(1,6,7,8,9)]
    ord.slope <- order(slope$DFslope, 
                       decreasing = TRUE)
    slope <- slope[ord.slope,]
    rownames(slope) <- NULL
    colnames(slope) <- c("Species removed", "Slope", "DFslope", "Change(%)", "Pval")
    
    sp.inter <-object$influential.species$influ.sp.intercept
    rows.inter <- match(sp.inter, object$influ.model.estimates$species)
    inter <- object$influ.model.estimates[rows.inter, c(1,2,3,4,5)]
    ord.inter <- order(inter$DFintercept, 
                       decreasing = TRUE)
    inter <- slope[ord.inter,]
    rownames(inter) <- NULL
    colnames(inter) <- c("Species removed", "Intercept", "DFintercept", "Change(%)", "Pval")
    sp.inter <-object$influential.species$influ.sp.intercept
    rows.inter <- match(sp.inter, object$influ.model.estimates$species)
    inter <- object$influ.model.estimates[rows.inter, c(1,6,7,8,9)]
    names(object$influ.model.estimates)
    
    res <- list("Influential species for the Slope" = sp.slope, "Estimates" = slope,
                "Influential species for the Intercept" = sp.inter, "Estimates" = inter)
    return(res)
    
}

### Summary method for class: sensiSamp:----------------------------------------

#' @export
summary.sensiSamp <- function(object, ...){
    simu <- nrow(object$samp.model.estimates)
    sig <- object$sign.analysis
    sig$perc.sign.intercept <- sig$perc.sign.intercept * 100
    sig$perc.sign.slope <- sig$perc.sign.slope * 100
    names(sig) <- c("% of Species Removed", 
                    "% of Significant Intercept",
                    "% of Significant Slope")
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