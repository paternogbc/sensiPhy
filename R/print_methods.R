### Summary method for class: data.phy:-----------------------------------------
#' @export
print.data.phy <- function (x, ...) 
{
    cat("Comparative dataset of", nrow(x$data), "taxa with", 
        ncol(x$data), "variables", "\n")
    if(class(x$phy) == "phylo"){
        cat("-Phylogeny:", "\n")
        cat("   ", length(x$phy$tip.label), " tips, ", x$phy$Nnode, 
            " internal nodes\n  ", sep = "")
        str(x$phy$tip.label)
    }
    if(class(x$phy) == "multiPhylo"){
        cat("-Multiphylo: \n")
        print(x$phy)
    }
    cat("-Dataset:", "\n")
    str(x$data)
    if (length(x$dropped) > 0) {
        cat("Dropped taxa: \n")
        cat(length(x$dropped), "taxa were dropped from phylogeny or datset \n")
    }
}  

### Print method for class: sensiIntra:-----------------------------------------
#' @export
print.sensiIntra <- function (x, ...) 
{
    cat("Sensitivity analysis for intraspecific variation\n")
    cat("Call:\n")
    print(x$call)
    cat("\n")
    cat("Transformation: \n")
    cat("X axis: ") 
    print(x$y.transf)
    cat("Y axis: ")
    print(x$x.transf)
    cat("Number of observations: ", nrow(x$data), "\n")
    cat("Number of simulations: ", nrow(x$model_results), "\n")
    cat(message("use summary(x) and sensi_plot(x) to check results"))
    
}

### Print method for class: sensiSamp:----------------------------------------
#' @export
print.sensiSamp <- function (x, ...) 
{
    cat("Sensitivity analysis for sampling uncertainty\n")
    cat("Call:\n")
    print(x$call)
    cat("\n")
    cat("Breaks: ", x$sign.analysis$percent_sp_removed, "\n")
    cat("Number of simulations: ", nrow(x$samp.model.estimates), "\n")
    cat(message("use summary(x) and sensi_plot(x) to check results"))
}

### Print method for class: sensiClade:----------------------------------------
#' @export
print.sensiClade <- function (x, ...) 
{
    cat("Sensitivity analysis of influential clades \n")
    cat("Call:\n")
    print(x$call)
    cat("\n")
    cat("Analyzed clades: ")
    cat(x$clade.model.estimates$clade, sep = "; ",  "\n")
    cat(message("use summary(x) and sensi_plot(x) to check results"))
}

### Print method for class: sensiTree:-----------------------------------------
#' @export
print.sensiTree <- function (x, ...) 
{
    cat("Sensitivity analysis for phylogenetic uncertainty \n")
    cat("Call:\n")
    print(x$call)
    cat("\n")
    cat("Number of observations: ", nrow(x$data), "\n")
    cat("Number of simulations: ", nrow(x$model_results), "\n")
    cat(message("use summary(x) and sensi_plot(x) to check results"))
    
}

### Print method for class: sensiInflu:-----------------------------------------
#' @export
print.sensiInflu <- function (x, ...) 
{
    cat("Sensitivity analysis for influential species \n")
    cat("Call:\n")
    print(x$call)
    cat("\n")
    cat("Number of simulations:", nrow(influ$influ.model.estimates), "\n")
    cat(message("use summary(x) and sensi_plot(x) to check results"))
    
}
