### 0. Print method for class: data.phy:-----------------------------------------
#' @importFrom utils str
#' @export

print.data.phy <- function(x, ...) 
{
    cat("Comparative dataset of", nrow(x$data), "taxa with", 
        ncol(x$data), "variables", "\n")
    if (class(x$phy) == "phylo") {
        cat("-Phylogeny:", "\n")
        cat("   ", length(x$phy$tip.label), " tips, ", x$phy$Nnode, 
            " internal nodes\n  ", sep = "")
        utils::str(x$phy$tip.label)
    }
    if (class(x$phy) == "multiPhylo") {
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

### 1. PRINT FOR PGLS--------------------------------------------------------------
### sensiIntra:-----------------------------------------
#' @export
print.sensiIntra <- function(x, ...) 
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

### sensiSamp:----------------------------------------
#' @export
print.sensiSamp <- function(x, ...) 
{
    cat("Sensitivity analysis for sampling uncertainty\n")
    cat("Call:\n")
    print(x$call)
    cat("\n")
    cat("Breaks: ", x$sign.analysis$percent_sp_removed, "\n")
    cat("Number of simulations: ", nrow(x$samp.model.estimates), "\n")
    cat(message("use summary(x) and sensi_plot(x) to check results"))
}

### sensiClade:----------------------------------------
#' @export
print.sensiClade <- function(x, ...) 
{
    cat("Sensitivity analysis of influential clades \n")
    cat("Call:\n")
    print(x$call)
    cat("\n")
    cat("Analyzed clades: ")
    cat(x$sensi.estimates$clade, sep = "; ",  "\n")
    cat(message("use summary(x) and sensi_plot(x) to check results"))
}

### sensiTree:-----------------------------------------
#' @export
print.sensiTree <- function(x, ...) 
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
print.sensiInflu <- function(x, ...) 
{
    cat("Sensitivity analysis for influential species \n")
    cat("Call:\n")
    print(x$call)
    cat("\n")
    cat("Number of simulations:", nrow(x$influ.model.estimates), "\n")
    cat(message("use summary(x) and sensi_plot(x) to check results"))
    
}

### 2. PRINT PHYLOGENETIC SIGNAL------------------------------------------
### clade.physig:---------------------------------------
#' @export
print.clade.physig <- function(x, ...) 
{
    cat("Sensitivity analysis of influential clades for Phylogenetic signal \n")
    cat("\n")
    cat("Call:\n")
    print(x$call)
    cat("\n")
    cat("Number of randomizations: ")
    cat(x$call$n.sim)
    cat("\n")
    cat("Clades analysed: ")
    cat(x$sensi.estimates$clade)
    cat(message("use summary(x) and sensi_plot(x) to check results"))
    cat(message("use x$sensi.estimates to access sensitivity analysis data"))
}

### influ.physig:---------------------------------------
#' @export
print.influ.physig <- function(x, ...) 
{
  cat("Sensitivity analysis of influential species for Phylogenetic signal \n")
  cat("\n")
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("Number of species: ")
  cat(nrow(x$data))
  cat("\n")
  cat(message("use summary(x) and sensi_plot(x) to check results"))
  cat(message("use x$influ.physig.estimates to access sensitivity analysis data"))
}

### samp.physig:---------------------------------------
#' @export
print.samp.physig <- function(x, ...) 
{
  cat("Sensitivity analysis of sampling uncertainty for Phylogenetic Signal \n")
  cat("\n")
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("Number of randomizations: ")
  cat(x$call$nsim)
  cat("\n")
  cat(message("use summary(x) and sensi_plot(x) to check results"))
  cat(message("use x$influ.physig.estimates to access sensitivity analysis data"))
}

### tree.physig:---------------------------------------
#' @export
print.tree.physig <- function(x, ...) 
{
  cat("Sensitivity analysis of phylogenetic uncertainty for Phylogenetic signal \n")
  cat("\n")
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("Number of trees evaluated: ")
  cat(nrow(x$tree.physig.estimates))
  cat("\n")
  cat(message("use summary(x) and sensi_plot(x) to check results"))
  cat(message("use x$tree.physig.estimates to access sensitivity analysis data"))
}

### intra.physig:---------------------------------------
#' @export
print.intra.physig <- function(x, ...) 
{
  cat("Sensitivity analysis of intraspecific variability for Phylogenetic signal \n")
  cat("\n")
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("Number of simulations: ")
  cat(nrow(x$intra.physig.estimates))
  cat("\n")
  cat(message("use summary(x) and sensi_plot(x) to check results"))
  cat(message("use x$tree.intra.estimates to access sensitivity analysis data"))
}


### 3. PRINT interaction PGLS------------------------------------------
#' @export
print.sensiTree_Clade <- function(x, ...) 
{
  cat("Sensitivity analysis for interaction between tree:clade \n")
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("Analyzed clades: ")
  cat(unique(x$sensi.estimates$clade), sep = "; ",  "\n")
  cat(message("use summary(x) and sensi_plot(x) to check results"))
}

#' @export
print.sensiTree_Influ <- function(x, ...) 
{
  cat("Sensitivity analysis for interaction between tree:influ \n")
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat(message("use summary(x) and sensi_plot(x) to check results"))
}
