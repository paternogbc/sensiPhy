#' Graphical sensitivity analysis for comparative methods
#'
#' Generic function for plotting results from
#' any sensitivity analysis performed with 'sensiPhy'
#' @aliases sensi_plot
#' @param x any output from the sensiPhy package. 
#' @param ... further arguments to methods 
#' @details sensi_plot recognize and print different sets of graphs depending 
#' on the function that generated 'x'. See the links below for details about
#' the graphs generated for each sensiPhy function:
#' \itemize{
#' \strong{PGLS regressions (single uncertainty):}
#'  \item{clade_phylm: \code{\link[sensiPhy]{sensi_plot.sensiClade}}}
#'  \item{influ_phylm: \code{\link[sensiPhy]{sensi_plot.sensiInflu}}}
#'  \item{samp_phylm: \code{\link[sensiPhy]{sensi_plot.sensiSamp}}}
#'  \item{intra_phylm: \code{\link[sensiPhy]{sensi_plot.sensiIntra}}}
#'  \item{tree_phylm: \code{\link[sensiPhy]{sensi_plot.sensiTree}}}
#'   }
#'  \strong{PGLS regressions (interacting uncertainties):}
#'  \itemize{
#'  \item{tree_intra_phylm: \code{\link[sensiPhy]{sensi_plot.sensiTree_Intra}}}
#'  \item{intra_clade_phylm: \code{\link[sensiPhy]{sensi_plot.sensiIntra_Clade}}}
#'  \item{intra_influ_phylm: \code{\link[sensiPhy]{sensi_plot.sensiIntra_Influ}}}
#'  \item{intra_samp_phylm: \code{\link[sensiPhy]{sensi_plot.sensiIntra_Samp}}}
#'  \item{tree_clade_phylm: \code{\link[sensiPhy]{sensi_plot.sensiTree_Clade}}}
#'  \item{tree_influ_phylm: \code{\link[sensiPhy]{sensi_plot.sensiTree_Influ}}}
#'  \item{tree_samp_phylm: \code{\link[sensiPhy]{sensi_plot.sensiTree_Samp}}}
#'  }
#'  \itemize{
#'  \strong{Phylogenetic signal:}
#'  \item{clade_physig: \code{\link[sensiPhy]{sensi_plot.clade.physig}}}
#'  \item{influ_physig: \code{\link[sensiPhy]{sensi_plot.influ.physig}}}
#'  \item{samp_physig:  \code{\link[sensiPhy]{sensi_plot.samp.physig}}}
#'  \item{tree_physig:  \code{\link[sensiPhy]{sensi_plot.tree.physig}}}
#'  \item{intra_physig:  \code{\link[sensiPhy]{sensi_plot.intra.physig}}}
#'  }
#' 
#' @author Gustavo Paterno
#' @importFrom grid grid.newpage pushViewport grid.layout viewport
#' @references The function `multiplot`, developped by Winston Chang, is used inside sensi_plot
#' to print multiple graphs in one frame. 
#' The source code is available here:
#' \url{http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/}
#' @export

### Start:
sensi_plot <- function(x, ...){
  UseMethod("sensi_plot")
}