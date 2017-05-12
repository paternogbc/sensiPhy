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
#'  \item{\code{\link[sensiPhy]{sensi_plot.sensiClade}}}
#'  \item{\code{\link[sensiPhy]{sensi_plot.sensiInflu}}}
#'  \item{\code{\link[sensiPhy]{sensi_plot.sensiSamp}}}
#'  \item{\code{\link[sensiPhy]{sensi_plot.sensiIntra}}}
#'  \item{\code{\link[sensiPhy]{sensi_plot.sensiTree}}}
#'  \item{\code{\link[sensiPhy]{sensi_plot.sensiIntra_Tree}}}
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
