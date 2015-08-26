#' Graphical sensitivity analysis for comparative methods
#'
#' Generic function for plotting results from
#' any sensitivity analysis performed with 'sensiPhy'
#' @aliases sensi_plot
#' @param x any output from the sensiPhy package. 
#' @param ... further arguments to methods 
#' @details sensi_plot recognize and print different sets of graphs depending 
#' on the function that has generated 'x'. See the links bellow for details about
#' the graphs generated for each sensiPhy function:
#' \itemize{
#'  \item{\code{\link[sensiPhy]{sensi_plot.sensiClade}}}
#'  \item{\code{\link[sensiPhy]{sensi_plot.sensiInflu}}}
#'  \item{\code{\link[sensiPhy]{sensi_plot.sensiSamp}}}
#'  \item{\code{\link[sensiPhy]{sensi_plot.sensiIntra}}}
#'  \item{\code{\link[sensiPhy]{sensi_plot.sensiTree}}}
#'  }
#' 
#' @author Gustavo Paterno
#' @importFrom grid grid.newpage pushViewport grid.layout viewport
#' @export

### Start:
sensi_plot <- function(x, ...){
    UseMethod("sensi_plot")
}
