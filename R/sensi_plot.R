#' Graphical sensitivity analysis for comparative methods
#'
#' Generic function for plotting results from
#' any sensitivity analysis performed with 'sensiPhy' (see \code{\link[sensiPhy]{samp_phylolm}},
#' \code{\link[sensiPhy]{samp_phyloglm}}, \code{\link[sensiPhy]{influ_phylolm}},
#' \code{\link[sensiPhy]{influ_phyloglm}}) 
#' @aliases sensi_plot
#' @param x output from \code{samp_phylolm}, \code{samp_phyloglm},
#'  \code{influ_phylolm} or \code{influ_phyloglm}
#' @param graphs choose which graph should be printed ("all", 1,2,3 or 4, see 
#' details). Default is \code{"all"}.
#' @param param choose which model parameter should be ploted  ("intercept" or 
#' "slope"). Default is "slope".
#' @details sensi_plot prints differents sets of graphs depending on the function
#' that has generated 'x'.
#' 
#' @section For 'x' from samp_phylolm or samp_phyloglm:
#' 
#' Graph 1: Estimated slopes or intercepts for each simution across  
#' percentage of species removed. Colours represent percentage 
#' of change in comparison with the full model (blue = lower than 5, orange = 
#' between 5 and 10 and red = higher than 10).
#' The red horizontal line represents the original slope or 
#' intercept from the full model (with all species). 
#' 
#' Graph 2: The proportion of estimated slopes and intercepts in each category 
#' across the percentage of species removed.
#' 
#' Graph 3: Estimated phylogenetic model parameter for each simulation across
#' the percentage of species removed.
#' 
#' Graph 4: The percentage of significant slopes or intercepts across the 
#' percentage of species removed.  
#' 
#' @section  For 'x' from influ_phylolm or influ_phyloglm:
#' Graph 1: Distribution of estimated slopes or intercepts for each 
#' simulation (leave-one-out deletion). Red vertical line represents the original
#' slope or intercept from the full model (with all species). 
#' 
#' Graph 2: Original regression plot (\eqn{trait~predictor}). Standardized 
#' difference in slope or intercept is represented by a continous colour scale, 
#' ranging from black (low \code{sDFintercept} or \code{sDFslope} values) to red
#' (high \code{sDFintercept} or \code{sDFslope} values).
#' 
#' Graph 3: Distribution of Standardized difference in slope or intercept. Red 
#' colour indicates inbfluential species (with a standardised difference above 
#' the value of \code{cutoff}).
#' 
#' Graph 4: Ditribution of the percentage of change in slope or intercept.
#' @examples
#' \dontrun{
#' library(sensiPhy)
#'
#' # Generating tree and traits:
#' set.seed(2468)
#' tree <- rtree(100)
#' pred<- rTraitCont(tree,root.value=0,sigma=1,model="BM")
#' cont_trait <- pred + rnorm(100,1,4)
#' dat<-data.frame(pred,cont_trait)
#' 
#' # Running sensitivity analysis:
#' fit1<-influ_phylolm(cont_trait ~ pred,data = dat,phy = tree)
#' fit2<-samp_phylolm(cont_trait ~ pred,data = dat,phy = tree)
#' 
#' # Plotting the results:
#' sensi_plot(fit1)
#' sensi_plot(fit2)
#' sensi_plot(fit2, graphs = 1, param = "intercept")
#' sensi_plot(fit2, graphs = 4, param = "slope")
#' }
#' @author Gustavo Paterno
#' @seealso \code{\link[sensiPhy]{samp_phylolm}},
#' \code{\link[sensiPhy]{samp_phyloglm}}, \code{\link[sensiPhy]{influ_phylolm}},
#' \code{\link[sensiPhy]{influ_phyloglm}}
#' @importFrom grid grid.newpage pushViewport grid.layout viewport
#' @export

### Start:
sensi_plot <- function(x,graphs="all",param="slope"){

### Basic error checking:
if (x[[1]] != "samp_phylolm" & x[[1]] != "samp_phyloglm" &
            x[[1]]!= "influ_phylolm" & x[[1]] != "influ_phyloglm")
        stop("x must be an output from one of these functions: samp_phyloglm,
             samp_phylolm, influ_phylolm and influ_phyloglm")

### samp_phylolm or samp_phylolm output:
if (x[[1]] == "samp_phylolm" | x[[1]] == "samp_phyloglm")
    plot_samp_phylolm(x,graphs,param)

### influ_phylolm or influ_phylolm output:
if (x[[1]] == "influ_phylolm" | x[[1]] == "influ_phyloglm")
    plot_influ_phylolm(x,graphs,param)
}
