#' Graphical sensitive analysis for comparative methods
#'
#' \code{sensi_plot} Plot results from \code{samp_influ_phyloglm},
#' \code{influ_pgls}
#' @aliases sensi_plot
#' @param x output from \code{samp_gls}, \code{influ_gls}
#' @param graphs choose which graphs should be printed on the output ("all", 1,2,3
#' or 4)
#' @param param choose which model parameter should be ploted  ("intercept" or "slope")
#' @export

### Start:
sensi_plot <- function(x,graphs="all",param="slope"){

### Basic error checking:
if (x[[1]] != "samp_phylolm" & x[[1]] != "samp_phyloglm" &
            x[[1]]!= "influ_phylolm" & x[[1]] != "influ_phyloglm")
        stop("x must be an output from one off these functions: samp_phyloglm,
             samp_phylolm, influ_phylolm and influ_phyloglm")

### samp_phylolm or samp_phylolm output:
        if (x[[1]] == "samp_phylolm" | x[[1]] == "samp_phyloglm")
                plot_samp_phylolm(x,graphs,param)

### influ_phylolm or influ_phylolm output:
if (x[[1]] == "influ_phylolm" | x[[1]] == "influ_phyloglm")
        plot_influ_phylolm(x,graphs,param)
}

