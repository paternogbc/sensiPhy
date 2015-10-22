#' Graphical diagnostics for class 'sensiIntra'
#'
#' \code{plot_tree.intra_phylm} Plot results from \code{tree_phylm},
#' \code{intra_phylm} and \code{intra_phyglm}
#' @param x output from \code{tree_phylm}, \code{tree_phyglm},
#' \code{intra_phylm} or \code{intra_phyglm}
#' @param graphs choose which graph should be printed in the output ("all", 1, 2 or 3)
#' @param ... further arguments to methods
#' @importFrom ggplot2 scale_color_manual geom_histogram geom_abline geom_density 
#' geom_vline xlab geom_point theme
#' @author Caterina Penone and Gustavo Paterno
#' @seealso \code{\link[ggplot2]{ggplot}}, \code{\link[sensiPhy]{tree_phylm}}
#' \code{\link[sensiPhy]{intra_phylm}}
#' @details For 'x' from \code{tree_phylm}, \code{tree_phyglm},
#' \code{intra_phylm} or \code{intra_phyglm}:
#' 
#' Graphs 1 and 2: Distribution of estimated slopes and intercepts for each tree (for \code{tree_phylm}) or 
#' value generated within a given interval (\code{intra_phylm})
#' Red vertical line represents the mean slope or intercept for all models. 
#'
#' Graph 3: Scatterplot with mean regression (black line) and standard deviation of the regression (blue dotted lines).
#' 
#' @examples
#' \dontrun{
#' library(sensiPhy)
#' 
#' }
#' 
#' @importFrom grid unit 
#' @export

sensi_plot.sensiIntra <- function(x, graphs="all", ...){
  
  
  # nulling variables
  sd <- formula <- slope <- ..density.. <- intercept <- NULL
  predictor <- response <-  s3 <- yy <- linety <- NULL

    resp<-all.vars(x$formula)[1]
    pred<-all.vars(x$formula)[2]
    dat<-data.frame("response"=x$datas[,resp],"predictor"=x$datas[,pred])
    result <- x$model_results
    statm<- x$stats
    slope.0 <-  as.numeric(statm[4,3])
    intercept.0 <-  as.numeric(statm[1,3])
    model_results<-x$model_results
    
    
    xf <- dat[, 2]
    yf <- plogis(intercept.0 + slope.0 * xf)
    yp <-plogis((intercept.0+statm[1,4]) + (slope.0 * xf + statm[4,4])) 
    ym <-plogis((intercept.0-statm[1,4]) + (slope.0 * xf - statm[4,4]))
    
    plot_data <- data.frame("xf" = c(xf,xf,xf),
                            "yy" = c(yf, yp, ym),
                            linety = rep(c("Mean","SD1","SD2"),each = length(yf)))

    #Distribution of estimated slopes:
    s1 <- ggplot2::ggplot(model_results,aes(x=slope,y=..density..),environment=environment())+
      geom_histogram(fill="lightyellow", alpha=.9,colour="grey60", size=.2) +
      geom_density(size=.2) +
      geom_vline(xintercept = slope.0,color="red",linetype=2,size=.7)+
      xlab("Estimated Slopes")

    #Distribution of estimated intercepts:
    i1 <- ggplot2::ggplot(model_results,aes(x=intercept,y=..density..),environment=environment())+
      geom_histogram(fill="lightyellow", alpha=.9,colour="grey60", size=.2) +
      geom_density(size=.2) +
      geom_vline(xintercept = intercept.0,color="red",linetype=2,size=.7)+
      xlab("Estimated Intercepts")
 
    #third plot: data visualisation
    s2 <- ggplot2::ggplot(dat,aes(x=predictor,y=response))+
      geom_point(size=3,alpha=.8)+
      theme(panel.background=element_rect(fill="white",colour="black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position="none")

    if(length(class(x)) == 1){
      s2.out <- s2 +
        geom_abline(intercept = intercept.0, slope=slope.0, aes(colour="mean"),size=1)+
        geom_abline(intercept = intercept.0+statm[1,4], slope=slope.0+statm[4,4], aes(color="sd Tree Uncert"),linetype=2,size=1,show_guide = T)+
        geom_abline(intercept = intercept.0-statm[1,4], slope=slope.0-statm[4,4], aes(color="sd Tree Uncert"),linetype=2,size=1,show_guide = T)+
        scale_color_manual("",values = c("black","blue"),guide=F)
    }
    
    
    if(length(class(x)) == 2){
      s2.out <- s2 + geom_line(data = plot_data, aes(x = xf, y = yy, linetype=linety, col=linety),size=1)+
        scale_linetype_manual(values = c(1,2,2)) +
        scale_color_manual("",values = c("black","blue","blue"),guide=F)
                   
    }
    
    
    
    ### Plotting:
    if (graphs=="all")
      suppressMessages(print(multiplot(s1,s2.out,i1,cols=2)))
    if (graphs==1)
      suppressMessages(print(s1))
    if (graphs==2)
      suppressMessages(print(i1))
    if (graphs==3)
      suppressMessages(print(s2.out))

}

#' Graphical diagnostics for class 'sensiTree'
#'
#' \code{plot_tree.intra_phylm} Plot results from \code{tree_phylm},
#' \code{intra_phylm} and \code{intra_phyglm}
#' @param x output from \code{tree_phylm}, \code{tree_phyglm},
#' \code{intra_phylm} or \code{intra_phyglm}
#' @param graphs choose which graph should be printed in the output ("all", 1, 2 or 3)
#' @param ... further arguments to methods
#' @importFrom ggplot2 scale_color_manual geom_histogram geom_abline geom_density 
#' geom_vline xlab geom_point theme
#' @author Caterina Penone and Gustavo Paterno
#' @seealso \code{\link[ggplot2]{ggplot}}, \code{\link[sensiPhy]{tree_phylm}}
#' \code{\link[sensiPhy]{intra_phylm}}
#' @details For 'x' from \code{tree_phylm}, \code{tree_phyglm},
#' \code{intra_phylm} or \code{intra_phyglm}:
#' 
#' Graphs 1 and 2: Distribution of estimated slopes and intercepts for each tree (for \code{tree_phylm}) or 
#' value generated within a given interval (\code{intra_phylm})
#' Red vertical line represents the mean slope or intercept for all models. 
#'
#' Graph 3: Scatterplot with mean regression (black line) and standard deviation of the regression (blue dotted lines).
#' 
#' @importFrom grid unit 
#' @export
sensi_plot.sensiTree <- function(x, graphs = "all", ...){
    sensi_plot.sensiIntra(x, graphs = "all", ...)
}
