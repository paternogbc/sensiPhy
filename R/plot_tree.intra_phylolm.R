#' Graphical diagnostics for \code{tree_phylolm} and \code{intra_phylolm}
#'
#' \code{plot_tree.intra_phylolm} Plot results from \code{tree_phylolm} and
#' \code{intra_phylolm}
#' @param x output from \code{tree_phylolm} or \code{intra_phylolm}
#' @param param choose which parameter ("intercept" or "slope" should be printed)
#' @param graphs choose which graph should be printed on the output ("all", 1, 2 or 3)
#' @importFrom ggplot2 scale_color_manual geom_histogram geom_abline geom_density 
#' geom_vline xlab geom_point
#' @author Caterina Penone and Gustavo Paterno
#' @seealso \code{\link[ggplot2]{ggplot}}, \code{\link[sensiPhy]{tree_phylolm}}
#' \code{\link[sensiPhy]{intra_phylolm}}
#' @details For 'x' from tree_phylolm or intra_phylolm:
#' 
#' Graphs 1 and 2: Distribution of estimated slopes and intercepts for each tree (for \code{tree_phylolm}) or 
#' value generated within a given interval (\code{intra_phylolm})
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

plot_tree.intra_phylolm <- function(x,graphs="all",param=NULL){
  
  
  # nulling variables
  sd <- formula <- slope <- ..density.. <- intercept <- NULL
  predictor <- response <-  s3 <- NULL

  
  ### Error check:
  if (x[[1]] != "tree_phylolm" & x[[1]] != "intra_phylolm")
    stop("x must be an output from tree_phylolm or intra_phylolm!")
  else
    
    resp<-all.vars(formula)[1]
    pred<-all.vars(formula)[2]
    dat<-data.frame("response"=x$datas$resp,"predictor"=x$datas$pred)
    result <- x$model_results
    statm<- x$stats
    slope.0 <-  as.numeric(statm[4,3])
    intercept.0 <-  as.numeric(statm[1,3])
    model_results<-x$model_results

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
      geom_abline(intercept = intercept.0, slope=slope.0, aes(colour="mean"),size=1)+
      geom_abline(intercept = intercept.0+statm[1,4], slope=slope.0+statm[4,4], aes(color="sd Tree Uncert"),linetype=2,size=1,show_guide = T)+
      geom_abline(intercept = intercept.0-statm[1,4], slope=slope.0-statm[4,4], aes(color="sd Tree Uncert"),linetype=2,size=1,show_guide = T)+
      scale_color_manual("",values = c("black","blue"),guide=F)

    ### Plotting:
    if (graphs=="all")
      suppressMessages(print(multiplot(s1,s2,i1,cols=2)))
    if (param == "slope" & graphs==1)
      suppressMessages(print(s1))
    if (param == "slope" & graphs==2)
      suppressMessages(print(s3))
    if (param == "slope" & graphs==3)
      suppressMessages(print(i1))

}
