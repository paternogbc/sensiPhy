#' Graphical diagnostics for class 'sensiInflu'
#'
#' \code{plot_influ_phylm} Plot results from \code{influ_phylm} and 
#' \code{influ_phyglm}
#' @param x output from \code{influ_phylm}
#' @param graphs choose which graph should be printed on the output ("all", 1,2,3 or 4)
#' @param param choose which parameter ("intercept" or "estimate" should be printed)
#' @param ... further arguments to methods
#' @importFrom ggplot2 aes geom_histogram geom_density geom_vline 
#' xlab theme element_text geom_point scale_colour_gradient element_rect ylab xlab
#' ggtitle element_blank
#' @author Gustavo Paterno
#' @seealso \code{\link[ggplot2]{ggplot}}
#' @details For 'x' from influ_phylm or influ_phyglm:
#' 
#' \strong{Graph 1:} Distribution of estimated slopes (estimates) or intercepts for each 
#' simulation (leave-one-out deletion). Red vertical line represents the original
#' slope or intercept from the full model (with all species). 
#' 
#' \strong{Graph 2:} Original regression plot (\eqn{trait~predictor}). Standardized 
#' difference in slope or intercept is represented by a continuous size scale. 
#' The names of the most influential species (sDF > cutoff) are ploted in the
#' graph. 
#' 
#' \strong{Graph 3:} Distribution of standardized difference in slope or intercept. Red 
#' colour indicates inbfluential species (with a standardised difference above 
#' the value of \code{cutoff}).
#' 
#' \strong{Graph 4:} Distribution of the percentage of change in slope or intercept.
#' @importFrom grid unit 
#' @export

### Start:
sensi_plot.sensiInflu <- function(x, graphs="all", param="estimate", ...){

# nulling variables:------------------------------------------------------------
estimate <- ..density.. <- intercept <- sDIFestimate <- estimate.perc <- NULL 
intercept.perc <- sDIFintercept <- species <-  NULL
        
        ### Organizing values:
        result <- x$sensi.estimates
        mappx <- x$formula[[3]]
        mappy <- x$formula[[2]]
        vars <- all.vars(x$formula)
        intercept.0 <-  as.numeric(x$full.model.estimates$coef[1])
        estimate.0     <-  as.numeric(x$full.model.estimates$coef[2])
        cutoff      <-  x$cutoff

        ### Removing species with error:
        if (isTRUE(class(x$errors) != "character" )){
                x$data <- x$data[-as.numeric(x$errors),]
        }
        result.tab <- data.frame(x$sensi.estimates,x$data[all.vars(x$formula)])

        ### Plots:
        # Distribution of estimated slopes:
        s1 <- ggplot2::ggplot(result,aes(x=estimate))+
                geom_histogram(fill="yellow",colour="black", size=.2,
                               alpha = .3) +
                geom_vline(xintercept = estimate.0,color="red",linetype=2,size=.7)+
                xlab("Estimates")+
                ylab("Frequency")+
                theme(axis.title=element_text(size=12),
                    axis.text = element_text(size=12),
                    panel.background = element_rect(fill="white",
                                                  colour="black"))
        
        # Distribution of estimated intercepts:
        i1 <- ggplot2::ggplot(result,aes(x=intercept))+
                geom_histogram(fill="yellow",colour="black", size=.2,
                               alpha = .3) +
                geom_vline(xintercept = intercept.0,color="red",linetype=2,size=.7)+
                xlab("Intercepts")+
                ylab("Frequency")+
                theme(axis.title=element_text(size=12),
                    axis.text = element_text(size=12),
                    panel.background = element_rect(fill="white",
                                                  colour="black"))
        
        # Original plot with Standardized DIFestimate as colour gradient
        s2 <- ggplot2::ggplot(result.tab, aes_string(y = mappy, x = mappx),
                              environment = environment())+
            geom_point(data = result.tab,
                       aes(size = abs(sDIFestimate)), alpha = .8)+
            ggplot2::scale_size_continuous(name = "|sDF|", range = c(1, 6))+
            ggplot2::geom_text(aes(label =  ifelse(abs(sDIFestimate) > cutoff, 
                                  as.character(species), ""),
                        vjust = 0, hjust = 0,
                        color = "red", size = .7), show.legend = F,
                        fontface = "bold") + 
            theme(legend.key.width = unit(.2,"cm"),
                  panel.background=element_rect(fill="white", colour = "black"),
                  legend.text = element_text(size = 12),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank())+
            scale_x_continuous(expand = c(.2, .2)) +
            ggtitle("Standardized Difference in estimate")+
            theme(axis.title = element_text(size = 12),
                  axis.text = element_text(size = 12),
                  panel.background = element_rect(fill = "white",
                                                  colour = "black"))
             
        # Original plot with Standardized DIFintercept as size gradient
        i2<-ggplot2::ggplot(result.tab,aes_string(y = mappy, x = mappx),
                            environment = environment())+
            geom_point(data = result.tab,
                       aes(size = abs(sDIFintercept)), alpha = .8)+
            ggplot2::scale_size_continuous(name = "sDF", range = c(1, 6))+
            ggplot2::geom_text(aes(label = ifelse(abs(sDIFintercept) > cutoff, 
                                                as.character(species), ""), 
                                 vjust = 0, hjust = 0, color = "red",
                                 size = .7), show.legend = F,  fontface = "bold") +
            theme(legend.key.width = unit(.2,"cm"),
                  panel.background=element_rect(fill="white",colour="black"),
                  legend.text = element_text(size=12),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank())+
          scale_x_continuous(expand = c(.2, .2)) +
            ggtitle("Standardized Difference in Intercept")+
            theme(axis.title=element_text(size=12),
                  axis.text = element_text(size=12),
                  panel.background = element_rect(fill="white",
                                                  colour="black"))

        # Influential points for slope estimate
        s3 <- ggplot2::ggplot(result,aes(x=sDIFestimate))+
                geom_histogram(fill="red",color="black",binwidth=.5) +
                xlab("Standardized Difference in estimate")+
                ylab("Frequency")+
                geom_histogram(data=subset(result,sDIFestimate<cutoff&sDIFestimate>-cutoff),
                               colour="black", fill="white",binwidth=.5)+
                geom_vline(xintercept = -cutoff,color="red",linetype=2,size=.7)+
                geom_vline(xintercept = cutoff,color="red",linetype=2,size=.7)+
                theme(axis.title=element_text(size=12),
                    axis.text = element_text(size=12),
                    panel.background = element_rect(fill="white",
                                                  colour="black"))

        # Influential points for intercept estimate
        i3 <- ggplot2::ggplot(result,aes(x=sDIFintercept))+
                geom_histogram(fill="red",color="black",binwidth=.5) +
                xlab("Standardized Difference in Intercept")+
                ylab("Frequency")+
                geom_histogram(data=subset(result,sDIFestimate<cutoff&sDIFestimate>-cutoff),
                               colour="black", fill="white",binwidth=.5)+
                geom_vline(xintercept = -cutoff,color="red",linetype=2,size=.7)+
                geom_vline(xintercept = cutoff,color="red",linetype=2,size=.7)+
                theme(axis.title=element_text(size=12),
                    axis.text = element_text(size=12),
                    panel.background = element_rect(fill="white",
                                              colour="black"))                

        # Distribution of slope.perc:

        s4 <- ggplot2::ggplot(result,aes(x=estimate.perc,y=..density..))+
                geom_histogram(data = result,
                               colour="black", fill="yellow",
                               alpha = .3)+
                xlab("% of change in estimate")+
                ylab("Frequency") +
                theme(axis.title=element_text(size=12),
                       axis.text = element_text(size=12),
                       panel.background = element_rect(fill="white",
                                                       colour="black"))

        # Distribution of slope.perc:
        i4 <- ggplot2::ggplot(result,aes(x=intercept.perc,y=..density..))+
                geom_histogram(
                              colour="black", fill="yellow",
                              alpha = .3)+
                xlab("% of change in Intercept")+
                ylab("Frequency")+
                theme(axis.title=element_text(size=12),
                    axis.text = element_text(size=12),
                    panel.background = element_rect(fill="white",
                                                  colour="black"))

        ### Plotting:
        if (param == "estimate" & graphs == "all")
            suppressMessages(return(multiplot(s1, s3, s2, s4, cols = 2)))
        if (param == "estimate" & graphs == 1)
            suppressMessages(return(s1))
        if (param == "estimate" & graphs == 2)
            suppressMessages(return(s2))
        if (param == "estimate" & graphs == 3)
            suppressMessages(return(s3))
        if (param == "estimate" & graphs == 4)
            suppressMessages(return(s4))
        if (param == "intercept" & graphs == "all")
            suppressMessages(return(multiplot(i1, i3, i2, i4, cols = 2)))
        if (param == "intercept" & graphs == 1)
            suppressMessages(return(i1))
        if (param == "intercept" & graphs == 2)
            suppressMessages(return(i2))
        if (param == "intercept" & graphs == 3)
            suppressMessages(return(i3))
        if (param == "intercept" & graphs == 4)
            suppressMessages(return(i4))

        ### Warnings
        if (isTRUE(class(x$errors) != "character" ))
                warnings("Deletion of some species caused error. These species were not ploted")

}



#' Graphical diagnostics for class 'sensiInflu.TraitEvol'
#'
#' \code{sensi_plot.sensiTree.TraitEvol} Plot results from \code{influ_discrete} and \code{influ_continuous}.
#' @param x output from \code{influ_discrete} or \code{influ_continuous}
#' @param graphs choose which graph should be printed in the output ("all", "q12", "q21", "aic" or" "optpar")
#' @param ... further arguments to methods
#' @importFrom ggplot2 geom_histogram geom_density geom_vline xlab theme
#' @author Gijsbert Werner
#' @seealso \code{\link[ggplot2]{ggplot}}, \code{\link[sensiPhy]{influ_discrete}}
#' \code{\link[sensiPhy]{influ_continuous}}
#' 
#'The following graphs are printed. 
#'
#' \strong{Graph aicc:} Distribution of estimated AICc-values across all single-species deletions. 
#' Red vertical line represents the mean signal among all estimates. 
#' Blue vertical line represents the median signal among all estimates. 
#' 
#' \strong{Graph optpar:} Distribution of estimated values for optimisation parameter specified using 'transform' (if applicable)
#' Red vertical line represents the mean signal among all estimates. 
#' Blue vertical line represents the median signal among all estimates. 
#' 
#' Additionally, only for \code{tree_discrete} the function creates the following graphs. 
#' 
#' \strong{Graph q12:} Distribution of estimated parameter values for transition rates q12 across all single-species deletions. 
#' Red vertical line represents the mean signal among all estimates. 
#' Blue vertical line represents the median signal among all estimates. 
#' 
#' \strong{Graph q21:} Distribution of estimated parameter values for transition rates q21. 
#' Red vertical line represents the mean signal among all estimates. 
#' Blue vertical line represents the median signal among all estimates. 
#' 
#' While only for \code{tree_continuous} the function creates the following graphs. 
#' 
#' \strong{Graph sigsq:} Distribution of estimated parameter values for rate of evolution sigsq across all single-species deletions. .
#' Red vertical line represents the mean signal among all estimates. 
#' Blue vertical line represents the median signal among all estimates. 
#' 
#' \strong{Graph z0:} Distribution of estimated parameter values for z0.
#' Red vertical line represents the mean signal among all estimates. 
#' Blue vertical line represents the median signal among all estimates. 
#' @importFrom grid unit 
#' @importFrom stats plogis
#' @importFrom stats reorder
#' @export
sensi_plot.sensiInflu.TraitEvol <- function(x, graphs="all", ...){
  if(as.character(x$call[[1]])=="influ_discrete"){
    q12_fig<-ggplot2::ggplot()+
      geom_histogram(aes(x$sensi.estimates$q12),
                     fill = "yellow",colour = "black", size = .2,
                     alpha = .3)+
      geom_vline(aes(xintercept=mean(x$sensi.estimates$q12)),colour="red")+
      geom_vline(aes(xintercept=median(x$sensi.estimates$q12)),colour="blue")+
      xlab("Estimated q12") +
      ylab("Frequency") +
      theme(axis.title = element_text(size = 12),
            axis.text = element_text(size = 12),
            panel.background = element_rect(fill = "white",
                                            colour = "black"))
    q21_fig<-ggplot2::ggplot()+
      geom_histogram(aes(x$sensi.estimates$q21),
                     fill = "yellow",colour = "black", size = .2,
                     alpha = .3)+
      geom_vline(aes(xintercept=mean(x$sensi.estimates$q21)),colour="red")+
      geom_vline(aes(xintercept=median(x$sensi.estimates$q21)),colour="blue")+
      xlab("Estimated q21") +
      ylab("Frequency") +
      theme(axis.title = element_text(size = 12),
            axis.text = element_text(size = 12),
            panel.background = element_rect(fill = "white",
                                            colour = "black"))
    aicc_fig<-ggplot2::ggplot()+
      geom_histogram(aes(x$sensi.estimates$aicc),
                     fill = "yellow",colour = "black", size = .2,
                     alpha = .3)+
      geom_vline(aes(xintercept=mean(x$sensi.estimates$aicc)),colour="red")+
      geom_vline(aes(xintercept=median(x$sensi.estimates$aicc)),colour="blue")+
      xlab("Estimated AICc") +
      ylab("Frequency") +
      theme(axis.title = element_text(size = 12),
            axis.text = element_text(size = 12),
            panel.background = element_rect(fill = "white",
                                            colour = "black"))
    optpar_fig<-ggplot2::ggplot()+
      geom_histogram(aes(x$sensi.estimates$optpar),
                     fill = "yellow",colour = "black", size = .2,
                     alpha = .3)+
      geom_vline(aes(xintercept=mean(x$sensi.estimates$optpar)),colour="red")+
      geom_vline(aes(xintercept=median(x$sensi.estimates$optpar)),colour="blue")+
      xlab(paste("Estimated",x$optpar,"parameter",sep=" ")) +
      ylab("Frequency") +
      theme(axis.title = element_text(size = 12),
            axis.text = element_text(size = 12),
            panel.background = element_rect(fill = "white",
                                            colour = "black"))
    
    if (graphs=="all"){
      if(x$optpar=="none"){
        suppressMessages(return(multiplot(q12_fig,q21_fig,aicc_fig, cols=2))) #When transformation = "none" was used, don't print the visualisation for the optimisation parameter
      } else
        suppressMessages(return(multiplot(q12_fig,q21_fig,aicc_fig,optpar_fig, cols=2)))
    }
    if (graphs=="q12")
      suppressMessages(return(q12_fig))
    if (graphs=="q21")
      suppressMessages(return(q21_fig))
    if (graphs=="aic")
      suppressMessages(return(aicc_fig))
    if (graphs=="optpar")
      suppressMessages(return(optpar_fig))
  } 
  
  if(as.character(x$call[[1]])=="influ_continuous"){
    sigsq_fig<-ggplot2::ggplot()+
      geom_histogram(aes(x$sensi.estimates$sigsq),
                     fill = "yellow",colour = "black", size = .2,
                     alpha = .3)+
      geom_vline(aes(xintercept=mean(x$sensi.estimates$sigsq)),colour="red")+
      geom_vline(aes(xintercept=median(x$sensi.estimates$sigsq)),colour="blue")+
      xlab("Estimated sigsq") +
      ylab("Frequency") +
      theme(axis.title = element_text(size = 12),
            axis.text = element_text(size = 12),
            panel.background = element_rect(fill = "white",
                                            colour = "black"))
    z0_fig<-ggplot2::ggplot()+
      geom_histogram(aes(x$sensi.estimates$z0),
                     fill = "yellow",colour = "black", size = .2,
                     alpha = .3)+
      geom_vline(aes(xintercept=mean(x$sensi.estimates$z0)),colour="red")+
      geom_vline(aes(xintercept=median(x$sensi.estimates$z0)),colour="blue")+
      xlab("Estimated z0") +
      ylab("Frequency") +
      theme(axis.title = element_text(size = 12),
            axis.text = element_text(size = 12),
            panel.background = element_rect(fill = "white",
                                            colour = "black"))
    aicc_fig<-ggplot2::ggplot()+
      geom_histogram(aes(x$sensi.estimates$aicc),
                     fill = "yellow",colour = "black", size = .2,
                     alpha = .3)+
      geom_vline(aes(xintercept=mean(x$sensi.estimates$aicc)),colour="red")+
      geom_vline(aes(xintercept=median(x$sensi.estimates$aicc)),colour="blue")+
      xlab("Estimated AICc") +
      ylab("Frequency") +
      theme(axis.title = element_text(size = 12),
            axis.text = element_text(size = 12),
            panel.background = element_rect(fill = "white",
                                            colour = "black"))
    optpar_fig<-ggplot2::ggplot()+
      geom_histogram(aes(x$sensi.estimates$optpar),
                     fill = "yellow",colour = "black", size = .2,
                     alpha = .3)+
      geom_vline(aes(xintercept=mean(x$sensi.estimates$optpar)),colour="red")+
      geom_vline(aes(xintercept=median(x$sensi.estimates$optpar)),colour="blue")+
      xlab(paste("Estimated",x$optpar,"parameter",sep=" ")) +
      ylab("Frequency") +
      theme(axis.title = element_text(size = 12),
            axis.text = element_text(size = 12),
            panel.background = element_rect(fill = "white",
                                            colour = "black"))
    
    if (graphs=="all"){
      if(x$optpar=="none"){
        suppressMessages(return(multiplot(sigsq_fig,z0_fig,aicc_fig, cols=2))) #When transformation = "none" was used, don't print the visualisation for the optimisation parameter
      } else
        suppressMessages(return(multiplot(sigsq_fig,z0_fig,aicc_fig,optpar_fig, cols=2)))
    }
    if (graphs=="sigsq")
      suppressMessages(return(sigsq_fig))
    if (graphs=="z0")
      suppressMessages(return(z0_fig))
    if (graphs=="aic")
      suppressMessages(return(aicc_fig))
    if (graphs=="optpar")
      suppressMessages(return(optpar_fig))
  }
  
}

