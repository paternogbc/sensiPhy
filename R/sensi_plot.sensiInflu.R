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
#' Graph 1: Distribution of estimated slopes (estimates) or intercepts for each 
#' simulation (leave-one-out deletion). Red vertical line represents the original
#' slope or intercept from the full model (with all species). 
#' 
#' Graph 2: Original regression plot (\eqn{trait~predictor}). Standardized 
#' difference in slope or intercept is represented by a continous size scale. 
#' The names of the most influential species (sDF > cutoff) are ploted in the
#' graph. 
#' 
#' Graph 3: Distribution of standardized difference in slope or intercept. Red 
#' colour indicates inbfluential species (with a standardised difference above 
#' the value of \code{cutoff}).
#' 
#' Graph 4: Distribution of the percentage of change in slope or intercept.
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

