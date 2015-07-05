#' Graphical sensitive analysis for comparative methods
#'
#' \code{sensi_plot} Plot results from \code{samp_pgls},
#' \code{influ_pgls}
#' @aliases sensi_plot
#' @param x output from \code{samp_gls}, \code{influ_gls}
#' @param graphs choose which graphs should be printed on the output ("all", 1,2,3 or 4)
#' @export

### Start:
plot_influ_phylolm <- function(x,graphs="all",param="slope"){

        ########### samp_pgls graphs ##################
        if (x[[1]] != "influ_phylolm")
                stop("x must be an output from influ_phylolm!")
        else
        ### Organizing values:
        result <- x$influ.model.estimates
        vars   <- all.vars(x$formula)
        vars2  <- gsub("list","",attr(terms(x$formula),"variables"))[-1]
        intercept.0 <-  as.numeric(x$full.model.estimates$coef[1])
        slope.0     <-  as.numeric(x$full.model.estimates$coef[2])

        ### Removing species with error:
        if (isTRUE(class(x$errors) != "character" )){
                x$data <- x$data[-as.numeric(x$errors),]
        }
        result.tab <- data.frame(x$influ.model.estimates,x$data[vars])

        ### Plots:
        # Distribution of estimated betas:
        s1 <- ggplot2::ggplot(result,aes(x=slope,y=..density..),
                              environment = environment())+
                geom_histogram(fill="lightyellow", alpha=.9,colour="grey60", size=.2) +
                geom_density(size=.2) +
                geom_vline(xintercept = slope.0,color="red",linetype=2,size=.7)+
                xlab("Estimated slopes")+
                theme(axis.text = element_text(size=14),
                      axis.title = element_text(size=16))
        # Distribution of estimated intercepts:
        i1 <- ggplot2::ggplot(result,aes(x=intercept,y=..density..),
                              environment = environment())+
                geom_histogram(fill="lightyellow", alpha=.9,colour="grey60", size=.2) +
                geom_density(size=.2) +
                geom_vline(xintercept = intercept.0,color="red",linetype=2,size=.7)+
                xlab("Estimated Intercepts")+
                theme(axis.text = element_text(size=14),
                      axis.title = element_text(size=16))

        # Original plot with Standardized DFslope as colour gradient
        s2<-ggplot2::ggplot(result.tab,aes(eval(parse(text=vars2[2])),
                                           eval(parse(text=vars2[1])),
                                           colour=abs(sDFslope)),
                            environment = environment())+
                geom_point(size=3,alpha=.8)+
                scale_colour_gradient( low="black",high="red",name="")+
                theme(legend.key.width = unit(.2,"cm"),
                      panel.background=element_rect(fill="white",colour="black"),
                      legend.text = element_text(size=14),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank())+
                ylab(vars2[1])+
                xlab(vars2[2])+
                ggtitle("Standardized Difference in slope")+
                theme(axis.text = element_text(size=14,colour="black"),
                      axis.title = element_text(size=16))

        # Original plot with Standardized DFintercept as colour gradient
        i2<-ggplot2::ggplot(result.tab,aes(eval(parse(text=vars2[2])),
                                           eval(parse(text=vars2[1])),
                                           colour=abs(sDFintercept)),
                            environment = environment())+
                geom_point(size=3,alpha=.8)+
                scale_colour_gradient( low="black",high="red",name="")+
                theme(legend.key.width = unit(.2,"cm"),
                      panel.background=element_rect(fill="white",colour="black"),
                      legend.text = element_text(size=14),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank())+
                ylab(vars2[1])+
                xlab(vars2[2])+
                ggtitle("Standardized Difference in slope")+
                theme(axis.text = element_text(size=14,colour="black"),
                      axis.title = element_text(size=16))

        # Statistically influential points for slope estimate
        s3 <- ggplot2::ggplot(result,aes(x=sDFslope),
                              environment = environment())+
                geom_histogram(fill="red",color="black",binwidth=.5) +
                xlab("Standardized Difference in Slope")+
                geom_histogram(data=subset(result,sDFslope<2&sDFslope>-2),
                               colour="black", fill="white",binwidth=.5)+
                theme(axis.text = element_text(size=14),
                      axis.title = element_text(size=16))+
                geom_vline(xintercept = -2,color="red",linetype=2,size=.7)+
                geom_vline(xintercept = 2,color="red",linetype=2,size=.7)

        # Statistically influential points for intercept estimate
        i3 <- ggplot2::ggplot(result,aes(x=sDFslope),
                              environment = environment())+
                geom_histogram(fill="red",color="black",binwidth=.5) +
                xlab("Standardized Difference in Slope")+
                geom_histogram(data=subset(result,sDFslope<2&sDFslope>-2),
                               colour="black", fill="white",binwidth=.5)+
                theme(axis.text = element_text(size=14),
                      axis.title = element_text(size=16))+
                geom_vline(xintercept = -2,color="red",linetype=2,size=.7)+
                geom_vline(xintercept = 2,color="red",linetype=2,size=.7)

print(grid.arrange(s1,s2,s3,i1,i2,i3,ncol=3,nrow=2))
}


