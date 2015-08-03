#' Graphical diagnostics for \code{influ_phylolm}
#'
#' \code{plot_influ_phylolm} Plot results from \code{influ_phylolm},
#' @param x output from \code{influ_phylolm}
#' @param graphs choose which graph should be printed on the output ("all", 1,2,3 or 4)
#' @param param choose which parameter ("intercept" or "slope" should be printed)
#' @importFrom ggplot2 aes geom_histogram geom_density geom_vline 
#' xlab theme element_text geom_point scale_colour_gradient element_rect ylab xlab
#' ggtitle element_blank
#' @importFrom grid unit 

### Start:
plot_influ_phylolm <- function(x,graphs="all",param="slope"){

# nulling variables:------------------------------------------------------------
slope <- ..density.. <- intercept <- sDFslope <- slope.perc <- NULL
intercept.perc <- sDFintercept <- NULL
        ########### samp_pgls graphs ##################
        if (x[[1]] != "influ_phylolm" & x[[1]] != "influ_phyloglm")
                stop("x must be an output from influ_phylolm or influ_phyloglm!")
        else

        ### Organizing values:
        result <- x$influ.model.estimates
        vars   <- all.vars(x$formula)
        vars2  <- gsub("list","",attr(terms(x$formula),"variables"))[-1]
        intercept.0 <-  as.numeric(x$full.model.estimates$coef[1])
        slope.0     <-  as.numeric(x$full.model.estimates$coef[2])
        cutoff      <-  x$cutoff

        ### Removing species with error:
        if (isTRUE(class(x$errors) != "character" )){
                x$data <- x$data[-as.numeric(x$errors),]
        }
        result.tab <- data.frame(x$influ.model.estimates,x$data[vars])

        ### Plots:
        # Distribution of estimated slopes:
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
                ggtitle("Standardized Difference in Intercept")+
                theme(axis.text = element_text(size=14,colour="black"),
                      axis.title = element_text(size=16))

        # Influential points for slope estimate
        s3 <- ggplot2::ggplot(result,aes(x=sDFslope),
                              environment = environment())+
                geom_histogram(fill="red",color="black",binwidth=.5) +
                xlab("Standardized Difference in Slope")+
                geom_histogram(data=subset(result,sDFslope<cutoff&sDFslope>-cutoff),
                               colour="black", fill="white",binwidth=.5)+
                theme(axis.text = element_text(size=14),
                      axis.title = element_text(size=16))+
                geom_vline(xintercept = -cutoff,color="red",linetype=2,size=.7)+
                geom_vline(xintercept = cutoff,color="red",linetype=2,size=.7)

        # Influential points for intercept estimate
        i3 <- ggplot2::ggplot(result,aes(x=sDFslope),
                              environment = environment())+
                geom_histogram(fill="red",color="black",binwidth=.5) +
                xlab("Standardized Difference in Intercept")+
                geom_histogram(data=subset(result,sDFslope<cutoff&sDFslope>-cutoff),
                               colour="black", fill="white",binwidth=.5)+
                theme(axis.text = element_text(size=14),
                      axis.title = element_text(size=16))+
                geom_vline(xintercept = -cutoff,color="red",linetype=2,size=.7)+
                geom_vline(xintercept = cutoff,color="red",linetype=2,size=.7)

        # Distribution of slope.perc:

        s4 <- ggplot2::ggplot(result,aes(x=slope.perc,y=..density..),
                              environment = environment())+
                geom_histogram(data=subset(result,sDFslope<cutoff&sDFslope>-cutoff),
                               colour="black", fill="white")+
                xlab("% of change in Slope")+
                theme(axis.text = element_text(size=14),
                      axis.title = element_text(size=16))

        # Distribution of slope.perc:
        i4 <- ggplot2::ggplot(result,aes(x=intercept.perc,y=..density..),
                              environment = environment())+
                geom_histogram(data=subset(result,sDFslope<cutoff&sDFslope>-cutoff),
                               colour="black", fill="white")+
                xlab("% of change in Intercept")+
                theme(axis.text = element_text(size=14),
                      axis.title = element_text(size=16))

        ### Ploting:
        if (param == "slope" & graphs=="all")
            suppressMessages(print(multiplot(s1,s3,s2,s4,cols=2)))
        if (param == "slope" & graphs==1)
            suppressMessages(print(s1))
        if (param == "slope" & graphs==2)
            suppressMessages(print(s2))
        if (param == "slope" & graphs==3)
            suppressMessages(print(s3))
        if (param == "slope" & graphs==4)
            suppressMessages(print(s4))
        if (param == "intercept" & graphs=="all")
            suppressMessages(print(multiplot(i1,i3,i2,i4,cols=2)))
        if (param == "intercept" & graphs==1)
            suppressMessages(print(i1))
        if (param == "intercept" & graphs==2)
            suppressMessages(print(i2))
        if (param == "intercept" & graphs==3)
            suppressMessages(print(i3))
        if (param == "intercept" & graphs==4)
            suppressMessages(print(i4))

        ### Warnings
        if (isTRUE(class(x$errors) != "character" ))
                warnings("Deletion of some species caused error. These species were not ploted")

}

