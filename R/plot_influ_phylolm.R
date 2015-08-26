#' Graphical diagnostics for class 'sensiInflu'
#'
#' \code{plot_influ_phylolm} Plot results from \code{influ_phylolm} and 
#' \code{influ_phyloglm}
#' @param x output from \code{influ_phylolm}
#' @param graphs choose which graph should be printed on the output ("all", 1,2,3 or 4)
#' @param param choose which parameter ("intercept" or "slope" should be printed)
#' @param ... further arguments to methods
#' @importFrom ggplot2 aes geom_histogram geom_density geom_vline 
#' xlab theme element_text geom_point scale_colour_gradient element_rect ylab xlab
#' ggtitle element_blank
#' @author Gustavo Paterno
#' @seealso \code{\link[ggplot2]{ggplot}}
#' @details For 'x' from influ_phylolm or influ_phyloglm:
#' 
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
#' #Generate a random tree
#' set.seed(2468)
#' tree <- rtree(100)
#'
#' #Generate random predictor variable (pred), evolving according to a BM model.
#' pred<- rTraitCont(tree,root.value=0,sigma=1,model="BM")
#'
#' #Generate two continous traits, one evolving highly correlated with the
#' #predictor (trait 1), and one evolving more randomly (trait 2)
#' cont_trait <- pred + rTraitCont(tree,model="BM",sigma=0.1)
#'
#' #Generate two binary traits, one highly correlated to pred (trait 1), the other less.
#' bin_trait <-rbinTrait(n=1,tree,beta=c(-1,0.5),alpha=0.1,
#'                      X=cbind(rep(1,length(tree$tip.label)),pred))
#' dat<-data.frame(pred,cont_trait,bin_trait)
#'
#' #Determine influential species for both regressions.
#' fit1<-influ_phylolm(cont_trait~pred,data = dat,phy = tree)
#' fit2<-influ_phyloglm(bin_trait~pred,data = dat,phy = tree)
#' 
#' # Plot results:
#' sensi_plot(fit1)
#' sensi_plot(fit2)
#' # You can also choose which graph and parameter should be ploted:
#' sensi_plot(fit1, graphs = 1, param = "intercept")
#' sensi_plot(fit2, graphs = "all", param = "slope")
#' }
#' @importFrom grid unit 
#' @export

### Start:
sensi_plot.sensiInflu <- function(x, graphs="all", param="slope"){

# nulling variables:------------------------------------------------------------
slope <- ..density.. <- intercept <- sDFslope <- slope.perc <- NULL
intercept.perc <- sDFintercept <- NULL

        ### Organizing values:
        result <- x$influ.model.estimates
        vars  <- gsub("list","",attr(stats::terms(x$formula),"variables"))[-1]
        intercept.0 <-  as.numeric(x$full.model.estimates$coef[1])
        slope.0     <-  as.numeric(x$full.model.estimates$coef[2])
        cutoff      <-  x$cutoff

        ### Removing species with error:
        if (isTRUE(class(x$errors) != "character" )){
                x$data <- x$data[-as.numeric(x$errors),]
        }
        result.tab <- data.frame(x$influ.model.estimates,x$data[all.vars(x$formula)])

        ### Plots:
        # Distribution of estimated slopes:
        s1 <- ggplot2::ggplot(result,aes(x=slope,y=..density..))+
                geom_histogram(fill="lightyellow", alpha=.9,colour="grey60", size=.2) +
                geom_density(size=.2) +
                geom_vline(xintercept = slope.0,color="red",linetype=2,size=.7)+
                xlab("Estimated slopes")+
                theme(axis.text = element_text(size=14),
                      axis.title = element_text(size=16))
        # Distribution of estimated intercepts:
        i1 <- ggplot2::ggplot(result,aes(x=intercept,y=..density..))+
                geom_histogram(fill="lightyellow", alpha=.9,colour="grey60", size=.2) +
                geom_density(size=.2) +
                geom_vline(xintercept = intercept.0,color="red",linetype=2,size=.7)+
                xlab("Estimated Intercepts")+
                theme(axis.text = element_text(size=14),
                      axis.title = element_text(size=16))

        # Original plot with Standardized DFslope as colour gradient
        s2<-ggplot2::ggplot(result.tab,aes(eval(parse(text=vars[2])),
                                           eval(parse(text=vars[1])),
                                           colour=abs(sDFslope)),
                            environment = environment())+
            geom_point(size=3,alpha=.8)+
            scale_colour_gradient( low="black",high="red",name="")+
            theme(legend.key.width = unit(.2,"cm"),
                  panel.background=element_rect(fill="white",colour="black"),
                  legend.text = element_text(size=14),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank())+
            ylab(vars[1])+
            xlab(vars[2])+
            ggtitle("Standardized Difference in slope")+
            theme(axis.text = element_text(size=14,colour="black"),
                  axis.title = element_text(size=16))

        # Original plot with Standardized DFintercept as colour gradient
        i2<-ggplot2::ggplot(result.tab,aes(eval(parse(text=vars[2])),
                                           eval(parse(text=vars[1])),
                                           colour=abs(sDFintercept)),
                            environment = environment())+
            geom_point(size=3,alpha=.8)+
            scale_colour_gradient( low="black",high="red",name="")+
            theme(legend.key.width = unit(.2,"cm"),
                  panel.background=element_rect(fill="white",colour="black"),
                  legend.text = element_text(size=14),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank())+
            ylab(vars[1])+
            xlab(vars[2])+
            ggtitle("Standardized Difference in Intercept")+
            theme(axis.text = element_text(size=14,colour="black"),
                  axis.title = element_text(size=16))

        # Influential points for slope estimate
        s3 <- ggplot2::ggplot(result,aes(x=sDFslope))+
                geom_histogram(fill="red",color="black",binwidth=.5) +
                xlab("Standardized Difference in Slope")+
                geom_histogram(data=subset(result,sDFslope<cutoff&sDFslope>-cutoff),
                               colour="black", fill="white",binwidth=.5)+
                theme(axis.text = element_text(size=14),
                      axis.title = element_text(size=16))+
                geom_vline(xintercept = -cutoff,color="red",linetype=2,size=.7)+
                geom_vline(xintercept = cutoff,color="red",linetype=2,size=.7)

        # Influential points for intercept estimate
        i3 <- ggplot2::ggplot(result,aes(x=sDFslope))+
                geom_histogram(fill="red",color="black",binwidth=.5) +
                xlab("Standardized Difference in Intercept")+
                geom_histogram(data=subset(result,sDFslope<cutoff&sDFslope>-cutoff),
                               colour="black", fill="white",binwidth=.5)+
                theme(axis.text = element_text(size=14),
                      axis.title = element_text(size=16))+
                geom_vline(xintercept = -cutoff,color="red",linetype=2,size=.7)+
                geom_vline(xintercept = cutoff,color="red",linetype=2,size=.7)

        # Distribution of slope.perc:

        s4 <- ggplot2::ggplot(result,aes(x=slope.perc,y=..density..))+
                geom_histogram(data=subset(result,sDFslope<cutoff&sDFslope>-cutoff),
                               colour="black", fill="white")+
                xlab("% of change in Slope")+
                theme(axis.text = element_text(size=14),
                      axis.title = element_text(size=16))

        # Distribution of slope.perc:
        i4 <- ggplot2::ggplot(result,aes(x=intercept.perc,y=..density..))+
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

