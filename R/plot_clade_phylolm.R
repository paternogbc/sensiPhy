plot_clade_phylolm <- function(x,clade){
    
    clades.names <- x$clade.model.estimates$clade
    clade.n <- which(clade == clades.names)
    if (length(clade.n) == 0) stop("Clade name is not correct")
    
    ### Organizing values:
    result <- x$clade.model.estimates
    vars   <- all.vars(x$formula)
    vars2  <- gsub("list","",attr(terms(x$formula),"variables"))[-1]
    intercept.0 <-  as.numeric(x$full.model.estimates$coef[1])
    slope.0     <-  as.numeric(x$full.model.estimates$coef[2])
    cutoff      <-  x$cutoff
    ### Removing species with error:
    if (isTRUE(class(x$errors) != "character" )){
        x$data <- x$data[-as.numeric(x$errors),]
    }
    
    inter <- c(x$clade.model.estimates$intercept[clade.n ],
               x$full.model.estimates$coef[1,1])
    slo <-  c(x$clade.model.estimates$slope[clade.n ],
              x$full.model.estimates$coef[2,1])
    estimates <- data.frame(inter,slo, model=c("Without clade","Full model"))
    
    result.tab <- data.frame(x$data[vars])
    ggplot2::ggplot(x$data,aes(x=eval(parse(text=vars2[2])),
                               y=eval(parse(text=vars2[1]))),
                    environment = environment())+
        #geom_point(size=3,
        #           aes(colour="Full Model"),
        #           alpha=.5)+
        geom_point(data=x$data[x$data$fam!=clade,],alpha=.7,
                   size=4)+
        geom_point(data=x$data[x$data$fam==clade,],alpha=.5,
                   size=4,aes(shape="Removed species"),colour="red")+
        geom_abline(data=estimates,aes(intercept=inter, slope=slo,linetype=factor(model)),
                    size=.8, show_guide = T)+
        scale_shape_manual(name="",values=c("Removed species"=16))+
        guides(shape = guide_legend(override.aes = list(linetype = 0)))+
        scale_linetype_manual(name="model",values=c("Full model"="solid",
                                                    "Without clade"="dashed"))+
        theme(axis.text = element_text(size = 18),
              axis.title = element_text(size = 18),
              legend.text = element_text(size = 16),
              plot.title = element_text(size = 20))+
        ylab(vars2[1])+
        xlab(vars2[2])+
        ggtitle(paste("Clade removed: ", clade,sep=""))
}

