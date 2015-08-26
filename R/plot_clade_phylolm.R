#' Graphical diagnostics for class 'sensiClade'
#'
#' Plot results from \code{clade_phylolm} and \code{clade_phyloglm}
#' @param x output from \code{clade_phylolm} or \code{clade_phyloglm}
#' @param clade The name of the clade to be evaluated (see details)
#' @param ... further arguments to methods.
#' @importFrom ggplot2 aes theme element_text geom_point element_rect ylab xlab
#' ggtitle element_blank geom_abline scale_shape_manual scale_linetype_manual 
#' guide_legend element_rect
#' guides
#' 
#' @author Gustavo Paterno
#' @seealso \code{\link[sensiPhy]{clade_phylolm}} 
#' @details This function plots the original scaterplot \eqn{y = a + bx} (with the 
#' full dataset), plus a comparison between the regression lines of the full model
#' and the model without the selected clade (set by \code{clade}). For further
#' details about this method, see \code{\link[sensiPhy]{clade_phylolm}}.
#' 
#' Species from the selected clade are represented in red (removed species),
#' solid line represents the regression with the full model and dashed line represent
#' the regression of the model without the species from the selected clade.
#' To check the available clades to plot, see \code{x$clade.model.estimates$clade} 
#' in the objected returned from \code{clade_phylolm} or \code{clade_phyloglm}. 
#' @export

sensi_plot.sensiClade <- function(x, clade = NULL, ...){
    
    clades.names <- x$clade.model.estimates$clade
    if (is.null(clade) == T){
        clade <- clades.names[1]
        warning("Clade argument was not defined. Ploting results for clade: ",
                clade,"
                Use clade = 'clade name' to plot results for other clades")
    }
    clade.n <- which(clade == clades.names)
    if (length(clade.n) == 0) stop(paste(clade,"is not a valide clade name"))
    
    ### Organizing values:
    result <- x$clade.model.estimates
    vars   <- all.vars(x$formula)
    vars2  <- gsub("list","",attr(stats::terms(x$formula),"variables"))[-1]
    intercept.0 <-  as.numeric(x$full.model.estimates$coef[1])
    slope.0     <-  as.numeric(x$full.model.estimates$coef[2])
    cutoff      <-  x$cutoff
    ### Removing clades with error:
    if (isTRUE(class(x$errors) != "character" )){
        x$data <- x$data[-as.numeric(x$errors),]
    }
    
    inter <- c(x$clade.model.estimates$intercept[clade.n ],
               x$full.model.estimates$coef[1,1])
    slo <-  c(x$clade.model.estimates$slope[clade.n ],
              x$full.model.estimates$coef[2,1])
    model <- NULL
    estimates <- data.frame(inter,slo, model=c("Without clade","Full model"))
    
    result.tab <- data.frame(x$data[vars])
    ggplot2::ggplot(x$data,aes(x=eval(parse(text=vars2[2])),
                               y=eval(parse(text=vars2[1]))),
                    environment = environment())+
        geom_point(data=x$data[x$data$fam!=clade,],alpha=.7,
                   size=4)+
        geom_point(data=x$data[x$data$fam==clade,],alpha=.5,
                   size=4,aes(shape="Removed species"),colour="red")+
        geom_abline(data=estimates,aes(intercept=inter, slope=slo,
                                       linetype=factor(model)),
                    size=.8, show_guide = T)+
        scale_shape_manual(name="",values=c("Removed species"=16))+
        guides(shape = guide_legend(override.aes = list(linetype = 0)))+
        scale_linetype_manual(name="model",values=c("Full model"="solid",
                                                    "Without clade"="dashed"))+
        theme(axis.text = element_text(size = 18),
              axis.title = element_text(size = 18),
              legend.text = element_text(size = 16),
              plot.title = element_text(size = 20),
              panel.background=element_rect(fill="white",colour="black"))+
        ylab(vars2[1])+
        xlab(vars2[2])+
        ggtitle(paste("Clade removed: ", clade,sep=""))
}

