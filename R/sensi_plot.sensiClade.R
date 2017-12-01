#' Graphical diagnostics for class 'sensiClade'
#'
#' Plot results from \code{clade_phylm} and \code{clade_phyglm}
#' @param x output from \code{clade_phylm} or \code{clade_phyglm}
#' @param clade The name of the clade to be evaluated (see details)
#' @param ... further arguments to methods.
#' @importFrom ggplot2 aes theme element_text geom_point element_rect ylab xlab
#' ggtitle element_blank geom_abline scale_shape_manual scale_linetype_manual 
#' guide_legend element_rect
#' guides
#' 
#' @author Gustavo Paterno
#' @seealso \code{\link[sensiPhy]{clade_phylm}} 
#' @details For 'x' from clade_phylm or clade_phyglm:
#' 
#' \strong{Graph 1:} The original scatterplot \eqn{y = a + bx} (with the 
#' full dataset) and a comparison between the regression lines of the full dataset
#' and the rerun without the selected clade (set by \code{clade}). For further
#' details about this method, see \code{\link[sensiPhy]{clade_phylm}}.
#' 
#' Species from the selected clade are represented in red (removed species), black
#' solid line represents the regression with the full model and red dashed line represents
#' the regression of the model without the species from the selected clade.
#' To check the available clades to plot, see \code{x$sensi.estimates$clade} 
#' in the object returned from \code{clade_phylm} or \code{clade_phyglm}. 
#' 
#' \strong{Graph 2:} Distribution of the simulated slopes (Null distribution
#' for a given clade sample size).
#' The red dashed line represents the estimated slope for the reduced model 
#' (without the focal clade) and the black line represents the slope for the 
#' full model.
#'  
#' @importFrom ggplot2 aes_string
#' @importFrom stats model.frame qt plogis 
#' @export

sensi_plot.sensiClade <- function(x, clade = NULL, ...){
    yy <- NULL; estimate <- NULL
    # start:
    full.data <- x$data
    mappx <- x$formula[[3]]
    mappy <- x$formula[[2]]
    vars <- all.vars(x$formula)
    clade.col <- x$clade.col
    
    clades.names <- x$sensi.estimates$clade
    if (is.null(clade) == T){
        clade <- clades.names[1]
        warning("Clade argument was not defined. Plotting results for clade: ",
                clade,"
                Use clade = 'clade name' to plot results for other clades")
    }
    clade.n <- which(clade == clades.names)
    if (length(clade.n) == 0) stop(paste(clade,"is not a valid clade name"))
    
    ### Organizing values:
    result <- x$sensi.estimates
    intercept.0 <-  as.numeric(x$full.model.estimates$coef[1])
    estimate.0     <-  as.numeric(x$full.model.estimates$coef[2])

    inter <- c(x$sensi.estimates$intercept[clade.n ],
               intercept.0)
    slo <-  c(x$sensi.estimates$estimate[clade.n ],
              estimate.0)
    model <- NULL
    estimates <- data.frame(inter,slo, model=c("Without clade", "Full data"))
    
    xf <- model.frame(formula = x$formula, data = full.data)[,2]
    yf <- plogis(estimates[2,1] + estimates[2,2] * xf)
    yw <- plogis(estimates[1,1] + estimates[1,2] * xf)
    plot_data <- data.frame("xf" = c(xf,xf),
                            "yy" = c(yw, yf),
                            model = rep(c("Without clade","Full data"),
                                        each = length(yf)))
                            
    match.y <- which(full.data[, clade.col] == clade)
    match.n <- which(full.data[, clade.col] != clade)
    
    g1 <- ggplot2::ggplot(full.data, aes_string(y = mappy, x = mappx),
                    environment = parent.frame())+
        geom_point(data = full.data[match.n, ], alpha = .7,
                   size = 4)+
        geom_point(data = full.data[match.y, ],alpha = .5,
                   size = 4, aes(shape = "Removed species"), colour = "red")+
        
        scale_shape_manual(name = "", values = c("Removed species" = 16))+
        guides(shape = guide_legend(override.aes = list(linetype = 0)))+
        scale_linetype_manual(name = "Model", values = c("Full data" = "solid",
                                                    "Without clade" = "dashed"))+
        scale_color_manual(name = "Model", values = c("Full data" = "black",
                                                    "Without clade" = "red"))+
        theme(axis.text = element_text(size = 12),
              axis.title = element_text(size = 12),
              legend.text = element_text(size = 12),
              plot.title = element_text(size = 12),
              panel.background = element_rect(fill = "white", colour = "black"),
              legend.position="bottom", legend.box = "horizontal")+
        ggtitle(paste("Clade removed: ", clade, sep = ""))
    
    ### Permuation Test plot:
    nd <- x$null.dist
    ces <- x$sensi.estimates
    
    nes <- nd[nd$clade == clade, ]
    slob <- ces[ces$clade == clade ,]$estimate
    slfu <- x$full.model.estimates$coef[[2]]
    
    ### P.value permutation test:
    p.values <- summary(x)[[1]]
    P <- p.values[p.values$clade.removed == clade, ]$Pval.randomization
  
    g2 <- ggplot2::ggplot(nes ,aes(x=estimate))+
      geom_histogram(fill="yellow",colour="black", size=.2,
                     alpha = .3) +
      geom_vline(xintercept = slob, color="red",linetype=2,size=.7)+
      geom_vline(xintercept = slfu, color="black",linetype=1,size=.7)+
      xlab(paste("Simulated estimates | N.species = ", 
                 ces[ces$clade==clade, ]$N.species, "| N.sim = ", 
                 nrow(nes))) +
      ylab("Frequency")+
      theme(axis.title=element_text(size=12),
            axis.text = element_text(size=12),
            panel.background = element_rect(fill="white",
                                            colour="black"))+
      ggtitle(paste("Randomization test for", clade, " | P = ", 
                    sprintf("%.3f", P)))
     
    ### plot lines: linear or logistic depending on output class
    if(length(class(x)) == 1){
        g.out <- g1 + geom_abline(data = estimates, aes(intercept = inter, slope = slo,
                                      linetype = factor(model),color=factor(model)),
                size=.8)
    }
    if(length(class(x)) == 2){
        g.out <- g1 + geom_line(data = plot_data, aes(x = xf, y = yy, linetype = factor(model),color=factor(model)))
    }
    return(suppressMessages(multiplot(g.out, g2, cols=2)))
}



#' Graphical diagnostics for class 'sensiClade.TraitEvol'
#'
#' \code{plot_clade_physig} Plot results from \code{clade_physig}
#' @param x output from \code{influ_physig}
#' @param clade The name of the clade to be evaluated (see details)
#' @param ... further arguments to methods.
#' @importFrom ggplot2 aes geom_histogram geom_density geom_vline 
#' xlab theme element_text geom_point scale_colour_gradient element_rect ylab xlab
#' ggtitle element_blank
#' @author Gustavo Paterno
#' @seealso \code{\link[ggplot2]{ggplot}}, \code{\link[phytools]{phylosig}}
#' @details For 'x' from clade_physig:
#' 
#' \strong{Graph 1:} Distribution of the simulated signal estimates (Null distribution
#' for a given clade sample size).
#' The red dashed line represents the estimated signal for the reduced data 
#' (without the focal clade) and the black line represents the signal estimate
#'  for the full data.
#' 
#' @export

### Start:
sensi_plot.sensiClade.TraiEvol <- function(x, clade = NULL,graph="all",...) {
  
  if(as.character(x$call[[1]])=="clade_continuous"){ #Check what type of TraitEvolution is evaluated
    if(is.null(graph)) stop("Specify what graph to print (sigsq or optpar)")
    if(any(graph %in% c("q12","q21"))) stop("q12 and q12 are not allowed for clade_discrete, specify 'all', 'sigsq' or 'optpar'")
    else
  # check clade
  clades.names <- x$sensi.estimates$clade
  if (is.null(clade) == TRUE){
    clade <- clades.names[1]
    message("Clade argument was not defined. Plotting results for clade: ",
            clade,"
            Use clade = 'clade name' to plot results for other clades")
  }
  
  times <- x$call$n.sim
  if(is.null(x$call$n.sim)) times <- 1000
  
  #Makde the sigsq graph
  ces       <- x$sensi.estimates 
  N.species <- ces[ces$clade==clade, ]$N.species
  p.values <- summary(x)$sigsq
  P <- p.values[p.values$clade.removed == clade, ]$Pval.randomization   
  
  wcf <- x$sensi.estimates$clade %in% clade
  wcn <- x$null.dist$clade %in% clade
  if(sum(wcf) == 0) stop(paste(clade, "is not a valid clade name", sep = " "))
  
  # Full estimate
  e.0 <- x$full.model.estimates$sigsq 
  # Withou clade estimate
  c.e <- x$sensi.estimates[wcf, ]$sigsq
  nd <- x$null.dist[wcn, ] ### CLADE NULL DIST
  
  ## Estimates dataframe
  vl <- data.frame(model = c("Full data", "Without clade"), 
                   estimate = c(e.0,c.e))
  leg.title <- paste("Estimated sigsq")
  
  ### Graph 1
  p1 <- ggplot2::ggplot(nd, aes(x = nd$sigsq)) + 
    geom_histogram(fill="yellow",colour="black", size=.2, alpha = .3) +
    geom_vline(data = vl, aes(xintercept = estimate,
                              colour = model, linetype = model),
               size = 1.2) + 
    scale_color_manual(leg.title, values = c("black","red")) +
    scale_linetype_manual(leg.title, values = c(1,2)) +
    theme(axis.title=element_text(size=12),
          axis.text = element_text(size=12),
          panel.background = element_rect(fill="white",
                                          colour="black"),
          legend.key = element_rect(fill = "transparent"),
          legend.background = element_rect(fill = "transparent"),
          legend.position = c(.9,.9))+
    ylab("Frequency")+
    xlab(paste("Simulated sigsq | N.species = ", 
               N.species, "| N.sim = ", times)) +
    ggtitle(paste("Randomization test for", clade, " | P = ", 
                  sprintf("%.3f", P)))
  #Make the optpar graph
    ces       <- x$sensi.estimates 
    N.species <- ces[ces$clade==clade, ]$N.species
    p.values <- summary(x)$optpar
    P <- p.values[p.values$clade.removed == clade, ]$Pval.randomization   
    
    wcf <- x$sensi.estimates$clade %in% clade
    wcn <- x$null.dist$clade %in% clade
    if(sum(wcf) == 0) stop(paste(clade, "is not a valid clade name", sep = " "))
    
    # Full estimate
    e.0 <- x$full.model.estimates$optpar 
    # Withou clade estimate
    c.e <- x$sensi.estimates[wcf, ]$optpar
    nd <- x$null.dist[wcn, ] ### CLADE NULL DIST
    
    ## Estimates dataframe
    vl <- data.frame(model = c("Full data", "Without clade"), 
                     estimate = c(e.0,c.e))
    leg.title <- paste("Estimated optpar")
    
    ### Graph 1
    p2 <- ggplot2::ggplot(nd, aes(x = nd$optpar)) + 
      geom_histogram(fill="yellow",colour="black", size=.2, alpha = .3) +
      geom_vline(data = vl, aes(xintercept = estimate,
                                colour = model, linetype = model),
                 size = 1.2) + 
      scale_color_manual(leg.title, values = c("black","red")) +
      scale_linetype_manual(leg.title, values = c(1,2)) +
      theme(axis.title=element_text(size=12),
            axis.text = element_text(size=12),
            panel.background = element_rect(fill="white",
                                            colour="black"),
            legend.key = element_rect(fill = "transparent"),
            legend.background = element_rect(fill = "transparent"),
            legend.position = c(.9,.9))+
      ylab("Frequency")+
      xlab(paste("Simulated optpar | N.species = ", 
                 N.species, "| N.sim = ", times)) +
      ggtitle(paste("Randomization test for", clade, " | P = ", 
                    sprintf("%.3f", P)))
    
    if (graph=="all"){
    suppressMessages(return(multiplot(p1,p2, cols=2))) 
      } 
    if (graph=="sigsq")
      suppressMessages(return(p1))
    if (graph=="optpar")
      suppressMessages(return(p2))
  }

  if(as.character(x$call[[1]])=="clade_discrete"){ #Check what type of TraitEvolution is evaluated
    if(is.null(graph)) stop("Specify what graph to print (q12 or q21)")
    if(any(graph %in% c("sigsq","optpar"))) stop("sigsq and optpar are not allowed for clade_discrete, specify 'all', 'q12' or 'q21'")
    else
      # check clade
      clades.names <- x$sensi.estimates$clade
    if (is.null(clade) == TRUE){
      clade <- clades.names[1]
      message("Clade argument was not defined. Plotting results for clade: ",
              clade,"
              Use clade = 'clade name' to plot results for other clades")
    }
    
    times <- x$call$n.sim
    if(is.null(x$call$n.sim)) times <- 1000
    
    #Makde the q12 graph
    ces       <- x$sensi.estimates 
    N.species <- ces[ces$clade==clade, ]$N.species
    p.values <- summary(x)$q12
    P <- p.values[p.values$clade.removed == clade, ]$Pval.randomization   
    
    wcf <- x$sensi.estimates$clade %in% clade
    wcn <- x$null.dist$clade %in% clade
    if(sum(wcf) == 0) stop(paste(clade, "is not a valid clade name", sep = " "))
    
    # Full estimate
    e.0 <- x$full.model.estimates$q12 
    # Withou clade estimate
    c.e <- x$sensi.estimates[wcf, ]$q12
    nd <- x$null.dist[wcn, ] ### CLADE NULL DIST
    
    ## Estimates dataframe
    vl <- data.frame(model = c("Full data", "Without clade"), 
                     estimate = c(e.0,c.e))
    leg.title <- paste("Estimated q12")
    
    ### Graph 1
    p1 <- ggplot2::ggplot(nd, aes(x = nd$q12)) + 
      geom_histogram(fill="yellow",colour="black", size=.2, alpha = .3) +
      geom_vline(data = vl, aes(xintercept = estimate,
                                colour = model, linetype = model),
                 size = 1.2) + 
      scale_color_manual(leg.title, values = c("black","red")) +
      scale_linetype_manual(leg.title, values = c(1,2)) +
      theme(axis.title=element_text(size=12),
            axis.text = element_text(size=12),
            panel.background = element_rect(fill="white",
                                            colour="black"),
            legend.key = element_rect(fill = "transparent"),
            legend.background = element_rect(fill = "transparent"),
            legend.position = c(.9,.9))+
      ylab("Frequency")+
      xlab(paste("Simulated q12 | N.species = ", 
                 N.species, "| N.sim = ", times)) +
      ggtitle(paste("Randomization test for", clade, " | P = ", 
                    sprintf("%.3f", P)))
    #Make the q21 graph
    ces       <- x$sensi.estimates 
    N.species <- ces[ces$clade==clade, ]$N.species
    p.values <- summary(x)$q21
    P <- p.values[p.values$clade.removed == clade, ]$Pval.randomization   
    
    wcf <- x$sensi.estimates$clade %in% clade
    wcn <- x$null.dist$clade %in% clade
    if(sum(wcf) == 0) stop(paste(clade, "is not a valid clade name", sep = " "))
    
    # Full estimate
    e.0 <- x$full.model.estimates$q21 
    # Withou clade estimate
    c.e <- x$sensi.estimates[wcf, ]$q21
    nd <- x$null.dist[wcn, ] ### CLADE NULL DIST
    
    ## Estimates dataframe
    vl <- data.frame(model = c("Full data", "Without clade"), 
                     estimate = c(e.0,c.e))
    leg.title <- paste("Estimated q21")
    
    ### Graph 1
    p2 <- ggplot2::ggplot(nd, aes(x = nd$q21)) + 
      geom_histogram(fill="yellow",colour="black", size=.2, alpha = .3) +
      geom_vline(data = vl, aes(xintercept = estimate,
                                colour = model, linetype = model),
                 size = 1.2) + 
      scale_color_manual(leg.title, values = c("black","red")) +
      scale_linetype_manual(leg.title, values = c(1,2)) +
      theme(axis.title=element_text(size=12),
            axis.text = element_text(size=12),
            panel.background = element_rect(fill="white",
                                            colour="black"),
            legend.key = element_rect(fill = "transparent"),
            legend.background = element_rect(fill = "transparent"),
            legend.position = c(.9,.9))+
      ylab("Frequency")+
      xlab(paste("Simulated q21 | N.species = ", 
                 N.species, "| N.sim = ", times)) +
      ggtitle(paste("Randomization test for", clade, " | P = ", 
                    sprintf("%.3f", P)))
    
    if (graph=="all"){
      suppressMessages(return(multiplot(p1,p2, cols=2))) 
    } 
    if (graph=="q12")
      suppressMessages(return(p1))
    if (graph=="q21")
      suppressMessages(return(p2))
    }
  
  }




