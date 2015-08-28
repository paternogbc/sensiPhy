### decide which plot should be printed (for plot_samp):
which_plot <- function(param = "slope", graphs = "all",
                       s1 = s1, s2 = s2, opt = opt, s4 = s4,
                       i1 = i1, i2 = i2, i4 = i4, model = model){
    
### ERROR checking:-------------------------------------------------------------
if(param != "slope" & param != "intercept")
    stop("param = ", param, " is not valid. Valid options are: `slope` or",
    "`intecept`")
if(graphs != 1 & graphs != 2 & graphs != 3 & graphs != 4 & graphs != "all")
    stop("graphs = ", graphs, " is not valid. Valid options are: `all`, 1, 2, 3, 4.")
    
    
### if model == "BM" or "trend" | SLOPE: ----------------------------------------
if (model == "BM" | model == "trend" & param == "slope" & graphs == "all")
    return(multiplot(s1, "No `optpar` plot available for model == `BM` or `trend`", 
                     s2, s4, cols=2))
if (model == "BM" | model == "trend" & param == "slope" & graphs == 1)
    return(s1)
if (model == "BM" | model == "trend" & param == "slope" & graphs == 2)
    return(s2)
if (model == "BM" | model == "trend" & param == "slope" & graphs == 3)
    stop("No `optpar` plot available for model == `BM` or `trend`")
if (model == "BM" | model == "trend" & param == "slope" & graphs == 4)
    return(s4)
    
### if model == "BM" or "trend" | INTERCEPT: -----------------------------------
if (model == "BM" | model == "trend" & param == "intercept" & graphs == "all")
    return(multiplot(i1, "No `optpar` plot available for model == `BM`", 
                     i2, i4, cols=2))
if (model == "BM" | model == "trend" & param == "intercept" & graphs == 1)
    return(i1)
if (model == "BM" | model == "trend" & param == "intercept" & graphs == 2)
    return(i2)
if (model == "BM" | model == "trend" & param == "intercept" & graphs == 3)
    stop("No `optpar` plot available for model == `BM`")
if (model == "BM" | model == "trend" & param == "intercept" & graphs == 4)
    return(i4)

### if model != "BM" | SLOPE: ----------------------------------------------
if (param == "slope" & graphs=="all")
    suppressMessages(return(multiplot(s1,opt, s2,s4,cols=2)))
if (param == "slope" & graphs==1)
    suppressMessages(return(s1))
if (param == "slope" & graphs==2)
    suppressMessages(return(s2))
if (param == "slope" & graphs==3)
    suppressMessages(return(opt))
if (param == "slope" & graphs==4)
    suppressMessages(return(s4))
if (param == "intercept" & graphs=="all")
    suppressMessages(return(multiplot(i1,opt,i2,i4,cols=2)))
if (param == "intercept" & graphs==1)
    suppressMessages(return(i1))
if (param == "intercept" & graphs==2)
    suppressMessages(return(i2))
if (param == "intercept" & graphs==3)
    suppressMessages(return(opt))
if (param == "intercept" & graphs==4)
    suppressMessages(return(i4))
}

### Function to plot multiple ggplo2 graphs:------------------------------------
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    plots <- c(list(...), plotlist)
    numPlots = length(plots)
    if (is.null(layout)) {
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
    }
    if (numPlots==1) {
        print(plots[[1]])
    } else {
        grid::grid.newpage()
        grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow(layout), ncol(layout))))
        for (i in 1:numPlots) {
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            print(plots[[i]], vp = grid::viewport(layout.pos.row = matchidx$row,
                                                  layout.pos.col = matchidx$col))
        }
    }
}
