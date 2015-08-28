### Plot decision:
### Ploting:
which_plot <- function(param = "slope", graphs = "all",
                       s1 = s1, s2 = s2, opt = opt, s4 = s4,
                       i1 = i1, i2 = i2, i4 = i4, model = model){
    
### ERROR checking:-------------------------------------------------------------
if(param != "slope" & param != "intercept")
    stop("param = ", param, " is not valid. Valid options are: `slope` or",
    "`intecept`")
if(graphs != 1 & graphs != 2 & graphs != 3 & graphs != 4 & graphs != "all")
    stop("graphs = ", graphs, " is not valid. Valid options are: `all`, 1, 2, 3, 4.")
    
    
### if model == "BM" | SLOPE: --------------------------------------------------
if (model == "BM" & param == "slope" & graphs == "all")
    return(multiplot(s1, "No `optpar` plot available for model == `BM`", 
                     s2, s4, cols=2))
if (model == "BM" & param == "slope" & graphs == 1)
    return(s1)
if (model == "BM" & param == "slope" & graphs == 2)
    return(s2)
if (model == "BM" & param == "slope" & graphs == 3)
    stop("No `optpar` plot available for model == `BM`")
if (model == "BM" & param == "slope" & graphs == 4)
    return(s4)
    
### if model == "BM" | INTERCEPT: ----------------------------------------------
if (model == "BM" & param == "intercept" & graphs == "all")
    return(multiplot(i1, "No `optpar` plot available for model == `BM`", 
                     i2, i4, cols=2))
if (model == "BM" & param == "intercept" & graphs == 1)
    return(i1)
if (model == "BM" & param == "intercept" & graphs == 2)
    return(i2)
if (model == "BM" & param == "intercept" & graphs == 3)
    stop("No `optpar` plot available for model == `BM`")
if (model == "BM" & param == "intercept" & graphs == 4)
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