### Function to match data and phylogeny:
match_dataphy <- function(formula, data, phy){

# original data set:
data.0 <- data
# Croping data frame by formula variables:
mf <- stats::model.frame(formula = formula, data = data.0, na.action = stats::na.exclude)
if (nrow(data.0) > nrow(mf)) warning("NA's in response or predictor,", 
                                     " rows with NA's were removed")

#Match data and phylogeny in comparative.data style
if(inherits(phy, "multiPhylo")){  
  phy1 <- phy[[1]]}
else
  phy1<-phy

tiplabl <- phy1$tip.label
taxa.nam <- as.character(rownames(mf))

in.both <- intersect(taxa.nam, tiplabl)

if (length(in.both) == 0)
  stop("No tips are common to the dataset and phylogeny")

mismatch <- union(setdiff(tiplabl,taxa.nam),setdiff(taxa.nam,tiplabl))
if (length(mismatch) != 0)   warning("Some phylogeny tips do not match species in data,",
                                     "species were dropped from phylogeny",
                                     " or data")

#Drop species from tree
if(inherits(phy, "multiPhylo")){ 
    phy <- lapply(phy, ape::drop.tip,tip = mismatch)
    class(phy)<-"multiPhylo"
    tip.order <- match(phy[[1]]$tip.label, rownames(mf))
}
if(inherits(phy, "phylo")){ 
    phy <- ape::drop.tip(phy,tip = mismatch)
    class(phy)<-"phylo"
    tip.order <- match(phy$tip.label, rownames(mf))
}

if (any(is.na(tip.order)))
    stop("Problem with sorting data frame: mismatch between tip labels and data frame labels")
 data <- data[tip.order, , drop = FALSE]
 data.out <- data.0[rownames(data),]

message(paste("Final dataset with ",nrow(data.out)," species in data and phylogeny"))
return(list(data = data.out, phy = phy))
}

### Function to plot multiple ggplo2 graphs:
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
