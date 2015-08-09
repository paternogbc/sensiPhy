### Function to match data and phylogeny:
match_dataphy <- function(formula,data,phy){

# original data set:
data.0 <- data
# Croping data frame by formula variables:
mf <- stats::model.frame(formula = formula, data = data, na.action = stats::na.pass )
vars <- all.vars(formula)
  
#Remove NA's before matching data and tips
if (sum(is.na(mf[,1]))!=0 || sum(is.na(mf[,2]))!=0)
{data <- mf[!is.na(mf[,1]) & !is.na(mf[,2]),]
 warning("NA's in response or predictor, rows with NA's were removed")}

#Match data and phylogeny in comparative.data style
if(inherits(phy, "multiPhylo")){  
  phy1<-phy[[1]]}
else
  phy1<-phy

tiplabl<-phy1$tip.label
taxa.nam<-as.character(rownames(data))

in.both <- intersect(taxa.nam, tiplabl)

if (length(in.both) == 0)
  stop("No tips are common to the dataset and phylogeny")

mismatch<-union(setdiff(tiplabl,taxa.nam),setdiff(taxa.nam,tiplabl))
if (length(mismatch) != 0)   warning("Phylogeny tips do not match the species list,
species were dropped from phylogeny or species list")

#Drop species from tree
if(inherits(phy, "multiPhylo")){ 
    phy<-lapply(phy,ape::drop.tip,tip=mismatch)
    class(phy)<-"multiPhylo"
    tip.order <- match(phy[[1]]$tip.label, rownames(data))
}
if(inherits(phy, "phylo")){ 
    phy<- ape::drop.tip(phy,tip=mismatch)
    class(phy)<-"phylo"
    tip.order <- match(phy$tip.label, rownames(data))
}

if (any(is.na(tip.order)))
    stop("Problem with sorting data frame: mismatch between tip labels and data frame labels")
data <- data[tip.order, , drop = FALSE]
data.out <- data.0[rownames(data),]

return(list(data = data.out,phy = phy))
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
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
        
        for (i in 1:numPlots) {
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            
            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
        }
    }
}
