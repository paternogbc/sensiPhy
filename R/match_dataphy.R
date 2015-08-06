match_dataphy <- function(formula,data,phy){

  resp<-formula[[2]]
  pred<-formula[[3]]
  
#Remove NA's before matching data and tips
  if (sum(is.na(resp))!=0 || sum(is.na(pred))!=0)
  {data<-data[!is.na(data$resp),]
  data<-data[!is.na(data$pred),]
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
phy<-lapply(phy,ape::drop.tip,tip=mismatch)
class(phy)<-"multiPhylo"

#Reorder rows according to tip labels
if(inherits(phy, "multiPhylo")){  
  phy1<-phy[[1]]}
else
  phy1<-phy

rownames(data)<-data$taxa.nam
tip.order <- match(phy1$tip.label, rownames(data))
if (any(is.na(tip.order)))
  stop("Problem with sorting data frame: mismatch between tip labels and data frame labels")
data <- data[tip.order, , drop = FALSE]
rownames(data) <- data$taxa.nam
}