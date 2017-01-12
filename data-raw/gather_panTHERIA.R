library(sensiPhy); library(ggplot2); library(dplyr)

pantheria.data <- read.csv("data-raw/pantheria.csv")
pantheria.phy  <- read.nexus("data-raw/pantheria.tre")

### Include species in data rownames
rownames(pantheria.data) <- pantheria.data$MSW93_Binomial
names(pantheria.data)

### Filtering variables:
pantheria.d <- dplyr::select(pantheria.data, 
                            order         = MSW93_Order,
                            family        = MSW93_Family,
                            adultMass     = X5.1_AdultBodyMass_g,
                            sexMaturity   = X23.1_SexualMaturityAge_d,
                            homeRange         = X22.1_HomeRange_km2)
pantheria.d <- subset(pantheria.d, pantheria.d$order == "Primates")
pantheria.d <- pantheria.d[, -1]

### Match data and phylogeny:
tree.drop    <- drop.tip(pantheria.phy[[1]],rownames(pantheria.d))              
primates.phy <-lapply(pantheria.phy, drop.tip, tip = tree.drop$tip.label)
class(primates.phy)<-"multiPhylo"

### Create package datset:
primates <- list(data = primates.data,
              phy  = primates.phy)

### Match data and Phy and remove missing data:
primates <- sensiPhy::match_dataphy(adultMass ~ sexMaturity + homeRange,
                                    data = primates.data, phy = primates.phy)
primates.data <- primates$data
primates.phy  <- primates$phy

### save package dataset:
devtools::use_data(primates, primates.data, primates.phy, overwrite = TRUE)
