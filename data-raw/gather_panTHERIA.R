library(sensiPhy); library(ggplot2); library(dplyr)

pantheria.data <- read.csv("data-raw/pantheria.csv")
pantheria.phy  <- read.nexus("data-raw/pantheria.tre")

### Include species in data rownames
rownames(pantheria.data) <- pantheria.data$MSW93_Binomial
names(pantheria.data)
rownames(pantheria.data)
### Filtering variables:
pantheria.data <- dplyr::select(pantheria.data, 
                            order         = MSW93_Order,
                            family        = MSW93_Family,
                            adultMass     = X5.1_AdultBodyMass_g,
                            sexMaturity   = X23.1_SexualMaturityAge_d)
pantheria.data <- subset(pantheria.data, pantheria.data$order == "Primates")
pantheria.data <- pantheria.data[, -1]
### Match data and phylogeny:
primates <- match_dataphy(sexMaturity ~ adultMass, pantheria.data, pantheria.phy)

### Create package datset:
primates <- list(data = primates$data,
              phy  = primates$phy)
devtools::use_data(primates, overwrite = TRUE)
