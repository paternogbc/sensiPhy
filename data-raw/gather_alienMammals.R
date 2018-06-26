library(sensiPhy); library(ggplot2); library(dplyr)

alien.data <- read.table("data-raw/alienMammals.txt", h=T)
alien.phy <- read.nexus("data-raw/alienMammals.tre")

### Include species in data rownames
rownames(alien.data) <- alien.data$Species_name

### Filtering variables:
alien.data <- dplyr::select(alien.data, 
                     family = Taxonomic_family,
                     adultMass = Mean_adult_body_mass_g,
                     gestaLen = Mean_gestation_length_d,
                     homeRange = Mean_home_range_size_km2,
                     SD_mass = SD_adult_body_mass_g,
                     SD_gesta = SD_gestation_length_d,
                     SD_range = SD_home_range_size_km2)

alien <- list(data = alien.data,
              phy = alien.phy)
devtools::use_data(alien, alien.data, alien.phy, overwrite = TRUE)


### Match data and Phy and remove missing data:
alien      <- sensiPhy::match_dataphy(family ~ 1, data = alien.data, phy = alien.phy)
alien.data <- alien$data
alien.phy  <- alien$phy

### save package dataset:
devtools::use_data(alien, alien.data, alien.phy, overwrite = TRUE)
