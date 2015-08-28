library(sensiPhy); library(ggplot2); library(dplyr)

alien.data <- read.table("data-raw/alienMammals.txt", h=T)
alien.phy <- read.nexus("data-raw/alienMammals.tre")

### Include species in data rownames
rownames(alien.data) <- alien.data$Species_name
names(alien.data)

### Filtering variables:
alien.data <- dplyr::select(alien.data, 
                     family = Taxonomic_family,
                     Mass = Mean_adult_body_mass_g,
                     gesta = Mean_gestation_length_d,
                     range = Mean_home_range_size_km2,
                     SE_mass = SE_adult_body_mass_g,
                     SE_gesta = SE_gestation_length_d,
                     SE_range = SE_home_range_size_km2)
alien <- list(data = alien.data,
              phy = alien.phy)
devtools::use_data(alien, overwrite = T)

