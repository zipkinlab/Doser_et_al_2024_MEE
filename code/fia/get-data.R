# get-data.R: this script generates the FIA data into the format necessary
#             for fitting models in spAbundance. This script will not run, 
#             as the raw FIA data are not included on GitHub. Note that the 
#             FIA data here are the fuzzed locations.
# Author: Jeffrey W. Doser
rm(list = ls())
library(tidyverse)
library(sf)
library(stars)

# Load data ---------------------------------------------------------------
# MAKE SURE YOU ARE USING THE FUZZED DATA
load("data/fuzzed_forest_plot_spp_data.RData")
# Spatial coordinates
load("data/fuzzed_forest_plot_coords.RData")
# Covariates
load("data/X_vars.rda")

coords <- fuzzed_forest_plot_coords %>%
  select(x = LON, y = LAT) %>%
  as.matrix()
# Convert to matrix
# Number of sites
J <- nrow(coords)
# Number of species
N <- n_distinct(fuzzed_forest_plot_spp_data$SPCD)
# Species ids
sp.id <- unique(fuzzed_forest_plot_spp_data$SPCD)
# Get data in species x site matrix.
y <- fuzzed_forest_plot_spp_data %>%
  pull(SP_PRESENT)
y <- matrix(y, N, J)
y.bio <- fuzzed_forest_plot_spp_data %>%
  pull(BIO_ACRE)
y.bio <- matrix(y.bio, N, J)
# Use total biomass at each site
y <- apply(y.bio, 2, sum)


# Get tree canopy cover from NLCD -----------------------------------------
my.proj <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs"
coords.sf <- st_as_sf(as.data.frame(coords),
		      coords = c('x', 'y'),
		      crs = my.proj)
tcc.us <- read_stars("/mnt/disk1/data/tcc-nlcd/nlcd_tcc_conus_2021_v2021-4.tif")
coords.tcc <- coords.sf %>%
  st_transform(crs = st_crs(tcc.us))
tcc.fia.sf <- st_extract(tcc.us, at = coords.tcc)
tcc.cov <- tcc.fia.sf$`nlcd_tcc_conus_2021_v2021-4.tif`
# Remove sites where there is no TCC
bad.sites <- which(tcc.cov == 254)
y <- y[-bad.sites]
X <- X[-bad.sites, ]
tcc.cov <- tcc.cov[-bad.sites]
coords <- coords[-bad.sites, ]
coords.sf <- coords.sf[-bad.sites, ]

# Get ecoregion associated with each data point ---------------------------
# Get ecoregion as strata for comparison to SVC model ---------------------
ecr <- st_read(dsn = "~/Dropbox/data/usgs-ecoregions/", layer = "us_eco_l3")

# Get the state that each
ecrs.albers <- ecr %>%
  st_transform(crs = my.proj)

# Join points toBCRs
tmp <- st_join(coords.sf, ecrs.albers, join = st_intersects)
ecoregionL3 <- as.numeric(tmp$US_L3CODE)
bad.sites <- which(is.na(ecoregionL3))
y <- y[-bad.sites]
X <- X[-bad.sites, ]
ecoregionL3 <- ecoregionL3[-bad.sites]
tcc.cov <- tcc.cov[-bad.sites]
coords <- coords[-bad.sites, ]
coords.sf <- coords.sf[-bad.sites, ]

# Format for spAbundance --------------------------------------------------
covs <- data.frame(elev = X[, 'elev'], 
		   ppt = X[, 'ppt'], 
                   tmax = X[, 'tmax'], 
                   tmin = X[, 'tmin'], 
                   ecoregion = ecoregionL3, 
                   tcc = tcc.cov)
data.list <- list(y = y,
		  covs = covs, 
		  coords = coords)
save(data.list, file = "data/fia-data.rda")


