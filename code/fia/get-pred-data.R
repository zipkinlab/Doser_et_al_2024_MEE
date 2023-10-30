# get-pred-data.R: this script extracts the grid for prediction across
#                  the US. This script will not run as the raw FIA coordinates
#                  are not available on GitHub.
# Author: Jeffrey W. Doser
rm(list = ls())
library(tidyverse)
library(sf)
library(sp)
library(stars)
library(AOI)
library(climateR)
library(raster)
library(rasterVis)
library(elevatr)
library(cowplot)

# Get area of prediction --------------------------------------------------
my.proj <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs"
usa <- st_as_sf(maps::map("state", fill = TRUE, plot = FALSE))
usa <- usa %>%
  st_transform(crs = my.proj)
# The dx x dy indicates the resolution in terms of km
grid.pred <- st_as_stars(st_bbox(usa), dx = 5, dy = 5)
# Convert to data frame
coords.pred <- as.data.frame(grid.pred, center = TRUE)
# Convert coordinates to an sf object
coords.pred.sf <- st_as_sf(coords.pred, 
			   coords = c('x', 'y'), 
			   crs = my.proj)

# Intersect with region of interest
coords.pred.sf <- st_intersection(coords.pred.sf, st_make_valid(usa))
coords.0 <- as.data.frame(st_coordinates(coords.pred.sf))
coords.lat.long <- coords.pred.sf %>%
  st_transform(crs = '+proj=longlat +datum=WGS84')
# Climate normals from TerraClimate -
# Time period to extract
period <- "19812010"
# Variables to extract
cvars <- c("tmax")
nc.vars <- length(cvars)
J <- nrow(coords.0)
X.0 <- matrix(NA, J, nc.vars + 1)
# Extract the data one variable at a time
for (j in 1:nc.vars){
  cdat <- getTerraClimNormals(coords.lat.long, cvars[j], period, 1:12)[[1]]
  plt_dat <- extract(cdat, coords.lat.long)
  if (cvars[j] == "tmax"){
      val <- apply(plt_dat[,2:ncol(plt_dat)], 1, max)
  }else if (cvars[j] == "tmin"){
      val <- apply(plt_dat[,2:ncol(plt_dat)], 1, min)
  }else if (cvars[j] == "ppt"){
      val <- apply(plt_dat[,2:ncol(plt_dat)], 1, sum)
  }else if (cvars[j] %in% c("pet","aet","def")){
      val <- apply(plt_dat[,5:9], 1, sum)
  }else {
      val <- apply(plt_dat[,5:9], 1, mean)
  }
  X.0[, j] <- val
}
# Get elevation data from AWS Terrain Tiles -
elev <- get_elev_point(coords.lat.long, src = "aws")
tmax <- X.0[, 1]

# Get TCC across the prediction grid --------------------------------------
tcc.us <- read_stars("~/Dropbox/data/tcc-nlcd/nlcd_tcc_conus_2021_v2021-4.tif")
coords.tcc <- coords.pred.sf %>%
  st_transform(crs = st_crs(tcc.us))
tcc.fia.sf <- st_extract(tcc.us, at = coords.tcc)
tcc.cov <- tcc.fia.sf$`nlcd_tcc_conus_2021_v2021-4.tif`

# Get Ecoregion for each point in the grid --------------------------------
ecr <- st_read(dsn = "~/Dropbox/data/usgs-ecoregions/", layer = "us_eco_l3")

# Get the state that each
ecrs.albers <- ecr %>%
  st_transform(crs = my.proj)

# Join points toBCRs
tmp <- st_join(coords.pred.sf, ecrs.albers, join = st_intersects)
ecoregionL3 <- as.numeric(tmp$US_L3CODE)

# Put all covariates together ---------------------------------------------
pred.covs <- data.frame(elev = elev, 
                        tmax = tmax, 
                        ecoregion = ecoregionL3, 
                        tcc = tcc.cov)
bad.indx <- which(apply(pred.covs, 1, function(a) sum(is.na(a))) > 0)
coords.0 <- coords.0[-bad.indx, ]
pred.covs <- pred.covs[-bad.indx, ]

# Save data ---------------------------------------------------------------
save(pred.covs, coords.0, file = 'data/fia-pred-data.rda')
