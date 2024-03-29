# 1b-get-cov-data.R: script to calculate the covariates for use in the FIA 
#                    case study. This calculates a variety of potential 
#                    covariates. 
# Author: Jeffrey W. Doser (with some adapted code from Malcolm Itter from a different project)

library(AOI)
library(climateR)
library(sf)
library(raster)
library(rasterVis)
library(elevatr)
library(ggplot2)
library(cowplot)

# NOTE: this script will not run as the raw FIA data are too large for GitHub.
#       If raw data files are desired, please contact the first author
#       (doserjef@msu.edu)

# Load plot coordinates ---------------------------------------------------
load("data/fuzzed_forest_plot_coords.RData")
plotCoords <- as.data.frame(fuzzed_forest_plot_coords)
plotCoords <- st_as_sf(plotCoords, coords = c("LON","LAT"),
                      crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs")
plotCoords <- st_transform(plotCoords, crs = "+proj=longlat +datum=WGS84")
n <- nrow(plotCoords)
site_idx <- plotCoords$pltID

# Download TerraClimate data ----------------------------------------------
period <- "19812010"
cvars <- c("tmax","tmin","ppt","pet","aet","def","vpd")

nc.vars <- length(cvars)

X <- matrix(NA, n, nc.vars + 1)
for (i in 1:nc.vars){
  cdat <- getTerraClimNormals(plotCoords, cvars[i], period, 1:12)[[1]]
  plt_dat <- extract(cdat, plotCoords)
  
  if (cvars[i] == "tmax"){
    val <- apply(plt_dat[,2:ncol(plt_dat)], 1, max)
  } else if (cvars[i] == "tmin"){
    val <- apply(plt_dat[,2:ncol(plt_dat)], 1, min)
  } else if (cvars[i] == "ppt"){
    val <- apply(plt_dat[,2:ncol(plt_dat)], 1, sum)
  } else if (cvars[i] %in% c("pet","aet","def")){
    val <- apply(plt_dat[,5:9], 1, sum)
  } else {
    val <- apply(plt_dat[,5:9], 1, mean)
  }
  X[,i] <- val
}

# Get elevation data ------------------------------------------------------
elev <- get_elev_point(plotCoords, src = "aws")
all(elev$pltID == site_idx) # quick check that data is in correct order
X[,nc.vars + 1] <- elev$elevation

colnames(X) <- c(cvars, "elev")
X_info <- data.frame(var = c(cvars,"elev"),
                     units = c(rep("degC",2), rep("mm",4), "kPa", "m"),
                     period = c(rep("Jan-Dec",3), rep("May-Sep",4), NA),
                     statistic = c("max","min", rep("sum",4),"mean",NA))

save(X, X_info, site_idx, file = "data/X_vars.rda")
