# 1a-prep-raw-data.R: script to prepare the raw FIA data into a more usable format.
#                     This script will not run as the "spp_plot_list.csv" is too 
#                     large to place on Github. 
# Authors: Andrew O. Finley and Jeffrey W. Doser
rm(list=ls())
library(tidyverse)
library(sf)

# NOTE: this script will not run as the raw FIA data are not included on GitHub.

# Read in the raw FIA data
dat <- read_csv("data/spp_plot_list.csv", col_types = list(SPCD = col_integer()))
                
dat %>% glimpse()
dim(dat)

# Keep most current year inventory and single condition plots.
dat <- dat %>% group_by(pltID) %>% slice_max(YEAR) %>% filter(COND_PROP == 1 & FORESTED_COND == 1)
dim(dat)

# There are duplicate locations after reprojection (see below).
drop_plts <- read_csv("data/drop_plots.csv")
dim(dat)
dat <- dat[!dat$pltID %in% drop_plts$pltID,]
dim(dat)

# Check that there are no duplicated species in a given plot, max_n should be 1.
dat %>% group_by(pltID, SPCD) %>% summarize(n = n()) %>% ungroup() %>% summarize(max_n = max(n))

# Complete spp record for each plot.
spp <- unique(dat$SPCD)
spp <- spp[!is.na(spp)]
n.spp <- length(spp)
n.loc <- length(unique(dat$pltID))

n.spp
n.loc
n.spp*n.loc

dat.comp <- dat %>% ungroup() %>% mutate(SP_PRESENT = 1) %>% complete(pltID, SPCD, fill = list(SP_PRESENT = 0, BIO_ACRE = 0))
dim(dat.comp)

# Remove rows NA SPCD rows, they're there because non-forested plots had NA SPCD
dat.comp <- dat.comp %>% filter(!is.na(SPCD))
dim(dat.comp) # Sanity check, n.spp*n.loc should equal nrow(dat.comp)

# Select the columns we want. 
dat.out <- dat.comp %>% ungroup() %>% select(pltID, SPCD, SP_PRESENT, BIO_ACRE)

dat.out <- dat.out %>% arrange(pltID, SPCD)
dim(dat.out)
# Get coords set up.

coords <- dat %>% group_by(pltID) %>% summarize(LON = first(LON), LAT = first(LAT))
coords.sf <- st_as_sf(coords, 
                      coords = c('LON', 'LAT'), 
                      crs = "+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs")
coords.sf <- coords.sf %>%
  st_transform(crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs")
coords.vals <- st_coordinates(coords.sf)
coords <- data.frame(pltID = coords.sf$pltID, 
                     LON = coords.vals[, 1], 
                     LAT = coords.vals[, 2])

any(duplicated(coords[,c("LON","LAT")]))

# Location of some plots are duplicated due to the reprojection (use this at the top of the script). 
# drop_plts <- data.frame(pltID = coords[which(duplicated(coords[,2:3])),1])
# write_csv(drop_plts, "data/drop_plots.csv")

dat.final <- left_join(dat.out, coords, by="pltID")

any(duplicated(coords$LON))

dat.final <- dat.final[order(dat.final$LON),]

## Reorder just in case the order above messed up the spp order.
dat.final <- dat.final %>% group_by(pltID) %>% arrange(SPCD, .by_group = TRUE)
dat.final

coords.final <- dat.final %>% group_by(pltID) %>% summarize(LON = first(LON), LAT = first(LAT))

dim(coords.final)

fuzzed_forest_plot_coords <- coords.final
save(fuzzed_forest_plot_coords, file="data/fuzzed_forest_plot_coords.RData")

fuzzed_forest_plot_spp_data <- dat.final
save(fuzzed_forest_plot_spp_data, file="data/fuzzed_forest_plot_spp_data.RData")
