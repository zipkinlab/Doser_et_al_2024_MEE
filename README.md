# spAbundance: An R package for single-species and multi-species spatially-explicit abundance models 

### Methods in Ecology and Evolution 

### [Jeffrey W. Doser](https://www.jeffdoser.com/), [Andrew O. Finley](https://www.finley-lab.com/), [Marc K&eacute;ry](https://www.mbr-pwrc.usgs.gov/pubanalysis/keryroylebook/), [Elise F. Zipkin](https://zipkinlab.org/) 

### Package [Website](https://www.jeffdoser.com/files/spabundance-web/) and [Repository](https://github.com/doserjef/spAbundance/)

### Please contact the first author for questions about the code or data used in the empirical case studies: Jeffrey W. Doser (doserjef@msu.edu)

---------------------------------

## Abstract

1. Numerous modeling techniques exist to estimate abundance of plant and animal populations. These methods seek to estimate abundance while accounting for multiple complexities found in ecological data, such as observational biases, spatial autocorrelation, and species correlations. There is, however, a lack of user-friendly and computationally efficient software to implement the various models, particularly for large data sets. 
2. We developed the `spAbundance` `R` package for fitting spatially-explicit Bayesian single-species and multi-species hierarchical distance sampling models, N-mixture models, and generalized linear mixed models. The models within the package can account for spatial autocorrelation using Nearest Neighbor Gaussian Processes and accommodate species correlations in multi-species models using a latent factor approach, which enables model fitting for data sets with large numbers of sites and/or species.
3. We provide three vignettes and three case studies that highlight `spAbundance` functionality. We used spatially-explicit multi-species distance sampling models to estimate density of 16 bird species in Florida, USA, an N-mixture model to estimate Black-throated Blue Warbler (*Setophaga caerulescens*) abundance in New Hampshire, USA, and a spatial linear mixed model to estimate forest aboveground biomass across the continental USA. 
4. `spAbundance` provides a user-friendly, formula-based interface to fit a variety of univariate and multivariate spatially-explicit abundance models. The package serves as a useful tool for ecologists and conservation practitioners to generate improved inference and predictions on the spatial drivers of populations and communities.


## Repository Directory

All code and resulting model objects were created and saved using spAbundance v0.1.1.

### [code/neon](./code/neon)

Contains all code to fit models and summarize results for the bird hierarchical distance sampling case study using data from NEON in the Disney Wilderness Preserve. Note the `extract-samples.R` script will not run successfully, as the full model objects are too large to include on GitHub. However, these objects can be generated first by running the `main-*.R` files if desired, and the complete model results can still be summarized using the `summary.R` figure using smaller objects created in the `extract-samples.R` script.

+ `extract-samples.R`: extracts MCMC samples and Rhat values from the full multi-species model objects, which are too large to put on GitHub.
+ `main-lfMsDS.R`: script to fit a multi-species distance sampling model with species correlations to estimate density of the 16 bird species.
+ `main-msDS.R`: script to fit a multi-species distance sampling model with species correlations to estimate density of the 16 bird species.
+ `main-sfMsDS.R`: script to fit a spatial multi-species distance sampling model with species correlations to estimate density of the 16 bird species.
+ `predict-sfMsDS.R`: script to predict density across the preserve using results from the spatially-explicit hierarchical distance sampling model.
+ `summary.R`: script to summarize the results from the Disney Wilderness Preserve case study.

### [code/hbef](./code/hbef)

Contains all code to fit models and summarize results for the Black-throated Blue Warbler case study using data from the Hubbard Brook Experimental Forest.

+ `main-NMix.R`: script to fit a non-spatial Poisson N-mixture model.
+ `main-NMix-NB.R`: script to fit a non-spatial negative binomial N-mixture model.
+ `main-spNMix.R`: script to fit a spatial Poisson N-mixture model.
+ `main-spNMix-NB.R`: script to fit a spatial negative binomial N-mixture model.
+ `predict-NMix.R`: script to predict abundance across HBEF using a Poisson N-mixture model.
+ `summary.R`: script to summarize results from the Black-throated Blue Warbler case study. 

### [code/fia](./code/fia)

Contains all code to fit models and summarize results for the forest aboveground biomass case study using data from the US Forest Service Forest Inventory and Analysis Program. Note that the data included here and used in the manuscript contains the "fuzzed coordinates". These coordinates are not the exact locations of the FIA plots, as FIA adds a small amount of random noise to plot locations to protect ownership privacy and ensure ecological integrity. Note that only the `main-*.R` files and `summary.R` files can be run using the data and files included on GitHub, as the raw FIA data and the full model result files from `spAbundance` are too large to include on GitHub. However, this still allows all models to be run and all figures and results to be extracted using smaller objects that are available on GitHub.

+ `combine-chains.R`: script to combine chains from three chains run in parallel into one model output file.
+ `get-cov-data.R`: script to calculate the covariates for use in the FIA case study. This calculates a variety of potential covariates, which are not all used in the case study.
+ `get-data.R`: script that generates the FIA data into the format necessary for fitting models in `spAbundance`.  
+ `get-pred-data.R`: script to extract the grid for prediction across the US.
+ `get-waic-and-samples.R`: this script extracts the WAIC for each of the four candidate models. 
+ `main-nonspatial.R`: script to fit a linear mixed model for the FIA biomass data with a random ecoregion intercept.
+ `main-ecoregion.R`: script to fit a linear mixed model with a random slope of tree canopy cover by ecoregion and a random intercept of ecoregion.
+ `main-sp.R`: script to fit a spatial linear mixed model.
+ `main-sp-ecoregion.R`: script to fit a spatial linear mixed model with a random slope of tree canopy cover by ecoregion.
+ `pred-main.R`: code to predict forest AGB across the continental US using a spatial linear mixed model with a random slope of tree canopy cover.
+ `prep-raw-data.R`: script to prepare the raw FIA data into a more usable format.
+ `summary.R`: script that summarizes results from the FIA biomass case study.

### [data](./data)

Data objects for fitting models in `spAbundance`. Note the data for the Central Florida bird case study and the HBEF case study are included within the package itself.

+ `fia-data.rda`: the FIA data (the fuzzed coordinates) for the FIA biomass case study.
+ `neon-locations/`: shapefiles for use in generating the figures for the Central Florida bird case study.
+ `fia-pred-data.rda`: the prediction data for the FIA case study.
+ `hbef-spatial/`: shapefiles for use in generating the figures for the Hubbard Brook case study.

### [results](./results)

All result files (or smaller versions of result files that can be placed on GitHub) for the three case studies.

+ `neon-msDS-small-results.rda`: subset of the results from a non-spatial multi-species HDS model for the central Florida case study.
+ `neon-lfMsDS-small-results.rda`: subset of the results from a non-spatial multi-species HDS model with species correlations for the central Florida case study.
+ `neon-sfMsDS-small-results.rda`: subset of the resulst from a spatial multi-species HDS model with species correlations for the central Florida case study.
+ `neon-pred-sfMsDS-results.rda`: prediction results from a spatial multi-species HDS model with species correlations for the central Florida case study.
+ `neon-waic-results.rda`: WAIC results from the central Florida case study.
+ `hbef-NMix-poisson-fit.rda`: results from a non-spatial Poisson N-mixture model for the HBEF case study.
+ `hbef-NMix-NB-fit.rda`: results from a non-spatial negative binomial N-mixture model for the HBEF case study.
+ `hbef-spNMix-poisson.rda`: results from a spatial Poisson N-mixture model for the HBEF case study.
+ `hbef-spNMix-nb.rda`: results from a spatial NB N-mixture model for the HBEF case study.
+ `hbef-NMix-poisson-pred-results.rda`: prediction results from a nonspatial Poisson N-mixture model for the HBEF case study.
+ `fia-top-model-results.rda`: a subset of the results that can fit on GitHub for the top performing model in the FIA case study (the spatial linear mixed model with a random ecoregion slope of tree canopy cover).
+ `fia-pred-results.rda`: prediction results from the top performing model for the FIA case study.
+ `fia-waic-results.rda`: WAIC results from all four candidate models for the FIA case study. 

