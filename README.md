# README

This repository contains the R code supporting the methodological analysis described in our paper "*Spatially explicit prediction of Nepal’s forest biomass stocks, a data-driven bioregionalisation and machine learning approach*" . The framework estimates the potential above-ground biomass (AGB) of forests and models deviations from that potential. The analysis proceeds in two main steps:

1. **Clustering and Prediction of Potential Forest AGB**  
   Locations are clustered into homogeneous groups based on tree species composition and topoclimatic variables. Using these clusters and plot-level estimates of forest AGB from the national forest inventory, we predict potential forest AGB.

2. **Modeling Deviations from Predicted Potential AGB**  
   Finer-scale predictors, including forest structure, disturbance likelihood, and elevation zones, are incorporated to explain deviations from predicted AGB.

## Repository Contents

This repository includes the following R scripts used for data preparation and modeling:

- `01_prepare_raster_covariates_and_polynomials.R`: Prepares input data by reading, transforming, and extracting predictor variables.  
- `02_fit_RCP_mod.R`: Fits the Region of Common Profile (RCP) model for clustering analysis.  
- `03_get_CI_SE_for_RCPs.R`: Obtains confidence intervals (CI) and standard errors (SE) for the fitted RCP models.  
- `04_predict_potential_forest_AGB.R`: Estimates potential forest AGB using repeated sampling of plots within RCPs and weighted averaging.  
- `05_predict_deviation_from_potential.R`: Models deviations from predicted potential AGB using finer-scale predictors (e.g., forest structure and disturbance proxies).  
 

## Contact

For questions or further information, please contact:  

**Shiva Khanal**  
📧 1khanalshiva@gmail.com  

