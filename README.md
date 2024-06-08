# README
This repository contains R codes used for the methodological analysis described in our paper. The analysis consists of two main steps:
1. **Clustering and Prediction of Potential Forest AGB**: In this step, we apply a rigorous clustering process to identify homogeneous groups of locations based on tree species and topoclimatic variables. Using these clusters and available plot-level estimates of forest AGB from national forest inventoy, we predict the potential above-ground biomass (AGB) of forests.
  
2. **Modeling Deviations from Predicted Potential AGB**: The second step involves incorporating finer-scale variables, such as forest structure, disturbance likelihood, and elevation zones, to model deviations from the predicted potential AGB.

## Repository Contents
Contains the R scripts used for data preparation and modeling.
  - `01_prepare_raster_covariates_and_polynomials.R`: Script for preparing input data, including reading, transforming, and extracting predictor variables.
  - `02_fit_RCP_mod.R`: Script for fitting Region of Common Profile (RCP)
  - `03_get_CI_SE_for_RCPs.R`:  Script for obtaining the CI and SE of the fitted RCP model. 
  - `04_predict_potential_forest_AGB.R`: Script for Applying repeated sampling of the AGB for plots belonging to different RCP and estimated potential AGB as a weighted estimated. 
  - `05_predict_deviation_from_potential.R`: Script for modeling deviations from the predicted potential AGB using finer-scale variables related to  proxies of  disturbance.

## Contact
For any questions or further information, please contact Shiva Khanal at 1khanalshiva@gmail.com.

