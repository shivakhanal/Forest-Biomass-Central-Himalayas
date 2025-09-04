### Load necessary libraries ###
library(raster)
library(sp)
library(corrplot)

### Read predictor rasters ###
raster_files <- list.files(full.names = TRUE, pattern = ".tif$")
stacked_rasters <- stack(lapply(raster_files, raster))

### Read species presence/absence data ###
# biodat0 <- read.csv("presabs_data_poly_with_bio4.csv")
biodat0 <- read.csv("presabs_data_raw.csv")

### Extract and examine raster data at species data points ###
species_points <- SpatialPoints(biodat0[, c("lon", "lat")], 
                                proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
transformed_points <- spTransform(species_points, projection(stacked_rasters))
extracted_data <- extract(stacked_rasters, transformed_points)

### Handle missing values ###
extracted_df <- as.data.frame(extracted_data)
extracted_df[is.na(extracted_df$swe), "swe"] <- 0
extracted_df[is.na(extracted_df$bio15), "bio15"] <- mean(extracted_df$bio15, na.rm = TRUE)
extracted_df[is.na(extracted_df$bio4), "bio4"] <- mean(extracted_df$bio4, na.rm = TRUE)
extracted_df[is.na(extracted_df$twi), "twi"] <- mean(extracted_df$twi, na.rm = TRUE)

### Combine species data with extracted raster data ###
raw_data <- cbind(biodat0, extracted_df)

# Uncomment to save raw data
# write.csv(raw_data, "presabs_data_raw.csv", row.names = FALSE)
raw_data <- read.csv("presabs_data_raw.csv")

### Convert selected variables to polynomials ###
rdat <- raw_data
req_poly_vars <- rdat[, names(extracted_df)]

# Generate polynomial terms for covariates
covar_data <- data.frame(poly(req_poly_vars$bio12, 1),
                         poly(req_poly_vars$bio15, 1),
                         poly(req_poly_vars$bio4, 1),
                         poly(req_poly_vars$daily_PISR, 1),
                         poly(req_poly_vars$gdd0, 1),
                         poly(req_poly_vars$pet, 1),
                         poly(req_poly_vars$slope_fill, 1),
                         poly(req_poly_vars$swe, 1),
                         poly(req_poly_vars$twi, 1))
names(covar_data) <- names(req_poly_vars)

### Combine polynomial data with categorical variables ###
fdat <- cbind(raw_data, covar_data)
write.csv(fdat, "presabs_data_poly.csv")

### Prediction on polynomial predictors ###
stacked_rasters <- stack(lapply(list.files(pattern = ".tif$", full.names = TRUE), raster))
covars <- as.data.frame(stacked_rasters)
stksub <- stacked_rasters[names(req_poly_vars)]

model_covars <- raw_data
stksub <- stacked_rasters

covarspoly <- data.frame(
  predict(poly(model_covars$bio12, 1), newdata = as.data.frame(stksub[[1]])$bio12),
  predict(poly(model_covars$bio15, 1), newdata = as.data.frame(stksub[[2]])$bio15),
  predict(poly(model_covars$bio4, 1), newdata = as.data.frame(stksub[[3]])$bio4),
  predict(poly(model_covars$daily_PISR, 1), newdata = as.data.frame(stksub[[4]])$daily_PISR),
  predict(poly(model_covars$gdd0, 1), newdata = as.data.frame(stksub[[5]])$gdd0),
  predict(poly(model_covars$pet, 1), newdata = as.data.frame(stksub[[6]])$pet),
  predict(poly(model_covars$slope_fill, 1), newdata = as.data.frame(stksub[[7]])$slope_fill),
  predict(poly(model_covars$swe, 1), newdata = as.data.frame(stksub[[8]])$swe),
  predict(poly(model_covars$twi, 1), newdata = as.data.frame(stksub[[9]])$twi)
)
names(covarspoly) <- names(stksub)[1:9]

# Beware: This will create a CSV larger than 10GB
write.csv(covars, "predictors_250m_poly.csv")

### Correlation structure of predictors ###
rc <- cor(covar_data)
dir.create("outputs")
pdf("outputs/correlation_matrix_preds.pdf", onefile = TRUE)
corrplot(rc, type = "upper", order = "hclust", tl.col = "black", tl.srt = 45)
dev.off()

### Relationship with total AGB ###
agb_data <- read.csv("data/All_trees_calculated_csv.csv")
names(agb_data)[11] <- "database_value"
biomass <- aggregate(total_mass_air_dry ~ plot_id, data = agb_data, FUN = "sum")
names(biomass) <- c("SiteNo", "agb")
biomass$agb <- biomass$agb / 1000

# Combine dataset with biomass data
final_data <- merge(fdat, biomass, by = "SiteNo")

