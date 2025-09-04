###### 1. Load required libraries and set up initial parameters ######
library(raster)
library(gtools)
library(doParallel)
library(ggplot2)
library(ggthemes)
library(GSIF)
library(gmodels)
library(caret)
library(Metrics)

set.seed(123)
###### 2. Read plot data ######
alltest <- read.csv('data/test_dat.csv')[, -1]
alltrain <- read.csv('data/train_dat.csv')[, -1]

plottrain <- SpatialPointsDataFrame(coords = alltrain[, c("lon", "lat")], 
                                    data = alltrain, 
                                    proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

plottest <- SpatialPointsDataFrame(coords = alltest[, c("lon", "lat")], 
                                   data = alltest, 
                                   proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

###### 3. Read predictors ######
alltrain <- plottrain[, c(1:4)]
alltest <- plottest[, c(1:4)]

###### 4. Extract the predictors ######
lst <- list.files("data/disturbance/", full.names = T, pattern = ".tif$")[-15]
stk <- stack(lapply(lst, raster))

plottrain <- spTransform(plottrain, projection(stk))
vari <- extract(stk, plottrain)
plottrain@data <- cbind(plottrain@data, vari)

plottest <- spTransform(plottest, projection(stk))
vari <- extract(stk, plottest)
plottest@data <- cbind(plottest@data, vari)

###### 5. Read potential AGB ######
path <- paste0("outputs/PotentialAGB/")
subDir <- "Deviation"
dir.create(file.path(path, subDir))
outpath <- paste0(path, "Deviation")

weightmax <- raster(paste0(path, "weighted_AGB_95per_allplots.tif"))
bootmax <- raster(paste0(path, "mean_boot_potential.tif"))

###### 6. Hard cluster and topographic cluster ######
rcp <- raster("outputs/optimumRCP/5/hard_class_5_RCPs.tif")
plottrain$rcp <- extract(rcp, spTransform(plottrain, projection(rcp)))
plottest$rcp <- extract(rcp, spTransform(plottest, projection(rcp)))

clust <- raster("data/Chap4/data/disturbance/cluster_no_gap.tif")
plottrain$clust <- extract(clust, plottrain)
plottest$clust <- extract(clust, plottest)

###### 7. Model the deviation ######
test <- plottest@data
training <- plottrain@data
training <- training[, names(test)]

maxname <- c("quantmax")

cl <- makePSOCKcluster(6)
registerDoParallel(cl)

for (i in 1:length(maxname)) { # In case there are more than one realization of potential AGB
  typename <- maxname[i]
  baha <- maxtype[[i]] # Raster(paste0(path, "/fit.tif"))
  projection(baha) <- "+proj=utm +zone=45 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
  
  plottrain$pred <- extract(baha, spTransform(plottrain, projection(baha)))
  plottest$pred <- extract(baha, spTransform(plottest, projection(baha)))
  
  test <- plottest@data
  training <- plottrain@data
  test$diff <- (test$pred - test$agb) # Absolute difference
  test <- test[test$diff >= 0,] # Excluding if potential is smaller
  
  training$diff <- (training$pred - training$agb)
  
  # Model fitting and prediction
  test[c("fragIndex", "ElevZone")] <- lapply(test[c("fragIndex", "ElevZone")], factor)
  training[c("fragIndex", "ElevZone")] <- lapply(training[c("fragIndex", "ElevZone")], factor)
  
  train_control <- trainControl(method = "cv", number = 10)
  vars <- c("percentTreeCover", "height", "fragIndex", "ElevZone", "costSurf", "diff")
  training <- training[!is.na(training$diff), ]
  
  model <- train(diff ~ ., data = training[, vars], trControl = train_control, method = "rf", importance = TRUE)
  importance <- varImp(model, scale = FALSE)$importance
  varImportance <- data.frame(Variables = row.names(importance), Importance = round(importance[, 'Overall'], 2))
  
  dir.create(file.path(outpath, typename))
  
  # Visualize the relative importance of variables
  vimplot <- ggplot(varImportance, aes(x = reorder(Variables, Importance), y = Importance)) +
    geom_bar(stat = 'identity') + 
    labs(x = 'Variables') +
    coord_flip() + theme(legend.position = "none") +
    theme_bw()
  
  ggsave(filename = paste0(outpath, "/", typename, "/boot_variable_imp_plots_", typename, ".png"), 
         plot = vimplot, width = 4, height = 5)
  
  predictions <- predict(model, test)
  result <- data.frame(Actual = test$diff, Predicted = predictions)
  result$Difference <- (result$Actual - result$Predicted)
  
  test$predDiff <- predictions
  
  # Plotting the results
  rms <- round(Metrics::rmse(result$Actual, result$Predicted), 2)
  rsq <- round(summary(lm(Actual ~ Predicted, data = result))$r.squared, 2)
  anno2 <- paste("RMSE ==", rms, "~R^2==~", rsq)
  
  predictionsFigure <- ggplot(result, aes(Predicted, Actual)) +
    geom_point(colour = "black", alpha = 0.3) +
    labs(x = "Predicted deviation of AGB (t/ha)", y = "Observed deviation of AGB (t/ha)") +
    geom_abline(slope = 1, intercept = 0, col = "blue") +
    theme_bw() +
    theme(aspect.ratio = 1) +
    scale_x_continuous(breaks = seq(0, 1000, by = 100)) +
    scale_y_continuous(breaks = seq(0, 1000, by = 100)) + 
    theme(plot.title = element_text(size = 20, vjust = 0),
          axis.text.x = element_text(angle = 50, size = 10, vjust = 1, hjust = 1),
          axis.text.y = element_text(size = 10)) +
    annotate(geom = "text", x = 350, y = 900, label = anno2, parse = TRUE) +
    expand_limits(x = 0, y = 0)
  
  ggsave(filename = paste0(outpath, "/", typename, "/Deviation_test_plots_", typename, ".png"), 
         plot = predictionsFigure, device = "png", width = 5, height = 5)
  
  # Deviation from maximum AGB colored by RCP
  col_vector <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', 
                  '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#e6195b', '#f038e6') 
  pal20 <- c(RColorBrewer::brewer.pal(n = 10, name = "Paired"))
  col_vector <- pal20
  
  resultdev <- cbind(test, result)
  
  predictionsFigureCol <- ggplot(resultdev, aes(Predicted, Actual, color = rcp)) +
    geom_point() +
    labs(x = expression(paste("Predicted deviation of AGB ", (tha^-1))), 
         y = expression(paste("Observed deviation of AGB ", (tha^-1)))) +
    geom_abline(slope = 1, intercept = 0, col = "blue") +
    theme_bw() +  
    theme(aspect.ratio = 1) +
    scale_x_continuous(breaks = seq(0, 1200, by = 100)) +
    scale_y_continuous(breaks = seq(0, 1200, by = 100)) +
    theme(plot.title = element_text(size = 20, vjust = 0),
          axis.text.x = element_text(angle = 50, size = 10, vjust = 1, hjust = 1),
          axis.text.y = element_text(size = 10)) +
    scale_color_manual(values = col_vector) +
    annotate(geom = "text", x = 350, y = 1000, label = paste0("RMSE = ", rms, " Rsq = ", rsq), color = "red") +
    expand_limits(x = 0, y = 0)
  
  ggsave(filename = paste0(outpath, "/", typename, "/Deviation_test_plots_color_", typename, ".png"), 
         plot = predictionsFigureCol, device = "png", width = 5, height = 5)
  
  plotmax <- baha - bootmax
  writeRaster(plotmax, paste0(outpath, "/Deviation_from_max_", typename, ".tif"), overwrite = TRUE)
  
  plot1 <- rasterToPoints(plotmax)
  df <- data.frame(plot1)
  colnames(df) <- c('lon', 'lat', 'AGB')
  
  world <- map_data('world')
  
  figure <- ggplot(df, aes(lon, lat)) +
    geom_raster(aes(fill = AGB)) + 
    scale_fill_gradient2(low = 'red', mid = 'white', high = 'blue') +
    labs(x = "Longitude", y = "Latitude", title = paste0('Deviation from Maximum of AGB', typename)) +
    theme_classic() + 
    coord_fixed() + 
    theme(legend.position = 'bottom') +
    geom_polygon(data = world, aes(x = long, y = lat, group = group), color = "black", fill = NA, size = 0.1)
  
  ggsave(filename = paste0(outpath, "/", typename, "/spatial_map_deviation_AGB_", typename, ".png"), 
         plot = figure, width = 9, height = 5)
}

stopCluster(cl)
