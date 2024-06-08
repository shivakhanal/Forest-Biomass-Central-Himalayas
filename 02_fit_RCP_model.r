################################################################################
# This code runs the regimix model in a loop for selected RCPs.                #
# It selects the best model based on the lowest BIC and then bootstraps        #
# the best model.                                                              #
################################################################################

# Set working directory
setwd("/data1/mounts/s01238ss-wsu-u0217-hie/Khanal_S/Chap4")

# Load necessary libraries
library(RCPmod)
library(raster)
library(gtools)
library(data.table)

# Read species presence/absence data
biodat0 <- read.csv("data/presabs_data_poly.csv")[, -1]
su <- colSums(biodat0[, 4:475])
su <- su[su <= 3]
biodat0 <- biodat0[, !names(biodat0) %in% names(su)]

# Define model covariates
model_covariates_vector <- c("bio12", "bio15", "gdd0", "bio4", "pet", "twi", "daily_PISR")
model_covariates_string <- paste(model_covariates_vector, collapse = " + ")

# Create the formula for the RCP model
species_vars <- paste(names(biodat0)[4:242], collapse = ",")
my_form_RCP <- paste0("cbind(", species_vars, ") ~ ", model_covariates_string)
my_form_spp <- NULL

print(my_form_RCP)

# Load predictor space
cat("Loading predictors...\n")
covars <- fread("data/predictors_250m_poly.csv", header = TRUE)
covars <- as.data.frame(covars)[, -1]

# Subset selected predictors
model_covariates_vector <- c("x", "y", model_covariates_vector)
covars <- covars[, model_covariates_vector]
gc()

# Filter complete cases and extract coordinates
covars <- covars[complete.cases(covars),]
locs <- covars[, 1:2]

# Create indices for data splitting
idx <- round(nrow(covars) / 200)
idxx <- c(seq(1, nrow(covars), idx), nrow(covars))

cat("Location data loaded\n")

# Define possible RCPs and set seed for reproducibility
possible_RCPs <- c(5, 10)
set.seed(123)
nstarts <- 1000  # Recommended to use 1000 or more

# Loop over possible RCPs
for (i in possible_RCPs) {
    # Record start time
    multifit <- Sys.time()

    # Create output directory
    output_dir <- file.path(getwd(), "outputs/noForest", i)
    dir.create(output_dir, showWarnings = FALSE)

    cat("Processing RCP", i, "\n")

    # Fit regimix model with multiple initializations
    nRCPs_samp <- regimix.multifit(form.RCP = my_form_RCP, form.spp = my_form_spp, data = biodat0, 
                                   nRCP = i, inits = "random2", nstart = nstarts, dist = "Bernoulli", 
                                   mc.cores = 25, titbits = FALSE)
    save(nRCPs_samp, file = file.path(output_dir, paste0("regimix.multifit_", i, "RCP.RData")))

    # Evaluate posterior probabilities
    RCPsamp_minPosteriorSites <- sapply(nRCPs_samp, function(x) min(colSums(x$postProbs)))
    print(range(RCPsamp_minPosteriorSites))
    print(length(RCPsamp_minPosteriorSites))

    # Extract BIC values and evaluate minimum posterior probabilities
    dat <- sapply(nRCPs_samp, function(y) y$BIC)
    RCPsamp_minPost <- sapply(nRCPs_samp, function(x) min(colSums(x$postProbs)))
    mdat <- data.frame(dat, RCPsamp_minPost)
    colnames(mdat) <- c(i, "postp")

    # Filter out models with bad posterior probabilities
    RCPsamp_ObviouslyBad <- RCPsamp_minPosteriorSites <= 12
    print(table(RCPsamp_ObviouslyBad))
    dat[RCPsamp_ObviouslyBad] <- NA
    fm.clean <- nRCPs_samp[!RCPsamp_ObviouslyBad]
    print(length(fm.clean))

    # Select the best model based on BIC
    goodUn <- which.min(sapply(fm.clean, BIC))
    fm.final <- regimix(form.RCP = my_form_RCP, form.spp = my_form_spp, data = biodat0,
                        dist = "Bernoulli", nRCP = i, inits = unlist(fm.clean[[goodUn]]$coef),
                        control = list(optimise = FALSE))
    save(fm.final, file = file.path(output_dir, paste0("2_best_model_", i, "RCP.RData")))

    # Plot the best model
    pdf(file = file.path(output_dir, paste0("best_model_plot_", i, ".pdf")))
    plot(fm.final)
    dev.off()

    # Record bootstrap start time
    bootstr <- Sys.time()

    # Perform bootstrap
    fitregi <- regiboot(fm.final, nboot = 100, type = "BayesBoot", mc.cores = 30, quiet = FALSE)
    save(fitregi, file = file.path(output_dir, paste0("regiboot_1000_", i, "RCP.RData")))

    # Record CSV write start time
    csvwrite <- Sys.time()

    # Predict new data
    for (j in 1:(length(idxx) - 1)) {
        cat("Prediction batch", j, "\n")
        predicted.sp.list <- predict.regimix(object = fm.final, object2 = fitregi,
                                             newdata = covars[idxx[j]:(idxx[j + 1]),], mc.cores = 30)
        write.csv(cbind(covars[idxx[j]:(idxx[j + 1]), 1:2], predicted.sp.list$bootPreds),
                  file = file.path(output_dir, paste0(j, ".csv")))
        gc()
    }

    # Merge prediction results
    lst <- list.files(output_dir, pattern = ".csv", full.names = TRUE)
    lst <- mixedsort(lst)
    csvs <- lapply(lst, function(x) read.csv(x)[, -1])
    findat <- do.call(rbind, csvs)
    dat <- findat

    # Record raster write start time
    rasterwr <- Sys.time()

    # Plot prediction results
    pdf(file = file.path(output_dir, paste0(i, "_RCP_preds.pdf")), onefile = TRUE)
    colour <- c("#dddddd", "#fff5f0", "#fee0d2", "#fcbba1", "#fc9272", "#fb6a4a", "#ef3b2c", 
                "#FF0000", "#cb181d", "#a50f15")
    breaks <- seq(0, 1, by = 0.1)
    beginCluster()
    for (k in 3:(i + 2)) {
        print(k)
        plot(raster::rasterFromXYZ(dat[, c(1, 2, k)]), breaks = breaks, col = colour, main = names(dat)[k])
        raster::writeRaster(raster::rasterFromXYZ(dat[, c(1, 2, k)]), 
                            filename = file.path(output_dir, names(dat)[k]), overwrite = TRUE)
    }
    endCluster()
    dev.off()

    # Record completion time
    completpro <- Sys.time()
    a <- c(multifit, bootstr, csvwrite, rasterwr, completpro)
    write.csv(a, file = file.path(output_dir, paste0("process_time_", i, ".csv")))
}
