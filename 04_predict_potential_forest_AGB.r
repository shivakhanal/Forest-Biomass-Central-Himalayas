# Load required libraries
library(raster)
library(gtools)
library(gmodels)

# Read fixed coordinates
ptrain <- shapefile("/data1/mounts/s01238ss-wsu-u0217-hie/Khanal_S/Chap4/data/ptrain.shp")
ptest <- shapefile("/data1/mounts/s01238ss-wsu-u0217-hie/Khanal_S/Chap4/data/ptest.shp")

# Set paths and variables
path <- "/outputs/optimumRCP/5/"
nRCP <- 5
subDir <- "PotentialAGB"
dir.create(file.path(subDir))
outpath <- subDir

# Load hard cluster and RCPs
hard <- raster(list.files(path, pattern="hard_class_.*\\.tif$", full.names=TRUE))
rcps <- mixedsort(list.files(path, pattern=".grd$", full.names=TRUE))
rcp <- stack(lapply(rcps, raster))

# Make data frame for prediction
rcpdat <- as.data.frame(rcp, xy=TRUE)
roudat <- rcpdat[!rowSums(is.na(rcpdat)) == nRCP,]

# Transform and extract values for train and test plots
plottrain <- spTransform(ptrain, projection(hard))
plottrain$rcp <- extract(hard, plottrain)
plottrain$pmax <- apply(extract(rcp, plottrain), 1, max)

plottest <- spTransform(ptest, projection(hard))
plottest$rcp <- extract(hard, plottest)
plottest$pmax <- apply(extract(rcp, plottest), 1, max)

# Predict potential AGB based on bootstrap
fintrain <- plottrain@data
plotsp <- SpatialPointsDataFrame(fintrain[,2:3], fintrain, proj4string=crs("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
plotsp <- spTransform(plotsp, projection(hard))
plotsp$rcp <- extract(hard, plotsp)

baha <- aggregate(plotsp@data$agb, by=list(plotsp@data$rcp), FUN=function(x) quantile(x, probs=0.95))
names(baha) <- c("rcp", "peragb")

# Select plots that have more than 95 percentile AGB
traindata <- merge(plotsp, baha, by="rcp")
traindata <- traindata@data

dat <- lapply(unique(traindata$rcp), function(i) {
  xx <- traindata[traindata$rcp == i,]
  xx[xx$agb >= xx$peragb[1],]
})
dat <- do.call("rbind", dat)

source("data/Chap4/functions/stratified.R")
source("/code/stratified.R")

# Bootstrap sampling and weighted AGB calculation
set.seed(123)
beginCluster()
raslst <- list()

for (p in 1:100) {
  cat("iteration =", p, "of 100\n")
  baha <- stratified(dat, group="rcp", size=1)
  baha$rcp <- paste0("RCP_", baha$rcp)
  rcp_selected <- rcp[[baha$rcp]]
  
  rasl <- lapply(1:nlayers(rcp_selected), function(i) {
    mras <- rcp_selected[[i]]
    if (names(mras) %in% unique(baha$rcp)) {
      mult <- baha[baha$rcp == names(mras),]$agb
      mras * mult
    }
  })
  raslst[[p]] <- calc(stack(rasl), sum, na.rm=TRUE)
}
endCluster()

newstk <- stack(raslst)
projection(newstk) <- "+proj=utm +zone=45 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

# Write percentile outputs
dir.create(file.path(outpath, "bootWeighted/100/"), recursive=TRUE)

beginCluster()
for (i in 1:100) {
  cat("iteration =", i, "of 100\n")
  writeRaster(raslst[[i]], paste0(outpath, "/bootWeighted/100/", i, ".tif"))
}
endCluster()

# Load and process bootstrap results
lst <- list.files(file.path(path, outpath, "bootWeighted/100"), full.names=TRUE)
stk <- stack(lapply(lst, raster))

kfun <- function(x) ci(x, confidence=0.95, na.rm=TRUE)
fun <- function(x) { x[x <= 0] <- NA; return(x) }

beginCluster()
stk <- clusterR(stk, calc, args=list(fun=fun))
endCluster()

beginCluster()
ov <- clusterR(stk, calc, args=list(fun=kfun))
names(ov) <- c("Estimate", "CI lower", "CI upper", "Std. Error")
endCluster()

writeRaster(ov[[1]], file.path(path, outpath, "mean_boot_potential.tif"))
writeRaster(ov[[2]], file.path(path, outpath, "LCI_potential.tif"))
writeRaster(ov[[3]], file.path(path, outpath, "UCI_potential.tif"))
writeRaster(ov[[4]], file.path(path, outpath, "SE_potential.tif"))
