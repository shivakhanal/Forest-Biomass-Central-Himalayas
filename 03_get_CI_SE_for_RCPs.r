# Load required libraries
library(RCPmod)
library(raster)
library(gtools)
library(data.table)
library(rasterVis)
library(RColorBrewer)

# Read input data
biodat0 <- read.csv("data/presabs_data_poly.csv")[,-1]
su <- colSums(biodat0[,4:236])
su <- su[su <= 3] 
biodat0 <- biodat0[, !names(biodat0) %in% names(su)]

model.covariates.vector <- c("bio12", "bio15", "gdd0", "bio4", "pet", "twi", "daily_PISR") 
model.covariates.string <- paste(model.covariates.vector, collapse=" + ")
st <- paste(names(biodat0)[4:242], collapse=",")
my.form.RCP <- paste('cbind(', st, ') ~ ', model.covariates.string, sep='')
my.form.spp <- NULL

# Load predictors
covars <- fread("data/predictors_250m_poly.csv", header=TRUE)
covars <- as.data.frame(covars)[,-1]
covars <- covars[, c("x", "y", model.covariates.vector)]
covars <- covars[complete.cases(covars),]
locs <- covars[, 1:2]

# Set up tiling for prediction
idx <- round(nrow(covars) / 200)
idxx <- c(seq(1, nrow(covars), idx), nrow(covars))

# Read the fitted model objects and export results
outpath <- 'data/noForest/5/5_1000_iter'
lst <- list.files(outpath, pattern=".RData", full.names=TRUE)
load(lst[1])
load(lst[2])

predicted.sp.list <- list()
for (j in 1:201) {
  predicted.sp.list[[j]] <- predict.regimix(object=fm.final, object2=fitregi, 
                                            newdata=covars[idxx[j]:(idxx[j+1]-1),], mc.cores=25)
}

dat <- do.call(rbind.data.frame, lapply(predicted.sp.list, function(i) cbind(i$bootPreds[,])))
dat.low <- do.call(rbind.data.frame, lapply(predicted.sp.list, function(i) cbind(i$bootCIs[,,1])))
dat.up <- do.call(rbind.data.frame, lapply(predicted.sp.list, function(i) cbind(i$bootCIs[,,2])))
dat <- cbind(locs[1:nrow(dat),], dat)
dat.low <- cbind(locs[1:nrow(dat.low),], dat.low)
dat.up <- cbind(locs[1:nrow(dat.low),], dat.up)

# Save results
dir.create(file.path(outpath, "CI"))
cipath <- file.path(outpath, "CI")
write.csv(dat, paste0(cipath, "/", "bootPreds.csv"))
write.csv(dat.low, paste0(cipath, "/", "dat.low.csv"))
write.csv(dat.up, paste0(cipath, "/", "dat.up.csv"))

beginCluster()
for (k in 3:(5+2)) { # The first two are lat, lon
  raster::writeRaster(raster::rasterFromXYZ(dat.low[,c(1,2,k)]), 
                      filename=paste0(cipath, "/", names(dat.low)[k], "_lo"), overwrite=TRUE)
  raster::writeRaster(raster::rasterFromXYZ(dat.up[,c(1,2,k)]), 
                      filename=paste0(cipath, "/", names(dat.up)[k], "_up"), overwrite=TRUE)
  raster::writeRaster(raster::rasterFromXYZ(dat[,c(1,2,k)]), 
                      filename=paste0(cipath, "/", names(dat)[k], "_pred"), overwrite=TRUE)  
}
endCluster()

# Plot the raster along with CI layers
myfl <- list.files(cipath, pattern=".grd$", full.names=TRUE)
myfl <- mixedsort(myfl)
myfl <- myfl[!grepl("_se", myfl)] 
nRCP <- 5
myfl <- myfl[c(seq(1, nRCP*3, 3), seq(2, nRCP*3, 3), seq(3, nRCP*3, 3))] 

colour <- c("#dddddd", "#fff5f0", "#fee0d2", "#fcbba1", "#fc9272", "#fb6a4a", "#ef3b2c", "#FF0000", "#cb181d", "#a50f15")
breaks <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)  
rasl <- lapply(myfl, raster) 
s <- stack(rasl)
aa <- lapply(1:nRCP, function(i) seq(i, nRCP*3, nRCP))
s <- s[[unlist(aa)]] 
nl <- nlayers(s)
m <- matrix(1:nl, ncol=3)
my.at <- breaks
myColorkey <- list(at=my.at, height=.7, labels=list(at=my.at))
rainbTheme10 <- rasterTheme(region=colour)

for (i in 1:nl) {
  p <- levelplot(s, layers=i, at=breaks, colorkey=myColorkey, par.settings=rainbTheme10, main=names(s)[i], 
                 between=list(x=0.0, y=0.0), margin=FALSE)
  print(p, split=c(col(m)[i], row(m)[i], ncol(m), nrow(m)), more=(i < nl))
}
dev.off()

mypalette <- c("#dddddd", brewer.pal(9, "Reds")[-1])  
rainbTheme10 <- rasterTheme(region=mypalette)

pdf(paste0(cipath, "/", "5_RCP_CI.pdf"))
png(paste0(cipath, "/", "5_RCP_CI.png"), width=1280, height=1880, units="px")
col.titles <- c('Lower CI', 'Membership Prob', 'Upper CI', rep("", (nRCP*3-3))) 
row.titles <- paste0('RCP-', nRCP:1)
levelplot(s, layout=c(3, nRCP), at=breaks, colorkey=myColorkey, par.settings=rainbTheme10, maxpixels=2e+06, 
          names.attr=col.titles, ylab=row.titles)
grid::grid.text('Probability', rot=90, y=unit(0.5, "npc"), x=unit(0.90, "npc"))
dev.off()
