library(terra)
library(proniche)

methods <- c("bioclim", "domain", "convexhull",
             "mahalanobis", "kernel", "mvnormal")

methods[3]

tmp <- terra::rast(c("GIS/wc2.0_bio_5m_01.tif", "GIS/wc2.0_bio_5m_06.tif"))
prc <- terra::rast(c("GIS/wc2.0_bio_5m_12.tif", "GIS/wc2.0_bio_5m_14.tif"))
vars <- c(tmp, prc)
chilus <- read.csv("GIS/chilus.csv", sep = ";")
vals <- terra::extract(vars, chilus, ID=FALSE)
head(vals)

# test with only 2 variables:
# vals <- vals[,1:2]
# vars <- vars[[1:2]]

# test with only 1 variable:
# vals <- vals[,1,drop=F]
# vars <- vars[[1]]

terra::plot(vars)
names(vars)

# test with NAs:
vals[1:3, 1:3] <- NA
head(vals)

# test duplicates:
vals[nrow(vals) + 1, ] <- vals[nrow(vals), ]
tail(vals)

bc <- promodel(vals, vars, method = "bioclim", dup.rm=T)
ch <- promodel(vals, vars, method = "convexhull", dup.rm=T)
dm <- promodel(vals, vars, method = "domain", dup.rm=T)
mm <- promodel(vals, vars, method = "mahalanobis", dup.rm=T)
km <- promodel(vals, vars, method = "kernel", dup.rm=T)
mv <- promodel(vals, vars, method = "mvnormal", dup.rm=T)

plot(bc[[1]], main = "bioclim")
plot(ch[[1]], main = "convexhull")
plot(dm[[1]], main = "domain")
plot(mm[[1]], main = "mahalanobis")
plot(km[[1]], main = "kernel")
plot(mv[[1]], main = "mvnormal")

plot(quantReclass(bc[[1]]), type = "continuous", main = "bioclim")
plot(quantReclass(ch[[1]]), type = "continuous", main = "convexhull")
plot(quantReclass(dm[[1]]), type = "continuous", main = "domain")
plot(quantReclass(mm[[1]]), type = "continuous", main = "mahalanobis")
plot(quantReclass(km[[1]]), type = "continuous", main = "kernel")
plot(quantReclass(mv[[1]]), type = "continuous", main = "mvnormal")


# test SpatRaster vs. df input ####

vals <- subset(dat, presence == 1, select = vars_sel)

method <- methods[1]
method <- methods[2]
method <- methods[3]
method <- methods[4]
method <- methods[5]
method <- methods[6]

m1 <- promodel(vals = vals, method = method, vars = layers_cut)
m2 <- promodel(vals = vals, method = method, vars = dat)

summary(m1[[1]])
as.matrix(summary(m2[[1]]))

plot(m1[[1]])
# points(dat[dat$presence == 1, c("x", "y")], pch = 20, col = "magenta", cex = 0.1)

dat$proniche <- m2[[1]]
dat_sv <- vect(dat, geom = c("x","y"))
# plot(dat_sv, "proniche", cex = 0.2)
plot(crop(dat_sv, layers_cut), "proniche", cex = 0.4, legend = F)
# points(dat[dat$presence == 1, c("x", "y")], pch = 20, col = "magenta", cex = 0.1)



# import example data from the web ####

occ_raw <- geodata::sp_occurrence(genus = "Chioglossa",
                                  species = "lusitanica",
                                  args = c("year=2022,2024"),
                                  fixnames = FALSE)

occ_clean <- fuzzySim::cleanCoords(data = occ_raw,
                                   coord.cols = c("decimalLongitude", "decimalLatitude"),
                                   uncert.col = "coordinateUncertaintyInMeters",
                                   uncert.limit = 10000,
                                   abs.col = "occurrenceStatus",
                                   plot = FALSE)

occ_coords <- occ_clean[ , c("decimalLongitude", "decimalLatitude")]



vars <- geodata::worldclim_global(var = "bio", res = 5, path = ".")


countries <- geodata::world(path = ".")
ib <- subset(countries, countries$NAME_0 %in% c("Portugal", "Spain"))
ib <- terra::disagg(ib)
largest <- (order(terra::expanse(ib), decreasing = TRUE))[1:2]
ib <- ib[largest, ]

vars <- terra::crop(vars, ib)
terra::plot(vars)


vals <- terra::extract(vars, occ_coords, ID=FALSE)
head(vals)

terra::plot(vars[[1]] * 0, col = "tan", background = "lightblue",
            legend = FALSE, main = "Presences")
points(occ_coords, pch = 20, cex = 0.2)



# test bioclim with different numbers of variables ####

par(mfrow = c(5, 4))
for (i in 1:terra::nlyr(vars)){
  bc <- promodel(vals[ , 1:i, drop = FALSE], vars, method = "bioclim")
  terra::plot(bc[[1]], type = "continuous", main=paste("Bioclim", i), plg = list(cex = 0.5, size = c(1, 0.3), tic = "none", border = NA), axes = FALSE, mar = c(0, 0.1, 1.5, 2.5))
  message(i, " range: ", range(terra::values(bc[[1]]), na.rm = TRUE))
}


# compare bioclim proniche with dismo / predicts ####

bc <- proniche::promodel(vals, vars, method = "bioclim")
bcd <- dismo::bioclim(raster::stack(vars), occ_coords)
bcp <- predicts::envelope(vars, occ_coords)


par(mfrow = c(3, 1))

plot(predict(vars, bcd), main = "dismo")
plot(predict(vars, bcp), main = "predicts")
plot(bc[[1]], type = "continuous", main = "proniche")

plot(quantReclass(predict(vars, bcd)), main = "dismo reclass")
plot(quantReclass(predict(vars, bcp)), main = "predicts reclass")
plot(quantReclass(bc[[1]]), type = "continuous", main = "proniche reclass")


# test convex hull with different numbers of variables ####

ch <- promodel(vals[,c(1, 5, 12)], vars[,c(1, 5, 12)], method = "convexhull")



# frequency plots ####

par(mfrow = c(5, 4), mar = c(2, 2, 2, 0.5))

freqPlot(vals)

freqPlot(vals, vars, cex.main = 0.8, cex.axis = 0.6)

freqPlot(vals, vars, type = "boxplot", cex.main = 0.8, cex.axis = 0.6)

freqPlot(vals, type = "boxplot", cex.main = 0.8, cex.axis = 0.6)


library(terra)
elevation <- terra::rast(system.file("ex/elev.tif", package = "terra"))
slope <- terra::terrain(elevation, v = "slope")
roughness <- terra::terrain(elevation, v = "roughness")
vars <- c(elevation, slope, roughness)
plot(vars)

# Get some hypothetical (random) presence points and their variables values:
vars_coords <- terra::crds(vars)
set.seed(123)
pres_coords <- vars_coords[sample(1:nrow(vars_coords), 20, replace = FALSE), ]
vals <- terra::extract(vars, pres_coords)
head(vals)

# Produce frequency plots:
freqPlot(vals)
freqPlot(vals, vars)
freqPlot(vals, vars, type = "boxplot")
