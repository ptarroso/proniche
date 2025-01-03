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

bc <- promodel(vals, vars, method = "bioclim")
ch <- promodel(vals, vars, method = "convexhull")
dm <- promodel(vals, vars, method = "domain")
mm <- promodel(vals, vars, method = "mahalanobis")
km <- promodel(vals, vars, method = "kernel")
mv <- promodel(vals, vars, method = "mvnormal")

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
