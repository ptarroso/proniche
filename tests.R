methods <- c("bioclim", "domain", "convexhull",
             "mahalanobis", "kernel", "mvnormal")

methods[3]

tmp <- terra::rast(c("GIS/wc2.0_bio_5m_01.tif", "GIS/wc2.0_bio_5m_06.tif"))
prc <- terra::rast(c("GIS/wc2.0_bio_5m_12.tif", "GIS/wc2.0_bio_5m_14.tif"))
vars <- c(tmp, prc)
vals <- terra::extract(vars, chilus, ID=FALSE)
# vals <- vals[,1:2]
# vars <- vars[[1:2]]
vals <- vals[,1,drop=F]
vars <- vars[[1]]
head(vals)
plot(vars)
names(vars)

plot(model(vals, vars, method = methods[1])[[1]])
plot(model(vals, vars, method = methods[2])[[1]])
plot(model(vals, vars, method = methods[3])[[1]])
plot(model(vals, vars, method = methods[4])[[1]])
plot(model(vals, vars, method = methods[5])[[1]])
plot(model(vals, vars, method = methods[6])[[1]])



# test SpatRaster vs. df input ####

vals <- subset(dat, presence == 1, select = vars_sel)

method <- methods[1]
method <- methods[2]
method <- methods[3]
method <- methods[4]
method <- methods[5]
method <- methods[6]

m1 <- model(vals = vals, method = method, vars = layers_cut)
m2 <- model(vals = vals, method = method, vars = dat)

summary(m1[[1]])
as.matrix(summary(m2[[1]]))

plot(m1[[1]])
# points(dat[dat$presence == 1, c("x", "y")], pch = 20, col = "magenta", cex = 0.1)

dat$proniche <- m2[[1]]
dat_sv <- vect(dat, geom = c("x","y"))
# plot(dat_sv, "proniche", cex = 0.2)
plot(crop(dat_sv, layers_cut), "proniche", cex = 0.4, legend = F)
# points(dat[dat$presence == 1, c("x", "y")], pch = 20, col = "magenta", cex = 0.1)
