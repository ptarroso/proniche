---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->



# proniche: PResence-Only NICHE modelling

R package for computing niche models using a variety of **truly presence-only methods** (i.e., methods that use only environmental info from the presence sites, with no background or pseudo-absence data):

- Bioclim (rectangular environmental envelope)
- Convex hull (in environmental space)
- Domain
- Mahalanobis distance
- Kernel density estimate
- Multivariate normal distribution


## User guide

### (Install and) load package


``` r
# devtools::install_github("https://github.com/ptarroso/proniche")

library(proniche)
```

### Import some example data


``` r
tmp <- terra::rast(c("GIS/wc2.0_bio_5m_01.tif", "GIS/wc2.0_bio_5m_06.tif"))
prc <- terra::rast(c("GIS/wc2.0_bio_5m_12.tif", "GIS/wc2.0_bio_5m_14.tif"))
vars <- c(tmp, prc)
terra::plot(vars)
```

<div class="figure">
<img src="man/figures/README-usage-1.png" alt="plot of chunk usage" width="100%" />
<p class="caption">plot of chunk usage</p>
</div>

``` r

chilus <- read.csv("GIS/chilus.csv", sep=";")
vals <- terra::extract(vars, chilus, ID=FALSE)
terra::plot(vars[[1]] * 0, col = "tan", background = "lightblue",
            legend = FALSE, main = "Presences")
points(chilus, pch = 20, cex = 0.2)
```

<div class="figure">
<img src="man/figures/README-usage-2.png" alt="plot of chunk usage" width="100%" />
<p class="caption">plot of chunk usage</p>
</div>

### Plot frequency distributions


``` r
par(mfrow = c(2, 2), mar = c(2, 2, 2, 1))
freqPlot(vals, vars)
```

<div class="figure">
<img src="man/figures/README-freqPlot-1.png" alt="plot of chunk freqPlot" width="100%" />
<p class="caption">plot of chunk freqPlot</p>
</div>



### Fit models


``` r
bc_fit <- promodel(vals, method = "bioclim")
ch_fit <- promodel(vals, method = "convexhull")
dm_fit <- promodel(vals, method = "domain")
mm_fit <- promodel(vals, method = "mahalanobis")
km_fit <- promodel(vals, method = "kernel")
mv_fit <- promodel(vals, method = "mvnormal")
```

### Check built models

``` r
par(mfrow = c(2, 3))
plot(bc_fit, main = "Bioclim")
plot(ch_fit, main="Convex Hull")
plot(dm_fit, main="Domain")
plot(mm_fit, main="Mahalanobis")
plot(km_fit, main="Kernel")
plot(mv_fit, main="Multivariate Normal")
```

<div class="figure">
<img src="man/figures/README-checkmodel-1.png" alt="plot of chunk checkmodel" width="100%" />
<p class="caption">plot of chunk checkmodel</p>
</div>

### Predict with models


``` r
bc <- predict(bc_fit, vars)
ch <- predict(ch_fit, vars)
dm <- predict(dm_fit, vars)
mm <- predict(mm_fit, vars)
km <- predict(km_fit, vars)
mv <- predict(mv_fit, vars)

par(mfrow = c(2, 3))
terra::plot(bc, type = "continuous", main="Bioclim")
terra::plot(ch, type = "continuous", main="Convex Hull")
terra::plot(dm, type = "continuous", main="Domain")
terra::plot(mm, type = "continuous", main="Mahalanobis")
terra::plot(km, type = "continuous", main="Kernel")
terra::plot(mv, type = "continuous", main="Multivariate Normal")
```

<div class="figure">
<img src="man/figures/README-predictmodel-1.png" alt="plot of chunk predictmodel" width="100%" />
<p class="caption">plot of chunk predictmodel</p>
</div>


### Reclassify predictions into comparable scale


``` r
bc_rcl <- quantReclass(bc[[1]])
ch_rcl <- quantReclass(ch[[1]])
dm_rcl <- quantReclass(dm[[1]])
mm_rcl <- quantReclass(mm[[1]])
km_rcl <- quantReclass(km[[1]])
mv_rcl <- quantReclass(mv[[1]])

par(mfrow = c(2, 3))
terra::plot(bc_rcl, range = c(0, 1), type = "continuous", main="Bioclim")
terra::plot(ch_rcl, range = c(0, 1), type = "continuous", main="Convex Hull")
terra::plot(dm_rcl, range = c(0, 1), type = "continuous", main="Domain")
terra::plot(mm_rcl, range = c(0, 1), type = "continuous", main="Mahalanobis")
terra::plot(km_rcl, range = c(0, 1), type = "continuous", main="Kernel")
terra::plot(mv_rcl, range = c(0, 1), type = "continuous", main="Multivariate Normal")
```

<div class="figure">
<img src="man/figures/README-reclass-1.png" alt="plot of chunk reclass" width="100%" />
<p class="caption">plot of chunk reclass</p>
</div>

### Ensemble predictions


``` r
preds <- c(bc_rcl, ch_rcl, dm_rcl, mm_rcl, km_rcl, mv_rcl)
ens_mean <- terra::app(preds, "mean")
ens_var <- terra::app(preds, "var")
par(mfrow = c(1, 2))
terra::plot(ens_mean, range = c(0, 1), main = "Ensemble mean")
terra::plot(ens_var, range = c(0, 1), main = "Ensemble variance")
```

<div class="figure">
<img src="man/figures/README-ensemble-1.png" alt="plot of chunk ensemble" width="100%" />
<p class="caption">plot of chunk ensemble</p>
</div>


## TO DO

- ~~Add multivariate normal distribution~~
- ~~Maybe add resampling to Mahalanobis (needs to estimate covariance)~~
- ~~Create an R package~~
- ~~Add documentation and references to original model publications~~
- ~~Maybe some accessory functions to clean data set? (Remove duplicates, NAs, etc)~~
- ~~Instead of receiving coordinates as input, perhaps modelling functions should get already the environmental values at presences?~~
- Function to plot models representation in environmental space?

- ~~Implement for data frames (not just SpatRasters)~~
- ~~Add function to reclassify predictions into comparable scale~~
