These files are derived from:

Fick S.E. & Hijmans R.J. (2017) WorldClim 2: new 1-km spatial resolution climate surfaces for global land areas. International Journal of Climatology 37:4302-4315. https://www.worldclim.org/

GBIF.org (26 January 2026) GBIF Occurrence Download https://doi.org/10.15468/dl


The data were downloaded and processed for package examples using the following code:

vars <- geodata::worldclim_global(var = "bio", res = 5, path = tempdir())
names(vars) <- sub("wc2.1_5m_", "", names(vars))  # simplifies layer names
vars <- vars[[c("bio_1", "bio_6", "bio_12", "bio_14")]]
vars <- terra::crop(vars, terra::ext(-10, 4, 36, 44))  # xmin, xmax, ymin, ymax
writeRaster(vars, "vars.tif", gdal = c("COMPRESS=DEFLATE"))

occs <- geodata::sp_occurrence(genus = "Tuta", species = "absoluta", ext = vars, args = "occurrenceStatus=Present")  # downloaded 26 Jan 2026
occs <- subset(occs, select = c("lon", "lat"))
occs <- na.omit(occs)
write.csv(occs, "occs.csv", row.names = FALSE)

