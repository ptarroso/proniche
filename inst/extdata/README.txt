These files are derived from:
- WorldClim v2 (Fick & Hijmans, 2017)
- GBIF occurrence data (GBIF.org, 2026)

Data were downloaded and processed for package examples using the following code:

occs <- geodata::sp_occurrence(genus = "Tuta", species = "absoluta", ext = vars, args = "occurrenceStatus=Present")  # downloaded 26 Jan 2026
occs <- subset(occs, select = c("lon", "lat"))
occs <- na.omit(occs)

vars <- geodata::worldclim_global(var = "bio", res = 5, path = "outputs")
names(vars) <- sub("wc2.1_5m_", "", names(vars))
vars <- vars[[c("bio_1", "bio_6", "bio_12", "bio_14")]]
vars <- terra::crop(vars, terra::ext(-10, 4, 36, 44))  # xmin, xmax, ymin, ymax


References

Fick S.E. & Hijmans R.J. (2017) WorldClim 2: new 1-km spatial resolution climate surfaces for global land areas. International Journal of Climatology 37:4302-4315
GBIF.org (26 January 2026) GBIF Occurrence Download https://doi.org/10.15468/dl

