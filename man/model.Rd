\name{model}
\alias{model}
\title{Compute model}
\description{
  Compute presence-only model predictions based on environmental variables and a specified modeling method.
}
\usage{
  model(vals, vars, method = "bioclim")
}
\arguments{
  \item{vals}{numeric values of the variables at the presence points.}

  \item{vars}{SpatRaster object with the variables onto which to map the model predictions.}

  \item{method}{character value indicating the method for computing the model. Available options are "bioclim" (the default), "convexhull", "domain", "mahalanobis", "kernel" or "mvnormal". Partial argument matching is used, i.e. you can use e.g. method="mahal".}
}
\value{
  A list containing the prediction and model-specific parameters, depending on the chosen method.
}
\examples{
\dontrun{
  # Get some predictor variables:
  library(terra)
  elevation <- terra::rast(system.file("ex/elev.tif", package="terra"))
  slope <- terra::terrain(elev, v = "slope")
  roughness <- terra::terrain(elev, v = "roughness")
  vars <- c(elevation, slope, roughness)
  plot(vars)

  # Get some random presence points and their variables values:
  vars_coords <- terra::crds(vars)
  set.seed(123)
  pres_coords <- vars_coords[sample(1:nrow(vars_coords), 20, replace = FALSE), ]
  vals <- terra::extract(vars, pres_coords)
  head(vals)

  # Compute a model:
  mod <- model(vals, vars, method = "mahalanobis")
  terra::plot(mod[[1]])
}
}
\seealso{
\code{\link[dismo:bioclim]{dismo::bioclim()}}, \code{\link[predicts:envelope]{predicts::envelope()}}, \code{\link[geometry:convhulln]{geometry::convhulln()}}, \code{\link[ks:kde]{ks::kde()}}
}
\references{
}
\keyword{models}
\keyword{spatial}
\keyword{multivariate}
\keyword{prediction}
\author{
Pedro Tarroso and A. Marcia Barbosa
}
