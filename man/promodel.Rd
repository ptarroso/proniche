\name{promodel}
\alias{promodel}
\title{Compute presence-only model}
\description{
  Compute truly presence-only model predictions based on environmental variables at species occurrence sites, using a specified modeling method.}
\usage{
  promodel(vals, vars, method = "bioclim", na.rm = TRUE, dup.rm = FALSE,
  verbosity = 2)
}
\arguments{
  \item{vals}{data frame with the values of the variables (columns) at the presence localities (rows). Column names must match the names of the corresponding variables in 'vars'.}
  \item{vars}{SpatRaster or data frame with the variables onto which to obtain the model predictions.}
  \item{method}{method for computing the model. The options are "bioclim" (the default), "convexhull", "domain", "mahalanobis", "kernel" or "mvnormal" (see Details). Case is ignored and partial argument matching is used (i.e., you can set e.g. method="Mahal").}
  \item{na.rm}{logical (default TRUE) indicating whether to remove rows with missing values in 'vals'. Otherwise, some methods may produce an error or empty results.}
  \item{dup.rm}{logical (default FALSE) indicating whether to remove rows with duplicate values for all variables in 'vals'.}
  \item{verbosity}{numeric value indicating the amount of messages to display. The default is 2, for the maximum number of messages available.}
}
\details{
Implemented methods are:
  \itemize{
    \item{Bioclim (Nix, 1986; Booth et al., 2014)}
    \item{Convex hull in environmental space (using \code{\link[geometry:convhulln]{geometry::convhulln()}} and \code{\link[geometry:inhulln]{geometry::inhulln()}})}
    \item{Domain (Carpenter et al., 1993)}
    \item{Mahalanobis distance (Farber & Kadmon, 2003)}
    \item{Kernel density estimate (using \code{\link[ks:kde]{ks::kde()}}}
    \item{Multivariate normal distribution}
  }
}
\value{
  A list whose first element is the prediction (a SpatRaster if 'vars' is of this class, or a numeric vector otherwise). There may be more elements containing model-specific parameters, depending on the chosen 'method'.
}
\note{
  Some of these methods are quite sensitive to outliers and to sampling representation.
}
\examples{
\dontrun{
  # Get some predictor variables:
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

  # Compute a model:
  mod <- promodel(vals, vars, method = "mahalanobis")
  terra::plot(mod[[1]])
}
}
\seealso{
\code{\link{quantReclass}}; \code{\link[dismo:bioclim]{dismo::bioclim()}}, \code{\link[predicts:envelope]{predicts::envelope()}}, \code{\link[geometry:convhulln]{geometry::convhulln()}}, \code{\link[ks:kde]{ks::kde()}}
}
\references{
Booth T.H., Nix H.A., Busby J.R. & Hutchinson M.F. (2014) bioclim: the first species distribution modelling package, its early applications and relevance to most current MaxEnt studies. Diversity & Distributions, 20(1):1-9

Carpenter G., Gillison A. N. Gillison & Winter J. (1993) DOMAIN: a flexible modelling procedure for mapping potential distributions of plants and animals. Biodiversity and Conservation 2(6):667-680

Farber O. & Kadmon R. (2003) Assessment of alternative approaches for bioclimatic modeling with special emphasis on the Mahalanobis distance. Ecological Modelling
160 (1-2):115-130

Nix, H.A. (1986) A biogeographic analysis of Australian elapid snakes. Atlas of elapid snakes of Australia: Australian flora and fauna series 7 (ed. by R. Longmore), pp. 4-15. Bureau of Flora and Fauna, Canberra.
}
\keyword{models}
\keyword{spatial}
\keyword{prediction}
\author{
Pedro Tarroso and A. Marcia Barbosa
}
