\name{model}
\alias{model}
\title{Compute model}
\description{
  Compute presence-only model predictions based on environmental variables at presence sites and a specified modeling method. Implemented methods are
  \itemize{
    \item{Bioclim (Nix, 1986; Booth et al., 2014)}
    \item{Convex hull (in environmental space)}
    \item{Domain (Carpenter et al., 1993)}
    \item{Mahalanobis distance (Farber & Kadmon, 2003)}
    \item{Kernel density estimate}
    \item{Multivariate normal distribution}
  }
}
\usage{
  model(vals, vars, method = "bioclim")
}
\arguments{
  \item{vals}{numeric values of the variables at the presence points.}

  \item{vars}{SpatRaster object with the variables onto which to map the model predictions.}

  \item{method}{character value indicating the method for computing the model. Available options are "bioclim" (the default), "convexhull", "domain", "mahalanobis", "kernel" or "mvnormal" (see Description). Partial argument matching is used, i.e. you can use e.g. method="mahal".}
}
\value{
  A list containing the prediction and model-specific parameters, depending on the chosen method.
}
\note{
  Some of these methods are sensitive to sampling and to outliers.
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
Booth T.H., Nix H.A., Busby J.R. & Hutchinson M.F. (2014) bioclim: the first species distribution modelling package, its early applications and relevance to most current MaxEnt studies. Diversity & Distributions, 20(1):1-9

Carpenter G., Gillison A. N. Gillison & Winter J. (1993) DOMAIN: a flexible modelling procedure for mapping potential distributions of plants and animals. Biodiversity and Conservation 2(6):667-680

Farber O. & Kadmon R. (2003) Assessment of alternative approaches for bioclimatic modeling with special emphasis on the Mahalanobis distance. Ecological Modelling
160 (1-2):115-130

Nix, H.A. (1986) A biogeographic analysis of Australian elapid snakes. Atlas of elapid snakes of Australia: Australian flora and fauna series 7 (ed. by R. Longmore), pp. 4-15. Bureau of Flora and Fauna, Canberra.
}
\keyword{models}
\keyword{spatial}
\keyword{multivariate}
\keyword{prediction}
\author{
Pedro Tarroso and A. Marcia Barbosa
}
