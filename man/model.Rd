\name{model}
\alias{model}
\title{Compute model}
\description{
  Compute model predictions based on environmental variables and a selected modeling method.
}
\usage{
  model(vals, vars, method = "bioclim")
}
\arguments{
  \item{vals}{numeric values of the variables at the presence points.}

  \item{vars}{SpatRaster object with the variables onto which to obtain the model predictions.}

  \item{method}{character value indicating the method for computing the model. Available options are "bioclim" (the default), "convexhull", "kernel", "domain", or "mvnormal".}
}
\value{
  A list containing the prediction and model-specific parameters, depending on the chosen method.
}
\examples{
\dontrun{
  # Example data
  library(terra)
  vals <- matrix(runif(50), ncol = 2)
  vars <- rast(matrix(runif(100), ncol = 10))

  result <- model(vals, vars, method = "bioclim")
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
