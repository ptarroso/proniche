\name{freqPlot}
\alias{freqPlot}
\title{Frequency plot}
\description{
  This function plots (in green) the frequency distribution of each variable at the presence localities, and optionally (in grey) also the frequency across a given study region, using either density or box plots.
}
\usage{
  freqPlot(vals, vars = NULL, type = "density", na.rm = TRUE, dup.rm = FALSE,
  verbosity = 2, ...)
}
\arguments{
  \item{vals}{
    A data frame with the values of the variables (columns) at the presence localities (rows). Column names must match the names of the corresponding variables in \code{vars}.
  }
  \item{vars}{
    An optional \code{SpatRaster} or data frame with the variables to visually compare against the frequency of \code{vals}.
  }
  \item{type}{
    The type of plot to produce: \code{"density"} or \code{"boxplot"}. Case is ignored, and partial argument matching is allowed (e.g., \code{type = "Box"}).
  }
  \item{na.rm}{
    Logical (default \code{TRUE}) indicating whether to remove rows with missing values in \code{vals}. Otherwise, some methods may produce an error or empty results.
  }
  \item{dup.rm}{
    Logical (default \code{FALSE}) indicating whether to remove rows with duplicate values for all variables in \code{vals}.
  }
  \item{verbosity}{
    Numeric value indicating the level of messages to display. Default is \code{2}, for the maximum number of messages available.
  }
  \item{...}{
    Optional additional arguments to be passed to \code{plot()}, such as \code{cex.main} or \code{cex.axis}.
  }
}
\details{
  This function provides a simple way to visualize the frequency distribution of one or more variables, optionally comparing the distribution of \code{vals} with \code{vars}. It supports two types of plots: density plots and boxplots. When \code{vars} is provided, the function overlays the comparison dataset onto the same plot.
}
\value{
  A plot of the specified type.
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

  # Produce frequency plots:
  freqPlot(vals)
  freqPlot(vals, vars)
  freqPlot(vals, vars, type = "boxplot")
  }
}
