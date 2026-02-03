#' Compute presence-only model
#'
#' Compute true presence-only models based on environmental variables at
#' species occurrence sites, using a specified modeling method.
#'
#' @param vals data frame with the values of the variables (columns) at the
#'   presence localities (rows).
#' @param method method for computing the model. Options are `"bioclim"`
#'   (default), `"convexhull"`, `"domain"`, `"mahalanobis"`, `"kernel"` or
#'   `"mvnormal"`. Case is ignored and partial matching is allowed (e.g.
#'   `method = "Mahal"`).
#' @param na.rm logical (default `TRUE`) indicating whether to remove rows with
#'   missing values in `vals`. Otherwise, some methods may produce an error or
#'   empty results.
#' @param dup.rm logical (default `FALSE`) indicating whether to remove rows
#'   with duplicate values for all variables in `vals`.
#' @param verbosity numeric value indicating the amount of messages to display.
#'   Default is 2 (maximum verbosity).
#' @param ... Additional arguments passed to the underlying model functions,
#'   such as \code{\link{bioclim}}, \code{\link{convexhull}},
#'   \code{\link{kernel}}, or \code{\link{mvnormal}}.
#'
#' @details
#' Implemented methods include:
#' \itemize{
#'   \item Bioclim (Nix, 1986; Booth et al., 2014)
#'   \item Convex hull in environmental space (using
#'     \code{\link[geometry:convhulln]{geometry::convhulln()}} and
#'     \code{\link[geometry:inhulln]{geometry::inhulln()}})
#'   \item Domain (Carpenter et al., 1993)
#'   \item Mahalanobis distance (Farber & Kadmon, 2003)
#'   \item Kernel density estimate (using \code{\link[ks:kde]{ks::kde()}})
#'   \item Multivariate normal distribution
#' }
#'
#' @return An object of class `proniche`, containing the fitted model in the
#'   element `proniche`. The internal structure of this element depends on the
#'   chosen `method`.
#'
#' @examples
#' \dontrun{
#' library(terra)
#'
#' # Predictor variables
#' elevation <- terra::rast(system.file("ex/elev.tif", package = "terra"))
#' slope <- terra::terrain(elevation, v = "slope")
#' roughness <- terra::terrain(elevation, v = "roughness")
#' vars <- c(elevation, slope, roughness)
#' terra::plot(vars)
#'
#' # Random presence points
#' vars_coords <- terra::crds(vars)
#' set.seed(123)
#' pres_coords <- vars_coords[sample(1:nrow(vars_coords), 20), ]
#' vals <- terra::extract(vars, pres_coords)
#' head(vals)
#'
#' # Compute and plot a model
#' mod <- promodel(vals, method = "mahalanobis")
#' plot(mod)
#' pred <- terra::predict(vars, mod)
#' terra::plot(pred)
#' }
#'
#' @seealso
#'   \code{\link{quantReclass}};
#'   \code{\link[dismo:bioclim]{dismo::bioclim()}};
#'   \code{\link[predicts:envelope]{predicts::envelope()}};
#'   \code{\link[geometry:convhulln]{geometry::convhulln()}};
#'   \code{\link[ks:kde]{ks::kde()}}
#'
#' @references
#' Booth T.H., Nix H.A., Busby J.R. & Hutchinson M.F. (2014)
#'   bioclim: the first species distribution modelling package, its early
#'   applications and relevance to most current MaxEnt studies.
#'   *Diversity & Distributions*, 20(1):1–9.
#'
#' Carpenter G., Gillison A.N. & Winter J. (1993)
#'   DOMAIN: a flexible modelling procedure for mapping potential distributions
#'   of plants and animals. *Biodiversity and Conservation*, 2(6):667–680.
#'
#' Farber O. & Kadmon R. (2003)
#'   Assessment of alternative approaches for bioclimatic modeling with special
#'   emphasis on the Mahalanobis distance. *Ecological Modelling*,
#'   160(1–2):115–130.
#'
#' Nix H.A. (1986)
#'   A biogeographic analysis of Australian elapid snakes.
#'   *Atlas of elapid snakes of Australia: Australian flora and fauna series 7*
#'   (ed. R. Longmore), pp. 4–15. Bureau of Flora and Fauna, Canberra.
#'
#' @keywords models spatial prediction
#'
#' @author Pedro Tarroso and A. Marcia Barbosa
#'
#' @importFrom graphics plot
#' @importFrom stats predict
#' @export
promodel <- function(vals, method = "bioclim", na.rm = TRUE, dup.rm = FALSE,
                     verbosity = 2, ...) {
    # version 1.3 (3 Jan 2025)

    method <- tolower(method)
    method <- match.arg(method,
        choices = c(
            "bioclim",
            "domain",
            "convexhull",
            "mahalanobis",
            "kernel",
            "mvnormal"
        )
    )

    vals <- as.data.frame(vals)
    vals <- dataPrune(vals, na.rm = na.rm, dup.rm = dup.rm, verbosity = verbosity)

    model <- switch(method,
        bioclim = bioclim(vals, ...),
        domain = domain(vals),
        convexhull = convexhull(vals, ...),
        mahalanobis = mahalanobis(vals),
        kernel = kernel(vals, ...),
        mvnormal = mvnormal(vals)
    )
    model <- list(proniche = model)
    class(model) <- c("proniche", class(model))
    return(model)
}


#' Predict method for proniche models
#'
#' @param object A model object of class `proniche`.
#' @param newdata A matrix, data frame, or SpatRaster containing the predictor
#'   variables. If a SpatRaster is supplied, a SpatRaster is returned.
#'
#' @return A numeric vector or SpatRaster with predicted suitability values.
#' @export
#'
predict.proniche <- function(object, newdata = NULL) {
    if (inherits(newdata, "SpatRaster")) {
        if (is.null(object$proniche$varnames)) {
            data <- as.matrix(newdata)
        } else {
            data <- as.matrix(newdata[[object$proniche$varnames]])
        }
    } else {
        if (is.null(object$proniche$varnames)) {
            data <- as.matrix(newdata)
        } else {
            data <- newdata[, object$proniche$varnames, drop = FALSE]
        }
    }
    p <- predict(object$proniche, newdata = data)
    if (inherits(newdata, "SpatRaster")) {
        pred <- newdata[[1]]
        pred[][, 1] <- p
        names(pred) <- class(object$proniche)
        return(pred)
    }
    return(p)
}


#' Plot method for proniche models
#'
#' @param x An object of class `proniche`.
#' @param ... Additional graphical parameters passed to the underlying plot
#'   method of the fitted model.
#'
#' @return A plot of the fitted niche model in environmental space.
#' @export

plot.proniche <- function(x, ...) {
    plot(x[["proniche"]], ...)
}

#' @export
print.proniche <- function(x) {
    print(x[["proniche"]])
}
