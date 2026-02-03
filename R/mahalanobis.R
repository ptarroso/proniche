#' Mahalanobis model
#'
#' This function fits a model based on mahalanobis distance.
#' Distance values are converted to p-values using chi-squared distribution
#'
#' @param x A matrix of data of environmental values at observed locations
#'
#' @return An object of class "mahalanobis"
#' @export

mahalanobis <- function(x) {
    x <- na.exclude(as.matrix(x))
    if (nrow(x) == 0) stop ("No rows left after NA exclusion.") # (AMB added)
    u <- colMeans(x)
    sigma <- stats::cov(x) # Need always a few points to estimate covariance
    model <- list(
        model = list(u = u, sigma = sigma),
        x = x,
        nvars = ncol(x),
        varnames = colnames(x)
    )
    class(model) <- "mahalanobis"
    return(model)
}

#' Predict Method for Mahalanobis
#'
#' This function predicts values using a `mahalanobis` object.
#'
#' @param object An object of class `mahalanobis`.
#' @param newdata New data for predictions in a form of matrix or conversible to matrix.  If NULL (default) predictions are given to training data.

#' @return A vector of predictions.
#' @importFrom stats pchisq
#' @export
predict.mahalanobis <- function(object, newdata = NULL) {
    if (is.null(newdata)) {
        data <- as.matrix(object$x)
    } else {
        data <- as.matrix(newdata)
    }
    mask <- is.na(rowSums(data))
    M <- rep(NA, nrow(data))
    for (i in which(!mask)) {
        M[i] <- mah_dist(unlist(data[i, ]), object$model$u, object$model$sigma)
    }
    p <- 1 - pchisq(M, object$nvars) # (AMB edited)
    return(p)
}


#' Plot Method for mahalanobis
#'
#' This function plots the `mahalanobis` model.
#'
#' @param x An object of class `mahalanobis`.
#' @param cols the columns indices of the data (variables) to use to plot. Only 2 values are used.
#' @param contours A vector of contours to plot defined as standard deviations away from the mean value.
#' @param border The color of the polygon borders deppicting model extent,
#' @param pnt.col The color of the training data points (if add=FALSE).
#' @param add Boolean to add to current device or, if FALSE, to generate a new plot.
#' @param ... Other plotting parameters to be passed (lwd, lty,...).
#'
#' @export
plot.mahalanobis <- function(x, cols = 1:2,
                             contours = c(1, 1.64, 1.96, 2.33),
                             border = "red", pnt.col = "gray", add = FALSE,
                             ...) {
    if (!add) {
        plot(x$x[, cols], col = pnt.col, ...)
    }
    for (sdev in contours) {
        pnt <- ellipse_from_cov(x$model$sigma, x$model$u, sdev, cols)
        polygon(pnt[, 1], pnt[, 2], border = border, col = NA, ...)
    }
}

#' Print Method for mahalanobis
#'
#' Print simple information for object of class `mahalanobis`.
#'
#' @param x An object of class `mahalanobis`.
#' @export
print.mahalanobis <- function(x) {
    print(paste(class(x), "model with", x$proniche$nvars, "variables."))
}
