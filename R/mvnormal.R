#' Multivariate Normal Distribution model
#'
#' This function fits a Multivariate Normal Distribution model
#'
#' @param x A matrix of data of environmental values at observed locations
#'
#' @return An object of class "mvnormal"
#' @export
mvnormal <- function(x) {
    data <- na.exclude(as.matrix(x))
    # Stop if not enough data points
    if (nrow(data) < 5) {
        stop("Low number of samples. Consider resampling with envResample.")
    }
    if (ncol(data) < 2) {
        stop("input 'method' requires more than one variable.")
    }
    # estimate mean value and covmat (sig)
    avg <- colMeans(x, na.rm = T)
    sig <- stats::cov(x)

    model <- list(
        model = list(u = avg, sigma = sig),
        x = x,
        nvars = ncol(x),
        varnames = colnames(x)
    )
    class(model) <- "mvnormal"
    return(model)
}

#' Predict Method for Multivariate Normal Distribution model
#'
#' This function predicts values using a `mvnormal` object.
#'
#' @param object An object of class `mvnormal`.
#' @param newdata New data for predictions in a form of matrix or conversible to matrix.  If NULL (default) predictions are given to training data.

#' @return A vector of predictions.
#' @export
predict.mvnormal <- function(object, newdata = NULL) {
    if (is.null(newdata)) {
        data <- as.matrix(object$x)
    } else {
        data <- as.matrix(newdata)
    }
    avg <- object$model$u
    sig <- object$model$sigma

    p <- dmnorm(data, avg, sig)
    return(p)
}

#' Plot Method for Multivariate Normal Distribution model
#'
#' This function plots the `mvnormal` model.
#'
#' @param x An object of class `mvnormal`.
#' @param cols the columns indices of the data (variables) to use to plot. Only 2 values are used.
#' @param contours A vector of contours to plot defined as standard deviations away from the mean value.
#' @param border The color of the polygon borders deppicting model extent,
#' @param pnt.col The color of the training data points (if add=FALSE).
#' @param add Boolean to add to current device or, if FALSE, to generate a new plot.
#' @param ... Other plotting parameters to be passed (lwd, lty,...).
#'
#' @export
plot.mvnormal <- function(x, cols = 1:2,
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

#' Print Method for Multivariate Normal Distribution model
#'
#' Print simple information for object of class `mvnormal`.
#'
#' @param x An object of class `mvnormal`.
#' @export
print.mvnormal <- function(x) {
    print(paste(class(x), "model with", x$proniche$nvars, "variables."))
}
