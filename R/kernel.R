#' Kernel model
#'
#' This function fits a kernel model
#'
#' Kernel estimations can be computational intensive.
#'
#' @param x A matrix of data of environmental values at observed locations
#' @param ... Arguments to be passed to kernel (seel ?ks::kde)
#'
#' @return An object of class "kernel"
#' @export
kernel <- function(x, ...) {
    x <- na.exclude(as.matrix(x))
    if (ncol(x) > 6) {
        warning(paste0(
            "Kernel density estimate is only implemented for up to ",
            "6 variables\n(see ?ks::kde)"
        ))
    }
    k <- ks::kde(x, compute.cont = FALSE, approx.cont = TRUE, ...)
    model <- list(
        model = k,
        x = x,
        nvars = ncol(x),
        varnames = colnames(x)
    )
    class(model) <- "kernel"
    return(model)
}

#' Predict Method for Kernel
#'
#' This function predicts values using a `kernel` object.
#'
#' @param object An object of class `kernel`.
#' @param newdata New data for predictions in a form of matrix or conversible to matrix.  If NULL (default) predictions are given to training data.

#' @return A vector of predictions.
#' @export
predict.kernel <- function(object, newdata = NULL) {
    if (is.null(newdata)) {
        data <- as.matrix(object$x)
    } else {
        data <- as.matrix(newdata)
    }
    mask <- is.na(rowSums(data))
    pred <- rep(NA, nrow(data))
    pred[!mask] <- predict(object$model, x = data[!mask, ])
    return(pred)
}

#' Plot Method for Kernel
#'
#' This function plots the `kernel` model.
#'
#' @param x A model object of class `kernel`.
#' @param cols the columns indices of the data (variables) to use to plot. Only 2 values are used.
#' @param contours How many contours to plot from kernel model.
#' @param border The color of the polygon borders deppicting model extent,
#' @param pnt.col The color of the training data points (if add=FALSE).
#' @param add Boolean to add to current device or, if FALSE, to generate a new plot.
#' @param ... Other plotting parameters to be passed (lwd, lty,...).
#'
#' @export
plot.kernel <- function(x, cols = 1:2, contours = 10, border = "red",
                        pnt.col = "gray", add = FALSE, ...) {
    if (!add) {
        plot(x$x[, cols], col = pnt.col, ...)
    }
    # To plot multidimensional data, it cuts predictions to provide contours
    pred <- predict(x)
    p <- seq(min(pred), max(pred), length.out = contours + 1)
    for (i in 1:contours) {
        pnt <- x$x[which(pred >= p[i]), cols]
        ch <- chull(pnt)
        polygon(pnt[ch, ], border = border, col = NA, ...)
    }
}

#' Print Method for Kernel
#'
#' Print simple information for object of class `kernel`.
#'
#' @param x An object of class `kernel`.
#' @export
print.kernel <- function(x) {
    print(paste(class(x), "model with", x$proniche$nvars, "variables."))
}
