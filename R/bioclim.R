#' Bioclim model
#'
#' This function fits a bioclim envelope model
#'
#' @param x A matrix of environmental values at presence locations
#' @param nq Number of quantiles to estimate bioclim envelopes.
#'
#' @return An object of class "bioclim"
#' @importFrom stats quantile
#' @importFrom graphics plot rect
#' @export
#'
bioclim <- function(x, nq = 10) {
    qt <- seq(0, 0.5, length.out = nq)
    env_rect <- matrix(NA, nq, ncol(x) * 2 + 1) # (AMB edited)
    env_rect[, 1] <- qt
    colnames(env_rect) <- c(
        "quantile",
        paste0(
            rep(colnames(x), each = 2),
            c("_min", "_max")
        )
    )
    for (i in 1:nq) {
        rng <- apply(x, 2, quantile, probs = c(qt[i], 1 - qt[i]))
        env_rect[i, -1] <- as.vector(rng)
    }
    model <- list(
        model = env_rect,
        x = x,
        nq = nq,
        nvars = ncol(x),
        varnames = colnames(x)
    )
    class(model) <- "bioclim"
    return(model)
}

#' Predict Method for Bioclim
#'
#' This function predicts values using a `bioclim` model object.
#'
#' @param object An object of class `bioclim`.
#' @param newdata New data for predictions, as a matrix or coercible to matrix. If NULL (default), predictions are given on the training data.

#' @return A vector of predictions.
#' @export
predict.bioclim <- function(object, newdata = NULL) {
    if (is.null(newdata)) {
        data <- object$x
    } else {
        data <- as.matrix(newdata)
    }
    pred <- rep(0, nrow(data))
    nvars <- object$nvars
    nq <- object$nq
    for (i in 1:nq) {
        tmp <- pred * 0
        rng <- matrix(object$model[i, -1], 2, nvars)
        for (j in 1:nvars) {
            tmp <- tmp + (data[, j] >= rng[1, j] & data[, j] < rng[2, j])
        }
        pred <- pred + (tmp == nvars)
    }
    return(pred)
}


#' Plot Method for Bioclim
#'
#' This function plots the `bioclim` model.
#'
#' @param x An object of class `bioclim`.
#' @param cols the columns indices of the data (variables) to use to plot. Only 2 values are used.
#' @param border The color of the polygon borders deppicting model extent,
#' @param pnt.col The color of the training data points (if add=FALSE).
#' @param add Boolean to add to current device or, if FALSE, to generate a new plot.
#' @param ... Other plotting parameters to be passed (lwd, lty,...).
#'
#' @export
plot.bioclim <- function(x, cols = 1:2, border = "red",
                         pnt.col = "gray", add = FALSE, ...) {
    if (!add) {
        plot(x$x[, cols], col = pnt.col, ...)
    }
    for (i in 1:x$nq) {
        rng <- matrix(x$model[i, -1], 2, x$nvars)
        rect(rng[1, cols[1]],
            rng[1, cols[2]],
            rng[2, cols[1]],
            rng[2, cols[2]],
            border = border, ...
        )
    }
}

#' Print Method for Bioclim
#'
#' Print simple information for object of class `bioclim`.
#'
#' @param x An object of class `bioclim`.
#' @export
print.bioclim <- function(x) {
    print(paste(class(x), "model with", x$nq, "variables."))
}
