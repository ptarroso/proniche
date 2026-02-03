#' Bioclim model
#'
#' This function fits a bioclim model
#'
#' @param x A matrix of data of environmental values at observed locations
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
#' This function predicts values using a `bioclim` object.
#'
#' @param model An object of class `bioclim`.
#' @param newdata New data for predictions in a form of matrix or conversible to matrix.  If NULL (default) predictions are given to training data.

#' @return A vector of predictions.
#' @export
predict.bioclim <- function(model, newdata = NULL) {
    if (is.null(newdata)) {
        data <- model$x
    } else {
        data <- as.matrix(newdata)
    }
    pred <- rep(0, nrow(data))
    nvars <- model$nvars
    nq <- model$nq
    for (i in 1:nq) {
        tmp <- pred * 0
        rng <- matrix(model$model[i, -1], 2, nvars)
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
#' @param model An object of class `bioclim`.
#' @param cols the columns indices of the data (variables) to use to plot. Only 2 values are used.
#' @param border The color of the polygon borders deppicting model extent,
#' @param pnt.col The color of the training data points (if add=FALSE).
#' @param add Boolean to add to current device or, if FALSE, to generate a new plot.
#' @param ... Other plotting parameters to be passed (lwd, lty,...).
#'
#' @export
plot.bioclim <- function(model, cols = 1:2, border = "red",
                         pnt.col = "gray", add = FALSE, ...) {
    if (!add) {
        plot(model$x[, cols], col = pnt.col, ...)
    }
    for (i in 1:model$nq) {
        rng <- matrix(model$model[i, -1], 2, model$nvars)
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
#' @param model An object of class `bioclim`.
#' @export
print.bioclim <- function(model) {
    print(paste(class(model), "model with", model$nq, "variables."))
}
