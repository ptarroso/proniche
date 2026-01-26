#' Convex hull model
#'
#' This function fits a convex hull model
#'
#' @param x A matrix of data of environmental values at observed locations
#'
#' @return An object of class "convexhull"
#' @export
convexhull <- function(x) {
    x <- na.exclude(as.matrix(x))
    original <- x
    nvars <- ncol(x)
    if (nvars < 2) stop("Convex hull requires more than one variable.")
    env_ch <- list()
    i <- 1
    id <- 1:nrow(x)
    while (nrow(x) > nvars) {
        ch <- geometry::convhulln(x, options = "QJ") # (AMB edited)
        # Have to redo the IDs in ch to map to the original matrix
        attr(ch, "pid") <- matrix(id[ch], nrow(ch), ncol(ch))
        rm <- unique(as.vector(ch))
        x <- x[-rm, , drop = FALSE] # (AMB edited)
        id <- id[-rm]
        env_ch[[i]] <- ch
        i <- i + 1
    }
    model <- list(
        model = env_ch,
        x = original,
        nch = length(env_ch),
        nvars = nvars,
        varnames = colnames(x)
    )
    class(model) <- "convexhull"
    return(model)
}

#' Predict Method for Convex Hull
#'
#' This function predicts values using a `convexhull` object.
#'
#' @param model An object of class `convexhull`.
#' @param newdata New data for predictions in a form of matrix or conversible to matrix.  If NULL (default) predictions are given to training data.

#' @return A vector of predictions.
#' @export
predict.convexhull <- function(model, newdata = NULL) {
    if (is.null(newdata)) {
        data <- as.matrix(model$x)
        mask <- NULL
    } else {
        data <- as.matrix(newdata)
        mask <- is.na(rowSums(data))
    }
    pred <- rep(0, nrow(data))
    for (i in 1:model$nch) {
        pred <- pred + geometry::inhulln(model$model[[i]], data)
    }
    pred[mask] <- NA
    return(pred)
}

#' Plot Method for Convex Hull
#'
#' This function plots the `convexhull` model.
#'
#' @param model An object of class `convexhull`.
#' @param cols the columns indices of the data (variables) to use to plot. Only 2 values are used.
#' @param border The color of the polygon borders deppicting model extent,
#' @param pnt.col The color of the training data points (if add=FALSE).
#' @param add Boolean to add to current device or, if FALSE, to generate a new plot.
#' @param ... Other plotting parameters to be passed (lwd, lty,...).
#'
#' @export
plot.convexhull <- function(model, cols = 1:2, border = "red",
                            pnt.col = "gray", add = FALSE, ...) {
    if (!add) {
        plot(model$x[, cols], col = pnt.col, ...)
    }
    # The representation of high dimensional chulls is difficult...
    # The trick used here is to check which points are within each
    # chull and generate a 2D chull of those on the user defined dims
    pred <- predict(model)
    for (i in 1:model$nch) {
        pnt <- model$x[which(pred >= i), cols]
        ch <- chull(pnt)
        polygon(pnt[ch, ], border = border, col = NA, ...)
    }
}

#' Print Method for Convex Hull
#'
#' Print simple information for object of class `convexhull`.
#'
#' @param model An object of class `convexhull`.
#' @export
print.convexhull <- function(model) {
    print(paste(class(model), "model with", model$nvars, "variables."))
}
