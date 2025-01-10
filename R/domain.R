#' Domain model
#'
#' This function fits a  domain model based on Gower's distance.
#'
#' Distances D are converted to similarity index by the complement (1-D).
#' Based on Carpenter et al. (1993). Biodiversity Conservation 2, 667-680
#'
#' @param x A matrix of data of environmental values at observed locations
#'
#' @return An object of class "domain"
#' @export

domain <- function(x) {
    x <- na.exclude(as.matrix(x))
    rng <- unlist(apply(x, 2, function(d) diff(range(d))))
    model <- list(
        model = rng,
        x = x,
        nvars = ncol(x),
        varnames = colnames(x)
    )
    class(model) <- "domain"
    return(model)
}

#' Predict Method for Domain
#'
#' This function predicts values using a `domain` object.
#'
#' @param model An object of class `domain`.
#' @param newdata New data for predictions in a form of matrix or conversible to matrix.  If NULL (default) predictions are given to training data.

#' @return A vector of predictions.
#' @export
predict.domain <- function(model, newdata = NULL) {
    if (is.null(newdata)) {
        data <- as.matrix(model$x)
    } else {
        data <- as.matrix(newdata)
    }
    D <- rep(NA, nrow(data))
    for (i in 1:nrow(data)) {
        if (!any(is.na(data[i, ]))) {
            G <- gower(unlist(data[i, ]), model$x, model$model)
            G[which(G > 1)] <- 1
            D[i] <- max(1 - G, na.rm = T)
        }
    }
    return(D)
}

#' Plot Method for Domain
#'
#' This function plots the `domain` model.
#'
#' @param model An object of class `domain`.
#' @param cols the columns indices of the data (variables) to use to plot. Only 2 values are used.
#' @param contouns Where to plot the extents on Gower's similarity values.
#' @param border The color of the polygon borders deppicting model extent,
#' @param pnt.col The color of the training data points (if add=FALSE).
#' @param add Boolean to add to current device or, if FALSE, to generate a new plot.
#' @param ... Other plotting parameters to be passed (lwd, lty,...).
#'
#' @export
plot.domain <- function(model, cols = 1:2, contours = seq(0.9, 1, 0.01),
                        border = "red", pnt.col = "gray", add = FALSE, ...) {
    if (!add) {
        plot(model$x[, cols], col = pnt.col, ...)
    }
    if (length(contours) == 1) {
        contours <- seq(0, 1, length.out = contours + 2)
        contours <- contours[-c(1, length(contours))]
    }
    for (d in contours) {
        # Using terra to retrive boundaries by aggregation
        bounds <- terra::vect()
        points <- igower(model, 1 - d)
        for (i in 1:length(points)) {
            rng1 <- range(points[[i]][, cols[1]])
            rng2 <- range(points[[i]][, cols[2]])
            e <- terra::ext(rng1[1], rng1[2], rng2[1], rng2[2])
            bounds <- rbind(bounds, terra::vect(e))
        }
        bounds <- terra::aggregate(bounds)
        terra::lines(bounds, col = border, ...)
    }
}

#' Print Method for Domain
#'
#' Print simple information for object of class `domain`.
#'
#' @param model An object of class `domain`.
#' @export
print.domain <- function(model) {
    print(paste(class(model), "model with", model$nq, "variables."))
}
