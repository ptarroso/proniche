#' @export
promodel <- function(x, method = "bioclim", na.rm = TRUE, dup.rm = FALSE,
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

    x <- as.data.frame(x)
    x <- dataPrune(x, na.rm = na.rm, dup.rm = dup.rm, verbosity = verbosity)

    model <- switch(method,
        bioclim = bioclim(x, ...),
        domain = domain(x),
        convexhull = convexhull(x, ...),
        mahalanobis = mahalanobis(x),
        kernel = kernel(x, ...),
        mvnormal = mvnormal(x)
    )
    model <- list(proniche = model)
    class(model) <- c("proniche", class(model))
    return(model)
}

#' @export
predict.proniche <- function(model, newdata = NULL) {
    if (inherits(newdata, "SpatRaster")) {
        if (is.null(model$proniche$varnames)) {
            data <- as.matrix(newdata)
        } else {
            data <- as.matrix(newdata[[model$proniche$varnames]])
        }
    } else {
        if (is.null(model$proniche$varnames)) {
            data <- as.matrix(newdata)
        } else {
            data <- newdata[, model$proniche$varnames, drop = FALSE]
        }
    }
    p <- predict(model$proniche, newdata = data)
    if (inherits(newdata, "SpatRaster")) {
        pred <- newdata[[1]]
        pred[][, 1] <- p
        names(pred) <- class(model$proniche)
        return(pred)
    }
    return(p)
}

#' @export
plot.proniche <- function(model, ...) {
    plot(model[["proniche"]], ...)
}

#' @export
print.proniche <- function(model) {
    print(model[["proniche"]])
}
