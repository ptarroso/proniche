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
        convexhull = convexhull(x),
        mahalanobis = mahalanobis(x),
        kernel = kernel(x, ...),
        mvnormal = mvnormal(x)
    )
    model <- list(proniche = model)
    class(model) <- c("proniche", class(model))
    return(model)
}

predict.proniche <- function(model, newdata = NULL) {
    if (inherits(newdata, "SpatRaster")) {
        data <- as.matrix(newdata[[names(model$proniche$x)]])
    } else {
        data <- newdata[, names(model$proniche$x), drop = FALSE]
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

plot.proniche <- function(model, ...) {
    plot(model[["proniche"]], ...)
}

print.proniche <- function(model) {
    print(model[["proniche"]])
}
