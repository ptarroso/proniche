# all functions (except the last one, dataPrune) by Pedro Tarroso
# some edited by A. Marcia Barbosa where indicated, to avoid errors or to accommodate dataframe (not only SpatRaster) vars
# not exported; called by wrapper function models()

bioclim_model <- function(x, nq = 10) {
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
        type = "bioclim", model = env_rect, x = x, nq = nq,
        nvars = ncol(x)
    )
    class(model) <- "proniche"
    return(model)
}

# Vars is a matrix. predict should be aware of format and convert
bioclim_predict <- function(model, newdata = NULL) {
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

bioclim_plot <- function(model, cols = 1:2, border = "red",
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

convexhull_model <- function(x) {
    x <- na.exclude(as.matrix(x))
    original <- x
    nvars <- ncol(x)
    env_ch <- list()
    i <- 1
    id <- 1:nrow(x)
    while (nrow(x) > nvars) {
        ch <- geometry::convhulln(x)
        # Have to redo the IDs in ch to map to the original matrix
        attr(ch, "pid") <- matrix(id[ch], nrow(ch), ncol(ch))
        rm <- unique(as.vector(ch))
        x <- x[-rm, , drop = FALSE] # (AMB edited)
        id <- id[-rm]
        env_ch[[i]] <- ch
        i <- i + 1
    }
    model <- list(
        type = "convex_hull",
        model = env_ch,
        x = original,
        nch = length(env_ch),
        nvars = nvars
    )
    class(model) <- "proniche"
    return(model)
}

convexhull_predict <- function(model, newdata = NULL) {
    if (is.null(newdata)) {
        data <- as.matrix(model$x)
    } else {
        data <- as.matrix(newdata)
    }
    pred <- rep(0, nrow(data))
    for (i in 1:model$nch) {
        pred <- pred + geometry::inhulln(model$model[[i]], data)
    }
    return(pred)
}

convexhull_plot <- function(model, cols = 1:2, border = "red",
                            pnt.col = "gray", add = FALSE, ...) {
    if (!add) {
        plot(model$x[, cols], col = pnt.col, ...)
    }
    # The representation of high dimensional chulls is difficult...
    # The trick used here is to check which points are within each
    # chull and generate a 2D chull of those on the user defined dims
    pred <- convexhull_predict(model)
    for (i in 1:model$nch) {
        pnt <- model$x[which(pred >= i), cols]
        ch <- chull(pnt)
        polygon(pnt[ch, ], border = border, col = NA, ...)
    }
}


# Probably Only works with max of 6 variables
kernel_model <- function(x, ...) {
    x <- na.exclude(as.matrix(x))
    k <- ks::kde(x, compute.cont = FALSE, approx.cont = TRUE, ...)
    model <- list(
        type = "kernel",
        model = k,
        x = x,
        nvars = ncol(x)
    )
    class(model) <- "proniche"
    return(model)
}

kernel_predict <- function(model, newdata = NULL) {
    if (is.null(newdata)) {
        data <- as.matrix(model$x)
    } else {
        data <- as.matrix(newdata)
    }
    data <- na.exclude(data)
    pred <- rep(0, nrow(data))
    pred <- predict(model$model, x = data)
    return(pred)
}

kernel_plot <- function(model, cols = 1:2, contours = 10, border = "red",
                        pnt.col = "gray", add = FALSE, ...) {
    if (!add) {
        plot(model$x[, cols], col = pnt.col, ...)
    }
    # To plot multidimensional data, it cuts predictions to provide contours
    pred <- kernel_predict(model)
    p <- seq(min(pred), max(pred), length.out = contours + 1)
    for (i in 1:contours) {
        pnt <- model$x[which(pred >= p[i]), cols]
        ch <- chull(pnt)
        polygon(pnt[ch, ], border = border, col = NA, ...)
    }
}

# Returns Gower's distance
gower <- function(x, train, rng = apply(train, 2, function(x) diff(range(x)))) {
    p <- ncol(train)
    G <- 0
    for (k in 1:p) {
        d <- abs(x[k] - train[, k]) / rng[k]
        G <- G + d
    }
    G / p
}

igower <- function(model, D) {
    p <- model$nvars
    rng <- model$model
    target <- p * D
    points <- list()

    for (i in 1:nrow(model$x)) {
        pnt <- model$x[i, ]

        points[[i]] <- expand.grid(
            lapply(1:p, function(k) {
                delta <- target * rng[k]
                c(pnt[k] - delta, pnt[k] + delta)
            })
        )
    }
    return(points)
}

# Based on Carpenter et al. (1993). Biodiversity Conservation 2, 667-680
domain_model <- function(x) {
    x <- na.exclude(as.matrix(x))
    rng <- unlist(apply(x, 2, function(d) diff(range(d))))
    model <- list(
        type = "domain",
        model = rng,
        x = x,
        nvars = ncol(x)
    )
    class(model) <- "proniche"
    return(model)
}

domain_predict <- function(model, newdata = NULL) {
    if (is.null(newdata)) {
        data <- as.matrix(model$x)
    } else {
        data <- na.exclude(as.matrix(newdata))
    }
    D <- rep(NA, nrow(data))
    for (i in 1:nrow(data)) {
        G <- gower(unlist(data[i, ]), model$x, model$model)
        G[which(G > 1)] <- 1
        D[i] <- max(1 - G, na.rm = T)
    }
    return(D)
}

domain_plot <- function(model, cols = 1:2, contours = seq(0.9, 1, 0.01),
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


# mahalanobis distance
mah_dist <- function(x, u, sigma) {
    sqrt(t(x - u) %*% solve(sigma) %*% (x - u))
}

ellipse_from_cov <- function(sigma, u, sdev_scale, dims = 1:2, n_points = 100) {
    eig <- eigen(sigma[dims, dims])
    vec <- eig$vectors
    val <- sqrt(eig$values)
    angles <- seq(0, 2 * pi, length.out = n_points)
    points <- cbind(cos(angles), sin(angles))
    points <- (points %*% (diag(val) * sdev_scale)) %*% t(vec)
    points <- t(t(points) + u[dims])
    return(points)
}


mahalanobis_model <- function(x) {
    x <- na.exclude(as.matrix(x))
    u <- colMeans(x)
    sigma <- stats::cov(x) # Need always a few points to estimate covaraince
    model <- list(
        type = "mahalanobis",
        model = list(u = u, sigma = sigma),
        x = x,
        nvars = ncol(x)
    )
    class(model) <- "proniche"
    return(model)
}

# Note: scales distances to zero-one range and inverts (near is 1, far is zero)
mahalanobis_predict <- function(model, newdata = NULL) {
    if (is.null(newdata)) {
        data <- as.matrix(model$x)
    } else {
        data <- na.exclude(as.matrix(newdata))
    }
    M <- rep(NA, nrow(data))
    for (i in 1:length(M)) {
        M[i] <- mah_dist(unlist(data[i, ]), model$model$u, model$model$sigma)
    }
    p <- 1 - (M / max(M)) # (AMB edited)
    return(p)
}


# contours are sdevs away from the mean
mahalanobis_plot <- function(model, cols = 1:2,
                             contours = c(1, 1.64, 1.96, 2.33),
                             border = "red", pnt.col = "gray", add = FALSE,
                             ...) {
    if (!add) {
        plot(model$x[, cols], col = pnt.col, ...)
    }
    for (sdev in contours) {
        pnt <- ellipse_from_cov(model$model$sigma, model$model$u, sdev, cols)
        polygon(pnt[, 1], pnt[, 2], border = border, col = NA, ...)
    }
}


# Multivariate normal distribution
dmnorm <- function(x, mu, sigma) {
    # mu is a vector of averages per variable
    # sigma is the covariance matrix
    if (is.null(dim(mu))) {
        mu <- matrix(mu)
    }
    n <- ncol(x)
    if (n != nrow(mu) | n != ncol(sigma)) {
        stop("Dimensions do not match!")
    }
    isigma <- solve(sigma)
    den <- ((2 * pi)**(n / 2) * sqrt(det(sigma)))
    centered <- sweep(x, 2, mu)
    dens <- apply(centered, 1, function(x) {
        exp(-(1 / 2) * t(x) %*% isigma %*% x)
    })
    return(dens / den)
}

mvnormal_model <- function(x) {
    data <- na.exclude(as.matrix(x))
    # Stop if not enough data points
    if (nrow(data) < 5) {
        stop("Low number of samples. Consider resampling with envResample.")
    }

    # estimate mean value and covmat (sig)
    avg <- colMeans(x, na.rm = T)
    sig <- stats::cov(x)

    model <- list(
        type = "mvnormal",
        model = list(u = avg, sigma = sig),
        x = x,
        nvars = ncol(x)
    )
    class(model) <- "proniche"
    return(model)
}


mvnormal_predict <- function(model, newdata = NULL) {
    if (is.null(newdata)) {
        data <- as.matrix(model$x)
    } else {
        data <- na.exclude(as.matrix(newdata))
    }
    avg <- model$model$u
    sig <- model$model$sigma

    p <- dmnorm(data, avg, sig)
    return(p)
}

# contours are sdevs away from the mean
mvnormal_plot <- function(model, cols = 1:2,
                          contours = c(1, 1.64, 1.96, 2.33),
                          border = "red", pnt.col = "gray", add = FALSE,
                          ...) {
    if (!add) {
        plot(model$x[, cols], col = pnt.col, ...)
    }
    for (sdev in contours) {
        pnt <- ellipse_from_cov(model$model$sigma, model$model$u, sdev, cols)
        polygon(pnt[, 1], pnt[, 2], border = border, col = NA, ...)
    }
}

# Adds samples to data by the nearest "maxsample" points to
# each original sample in the environmental space
# Note: might have non-unique environmental samples!
envResample <- function(vals, vars, maxsample = 50) {
    data <- t(data.frame(vars))
    ns <- round(maxsample / nrow(vals), 0)
    for (i in 1:nrow(vals)) {
        d <- colSums((data - unlist(vals[i, ]))**2)
        data.new <- t(data[, order(d)[1:ns]])
        colnames(data.new) <- colnames(vals)
        vals <- rbind(vals, data.new)
    }
    vals
}


dataPrune <- function(vals, na.rm = TRUE, dup.rm = FALSE, verbosity = 2) {
    if (isTRUE(na.rm)) {
        n_in <- nrow(vals)
        vals <- na.omit(vals)
        n_out <- nrow(vals)
        if (verbosity > 1 && n_in > n_out) {
            message(n_in - n_out, " rows removed due to missing data.")
        }
    }

    if (isTRUE(dup.rm)) {
        n_in <- nrow(vals)
        vals <- unique(vals)
        n_out <- nrow(vals)
        if (verbosity > 1 && n_in > n_out) {
            message(n_in - n_out, " rows removed due to duplicates.")
        }
    }

    vals
}
