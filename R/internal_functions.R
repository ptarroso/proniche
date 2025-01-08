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
                            pnt.col = "gray", add = TRUE, ...) {
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
kernelModel <- function(vals, vars, ...) {
    # k <- ks::kde(vals, compute.cont=FALSE, approx.cont=TRUE)
    k <- ks::kde(as.matrix(vals), compute.cont = FALSE, approx.cont = TRUE, ...) # (AMB edited)
    x <- as.data.frame(vars, na.rm = F)
    mask <- is.na(rowSums(x))
    pred <- vars[[1]]
    pred[!mask] <- predict(k, x = x[!mask, ])
    list(pred, k)
}

gower <- function(x, train, rng = apply(train, 2, function(x) diff(range(x)))) {
    p <- ncol(train)
    for (k in 1:p) {
        d <- abs(x[k] - train[, k]) / rng[k]
        if (k == 1) {
            G <- d
        } else {
            G <- G + d
        }
    }
    G / p
}

domainmodel <- function(vals, vars) {
    data <- as.data.frame(vars, na.rm = F)
    D <- rep(NA, nrow(data))
    mask <- !is.na(rowSums(data))
    rng <- unlist(apply(vals, 2, function(d) diff(range(d))))
    for (i in which(mask)) {
        G <- gower(unlist(data[i, ]), vals, rng)
        G[which(G > 1)] <- 1
        D[i] <- max(1 - G, na.rm = T)
    }
    pred <- vars[[1]] * NA
    if (inherits(pred, "SpatRaster")) { # (AMB added)
        terra::values(pred) <- D
    } else {
        pred <- D
    } # (AMB edited)
    list(pred)
}

# mahalanobis distance
mah.dist <- function(x, u, sigma) {
    sqrt(t(x - u) %*% solve(sigma) %*% (x - u))
}

# Note: scales distances to zero-one range and inverts (near is 1, far is zero)
mahalanobisModel <- function(vals, vars) {
    data <- as.data.frame(vars, na.rm = F)
    M <- rep(NA, nrow(data))
    mask <- !is.na(rowSums(data))
    u <- colMeans(vals)
    sigma <- stats::cov(vals) # Need always a few points to estimate covaraince
    for (i in which(mask)) {
        M[i] <- mah.dist(unlist(data[i, ]), u, sigma)
    }
    # values(pred) <- 1 - (M/max(M, na.rm=T))
    p <- 1 - (M / max(M, na.rm = T)) # (AMB edited)
    if (inherits(vars, "SpatRaster")) { # (AMB added)
        pred <- vars[[1]] * NA
        terra::values(pred) <- p # (AMB edited)
    } else {
        pred <- p
    } # (AMB edited)
    list(pred)
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

mvnormalModel <- function(vals, vars) {
    # Stop if not enough data points
    if (nrow(vals) < 5) {
        stop("Very low number of samples. Consider resampling first with envResample.")
    }
    data <- as.data.frame(vars, na.rm = F)
    mask <- !is.na(rowSums(data))

    # estimate mean value and covmat (sig)
    avg <- colMeans(vals, na.rm = T)
    sig <- stats::cov(vals)

    # values(pred)[mask,] <- dmnorm(data[mask,], avg, sig)
    p <- dmnorm(data[mask, ], avg, sig)
    if (inherits(vars, "SpatRaster")) { # (AMB added)
        pred <- vars[[1]] * NA
        terra::values(pred)[mask, ] <- p # (AMB edited)
    } else {
        pred <- p
    } # (AMB added)
    list(pred, list(avg, sig))
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
            message(n_in - n_out, " rows removed due to missing data in 'vals'.")
        }
    }

    if (isTRUE(dup.rm)) {
        n_in <- nrow(vals)
        vals <- unique(vals)
        n_out <- nrow(vals)
        if (verbosity > 1 && n_in > n_out) {
            message(n_in - n_out, " rows removed due to duplicate 'vals'.")
        }
    }

    vals
}
