# Returns Gower's distance

#' @importFrom stats na.omit

gower <- function(x, train, rng = apply(train, 2, function(x) diff(range(x)))) {
    p <- ncol(train)
    G <- 0
    for (k in 1:p) {
        d <- abs(x[k] - train[, k]) / rng[k]
        G <- G + d
    }
    G / p
}

# points at a Gower's distance value of training data
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
