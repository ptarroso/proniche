bioclim <- function(vals, vars, nq=10) {
    qt <- seq(0, 0.5, length.out=nq)
    env_rect <- matrix(NA, nq, nlyr(vars)*2+1)
    env_rect[,1] <- qt
    colnames(env_rect) <- c("quantile", paste0(rep(names(vars), each=2), c("_min", "_max")))
    pred <- vars[[1]]*0
    for (i in 1:nq) {
        rng <- apply(vals, 2, quantile, probs=c(qt[i], 1-qt[i]))
        env_rect[i,-1] <- as.vector(rng)
        tmp <- pred * 0
		for (j in 1:nlyr(vars)) {
		    tmp <- tmp + (vars[[j]] >= rng[1,j] & vars[[j]] < rng[2,j])
		}
		pred <- pred + (tmp == nlyr(vars))
	}
	list(pred, env_rect)
}


convexHullModel <- function(vals, vars) {
    vals <- as.matrix(vals)
    pred <- vars[[1]]*0
    data <- as.matrix(vars)
    env_ch <- list()
    i <- 1
    while (nrow(vals) > 4) {
        ch <- geometry::convhulln(vals)
        values(pred) <- values(pred) + geometry::inhulln(ch, data)
        vals <- vals[-unique(ch),]
        env_ch[[i]] <- ch
        i <- i + 1
    }
    list(pred, env_ch)
}

# Probably Only works with max of 6 variables
kernelModel <- function(vals, vars) {
    k <- ks::kde(vals, compute.cont=FALSE, approx.cont=TRUE)
    x <- as.data.frame(vars, na.rm=F)
    mask <- is.na(rowSums(x))
    pred <- vars[[1]]
    pred[!mask] <- predict(k, x=x[!mask,])
    list(pred, k)
}

gower <- function(x, train, rng = apply(train, 2, function(x) diff(range(x)))) {
    p <- ncol(train)
    for (k in 1:p) {
        d <- abs(x[k] - train[,k])/rng[k]
        if (k == 1)
            G <- d
        else
            G <- G + d
    }
    G/p
}

domainmodel <- function(vals, vars) {
    data <- as.data.frame(vars, na.rm=F)
    D <- rep(NA, nrow(data))
    mask <- !is.na(rowSums(data))
    rng <- unlist(apply(vals, 2, function(d) diff(range(d))))
    for (i in which(mask)) {
        G <- gower(unlist(data[i,]), vals, rng)
        G[which(G>1)] <- 1
        D[i] <- max(1-G, na.rm=T)
    }
    pred <- vars[[1]] * NA
    values(pred) <- D
    list(pred)
}

# mahalanobis distance
mah.dist <- function(x, u, sigma) {
    sqrt(t(x-u)%*%solve(sigma)%*%(x-u))
}

#Note: scales distances to zero-one range and inverts (near is 1, far is zero)
mahalanobisModel <- function(vals, vars) {
    data <- as.data.frame(vars, na.rm=F)
    M <- rep(NA, nrow(data))
    mask <- !is.na(rowSums(data))
    u <- colMeans(vals)
    sigma <- stats::cov(vals) # Need always a few points to estimate covaraince
    for (i in which(mask)) {
        M[i] <- mah.dist(unlist(data[i,]), u, sigma)
    }
    pred <- vars[[1]] * NA
    values(pred) <- 1 - (M/max(M, na.rm=T))
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
    den <- ((2*pi)**(n/2) * sqrt(det(sigma)))
    centered <- sweep(x, 2, mu)
    dens <- apply(centered, 1, function(x){exp(-(1/2) * t(x) %*% isigma %*% x)})
    return(dens/den)
}

mvnormalModel <- function(vals, vars) {
    # Stop if not enough data points
    if (nrow(vals) < 5) {
        stop("Very low number of samples. Consider resampling first with envResample.")
    }
    data <- as.data.frame(vars, na.rm=F)
    mask <- !is.na(rowSums(data))

    # estimate mean value and covmat (sig)
    avg <- colMeans(vals, na.rm=T)
    sig <- stats::cov(vals)

    pred <- vars[[1]] * NA
    values(pred)[mask,] <- dmnorm(data[mask,], avg, sig)
    list(pred, list(avg, sig))
}

# Adds samples to data by the nearest "maxsample" points to
# each original sample in the environmental space
# Note: might have non-unique environmental samples!
envResample <- function(vals, vars, maxsample=50) {
    data <- t(data.frame(vars))
    ns <- round(maxsample/nrow(vals), 0)
    for (i in 1:nrow(vals)) {
        d <- colSums((data - unlist(vals[i,]))**2)
        data.new <- t(data[,order(d)[1:ns]])
        colnames(data.new) <- colnames(vals)
        vals <- rbind(vals, data.new)
    }
    vals
}
