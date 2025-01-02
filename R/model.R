model <- function(vals, vars, method = "bioclim") {

  # version 1.1 (2 Jan 2025)

  method <- tolower(method)
  method <- match.arg(method, choices = c("bioclim", "domain", "convexhull", "mahalanobis", "kernel", "mvnormal"))

  vals <- as.data.frame(vals)  # else errors for some methods if only one variable
  if (inherits(vars, "SpatRaster")) vars <- vars[[names(vals)]]
  else vars <- vars[ , names(vals), drop = FALSE]

  if (method == "bioclim")
    bioclim(vals, vars)

  else if (method == "convexhull") {
    if (ncol(vals) < 2) stop("input 'method' requires more than one variable.")

    # otherwise:
    # Error in doTryCatch(return(expr), name, parentenv, handler) :
    #   Received error code 1 from qhull. Qhull error:
    #   QH6050 qhull error: dimension 1 must be > 1
    #
    # While executing:  | qhull Tv  Qt
    # Options selected for Qhull 2020.2.r 2020/08/31:
    #   run-id 1859778841  Tverify  Qtriangulate  _pre-merge  _zero-centrum
    # _maxoutside  0

    convexHullModel(vals, vars)
  }

  else if (method == "domain")
    domainmodel(vals, vars)

  else if (method == "mahalanobis")
    mahalanobisModel(vals, vars)

  else if (method == "kernel") {
    # if (ncol(vals) > 6)
    #   stop ("Kernel density estimate is only implemented for up to 6 variables.")
    if (ncol(vals) > 4) {
      message ("Setting binned=FALSE, as required by ks::kde() when >4 variables.\nThis is computationally intensive!")
      # otherwise:
      # Error in ks::kde(as.matrix(vals), compute.cont = FALSE, approx.cont = TRUE) : Binned estimation for d>4 not implemented. Set binned=FALSE for exact estimation.
      kernelModel(vals, vars, binned = FALSE)
    }
    else kernelModel(vals, vars)
  }

  else if (method == "mvnormal") {
    if (ncol(vals) < 2) stop("input 'method' requires more than one variable.")
    mvnormalModel(vals, vars)
  }
}
