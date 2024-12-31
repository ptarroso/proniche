model <- function(vals, vars, method = "bioclim") {

  # version 1.0 (31 Dec 2024)

  method <- match.arg(method, choices = c("bioclim", "domain", "convexhull", "mahalanobis", "kernel", "mvnormal"))

  vals <- as.data.frame(vals)  # else errors for some methods if only one variable

  if (method == "bioclim")
    bioclim(vals, vars)

  else if (method == "convexhull") {
    if (ncol(vals) < 2) stop("'convexhull' method requires at least 2 variables.")

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

  else if (method == "kernel")
    kernelModel(vals, vars)

  else if (method == "mvnormal")
    mvnormalModel(vals, vars)
}
