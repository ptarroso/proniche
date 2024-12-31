model <- function(vals, vars, method = "bioclim") {

  method <- match.arg(method, choices = c("bioclim", "convexhull", "kernel", "domain", "mahalanobis", "mvnormal"))

  if (method == "bioclim")
    bioclim(vals, vars)

  else if (method == "convexhull")
    convexHullModel(vals, vars)

  else if (method == "kernel")
    kernelModel(vals, vars)

  else if (method == "domain")
    domainmodel(vals, vars)

  else if (method == "mahalanobis")
    mahalanobisModel(vals, vars)

  else if (method == "mvnormal")
    mvnormalModel(vals, vars)
}
