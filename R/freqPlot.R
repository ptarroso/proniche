freqPlot <- function(vals, vars = NULL, type = "density", na.rm = TRUE, dup.rm = FALSE, verbosity = 2, ...) {

  type <- match.arg(tolower(type), choices = c("density", "boxplot"))

  vals <- as.data.frame(vals)

  vals <- dataPrune(vals, na.rm = na.rm, dup.rm = dup.rm, verbosity = verbosity)

  if (!is.null(vars)) {
    if (inherits(vars, "SpatRaster")) vars <- terra::as.data.frame(vars)
    vars <- vars[ , names(vals)]
  }

  for (i in 1:ncol(vals)) {
    if (is.null(vars)) d0 <- NULL
    else d0 <- density(vars[,i], na.rm = na.rm)
    d1 <- density(vals[,i], na.rm = na.rm)

    if (type == "density") {
      plot(NULL,
           xlim = range(c(d0$x, d1$x)),
           ylim = range(c(d0$y, d1$y)),
           main = names(vals)[i],
           ...)

      if (!is.null(vars)) polygon(c(d0$x, rev(d0$x)), c(d0$y, rep(0, length(d0$y))), col = rgb(0.5, 0.5, 0.5, 0.5), lwd = 0.3)
      polygon(c(d1$x, rev(d1$x)), c(d1$y, rep(0, length(d1$y))), col = rgb(0, 0.4, 0, 0.7), lwd = 0.3)

    } else if (type == "boxplot")
      if (is.null(vars))
        boxplot(vals[ , i],
                col = rgb(0, 0.4, 0, 0.7),
                notch = TRUE, lwd = 0.3,
                main = names(vals)[i],
                ...)
      else
        boxplot(vals[ , i], vars[ , i],
                col = c(rgb(0, 0.4, 0, 0.7), rgb(0.5, 0.5, 0.5, 0.5)),
                notch = TRUE, lwd = 0.3,
                main = names(vals)[i], names = c("vals", "vars"),
                ...)
  }
}
