#' Frequency plot
#'
#' This function plots (in green) the frequency distribution of each variable at the presence localities, and optionally (in grey) also the frequency across a given study region, using either density or box plots.
#'
#' This function provides a simple way to visualize the frequency distribution of one or more variables, optionally comparing the distribution of \code{vals} with \code{vars}. It supports two types of plots: density plots and boxplots. When \code{vars} is provided, the function overlays the comparison dataset onto the same plot.
#'
#' @param vals data frame with the values of the variables (columns) at the presence localities (rows). Column names must match the names of the corresponding variables in 'vars'.
#' @param vars optional SpatRaster or data frame with the variables with which to visually compare the frequency of 'vals'.
#' @param type type of plot to produce: "density" or "boxplot". Case is ignored and partial argument matching is used (i.e., you can set e.g. type="Box")
#' @param na.rm logical (default TRUE) indicating whether to remove rows with missing values in 'vals'. Otherwise, some methods may produce an error or empty results.
#' @param dup.rm logical (default FALSE) indicating whether to remove rows with duplicate values for all variables in 'vals'.
#' @param verbosity numeric value indicating the amount of messages to display. The default is 2, for the maximum number of messages available.
#' @param ... Optional additional arguments that can be passed to plot(), e.g. 'cex.main' or 'cex.axis'
#'
#' @returns A plot of the specified type.
#' @export
#'
#' @examples

freqPlot <- function(vals, vars = NULL, type = "density", na.rm = TRUE, dup.rm = FALSE, verbosity = 2, ...) {

    type <- match.arg(tolower(type), choices = c("density", "boxplot"))

    vals <- as.data.frame(vals)

    vals <- dataPrune(vals, na.rm = na.rm, dup.rm = dup.rm, verbosity = verbosity)

    if (!is.null(vars)) {
        if (inherits(vars, "SpatRaster")) vars <- terra::as.data.frame(vars)
        vars <- vars[, names(vals)]
    }

    for (i in 1:ncol(vals)) {
        if (is.null(vars)) {
            d0 <- NULL
        } else {
            d0 <- density(vars[, i], na.rm = na.rm)
        }
        d1 <- density(vals[, i], na.rm = na.rm)

        if (type == "density") {
            plot(NULL,
                xlim = range(c(d0$x, d1$x)),
                ylim = range(c(d0$y, d1$y)),
                main = names(vals)[i],
                ...
            )

            if (!is.null(vars)) polygon(c(d0$x, rev(d0$x)), c(d0$y, rep(0, length(d0$y))), col = rgb(0.5, 0.5, 0.5, 0.5), lwd = 0.3)
            polygon(c(d1$x, rev(d1$x)), c(d1$y, rep(0, length(d1$y))), col = rgb(0, 0.4, 0, 0.7), lwd = 0.3)
        } else if (type == "boxplot") {
            if (is.null(vars)) {
                boxplot(vals[, i],
                    col = rgb(0, 0.4, 0, 0.7),
                    notch = TRUE, lwd = 0.3,
                    main = names(vals)[i],
                    ...
                )
            } else {
                boxplot(vals[, i], vars[, i],
                    col = c(rgb(0, 0.4, 0, 0.7), rgb(0.5, 0.5, 0.5, 0.5)),
                    notch = TRUE, lwd = 0.3,
                    main = names(vals)[i], names = c("vals", "vars"),
                    ...
                )
            }
        }
    }
}
