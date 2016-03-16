
filled.contour2 = function (x = seq(0, 1, length.out = nrow(z)), y = seq(0, 1,
                                                                         length.out = ncol(z)), z, xlim = range(x, finite = TRUE),
                            ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE),
                            levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors,
                            col = color.palette(length(levels) - 1), plot.title, plot.axes,
                            key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1,
                            key.extend = FALSE,
                            axes = TRUE, frame.plot = axes, ...)
{
    if (missing(z)) {
        if (!missing(x)) {
            if (is.list(x)) {
                z <- x$z
                y <- x$y
                x <- x$x
            }
            else {
                z <- x
                x <- seq.int(0, 1, length.out = nrow(z))
            }
        }
        else stop("no 'z' matrix specified")
    }
    else if (is.list(x)) {
        y <- x$y
        x <- x$x
    }
    if (any(diff(x) <= 0) || any(diff(y) <= 0))
        stop("increasing 'x' and 'y' values expected")
    mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
    on.exit(par(par.orig))
    w <- (3 + mar.orig[2L]) * par("csi") * 2.54
    w <- lcm(w * ifelse(key.extend, 0.9, 0.9))
    layout(matrix(c(2, 1), ncol = 2L), widths = c(1, w))
    par(las = las)
    mar <- mar.orig
    mar[4L] <- mar[2L]
    mar[2L] <- 1
    par(mar = mar)
    plot.new()
    plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i",
                yaxs = "i")

    if (key.extend) {
        # expand levels by one step above and below
        dl <- diff(levels[1:2])   # level to level distance
        # draw key-color rectangles but skip the first and last level
        last <- length(levels)
        xi <- 0
        xa <- 1
        rect(xi, levels[1:(last-2)],
             xa, levels[2:(last-1)],
             col = col[1:(length(col)-1)])
        # allow drawing triangles into the margins
        apex <- 0.6   # apex height as factor of dl
        clipmax <- apex + (0.05*apex)  # add fudge factor 5%
        # to account for line width
        clip(xi,xa, levels[1]-(dl*clipmax), levels[last]+(dl*clipmax))
        # draw the range extension polygons
#         polygon(c(xi,xi,xa,xa,xa/2),
#                 c(levels[2]-(dl), levels[2], levels[2],
#                   levels[2]-(dl), levels[1]-(dl*apex)),
#                 col = col[1])
        polygon(c(xi,xi,xa,xa,xa/2),
                c(levels[last-1]+(dl), levels[last-1], levels[last-1],
                  levels[last-1]+(dl), levels[last]+(dl*apex)),
                col = col[length(col)])
    }
    else {
        rect(0, levels[-length(levels)], 1, levels[-1L], col = col)
    }

    if ((missing(key.axes) && axes)) {
        if (key.extend) {axis(4, lwd = 0, lwd.tick=1)}
            else {axis(4)}
        }
    else {
        key.axes
    }


    #####
     if (key.extend) {
        #clip(xi,xa, levels[1]-(dl*apex), levels[last]+(dl* apex))
        polygon(c(xi,xa/2,xa,xa,xa/2,xi),
                c(levels[2]-(dl),
                  levels[1]-(dl*1.5),
                  levels[2]-(dl),
                  levels[last-1]+(dl),
                  levels[last]+(dl*1.5),
                  levels[last-1]+(dl) ),
                lwd = 1.1 )
    }
    else {
        box()
    }
    if (!missing(key.title))
        key.title
    mar <- mar.orig
    mar[4L] <- 1
    par(mar = mar)
    plot.new()
    plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
    .filled.contour(x, y, z, levels, col)
    if (missing(plot.axes)) {
        if (axes) {
            title(main = "", xlab = "", ylab = "")
            Axis(x, side = 1)
            Axis(y, side = 2)
        }
    }
    else plot.axes
    if (frame.plot)
        box()
    if (missing(plot.title))
        title(...)
    else plot.title
    invisible()
}
