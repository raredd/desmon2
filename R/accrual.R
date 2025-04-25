### accrual utils
# accrual, accrual_months
###


#' Accrual plot
#' 
#' Plot (unitless but usually monthly) accrual numbers. \code{accrual_months}
#' is a helper function to convert a vector of dates into the proper format
#' for \code{actual}.
#' 
#' @param actual,expected integer vectors of actual and expected for each time
#'   unit; note these are expected to be sequential with no gaps, so units with
#'   0 accrued should be included
#' @param target target accrual
#' @param time0 accrual at time 0; default is 0
#' @param lag time units to shift the accrual line
#' @param, col,pch,type arguments controlling the accrual line, passed to
#'   \code{\link{points}}
#' @param start,end,full labels for first, last, and full accrual times; if
#'   \code{NA} or \code{NULL}, these will not be added
#' @param pos.full for the \code{full} label, \code{pos} value pased to
#'   \code{\link[text]}
#' @param pct,fill colored (\code{fill}) bands for accrual ranges based on
#'   expected accrual; use \code{NA} or \code{NULL} to suppress
#' @param pct.legend label for \code{pct} legend; use \code{NA} or \code{NULL}
#'   to suppress
#' @param args.pct a \emph{named} list of arguments passed to \code{\link{legend}}
#'   controlling the \code{pct.legend}
#' @param legend optional \emph{named} vector to annotate plot; if \code{TRUE},
#'   one will be made; use \code{NA} or \code{NULL} to suppress; see examples
#' @param xlim,ylim the x- and y-axis limits
#' @param xlab,ylab,ylab2 the x-, y-, and secondary y-axis labels
#' @param axes logical; if \code{TRUE}, plot axes and box are drawn
#' @param ... additional arguments passed to \code{plot}
#' @param dates a vector of \code{\link{Dates}} or a class to be coerced
#' 
#' @return
#' A list with the \code{x} and \code{y} coordinates of the accrual line.
#' Additionally, the coordinates for polygon(s) for \code{pct}, if given.
#' 
#' @examples
#' ## basic usage
#' set.seed(1)
#' ac <- rpois(12, 2)
#' ex <- rep(c(1, 3, 5), times = c(1, 3, 7))
#' accrual(ac, ex)
#' accrual(ac, ex, start = 'Jan 2000', legend = TRUE, fill = adjustcolor(1:3, 0.5))
#' accrual(ac, ex, start = 'Jan 2000', args.pct = list(x = 'topleft'),
#'         pct = c(1.25, 1, 0.75, 0.5, 0.25), fill = topo.colors(4))
#' 
#' ## minimal plot
#' accrual(ac, ex, pct = NA, start = NA, end = NA, full = NA, xlab = 'Time',
#'         axes = FALSE, panel.first = {box(); axis(1); axis(2)}, ylim = c(0, 50))
#' 
#' ## adding multiple accrual lines
#' accrual(ac, ex)
#' accrual(ac * 0.75, col = 'blue')
#' accrual(ac[1:4], col = 'blue', lag = 6)
#' points(c(0, seq_along(ac)), c(0, cumsum(ac)), cex = 2)
#' 
#' ## adding legends
#' st <- 'May 2025'
#' en <- 'Apr 2026'
#' lg <- c(
#'   'Start Date' = st,
#'   'Expected' = sum(ex), 'Exp. end date' = en, 'Exp. per month' = '2-4',
#'   'Total to date' = sprintf('%s (%.0f%%)', sum(ac), sum(ac) / sum(ex) * 100),
#'   'N per month' = round(mean(ac), 1)
#' )
#' op <- par(mar = c(5, 5, 2, 5), las = 1L)
#' accrual(ac, ex, legend = lg, start = st, end = en, fill = c(3, 7, 2),
#'         pos.full = 1L)
#' par(op)
#' 
#' ## using a vector of dates
#' x <- Sys.Date() + unlist(sapply(seq_along(ac), function(ii)
#'   30 * ii + rpois(ac[ii], 5)))
#' ## these are equivalent ways to call accrual with date vectors
#' a1 <- accrual(x, ex)
#' a2 <- accrual(accrual_months(x), ex)
#' identical(a1, a2)
#' 
#' ## the plot data is returned
#' plot(
#'   a1$x, a1$y, ylim = c(0, max(unlist(a1))),
#'   panel.first = matlines(a1$pct$x, a1$pct[, -1])
#' )
#' 
#' @export

accrual <-
  function(actual, expected, target = sum(expected),
           col = par('col'), pch = 18L, type = 'b',
           time0 = 0, lag = 0,
           start = 'Month 0', end = 'Last month', full = 'Full accrual expected',
           pos.full = 2L,
           pct = c(1, 0.75, 0.25, 0), fill = c('yellow', 'orange', 'red'),
           pct.legend = 'Percent of expected', args.pct = list(),
           legend = NULL, xlim = NULL, ylim = NULL,
           xlab = NULL, ylab = NULL, ylab2 = NULL, axes = TRUE, ...) {
    ok <- function(x) {
      !(is.null(x) || all(is.na(x)) || identical(x, FALSE))
    }
    if (inherits(actual, c('Date', 'POSIXct', 'POSIXlt', 'POSIXt')))
      actual <- accrual_months(actual)
    if (inherits(actual, 'accrual_months')) {
      if (missing(start))
        start <- attr(actual, 'start')
      if (missing(end))
        end <- attr(actual, 'end')
    }
    xa <- seq_along(actual)
    ya <- actual
    if (missing(expected)) {
      x <- c(0, xa) + lag
      y <- c(time0, cumsum(ya))
      points(x, y, xpd = NA, col = col, pch = pch, type = type)
      return(invisible(list(x = x, y = y, pct = NULL)))
    }
    xe <- seq_along(expected)
    ye <- expected
    
    xmax <- max(c(xa, xe))
    ymax <- max(target, sum(actual))
    plot(
      xa, ya, ..., axes = FALSE, type = 'n',
      xlim = xlim %||% c(0, max(xa)),
      ylim = ylim %||% c(-ymax / 5, ymax),
      xlab = xlab %||% 'Months since accrual start',
      ylab = ylab %||% 'Percent of target accrual'
    )
    
    p <- par('usr')
    if (axes) {
      box(bty = par('bty'))
      abline(h = 0, lwd = 0.5)
      axis(1L, seq(0, xmax), mgp = par('mgp') * -1.5, tcl = 0.2)
      if (ok(start))
        axis(1L, 0, start)
      if (ok(end))
        axis(1L, xmax, end)
      at <- seq(0, target, length.out = 5L)
      axis(2L, at, at / target * 100)
      axis(4L, pretty(c(0, target)))
      text(diff(p[1:2]) * 1.05, sum(p[3:4]) / 2,
           ylab2 %||% 'Number accrued', srt = 90, xpd = NA)
    }
    
    ## summary legend
    if (ok(legend)) {
      if (isTRUE(legend)) {
        legend <- c(
          'Start date' = ifelse(missing(start), NA, start),
          'End date' = ifelse(missing(end), NA, end),
          'Expected accrual' = target,
          'Expected per month' = round(mean(ye), 1L),
          'Actual to date' = sprintf('%s (%.0f%%)', sum(ya), sum(ya) / sum(ye) * 100),
          'Actual per month' = round(mean(ya), 1L)
        )
        legend <- na.omit(legend)
      }
      if (is.null(names(legend)))
        names(legend) <- legend
      pad <- max(strwidth(names(legend)))
      d1 <- diff(p[1:2]) / 50
      d2 <- max(strwidth(legend))
      p0 <- function(x) {
        paste0(x, collapse = '\n')
      }
      text(0, p[4L], p0(names(legend)), adj = c(0, 1.1), xpd = NA)
      text(pad + d1 + d2, p[4L], p0(legend), adj = c(1, 1.1), xpd = NA)
    }
    
    co <- NULL
    if (ok(pct)) {
      n <- c(0, ye)
      n <- cumsum(n)
      x <- seq(0, length(n) - 1L)
      pct <- sort(pct, decreasing = TRUE)
      fill <- rep_len(fill, length(pct))
      px <- c(x, rev(x))
      co <- data.frame(x = px)
      
      for (ii in seq_along(pct[-1L])) {
        y <- n * pct[ii]
        z <- n * pct[ii + 1L]
        py <- c(y, rev(z))
        polygon(px, py, col = fill[ii], border = NA)
        lines(x, y, lwd = 0.25)
        points(x, y, pch = 19L, cex = 0.1)
        co[, paste0('y', ii)] <- py
      }
      
      leg <- rev(pct) * 100
      leg <- rev(sprintf('%s - %s', c(0, leg[-length(leg)]), leg)[-1L])
      if (ok(pct.legend)) {
        args <- list(
          title = pct.legend, xpd = NA,
          x = 0, y = ymax / 2,
          # x = 'bottom', horiz = TRUE, inset = c(0, 0),
          legend = leg, bty = 'n', fill = fill, border = NA
        )
        do.call('legend', modifyList(args, args.pct))
      }
    }
    
    ## expected full accrual marker
    if (ok(full)) {
      # text(max(xe), p[3L], full, xpd = NA, adj = c(0.5, -0.5))
      # points(max(xe), p[3L], pch = 18L, xpd = NA, cex = 1.5)
      pad <- p[3L] * 2
      segments(max(xe), p[3L], max(xe), pad, xpd = NA, lty = 'dashed')
      text(max(xe), pad, full, xpd = NA, pos = pos.full)
    }
    
    ## actual accrual by month
    x <- c(0, xa) + lag
    y <- c(time0, cumsum(ya))
    points(x, y, xpd = NA, col = col, pch = pch, type = type)
    
    invisible(list(x = x, y = y, pct = co))
  }

#' @rdname accrual
#' @export
accrual_months <- function(dates) {
  x <- as.Date(sort(dates))
  f <- ftable(factor(format(x, '%Y')), factor(as.numeric(format(x, '%m')), 1:12))
  f <- t(as.matrix(f))
  d <- setNames(c(f), sprintf('%s %s', month.abb[row(f)], colnames(f)[col(f)]))
  r <- format(range(x), '%b %Y')
  i <- which(names(d) %in% r)
  
  structure(
    d[seq(i[1L], i[2L])],
    class = 'accrual_months',
    start = r[1L], end = r[2L]
  )
}
