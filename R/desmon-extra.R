### desmon extras
# bin1samp_power, bin1samp_sim, twostg_power, twostg_sim
###


#' Power for one-sample exact binomial designs
#' 
#' @description
#' Determines the power and significance level for a one-sided, one-sample
#' exact binomial test.
#' 
#' \code{bin1samp_power} accepts a single value or vector for \code{n} and/or
#' \code{r} and will return a matrix with results for each combination. If
#' only one value is given for each, a vector is returned.
#' 
#' @param p0,pa probability of success under the null and alternative
#'   hypotheses, respectively
#' @param n sample size
#' @param r a vector of critical values, typically the minimum number of
#'   successes required to reject \code{pa}
#' @param plot logical; if \code{TRUE}, a the sequence of \code{r} versus
#'   type-I and type-II errors is plotted
#' 
#' @return
#' \code{bin1samp_power} returns a vector or matrix with the following:
#' 
#' \item{\code{type1}}{the overall type-I error}
#' \item{\code{type2}}{the overall type-II error}
#' 
#' \code{bin1samp_sim} returns a data frame with the following columns:
#' 
#' \item{\code{r}}{critical values}
#' \item{\code{type1}}{the overall type-I errors}
#' \item{\code{type2}}{the overall type-II errors}
#' 
#' @seealso
#' \code{\link{twostg_power}}; \code{\link{twostg_sim}}
#' 
#' @examples
#' p0 <- 0.1
#' pa <- 0.3
#' des <- desmon2:::bin1samp(p0, pa)
#' bin1samp_power(p0, pa, des['n'], des['r'] + 1L)
#' 
#' ## compare
#' des[c('size', 'type2')]
#' 
#' bin1samp_power(p0, pa, des['n'], des['r'] + -2:2)
#' bin1samp_power(p0, pa, des['n'] + 0:1, des['r'] + -2:2)
#' 
#' 
#' ## simulate over critical values
#' bin1samp_sim(p0, pa, des['n'])
#' 
#' @export

bin1samp_power <- function(p0, pa, n, r) {
  stopifnot(
    length(p0) == 1L,
    length(pa) == 1L
  )
  
  exp <- expand.grid(r = r, n = n)
  
  bin1pow_ <- function(n, p, r) {
    sum(dbinom(seq(r, n), n, p))
  }
  
  res <- t(sapply(seq.int(nrow(exp)), function(ii)
    c(type1 = bin1pow_(exp$n[ii], p0, exp$r[ii]),
      type2 = 1 - bin1pow_(exp$n[ii], pa, exp$r[ii]))))
  rownames(res) <- sprintf('n=%s, r=%s', exp$n, exp$r)
  
  if (nrow(exp) == 1L)
    drop(res) else res
}

#' @rdname bin1samp_power
#' @export
bin1samp_sim <- function(p0, pa, n, r = seq.int(n), plot = TRUE) {
  res <- sapply(r, function(x) bin1samp_power(p0, pa, n, x))
  res <- data.frame(r = r, t(res))
  
  if (plot) {
    matplot(
      r, res[, -1L], type = 'l', las = 1L,
      xlab = sprintf('critical value, r of %s', n),
      ylab = 'Error probability'
    )
    legend(
      'right', col = 1:2, lty = 1:2, bty = 'n',
      legend = paste('Type', c('I', 'II'))
    )
  }
  
  ## alpha < 0.1 and power > 0.8
  res$` ` <- ifelse(res$type1 < 0.1 & res$type2 < 0.2, '*', '')
  
  res
}

#' Power for two-stage designs
#' 
#' Determines the operating characteristics of single-arm, two-stage designs.
#' 
#' @param p0,pa probability of success under the null and alternative
#'   hypotheses, respectively
#' @param n1,n2 sample size of first and second stage
#' @param r1,r2 maximum number of responses in first stage and overall where
#'   treatment is declared ineffective
#' @param plot logical; if \code{TRUE}, the type-I and type-II errors are
#'   plotted for each combination of (valid) \code{r1} and \code{r2} values
#' 
#' @return
#' \code{twostg_power} returns vector with the following elements:
#' 
#' \item{\code{Pr.stop1.H0}}{probability of stopping after the first stage
#'   if \code{p0} is true}
#' \item{\code{Pr.stop1.H1}}{probability of stopping after the first stage
#'   if \code{pa} is true}
#' \item{\code{type1}}{the overall type-I error}
#' \item{\code{type2}}{the overall type-II error}
#' \item{\code{E.tot.n.H0}}{expected total sample size if \code{p0} is true}
#' \item{\code{E.tot.n.H1}}{expected total sample size if \code{pa} is true}
#' 
#' \code{twostg_sim} returns a data frame with columns for each of the above
#' plus the following:
#' 
#' \item{\code{r1}}{critical values for the first stage}
#' \item{\code{r2}}{critical values for the second stage}
#' 
#' @seealso
#' \code{\link{bin1samp_power}}; \code{\link{bin1samp_sim}}
#' 
#' @examples
#' p0 <- 0.1
#' pa <- 0.3
#' des <- desmon2:::simon(p0, pa)$designs[1L, ]
#' twostg_power(p0, pa, des[['n1']], des[['n2']], des[['r1']], des[['r2']])
#' 
#' ## compare
#' des
#' 
#' 
#' ## simulate over critical values
#' twostg_sim(p0, pa, des[['n1']], des[['n2']])
#' 
#' \dontrun{
#' res <- twostg_sim(p0, pa, des[['n1']], des[['n2']])
#' with(res, {
#'   iplotr::iscatter(
#'     type1, type2, group = type1 < 0.1 & type2 < 0.2,
#'     labels = list(
#'       r1 = r1, r2 = r2,
#'       alpha = round(type1, 3), power = round(1 - type2, 3)
#'     )
#'   )
#' })
#' }
#' 
#' @export

twostg_power <- function(p0, pa, n1, n2, r1, r2) {
  stopifnot(
    length(p0) == 1L,
    length(pa) == 1L,
    length(n1) == 1L,
    length(n2) == 1L,
    length(r1) == 1L,
    length(r2) == 1L
  )
  
  nul <- twostg(n1, n2, p0, r1, r2)$prob
  alt <- twostg(n1, n2, pa, r1, r2)$prob
  
  c(Pr.stop1.H0 = nul[[2L]], Pr.stop1.H1 = alt[[2L]],
    type1 = 1 - nul[[1L]], type2 = alt[[1L]],
    E.tot.n.H0 = n1 + n2 * (1 - nul[[2L]]),
    E.tot.n.H1 = n1 + n2 * (1 - alt[[2L]]))
}

#' @rdname twostg_power
#' @export
twostg_sim <- function(p0, pa, n1, n2, r1 = seq.int(n1),
                       r2 = seq.int(n1 + n2), plot = TRUE) {
  exp <- expand.grid(r1 = r1, r2 = r2)
  exp <- exp[exp[, 2L] > exp[, 1L] & n1 > exp[, 1L] & (n1 + n2) > exp[, 2L], ]
  res <- Map(twostg_power, p0, pa, n1, n2, exp[, 1L], exp[, 2L])
  res <- data.frame(do.call('rbind', res), r1 = exp[, 1L], r2 = exp[, 2L])
  grp <- (res$type1 < 0.1 & res$type2 < 0.2) + 1L
  res$` ` <- ifelse(grp == 2L, '*', '')
  
  if (plot) {
    plot(type2 ~ type1, res, las = 1L, col = grp, pch = c(1L, 16L)[grp])
    legend(
      'topright', col = 2L, pch = 16L, bty = 'n',
      legend = expression(alpha < 0.1 ~ power > 0.8)
    )
  }
  
  res
}
