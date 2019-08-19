### desmon extras
# bin1samp_power, twostg_power
###


#' Power for one-sample exact binomial designs
#' 
#' Determines the power and significance level for a one-sided, one-sample
#' exact binomial test.
#' 
#' @param p0,pa probability of success under the null and alternative
#' hypotheses, respectively
#' @param n sample size
#' @param r critical value
#' 
#' @return
#' A vector with the following elements:
#' 
#' \item{\code{type1}}{the overall type-I error}
#' \item{\code{type2}}{the overall type-II error}
#' 
#' @seealso
#' \code{\link{twostg_power}}
#' 
#' @examples
#' p0 <- 0.1
#' pa <- 0.3
#' des <- desmon2:::bin1samp(p0, pa)
#' 
#' bin1samp_power(p0, pa, des['n'], des['r'] + 1L)
#' 
#' ## compare
#' des[c('size', 'type2')]
#' 
#' @export

bin1samp_power <- function(p0, pa, n, r) {
  stopifnot(
    length(p0) == 1L,
    length(pa) == 1L,
    length(n) == 1L,
    length(r) == 1L
  )
  
  bin1pow_ <- function(n, p, r) {
    sum(dbinom(seq(r, n), n, p))
  }
  
  c(type1 = bin1pow_(n, p0, r), type2 = 1 - bin1pow_(n, pa, r))
}

#' Power for two-stage designs
#' 
#' Determines the operating characteristics of single-arm, two-stage designs.
#' 
#' @param p0,pa probability of success under the null and alternative
#' hypotheses, respectively
#' @param n1,n2 sample size of first and second stage
#' @param r1,r2 maximum number of responses in first stage and overall where
#' treatment is declared ineffective
#' 
#' @return
#' A vector with the following elements:
#' 
#' \item{\code{Pr.stop1.H0}}{probability of stopping after the first stage
#' if \code{p0} is true}
#' \item{\code{Pr.stop1.H1}}{probability of stopping after the first stage
#' if \code{pa} is true}
#' \item{\code{type1}}{the overall type-I error}
#' \item{\code{type2}}{the overall type-II error}
#' \item{\code{E.tot.n.H0}}{expected total sample size if \code{p0} is true}
#' \item{\code{E.tot.n.H1}}{expected total sample size if \code{pa} is true}
#' 
#' @seealso
#' \code{\link{bin1samp_power}}
#' 
#' @examples
#' p0 <- 0.1
#' pa <- 0.3
#' des <- desmon2:::simon(p0, pa)$designs[1L, ]
#' 
#' twostg_power(p0, pa, des[['n1']], des[['n2']], des[['r1']], des[['r2']])
#' 
#' ## compare
#' des
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
