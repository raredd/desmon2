### desmon extras
# bin1pow
###


#' Power for one-sample exact binomial tests
#' 
#' Determines the power and significance level for a one-sided, one-sample
#' exact binomial test.
#' 
#' @param n sample size
#' @param p0,pa probability of success under the null and alternative
#' hypotheses, respectively
#' @param r critical value
#' 
#' @examples
#' p0 <- 0.1
#' pa <- 0.3
#' 
#' des <- desmon2:::bin1samp(p0, pa)
#' des[c('size', 'type2')]
#' 
#' bin1pow(des['n'], p0, pa, des['r'] + 1L)
#' 
#' @export

bin1pow <- function(n, p0, pa, r) {
  stopifnot(
    length(n) == 1L,
    length(pa) == 1L
  )
  
  bin1pow_ <- function(n, p, r) {
    sum(dbinom(seq(r, n), n, p))
  }
  
  c(type1 = bin1pow_(n, p0, r), type2 = 1 - bin1pow_(n, pa, r))
}
