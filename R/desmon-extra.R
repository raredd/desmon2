### desmon extras
# bin1pow
###


#' Power for one-sample exact binomial tests
#' 
#' Determines the power for a one-sided, one-sample exact binomial test.
#' 
#' @param n sample size
#' @param pa probability of success under alternative hypothesis
#' @param r critical value
#' 
#' @examples
#' p0 <- 0.1
#' pa <- 0.3
#' 
#' des <- desmon2:::bin1samp(p0, pa)
#' 1 - des[['type2']]
#' 
#' bin1pow(des['n'], pa, des['r'] + 1L)
#' 
#' @export

bin1pow <- function(n, pa, r) {
  stopifnot(
    length(n) == 1L,
    length(pa) == 1L
  )
  
  bin1pow_ <- function(n, pa, r) {
    sum(dbinom(seq(r, n), n, pa))
  }
  
  Vectorize(bin1pow_, 'r', TRUE, FALSE)(n, pa, r)
}
