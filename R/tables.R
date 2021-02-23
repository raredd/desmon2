### probability tables
# blt_table, dlt_table, pr_table
###


#' BLT table
#' 
#' Bayes (dose)-limiting toxicities table for designing phase I dose-
#' escalation studies.
#' 
#' @param r,n number of DLTs observed and total size
#' @param prior a vector of length two giving the distribution assumptions
#' for a beta prior (default is \code{beta(0.5, 0.5)})
#' @param alpha type-I error probability
#' @param cval for specifed proportions of DLTs, the probability that the rate
#' will be greater than \code{cval} is computed
#' @param digits number of digits to print
#' 
#' @return
#' A vector of five elements with the total size, number of DLTs, \code{alpha}
#' critical value, \code{1 - alpha} critical value, and calculated probability
#' that the DLT rate is greater than \code{cval}.
#' 
#' @seealso
#' \code{\link{dlt_table}}; \code{\link{pr_table}}
#' 
#' @examples
#' x <- seq(0.005, 0.995, length.out = 500)
#' plot(x, dbeta(x, 2, 3), type = 'l')
#' 
#' blt_table(3, 6, alpha = .0125)
#' t(Vectorize(blt_table)(0:8, 8))
#' 
#' @export

blt_table <- function(r, n, prior = c(1, 1) / 2, alpha = 0.025,
                      cval = 1 / 3, digits = 2L) {
  ciw <- round((1 - alpha * 2) * 100, 1L)
  prior <- rep_len(prior, 2L)
  p0 <- r + prior[1L]
  p1 <- n - r + prior[2L]
  
  l <- c(
    n, r,
    c(qbeta(alpha, p0, p1), qbeta(1 - alpha, p0, p1),
      pbeta(cval, p0, p1, lower.tail = FALSE))
  )
  
  setNames(
    as.list(l),
    c('Total size', 'No. DLTs', sprintf('L%s%%CI', ciw),
      sprintf('U%s%%CI', ciw), 'Pr(>|cval|)')
  )
}

#' DLT table
#' 
#' Create a standard 3+3 dose-limiting toxicity table with probabilities of
#' dose-escalation for true but unknown probabilities of experiencing
#' toxicity.
#' 
#' @param prob a vector of probabilities for event
#' @param digits number of digits past the decimal point to keep
#' 
#' @seealso
#' \code{\link{blt_table}}; \code{\link{pr_table}}
#' 
#' @examples
#' ## probabilities of experiencing event
#' prob <- c(1, 5, 1:5 * 10) / 100
#' 
#' dlt_table(prob)
#' t(dlt_table(prob, 3))
#' 
#' @export

dlt_table <- function(prob, digits = getOption('digits')) {
  res <- sapply(prob, function(pr) {
    dbinom(0L, 3L, pr) + dbinom(0L, 3L, pr) * dbinom(1L, 3L, pr)
  })
  
  res <- rbind(
    'Pr(DLT)' = prob,
    'Pr(Escalation)' = res
  )
  
  round(res, digits)
}

#' Probability table
#' 
#' A single-stage function for calculating the probabilities of \code{crit}
#' or fewer (greater if \code{greater = TRUE}) events for true but unknown
#' \code{prob}abilities of the event occurring.
#' 
#' @param prob a vector of probabilities for event
#' @param n,crit sample size and critical number of events, respectively; note
#' that the probabilities will always be calculated as weakly less or greater
#' than \code{crit}, i.e., \code{Pr(>= crit)} or \code{Pr(<= crit)} depending
#' if \code{greater} is \code{TRUE} or \code{FALSE}
#' @param greater \code{logical}; the direction of \code{crit}: if
#' \code{FALSE}, the probabilities are calculated for \code{crit}
#' \emph{or fewer} events; if \code{TRUE}, the probabilities are calculated
#' for \code{crit} \emph{or more} events
#' @param digits number of digits past the decimal point to keep
#' 
#' @seealso
#' \code{\link{blt_table}}; \code{\link{dlt_table}}
#' 
#' @examples
#' prob <- c(1, 5, 1:5 * 10) / 100
#' 
#' ## probabilities of crit or fewer events given true but unknown prob
#' pr_table(prob, n = 15, crit = 3, greater = FALSE)
#' 
#' ## probabilities of crit or more events given true but unknown prob
#' pr_table(prob, n = 15, crit = 3, greater = TRUE)
#' 
#' ## compare
#' 1 - as.numeric(pr_table(prob, n = 15, crit = 2, greater = FALSE)[2L, ])
#' 
#' @export

pr_table <- function(prob, n, crit, greater = FALSE, digits = getOption('digits')) {
  n <- as.integer(n)
  crit <- c(0L, sequence(max(as.integer(crit - greater))))
  
  res <- sapply(prob, function(pr) {
    x <- sum(mapply(function(x) dbinom(x, n, pr), crit))
    if (greater)
      1 - x else x
  })
  
  nn  <- c(
    'Pr(Event)',
    sprintf('Pr(%s%s)', c('<=', '>=')[greater + 1L], max(crit) + greater)
  )
  res <- structure(
    rbind(prob, res),
    dimnames = list(nn, NULL)
  )
  
  round(res, digits)
}
