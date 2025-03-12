### probability tables, rand tables
# blt_table, dlt_table, pr_table, ransch, ranschtbl, summary.ransch
#
# unexported:
# ransch_
###


#' BLT table
#' 
#' Bayes (dose)-limiting toxicities table for designing phase I dose-
#' escalation studies.
#' 
#' @param r,n number of DLTs observed and total size
#' @param prior a vector of length two giving the distribution assumptions
#'   for a beta prior (default is \code{beta(0.5, 0.5)})
#' @param alpha type-I error probability
#' @param cval for specifed proportions of DLTs, the probability that the rate
#'   will be greater than \code{cval} is computed
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
#'   that the probabilities will always be calculated as weakly less or greater
#'   than \code{crit}, i.e., \code{Pr(>= crit)} or \code{Pr(<= crit)} depending
#'   if \code{greater} is \code{TRUE} or \code{FALSE}
#' @param greater \code{logical}; direction of \code{crit}: if \code{FALSE},
#'   the probabilities are calculated for \code{crit} \emph{or fewer} events;
#'   if \code{TRUE}, the probabilities are calculated for \code{crit}
#'   \emph{or more} events
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

#' ransch
#'
#' Generate block randomization schedule tables with fixed or random block
#' sizes.
#'
#' @param n sample size of study or each stratum
#' @param block block size; note if \code{block} is not a factor of \code{n},
#'   \code{n} will be increased to accommodate a full block
#'
#'   for randomly-sized blocks, a vector of potential block sizes; note that
#'   a block size must be a multiple of \code{sum(r)}
#' @param arms names of the treatment arms
#' @param r randomization ratio; see examples
#' @param strata an optional named list of vectors for each stratum
#' @param write a file path to write a directory with csv files for each
#'   randomization table
#'
#' @examples
#' ## no strata, 2-3 treatments with varying randomization ratios
#' ransch(24, 4, 1:2) ## 1:1
#' ransch(24, 6, 1:3) ## 1:1:1
#' ransch(24, 8, 1:3, c(1, 2, 1)) ## 1:2:1
#'
#'
#' ## randomly-sized blocks
#' ransch(24, c(2, 4, 6), 1:2)
#'
#' set.seed(1)
#' r1 <- ransch(24, c(3, 6, 9), 1:3)
#' set.seed(1)
#' r2 <- ransch(24, 1:10, 1:3)
#'
#' ## note that these two are the same since only blocks sized 3, 6, 9
#' ## work for 1:1:1 randomization
#' identical(r1, r2)
#' addmargins(table(r1[[1]][, -1]))
#' summary(r1)
#'
#' 
#' ## one two-level stratum
#' ransch(24, 4, 1:2, strata = list(Age = c('<65', '>=65')))
#'
#' ## multiple strata
#' strata <- list(Site = LETTERS[1:3], Age = c('<65', '>=65'))
#' r3 <- ransch(24, 4, 1:2, strata = strata)
#' summary(r3)
#'
#'
#' ## tables for printing
#' ranschtbl(24, 4, c('Pbo', 'Trt'))
#' ranschtbl(24, 4, c('Pbo', 'Trt'), c(1, 3), strata)
#'
#' \dontrun{
#' ## use a path to write tables to file (one csv per stratum)
#' ranschtbl(24, 4, c('Pbo', 'Trt'), strata = strata, write = '~/desktop')
#' }
#'
#' @export

ransch <- function(n, block, arms, r = rep_len(1L, length(arms)),
                   strata = NULL) {
  if (!is.null(strata)) {
    if (is.null(names(strata)))
      names(strata) <- paste0('Stratum', seq_along(strata))
    
    strata <- Map(paste, names(strata), strata)
    strata <- apply(expand.grid(strata), 1L, toString)
  }
  
  res <- replicate(pmax(1L, length(strata)), simplify = FALSE, {
    ransch_(n, block, arms, rep_len(r, length(arms)))
  })
  class(res) <- c('ransch')
  
  setNames(res, if (is.null(strata)) 'ransch' else strata)
}

#' @rdname ransch
#' @export
ranschtbl <- function(n, block, arms, r = rep_len(1L, length(arms)),
                      strata = NULL, write = NULL) {
  res <- ransch(n, block, arms, r, strata)
  
  ## add columns to each table
  res[] <- lapply(seq_along(res), function(ii) {
    within(res[[ii]], {
      Stratum <- names(res)[ii]
      Name <- ID <- Date <- NA
    })
  })
  
  class(res) <- c('ransch')
  
  if (is.null(names(res)))
    names(res) <- 'ransch'
  
  if (is.character(write) && dir.exists(write)) {
    path <- sprintf('%s/ransch-%s', write, format(Sys.time(), '%Y%m%d%H%M'))
    dir.create(path, showWarnings = FALSE, recursive = TRUE)
    
    for (ii in seq_along(res))
      write.csv(res[[ii]], sprintf('%s/%s.csv', path, names(res)[ii]),
                na = '', row.names = FALSE)
    
    invisible(res)
  } else res
}

#' @rdname ransch
#' @export
summary.ransch <- function(x, totals = TRUE, verbose = TRUE, ...) {
  if (inherits(x, 'data.frame'))
    x <- list(ransch = x)
  
  if (verbose) {
    nstrata <- if (length(x) > 1)
      length(x) else 0
    nmstrata <- if (nstrata > 0)
      names(x)
    arms <- sapply(x, function(y) length(unique(y$Assignment)))
    n_strata <- sapply(x, nrow)
    n_total <- sum(n_strata)
    bl <- sapply(x, function(y)
      toString(unique(table(y$Block))))
    
    message('Total randomized: ', n_total)
    if (nstrata > 0) {
      message(
        sprintf('%s: %s (Block size: %s, Arms: %s)\n',
                nmstrata, n_strata, bl, arms)
      )
    } else {
      message('Arms: ', toString(arms))
      message('Block size: ', bl)
    }
  }
  
  lapply(x, function(x) {
    res <- table(x[, c('Block', 'Assignment')])
    if (totals)
      addmargins(res, 1:2) else res
  })
}

ransch_ <- function(n, block, arms, r) {
  ## table(ransch_(12, 6, c('Pbo', 'Trt'), c(1, 1))[, -1])
  ## table(ransch_(12, 1:4, c('Pbo', 'Trt'), c(1, 1))[, -1])
  stopifnot(length(arms) == length(r))
  
  sample <- function(x, ...) {
    x[sample.int(length(x), ...)]
  }
  rblock <- function(b, arms, r) {
    arms <- rep_len(rep(arms, r), b)
    sample(arms)
  }
  
  ok <- block %% sum(r) == 0L
  if (any(!ok)) {
    warning(
      sprintf('removing block sizes (%s) not compatible with r (%s)',
              toString(block[!ok]), paste0(r, collapse = ':')),
      call. = FALSE
    )
  }
  block <- block[ok]
  # block <- sample(block, n, replace = TRUE)
  block <- if (length(unique(block)) == 1L)
    rep_len(block, n) else sample(block, n, replace = TRUE)
  
  idx <- cumsum(block) < n
  block <- block[c(which(idx), sum(idx) + 1L)]
  
  res <- data.frame(
    Number = seq.int(sum(block)),
    Block = rep(seq_along(block), block),
    Assignment = unlist(lapply(block, rblock, arms = arms, r = r))
  )
  
  structure(res, class = c('ransch', class(res)))
}
