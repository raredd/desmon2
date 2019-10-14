### stat text
# bin1samp_text, mtd_text, simon_text, twostg_text
# 
# S3 methods
# print.stat_text
###


#' Stat text printing method
#' 
#' @param x an object of class \code{"stat_text"}
#' @param width a positive integer giving the target column for wrapping
#' lines in the output; see \code{\link{strwrap}}
#' @param ... ignored
#' 
#' @export

print.stat_text <- function(x, width = NULL, ...) {
  if (is.null(width))
    width <- getOption('width') * 0.8
  cat(strwrap(x, width = width), sep = '\n')
  cat('\n')
  
  if (!is.null(extra <- attr(x, 'extra'))) {
    cat(strwrap(extra, width = width), sep = '\n')
    cat('\n')
  }
  
  invisible(x)
}

#' One-sample design text
#' 
#' Output one-sided, one-sample, single-stage exact binomial design text.
#' 
#' @param p0,pa the null and alternative hypotheses
#' @param ... additional arguments passed to \code{\link[desmon]{bin1samp}}
#' @param conf confidence level for confidence intervals
#' @param outcome text string describing the outcome; if a string such as
#' \code{"long description (abbr)"} is given, then entire string is used in
#' the first instance, and only the text in parens is used subsequently
#' 
#' @family designs
#' 
#' @seealso
#' \code{\link[desmon]{bin1samp}}; \code{\link[desmon]{binci}}
#' 
#' @return
#' A character string describing the design with \code{attr(, "design")}
#' containing the design.
#' 
#' @examples
#' ## basic usage
#' bin1samp_text(.3, .5)
#' 
#' \dontrun{
#' ## use cat and/or strwrap for improved formatting and writing
#' cat(strwrap(bin1samp_text(.3, .5), width = 80), sep = '\n')
#' 
#' cat(bin1samp_text(.3, .5), file = '~/bin1samp_design.txt')
#' 
#' cat(bin1samp_text(.3, .5, outcome = 'overall response (OR) rate'),
#'     file = '~/bin1samp_design2.txt')
#' }
#' 
#' @export

bin1samp_text <- function(p0, pa, ..., conf = 0.95,
                          outcome = 'OUTCOME (OUT)') {
  if (!nzchar(outcome2 <- gsub('\\((.*?)\\)|.', '\\1', outcome)))
    outcome2 <- outcome
  
  call <- bin1samp(p0, pa, ...)
  psep <- '\n\n'
  
  args <- as.list(call)
  
  onewid <- sapply(seq.int(args$n), function(x) {
    ci <- binci(x, args$n, conf = conf)
    unname(diff(ci[1:2]) / 2)
  })
  
  txt <- paste(
    'The study will use a single-stage, exact binomial design. The', outcome,
    'of at least', sprintf('%s%%', pa * 100),
    'will be considered promising whereas', outcome2, 'of',
    paste0(p0 * 100, '%'), 'or less will be considered non-promising.',
    
    psep,
    
    args$n, 'eligible patients will be enrolled. If', args$r, 'or fewer',
    outcome2, 'are observed, the regimen will be considered non-promising,',
    'and the study will be unsuccessful. If at least', args$r + 1L,
    outcome2, 'are observed in', args$n,
    'patients, the study will be considered successful and the regimen',
    'worthy of further study.',
    
    psep,
    
    'The study has an overall power and type-I error of',
    sprintf('%.3f', 1 - args$type2), 'and', sprintf('%.3f', args$size),
    ', respectively. With a total size of', args$n,
    'patients, the single-stage exact 95% confidence interval for', outcome2,
    'will be no wider than +/-', sprintf('%s%%.', round(max(onewid) * 100)),
    
    psep,
    
    'If the true', outcome2, 'is', sprintf('%s%%,', p0 * 100),
    'the regimen will be considered non-promising with probability',
    sprintf('%.3f;', 1 - args$size), 'and promising with probability',
    sprintf('%.3f', args$size), '(type-I error).',
    
    psep,
    
    'If the true', outcome2, 'is', sprintf('%s%%,', pa * 100),
    'the regimen will be considered non-promising with probability',
    sprintf('%.3f', args$type2),
    '(type-II error) and promising with probability',
    sprintf('%.3f.', 1 - args$type2)
  )
  
  structure(
    trim(txt, TRUE), design = call, class = 'stat_text'
  )
}

#' MTD design text
#' 
#' Output standard 3+3 maximum tolerated dose design text.
#' 
#' @param prob a vector of toxicity event rates
#' @param ndose number of dose levels
#' @param expansion size of the expansion cohort
#' @param digits number of digits past the decimal point to keep
#' 
#' @family designs
#' 
#' @seealso
#' \code{\link{dlt_table}}; \code{\link{sim3p3}}
#' 
#' @return
#' A character string describing the design with \code{attr(, "table")}
#' containing the DLT table in markdown format.
#' 
#' @examples
#' ## basic usage
#' cat(mtd_text())
#' cat(mtd_text(1:10 / 20, 5, 15))
#' 
#' \dontrun{
#' cat(mtd_text(), file = '~/mtd_design.txt')
#' }
#' 
#' @export

mtd_text <- function(prob = 1:5 / 10, ndose = 3L, expansion = 10L,
                     digits = 2L) {
  ndose <- ndose[1L]
  
  bin <- Vectorize(binci)
  ciw <- 3 + expansion
  ciw <- apply(bin(seq.int(ciw), ciw), 2L, diff) / 2
  
  psep <- '\n\n'
  
  txt <- paste(
    'This study will utilize a standard 3+3 dose-escalation to establish the',
    'maximum tolerated dose (MTD).',
    
    psep,
    
    'Three patients will be treated at the starting dose. If no DLTs are',
    'observed in the first three patients, the dose is escalated, and three',
    'additional patients are treated. If all dose levels are exhausted, then',
    'the highest level with at least three treated patients with no DLT will',
    'be the MTD.',
    
    psep,
    
    'If one DLT is observed in the first three patients treated per level,',
    'an additional three patients are treated at the current level. If no',
    'additional DLTs occur, the dose will be escalated. If two or more DLTs',
    'are observed on a dose level, the previous level will be the MTD,',
    'and the dose-escalation phase of the study will stop. No re-escalations',
    'will occur.',
    
    psep,
    
    'If two or more DLTs are observed on a dose level with no lower or safe',
    'doses, then no dose will be considered tolerable, and the study will',
    'be stopped for safety.',
    
    psep,
    
    'Once the MTD has been established, an additional expansion cohort of',
    expansion, 'patients will be entered at the MTD to confirm safety.',
    
    psep,
    
    'The dose-finding phase of the study will enroll up to', ndose * 6,
    'patients over', ndose, 'dose levels. With at least 3 +', expansion,
    'treated at the MTD, the 95% exact binomial confidence interval for the',
    'observed rate of DLTs will be no wider than +/-',
    sprintf('%s%%.', round(max(ciw) * 100))
  )
  
  if (!is.null(prob)) {
    tbl <- dlt_table(prob, digits)
    extra <- t(c(sprintf('**%s**', rownames(tbl)[2L]), tbl[2L, ]))
    colnames(extra) <- c(rownames(tbl)[1L], tbl[1L, ])
    txt <- c(
      txt,
      psep,
      unclass(knitr::kable(extra, format = 'markdown')),
      psep
    )
    txt <- paste(txt, collapse = '\n')
  }
  
  structure(
    trim(txt, TRUE), table = tbl, class = 'stat_text'
  )
}

#' Simon design text
#' 
#' Output Simon two-stage design text.
#' 
#' @param p0,pa the null and alternative hypotheses
#' @param ... additional arguments passed to \code{\link[desmon]{simon}}
#' @param conf confidence level for single- and two-stage confidence intervals
#' @param which optional; an integer selecting the design to use if multiple
#' are found
#' @param outcome text string describing the outcome; if a string such as
#' \code{"long desc (abbr)"} is given, then entire string is used in the first
#' instance, and only the text in parens is used subsequently
#' 
#' @family designs
#' 
#' @seealso
#' \code{\link[desmon]{simon}}; \code{\link[desmon]{twostg}};
#' \code{\link[desmon]{twocon}}; \code{\link[desmon]{binci}}
#' 
#' @return
#' A character string describing the Simon design with
#' \code{attr(, "design")} containing the designs.
#' 
#' @examples
#' ## basic usage
#' p0 <- 0.3
#' pa <- 0.5
#' simon_text(p0, pa)
#' 
#' ## select a different design
#' simon_text(p0, pa, which = 2)
#' 
#' \dontrun{
#' ## use cat and/or strwrap for improved formatting and writing
#' cat(strwrap(simon_text(p0, pa), width = 80), sep = '\n')
#' 
#' cat(simon_text(p0, pa), file = '~/simon_design.txt')
#' 
#' cat(simon_text(p0, pa, outcome = 'overall response (OR) rate'),
#'     file = '~/simon_design2.txt')
#' }
#' 
#' @export

simon_text <- function(p0, pa, ..., conf = 0.95, which = 1L,
                       outcome = 'OUTCOME (OUT)') {
  if (!nzchar(outcome2 <- gsub('\\((.*?)\\)|.', '\\1', outcome)))
    outcome2 <- outcome
  
  call <- simon(p0, pa, ...)
  psep <- '\n\n'
  
  args <- call$designs
  if (which > nrow(args))
    which <- 1L
  args <- as.list(args[which, ])
  
  total <- args$n1 + args$n2
  
  onewid <- sapply(seq.int(args$n1), function(x) {
    ci <- binci(x, args$n1, conf = conf)
    unname(diff(ci[1:2]) / 2)
  })
  twowid <- sapply(seq(args$r1 + 1L, total), function(x) {
    ci <- twocon(args$n1, args$n2, args$r1, x, conf = conf)
    unname(diff(ci[1:2]) / 2)
  })
  
  p0stg <- twostg(args$n1, args$n2, p0, args$r1, args$r2)$prob
  pastg <- twostg(args$n1, args$n2, pa, args$r1, args$r2)$prob
  
  
  txt <- paste(
    'The study will use a Simon optimal two-stage design to allow for early',
    'termination for lack of efficacy. The', outcome, 'of at least',
    sprintf('%s%%', pa * 100), 'will be considered promising whereas',
    outcome2, 'of', paste0(p0 * 100, '%'),
    'or less will be considered non-promising.',
    
    psep,
    
    args$n1, 'eligible patients will be enrolled to stage one. If', args$r1,
    'or fewer', outcome2,
    'are observed, the regimen will be considered non-promising, and',
    'the study will stop early for lack of efficacy.',
    
    psep,
    
    'If', args$r1 + 1L, 'or more', outcome2, 'are observed within the first',
    args$n1, 'patients, an additional', args$n2,
    'patients will be entered to stage two for a total sample size of',
    total, 'patients.',
    
    'If at least', args$r2 + 1L, outcome2, 'are observed in', total,
    'patients, the study will be considered successful and the regimen',
    'worthy of further study. If the total number of responses observed is',
    args$r2, 'or fewer, the regimen will be considered non-promising.',
    
    psep,
    
    'The study has an overall power and type-I error of',
    sprintf('%.3f', 1 - args$type2), 'and', sprintf('%.3f', args$size),
    ', respectively. With a total size of', total,
    'patients, the two-stage exact 95% confidence interval for', outcome2,
    'will be no wider than +/-', sprintf('%s%%.', round(max(twowid) * 100)),
    
    psep,
    
    'If the true', outcome2, 'is', sprintf('%s%%,', p0 * 100),
    'the regimen will be considered non-promising with probability',
    sprintf('%.3f;', p0stg[1L]), 'and stop early with a probability of',
    sprintf('%.3f.', p0stg[2L]),
    
    'With the stage-one sample size of', args$n1, 'patients, the exact',
    sprintf('%s%%', conf * 100), 'confidence interval for', outcome2,
    'will be no wider than +/-', sprintf('%s%%.', round(max(onewid) * 100)),
    
    psep,
    
    'If the true', outcome2, 'is', sprintf('%s%%,', pa * 100),
    'the regimen will be considered non-promising with probability',
    sprintf('%.3f', pastg[1L]),
    '(type-II error) and stop early with a probability of',
    sprintf('%.3f.', pastg[2L])
  )
  
  structure(
    trim(txt, TRUE), design = call, class = 'stat_text'
  )
}

#' Two-stage design text
#' 
#' Output two-stage design text.
#' 
#' @param p0,pa the null and alternative hypotheses
#' @param n1,n2 sample size of first and second stage
#' @param r1,r2 maximum number of responses in first stage and overall where
#' treatment is declared ineffective
#' @param conf confidence level for single- and two-stage confidence intervals
#' @param outcome text string describing the outcome; if a string such as
#' \code{"long description (abbr)"} is given, then entire string is used in
#' the first instance, and only the text in parens is used subsequently
#' 
#' @family designs
#' 
#' @seealso
#' \code{\link{twostg_power}}; \code{\link[desmon]{simon}};
#' \code{\link[desmon]{twostg}}; \code{\link[desmon]{twocon}};
#' \code{\link[desmon]{binci}}
#' 
#' @return
#' A character string describing the two-stage design with
#' \code{attr(., "design")} containing the design.
#' 
#' @examples
#' ## basic usage (compare to a simon design)
#' p0 <- 0.3
#' pa <- 0.5
#' s <- as.list(desmon2:::simon(p0, pa)$designs[1L, ])
#' 
#' twostg_text(p0, pa, s$n1, s$n2, s$r1, s$r2)
#' 
#' \dontrun{
#' ## use cat and/or strwrap for improved formatting and writing
#' cat(strwrap(twostg_text(p0, pa, s$n1, s$n2, s$r1, s$r2), width = 80),
#'     sep = '\n')
#' 
#' cat(twostg_text(p0, pa, s$n1, s$n2, s$r1, s$r2),
#'     file = '~/simon_design.txt')
#' 
#' cat(twostg_text(p0, pa, s$n1, s$n2, s$r1, s$r2,
#'                 outcome = 'overall response (OR) rate'),
#'     file = '~/simon_design2.txt')
#' }
#' 
#' @export

twostg_text <- function(p0, pa, n1, n2, r1, r2, conf = 0.95,
                       outcome = 'OUTCOME (OUT)') {
  if (!nzchar(outcome2 <- gsub('\\((.*?)\\)|.', '\\1', outcome)))
    outcome2 <- outcome
  
  psep <- '\n\n'
  
  args <- twostg_power(p0, pa, n1, n2, r1, r2)
  
  p0stg <- abs(c(1, 0) - args[c('type1', 'Pr.stop1.H0')])
  pastg <- args[c('type2', 'Pr.stop1.H1')]
  
  args <- list(
    n1 = n1, r1 = r1, n2 = n2, r2 = r2,
    Pstop1.H0 = args[['Pr.stop1.H0']],
    size = args[['type1']], type2 = args[['type2']],
    E.tot.n.H0 = args[['E.tot.n.H0']]
  )
  
  total <- args$n1 + args$n2
  
  onewid <- sapply(seq.int(args$n1), function(x) {
    ci <- binci(x, args$n1, conf = conf)
    unname(diff(ci[1:2]) / 2)
  })
  twowid <- sapply(seq(args$r1 + 1L, total), function(x) {
    ci <- twocon(args$n1, args$n2, args$r1, x, conf = conf)
    unname(diff(ci[1:2]) / 2)
  })
  
  
  txt <- paste(
    'The study will use a two-stage design to allow for early',
    'termination for lack of efficacy. The', outcome, 'of at least',
    sprintf('%s%%', pa * 100), 'will be considered promising whereas',
    outcome2, 'of', paste0(p0 * 100, '%'),
    'or less will be considered non-promising.',
    
    psep,
    
    args$n1, 'eligible patients will be enrolled to stage one. If', args$r1,
    'or fewer', outcome2,
    'are observed, the regimen will be considered non-promising, and',
    'the study will stop early for lack of efficacy.',
    
    psep,
    
    'If', args$r1 + 1L, 'or more', outcome2, 'are observed within the first',
    args$n1, 'patients, an additional', args$n2,
    'patients will be entered to stage two for a total sample size of',
    total, 'patients.',
    
    'If at least', args$r2 + 1L, outcome2, 'are observed in', total,
    'patients, the study will be considered successful and the regimen',
    'worthy of further study. If the total number of responses observed is',
    args$r2, 'or fewer, the regimen will be considered non-promising.',
    
    psep,
    
    'The study has an overall power and type-I error of',
    sprintf('%.3f', 1 - args$type2), 'and', sprintf('%.3f', args$size),
    ', respectively. With a total size of', total,
    'patients, the two-stage exact 95% confidence interval for', outcome2,
    'will be no wider than +/-', sprintf('%s%%.', round(max(twowid) * 100)),
    
    psep,
    
    'If the true', outcome2, 'is', sprintf('%s%%,', p0 * 100),
    'the regimen will be considered non-promising with probability',
    sprintf('%.3f;', p0stg[1L]), 'and stop early with a probability of',
    sprintf('%.3f.', p0stg[2L]),
    
    'With the stage-one sample size of', args$n1, 'patients, the exact',
    sprintf('%s%%', conf * 100), 'confidence interval for', outcome2,
    'will be no wider than +/-', sprintf('%s%%.', round(max(onewid) * 100)),
    
    psep,
    
    'If the true', outcome2, 'is', sprintf('%s%%,', pa * 100),
    'the regimen will be considered non-promising with probability',
    sprintf('%.3f', pastg[1L]),
    '(type-II error) and stop early with a probability of',
    sprintf('%.3f.', pastg[2L])
  )
  
  structure(
    trim(txt, TRUE), design = args, class = 'stat_text'
  )
}

#' Two-proportion design text
#' 
#' Output two-proportion design text. The power or sample size (per group)
#' will be calculated depending on the arguments given (i.e., either
#' \code{n1} or \code{power} should be \code{NULL}).
#' 
#' @param p1,p2 the success rates for both groups (\code{p1 > p2})
#' @param n1,n2 sample size for both groups
#' @param power,alpha the power and one-sided significance level (type I error)
#' @param r proportion assigned to group 1
#' @param arms labels for the two groups
#' @param type the type of test used to calculate power, see
#' \code{\link[desmon]{b2p}}
#' @param cont.cor logical; if \code{TRUE} (default), the sample size will be
#' calculated for the continuity corrected statistic
#' 
#' @family designs
#' 
#' @seealso
#' \code{\link[desmon]{b2p}}; \code{\link[desmon]{b2n}}
#' 
#' @return
#' A character string describing the two-stage design with
#' \code{attr(., "design")} containing the design.
#' 
#' @examples
#' ## basic usage (compare to a simon design)
#' p1 <- 0.5
#' p2 <- 0.3
#' 
#' b2p_text(p1, p2, power = 0.8)
#' b2p_text(p1, p2, n1 = 103)
#' b2p_text(p1, p2, power = 0.8, r = 2/3)
#' 
#' \dontrun{
#' ## use cat and/or strwrap for improved formatting and writing
#' cat(strwrap(b2p_text(p1, p2, power = 0.8), width = 80), sep = '\n')
#' 
#' cat(b2p_text(p1, p2, power = 0.8), file = '~/simon_design.txt')
#' }
#' 
#' @export

b2p_text <- function(p1, p2, n1 = NULL, n2 = n1,
                     power = NULL, alpha = 0.025, r = 0.5,
                     arms = c('Experimental', 'Control'),
                     type = c('fisher', 'UMPU', 'approx.cor', 'approx.unc'),
                     cont.cor = TRUE) {
  if (is.null(n1) & is.null(power))
    stop('either \'n1\' and \'n2\' or \'power\' must be given')
  
  psep <- '\n\n'
  
  type <- match.arg(type)
  test <- switch(
    type,
    fisher = 'exact power for Fisher\'s exact test',
    UMPU = 'exact power for the uniformly most powerful unbiased test',
    approx.cor = 'normal approximation with continuity correction ',
    approx.unc = 'normal approximation without continuity correction'
  )
  
  r2 <- c(r, 1 - r)
  r2 <- r2 / min(r2)
  or <- p1 * (1 - p2) / (p2 * (1 - p1))
  
  args <- list(p1 = p1, p2 = p2, n1 = n1, n2 = n2, power = power, or = or,
               alpha = alpha, r = r, r2 = r2, type = type, test = test)
  
  if (is.null(power)) {
    args$power <- desmon::b2p(p1, p2, n1, n2, alpha, TRUE)[[type]]
  }
  if (is.null(n1)) {
    des <- desmon::b2n(p1, p2, power, r, alpha)[[(!cont.cor) + 1L]]
    args$n1 <- ceiling(r * des)
    args$n2 <- ceiling((1 - r) * des)
  }
  
  
  txt <- paste(
    'Patients will be randomly assigned', paste(r2, collapse = ':'),
    'to the', paste(arms, collapse = ' or '), 'arms.',
    'The success rates for are assumed to be', p1, 'and', p2,
    sprintf('(odds ratio of %.3f)', args$or), 'for the',
    paste(arms, collapse = ' and '), 'arms, respectively.',
    
    psep,
    
    'With', args$n1, 'and', args$n2, 'assigned to each arm, there will be',
    'at least', sprintf('%.3f', args$power), 'power',
    sprintf('(%s)', args$test), 'with a one-sided',
    sprintf('%.3f', args$alpha), 'significance level.',
    
    psep,
    
    'The overall sample size will be', args$n1 + args$n2, 'with', args$n1,
    'on the', arms[1L], 'arm and', args$n2, 'on the', arms[2L], 'arm.'
  )
  
  structure(
    trim(txt, TRUE), design = args, class = 'stat_text'
  )
}
