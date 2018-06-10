### stat text
# bin1samp_text, mtd_text, simon_text
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
#' \code{"long desc (abbr)"} is given, then entire string is used in the first
#' instance, and only the text in parens is used subsequently
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
    'The study will use a single-stage, exact binomial design. The',
    outcome,
    'of at least',
    sprintf('%s%%', pa * 100),
    'will be considered promising whereas',
    outcome2,
    'of',
    paste0(p0 * 100, '%'), 'or less will be considered non-promising.',
    
    psep,
    
    args$n,
    'eligible patients will be enrolled. If',
    args$r,
    'or fewer',
    outcome2,
    'are observed, the regimen will be considered non-promising, and',
    'the study will be unsuccessful. If at least',
    args$r + 1L,
    outcome2,
    'are observed in',
    args$n,
    'patients, the study will be considered successful and the regimen',
    'worthy of further study.',
    
    psep,
    
    'The study has an overall power and type-I error of',
    sprintf('%.3f', 1 - args$type2),
    'and',
    sprintf('%.3f', args$size),
    ', respectively. With a total size of',
    args$n,
    'patients, the single-stage exact 95% confidence interval for',
    outcome2,
    'will be no wider than +/-',
    sprintf('%s%%.', round(max(onewid) * 100)),
    
    psep,
    
    'If the true',
    outcome2,
    'is',
    sprintf('%s%%,', p0 * 100),
    'the regimen will be considered non-promising with probability',
    sprintf('%.3f;', 1 - args$size),
    'and promising with probability',
    sprintf('%.3f', args$size),
    '(type-I error).',
    
    psep,
    
    'If the true',
    outcome2,
    'is',
    sprintf('%s%%,', pa * 100),
    'the regimen will be considered non-promising with probability',
    sprintf('%.3f', args$type2),
    '(type-II error) and promising with probability',
    sprintf('%.3f.', 1 - args$type2)
  )
  
  
  structure(
    trim(txt),
    design = call,
    class = 'stat_text'
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
#' \code{\link{dlt_table}}
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
    
    'Three patients will be entered per dose level. If no DLTs occur in the',
    'first three patients, the dose is escalated to the next highest level',
    'for next cohort of three patients. If all dose levels are exhausted,',
    'then the highest dose will be considered the MTD.',
    
    psep,
    
    'If one DLT is experienced in the first three patients for a dose level,',
    'an additional three patients are entered at the same dose. If no additional',
    'DLTs occur at the current level, the dose will be escalated; if two or more',
    'DLTs occur in the six patients, the previous level is defined as the MTD,',
    'and the dose-escalation phase of the study will stop. No re-escalations',
    'will occur',
    
    psep,
    
    'If two or more DLTs occur at level 0, then the dose will be de-escalated to',
    'level -1 and proceed as above. Furthermore, if two DLTs are experienced at',
    'level -1, the study will be stopped, and no dose will be considered safe.',
    
    psep,
    
    'Once the MTD has been established, an additional expansion cohort of',
    expansion,
    'patients will be entered at this dose level to confirm safety.',
    
    psep,
    
    'The dose-finding phase of the study will enroll between 4 and',
    ndose * 6,
    'patients. With at least 3 + ',
    expansion,
    'treated at the MTD, the 95% exact binomial confidence interval for the',
    'observed rate of MTD will be no wider than +/-',
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
    trim(txt),
    table = tbl,
    class = 'stat_text'
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
#' A character string describing the Simon design with \code{attr(, "simon")}
#' containing the designs.
#' 
#' @examples
#' ## basic usage
#' simon_text(.3, .5)
#' 
#' ## select a different design
#' simon_text(.3, .5, which = 2)
#' 
#' \dontrun{
#' ## use cat and/or strwrap for improved formatting and writing
#' cat(strwrap(simon_text(.3, .5), width = 80), sep = '\n')
#' 
#' cat(simon_text(.3, .5), file = '~/simon_design.txt')
#' 
#' cat(simon_text(.3, .5, outcome = 'overall response (OR) rate'),
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
    'The study will use a Simon two-stage design to allow for early',
    'termination for lack of efficacy. The',
    outcome,
    'of at least',
    sprintf('%s%%', pa * 100),
    'will be considered promising whereas',
    outcome2,
    'of',
    paste0(p0 * 100, '%'), 'or less will be considered non-promising.',
    
    psep,
    
    args$n1,
    'eligible patients will be enrolled to stage one. If',
    args$r1,
    'or fewer',
    outcome2,
    'are observed, the regimen will be considered non-promising, and',
    'the study will stop early for lack of efficacy.',
    
    psep,
    
    'If',
    args$r1 + 1L,
    'or more',
    outcome2,
    'are observed within the first',
    args$n1,
    'patients, an additional',
    args$n2,
    'patients will be entered to stage two for a total sample size of',
    total,
    'patients.',
    'If at least',
    args$r2 + 1L,
    outcome2,
    'are observed in',
    total,
    'patients, the study will be considered successful and the regimen',
    'worthy of further study. If the total number of responses observed is',
    args$r2,
    'or fewer, the regimen will be considered non-promising.',
    
    psep,
    
    'The study has an overall power and type-I error of',
    sprintf('%.3f', 1 - args$type2),
    'and',
    sprintf('%.3f', args$size),
    ', respectively. With a total size of',
    total,
    'patients, the two-stage exact 95% confidence interval for',
    outcome2,
    'will be no wider than +/-',
    sprintf('%s%%.', round(max(twowid) * 100)),
    
    psep,
    
    'If the true',
    outcome2,
    'is',
    sprintf('%s%%,', p0 * 100),
    'the regimen will be considered non-promising with probability',
    sprintf('%.3f;', p0stg[1L]),
    'furthermore, the study will stop early with a probability of',
    sprintf('%.3f', p0stg[2L]),
    'if the true',
    outcome2,
    'is',
    sprintf('%s%%.', p0 * 100),
    'With the stage-one sample size of',
    args$n1,
    'patients, the exact',
    sprintf('%s%%', conf * 100),
    'confidence interval for',
    outcome2,
    'will be no wider than +/-',
    sprintf('%s%%.', round(max(onewid) * 100)),
    
    psep,
    
    'If the true',
    outcome2,
    'is',
    sprintf('%s%%,', pa * 100),
    'the regimen will be considered non-promising with probability',
    sprintf('%.3f', pastg[1L]),
    '(type-II error); furthermore, the study will stop early with a',
    'probability of',
    sprintf('%.3f', pastg[2L]),
    'if the true',
    outcome2,
    'is',
    sprintf('%s%%.', pa * 100)
  )
  
  
  structure(
    trim(txt),
    design = call,
    class = 'stat_text'
  )
}
