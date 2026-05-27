### ct.gv utils
# ctae, ctae_table, ctae_report
###


#' Clinical trials adverse events
#' 
#' Summarize and format adverse events for ct.gov/FDA reporting.
#' 
#' @param data a data frame with adverse events
#' @param id,toxdesc,toxcat,comment column names with patient IDs, toxicity
#'   descriptions/terms, toxicity categories/organ systems, and comments/additional
#'   details (e.g., clarifications for other/specify)
#' @param sae logical; if \code{TRUE}, result is formatted for reporting SAEs;
#'   otherwise, non-SAEs are assumed
#' @param atrisk number of patients at-risk (for each stratum)
#' @param strata optional column name defining groups of patients; note that
#'   for \code{ctae_report}, each unique level of \code{strata} should have
#'   an \code{atrisk} value
#' 
#' @examples
#' dd <- data.frame(
#'   casenum = rep(1:2, c(3, 5)),
#'   toxdesc = c('aa', 'aa', 'b', 'aa', 'aaa', 'b', 'bb', 'c'),
#'   toxcat = c('A', 'A', 'B', 'A', 'A', 'B', 'B', 'C')
#' )
#' ae <- ctae(dd)
#' ae
#' c1 <- ctae_table(ae, FALSE, 20)
#' c1
#' ## alternatively
#' c2 <- ctae_report(dd, sae = FALSE, atrisk = 20)
#' c2
#' identical(c1, c2)
#' 
#' ## stratify table to report by cohort
#' dd$group <- factor(dd$casenum, 1:2, c('Arm A', 'Arm B'))
#' ctae_report(dd, sae = FALSE, atrisk = c(10, 10), strata = 'group')
#' 
#' @export

ctae <- function(data, id = 'casenum', toxdesc = 'toxdesc', toxcat = 'toxcat',
                 comment = NULL) {
  p0 <- function(x) {
    x <- type.convert(trimws(x), as.is = TRUE, na.strings = c('', 'NA', '.'))
    t <- table(x)
    toString(sprintf('%s (n=%s)', names(t), t))
  }
  codes <- names(sort(table(data[, toxdesc]), decreasing = TRUE))
  
  dd <- lapply(codes, function(x) {
    d <- data[data[, toxdesc] %in% x, ]
    o <- if (!is.null(comment) && grepl('(?i)other|specify', x))
      p0(d[, comment]) else NA
    data.frame(
      term = x,
      system = unique(d[, toxcat]),
      npatient = length(unique(d[, id])),
      nevent = nrow(d),
      other = o
    )
  })
  dd <- do.call('rbind', dd)
  
  structure(dd, class = c('ctae', class(dd)))
}

#' @rdname ctae
#' @export
ctae_table <- function(x, sae, atrisk) {
  stopifnot(
    inherits(x, 'ctae'),
    is.logical(sae),
    length(sae) == 1L,
    is.numeric(atrisk),
    length(atrisk) == 1L,
    all(atrisk >= x$npatient)
  )
  
  data.frame(
    adverseEventType = if (sae) 'Serious' else 'Other',
    assessmentType = 'Systematic Assessment',
    additionalDescription = x$other,
    organSystemName = x$system,
    sourceVocabulary = NA,
    term = x$term,
    'trunc-arm-name{numEvents}' = x$nevent,
    'trunc-arm-name{numSubjectsAffected}' = x$npatient,
    'trunc-arm-name{numSubjectsAtRisk}' = atrisk,
    check.names = FALSE
  )
}

#' @rdname ctae
#' @export
ctae_report <- function(data, id = 'casenum', toxdesc = 'toxdesc', toxcat = 'toxcat',
                        comment = NULL, sae, atrisk, strata = NULL) {
  data$...strata <- if (is.null(strata))
    rep_len(factor('trunc-arm-name'), nrow(data)) else factor(data[, strata])
  spl <- split(data, data$...strata)
  stopifnot(
    'need number at-risk for each group' = length(atrisk) == length(spl)
  )
  res <- lapply(seq_along(spl), function(ii) {
    x <- spl[[ii]]
    n <- atrisk[ii]
    l <- unique(as.character(x$...strata))
    ct <- ctae(x, id, toxdesc, toxcat, comment)
    ct <- ctae_table(ct, sae, n)
    names(ct) <- gsub('.*(?=\\{)', l, names(ct), perl = TRUE)
    ct
  })
  
  ii <- grepl('\\{num', names(res[[1L]]))
  nn <- names(res[[1L]])[!ii]
  mm <- Reduce(function(xx, yy) merge(xx, yy, by = nn, all = TRUE), res)
  ii <- grepl('\\{num', names(mm))
  mm[, ii][is.na(mm[, ii])] <- 0
  mm
}
