### ct.gv utils
# ctae, ctae_table
###


#' Clinical trials adverse events
#' 
#' Summarize and format adverse events for ct.gov/FDA reporting.
#' 
#' @param data a data frame with adverse events
#' @param id,toxdesc,toxcat,comment column names with patient IDs, toxicity
#'   descriptions/terms, toxicity categories/organ systems, and comments/additional
#'   details (e.g., clarifications for other/specify)
#' @param sae logical; if \code{TRUE}, result is formatted for reporting SAEs
#' @param atrisk number of patients at-risk
#' 
#' @examples
#' dd <- data.frame(
#'   casenum = rep(1:2, c(3, 5)),
#'   toxdesc = c('aa', 'aa', 'b', 'aa', 'aaa', 'b', 'bb', 'c'),
#'   toxcat = c('A', 'A', 'B', 'A', 'A', 'B', 'B', 'C')
#' )
#' ae <- ctae(dd)
#' ae
#' ctae_table(ae, FALSE, 20)
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
    'Disease Progres{numEvents}' = x$nevent,
    'Disease Progres{numSubjectsAffected}' = x$npatient,
    'Disease Progres{numSubjectsAtRisk}' = atrisk,
    check.names = FALSE
  )
}
