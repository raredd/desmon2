### simulate phase 1 designs
# sim3p3, simcrm
# 
# unexported:
# sim3p3_
### 


#' Simulate 3+3 design
#' 
#' Simulate a standard 3+3 dose (de-)escalation study.
#' 
#' @param probs vector of true rates of DLTs
#' @param expansion (optional) integer giving size of expansion cohort
#' @param escalation logical; if \code{TRUE}, dose-escalation rules are used;
#'   otherwise, de-escalation rules are used with decreasing \code{probs}
#' @param n.max maximum number of patients available for dose-finding and
#'   expansion cohorts combined
#' 
#' @return
#' A list with the following components:
#' 
#' \item{\code{$data}}{dose by patient matrix: 0 - no DLT, 1 - DLT}
#' \item{\code{$summary}}{total number of DTLs, patients, and probability of
#'   DLT per dose level for \strong{esc}alation and \strong{tot}al}
#' \item{\code{$level}}{mtd level (i.e., row of \code{data} selected)}
#' \item{\code{$mtd}}{DLT rate of dose level selected}
#' \item{\code{$expansion}}{observed rate of DLT in expansion cohort}
#' \item{\code{$cohort}}{observed rate of DLT in all patients}
#' \item{\code{$enrolled}}{total number of patients enrolled on study}
#' 
#' @seealso
#' \code{\link{simcrm}}
#' 
#' @examples
#' ## standard 3+3 with 10 patient expansion cohort
#' pr <- 1:5 / 10
#' sim3p3(pr)
#' 
#' 
#' ## run 1000 3+3 studies
#' sim1 <- replicate(1000L, sim3p3(pr), simplify = FALSE)
#' 
#' op <- par(mfrow = c(1, 2), las = 1L)
#' out <- sapply(sim1, function(x) c(x$level, x$enrolled))
#' barplot(
#'   prop.table(table(out[1, ])),
#'   main = 'Dose selected',
#'   xlab = 'RP2D',
#'   ylab = 'Proportion selected'
#' )
#' boxplot(
#'   split(out[2, ], out[1, ]),
#'   ylim = c(0, length(pr) * 6 + 10),
#'   main = 'N used to select dose x',
#'   xlab = 'Dose selected\n(0 - no safe dose)',
#'   ylab = 'Total size'
#' )
#' par(op)
#' 
#' 
#' ## 3+3 de-escalation, starting at highest dose
#' ## run 1000 3+3 studies
#' pr2 <- rev(pr)
#' sim2 <- replicate(1000L, sim3p3(pr2, escalation = FALSE), simplify = FALSE)
#' 
#' op <- par(mfrow = c(1, 2), las = 1L)
#' out <- sapply(sim2, function(x) c(x$level, x$enrolled))
#' barplot(
#'   prop.table(table(out[1, ])),
#'   main = 'Dose selected',
#'   xlab = 'RP2D',
#'   ylab = 'Proportion selected'
#' )
#' boxplot(
#'   split(out[2, ], out[1, ]),
#'   ylim = c(0, length(pr) * 6 + 10),
#'   main = 'N used to select dose x',
#'   xlab = 'Dose selected\n(0 - no safe dose)',
#'   ylab = 'Total size'
#' )
#' par(op)
#' 
#' @export

sim3p3 <- function(probs, expansion = 10, escalation = TRUE, n.max = Inf) {
  if (escalation && !all(diff(order(probs)) == 1L)) {
    warning('if escalation = TRUE, \'probs\' should be increasing - sorting')
  } else if (!escalation && !all(diff(order(probs)) == -1L)) {
    warning('if escalation = FALSE, \'probs\' should be decreasing - sorting')
  }
  
  probs <- sort(probs, decreasing = !escalation)
  
  ## groups of 3 entered
  grp <- c(3L, 3L)
  npl <- sum(grp)
  
  ## initiate results matrix
  if (expansion < 0)
    expansion <- 0L
  dnn <- list(
    paste('dose', seq_along(probs)),
    paste(rep(c('grp1', 'grp2', 'exp'), c(grp, expansion)),
          sequence(c(grp, expansion)))
  )
  mat <- matrix(NA, length(probs), sum(c(npl, expansion)), dimnames = dnn)
  
  res <- sim3p3_(probs, !escalation, mat, 1L)
  mat <- res$mat
  idx <- res$idx
  
  f <- function(matrix, n.max = Inf) {
    if (is.infinite(n.max))
      return(matrix)
    ii <- t(matrix(cumsum(t(+!is.na(matrix))), nrow(matrix), byrow = TRUE))
    ii <- ii <= n.max & matrix(!duplicated(c(ii)), nrow(ii))
    matrix <- t(matrix)
    matrix[!ii] <- NA
    t(matrix)
  }
  
  ## update tox and dose level based on only n.max patients
  if (is.finite(n.max)) {
    mat[, 1:6] <- f(mat[, 1:6], n.max)
    idx <- rowSums(mat, na.rm = TRUE) < 2 & rowSums(!is.na(mat)) > 0
    idx <- max(cumsum(idx))
    expansion <- pmin(expansion, pmax(0, n.max - sum(!is.na(mat))))
  }
  
  pr_all <- pr_exp <- NA
  idx <- pmin(idx, length(probs))
  
  ## optional - add expansion cohort at safe dose
  if (idx > 0) {
    expansion <- rbinom(expansion, 1L, probs[idx])
    len <- ncol(mat) - length(expansion) - 6L
    mat[idx, -(1:6)] <- c(expansion, rep_len(NA, len))
    ## observed probs for all treated at dose and expansion cohort
    pr_all <- mean(mat[idx, ], na.rm = TRUE)
    if (length(expansion))
      pr_exp <- mean(expansion)
  }
  
  sum <- data.frame(
    level = seq_along(probs),
    dlt_esc = rowSums(mat[, 1:6], na.rm = TRUE),
    npt_esc = rowSums(!is.na(mat[, 1:6])),
    prob_esc = NA,
    dlt_tot = rowSums(mat, na.rm = TRUE),
    npt_tot = rowSums(!is.na(mat)),
    prob_tot = NA
  )
  sum <- within(sum, {
    prob_esc <- replace(dlt_esc / npt_esc, -idx, NA)
    prob_tot <- replace(dlt_tot / npt_tot, -idx, NA)
  })
  
  l <- function(x, label) {
    structure(x, label = label)
  }
  
  list(
    data = mat, summary = sum,
    level = l(idx, 'Level selected as MTD'),
    mtd = l(probs[idx], 'Probability of DLT in level selected'),
    expansion = l(pr_exp, 'DLT rate in expansion cohort'),
    cohort = l(pr_all, 'DLT rate in MTD cohort'),
    enrolled = l(length(sort(mat)), 'Total enrolled on study')
  )
}

sim3p3_ <- function(probs, d, mat, idx = 1L) {
  ## idx: starting level
  
  while ((length(probs) + 1L) > idx) {
    ## dlt indicators for all pts
    dlt <- rbinom(6L, 1L, probs[idx])
    
    ## rules based on first three entered
    lvl <- sum(dlt[1:3])
    if (lvl > 1) {
      ## 2 or 3 dlt
      ## d: de-escalate to next level
      ## e: stop - prev level is mtd
      mat[idx, 1:3] <- dlt[1:3]
      idx <- idx + if (d) 1L else -1L
      if (d)
        next else break
    } else if (lvl == 0) {
      ## 0 dlt in 3
      ## d: stop - current level is mtd
      ## e: escalate to next level
      mat[idx, 1:3] <- dlt[1:3]
      idx <- if (d)
        idx else idx + 1L
      if (d)
        break else next
    } else if (lvl == 1) {
      ## 1 dlt in 3 -- add 3
      mat[idx, 1:6] <- dlt
      if (sum(dlt) > 1) {
        ## 1 or more dlt in last 3
        #E d: de-escalate to next level
        ## e: stop - prev level is mtd
        idx <- idx + if (d) 1L else -1L
        if (d)
          next else break
      } else {
        ## no more dlt in last 3
        ## d: stop - current level is mtd
        ## e: escalate to next level
        idx <- if (d)
          idx else idx + 1L
        if (d)
          break else next
      }
    }
  }
  
  list(mat = mat, idx = idx)
}

#' Simulate CRM design
#' 
#' Simulate a \link[=simCrm]{CRM} dose-escalation study.
#' 
#' @param probs vector of true rates of DLTs
#' @param expansion (optional) integer giving size of expansion cohort
#' @param target target AE rate
#' @param nptmax max number of patients used on study (not including the
#'   expansion cohort)
#' @param ... additional arguments passed to \code{\link[UBCRM]{simCrm}}
#' 
#' @return
#' A list with following components:
#' 
#' \item{\code{$data}}{\code{UBCRM::simCrm()$data}}
#' \item{\code{$level}}{mtd level (i.e., row of \code{data} selected)}
#' \item{\code{$mtd}}{DLT rate of dose level selected}
#' \item{\code{$expansion}}{observed rate of DLT in expansion cohort}
#' \item{\code{$cohort}}{observed rate of DLT in all patients}
#' \item{\code{$enrolled}}{total number of patients enrolled on study}
#' 
#' @seealso
#' \code{\link{sim3p3}}
#' 
#' @examples
#' pr <- 1:5 / 10
#' simcrm(pr)
#' 
#' 
#' ## run 100 crm studies
#' sim1 <- replicate(100L, simcrm(pr), simplify = FALSE)
#' 
#' op <- par(mfrow = c(1, 2), las = 1L)
#' out <- sapply(sim1, function(x) c(x$level, x$enrolled))
#' barplot(
#'   prop.table(table(out[1, ])),
#'   main = 'Dose selected',
#'   xlab = 'RP2D',
#'   ylab = 'Proportion selected'
#' )
#' boxplot(
#'   split(out[2, ], out[1, ]),
#'   ylim = c(0, length(pr) * 6 + 10),
#'   main = 'N used to select dose x',
#'   xlab = 'Dose selected\n(0 - no safe dose)',
#'   ylab = 'Total size'
#' )
#' par(op)
#' 
#' @export

simcrm <- function(probs, expansion = 0L, target = 1 / 6,
                   nptmax = length(probs) * 6L, ...) {
  lvl <- paste('dose', seq_along(probs))
  crm <- UBCRM::simCrm(
    probs, firstdose = 1L, target = target, nptmax = nptmax, ...
  )
  if (!sum(crm$mtd %in% crm$dose))
    crm$mtd <- max(crm$dose)
  
  size <- expansion + sum(crm$data$npt)
  
  ## posterior dlt prob
  coh <- expansion + sum(crm$dose %in% crm$mtd) * 3
  coh <- rbinom(coh, 1L, crm$prob[crm$mtd])
  
  ## dlt rate in expansion cohort
  exp <- if (expansion > 0)
    mean(rbinom(expansion, 1L, crm$prob[crm$mtd])) else NA
  
  list(data = crm$data, mtd = probs[crm$mtd], level = crm$mtd,
       expansion = exp, cohort = mean(coh), enrolled = size)
}
