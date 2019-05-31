### simulate
# sim3p3, simcrm
### 


#' Simulate 3+3 design
#' 
#' Simulate a standard 3+3 dose-escalation study.
#' 
#' @param probs vector of true rates of DLTs
#' @param expansion (optional) integer giving size of expansion cohort
#' 
#' @return
#' A list with the following components:
#' 
#' \item{\code{$data}}{dose by patient matrix: 0 - no DLT, 1 - DLT}
#' \item{\code{$mtd}}{numeric dose level selected}
#' \item{\code{$level}}{integer mtd level (ie, row of \code{data} selected)}
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
#' x <- replicate(1000L, sim3p3(pr)$level)
#' barplot(table(x), xlab = 'RP2D')
#' 
#' @export

sim3p3 <- function(probs, expansion = 10) {
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
  
  idx <- 1L
  while ((length(probs) + 1L) > idx) {
    ## dlt indicators for all pts
    dlt <- rbinom(npl, 1L, probs[idx])
    
    ## rules based on first three entered
    lvl <- sum(dlt[1:3])
    if (lvl > 1) {
      ## 2 or 3 dlt -- stop, prev level is mtd
      mat[idx, 1:3] <- dlt[1:3]
      idx <- idx - 1L
      break
    } else if (lvl == 0) {
      ## 0 dlt in 3 -- escalate
      mat[idx, 1:3] <- dlt[1:3]
      idx <- idx + 1L
      next
    } else if (lvl == 1) {
      ## 1 dlt in 3 -- add 3
      mat[idx, 1:6] <- dlt
      if (sum(dlt) > 1) {
        ## 1 or more dlt in last 3 -- prev level is mtd
        idx <- idx - 1L
        break
      } else {
        ## no more dlt in last 3 -- escalate
        idx <- idx + 1L
        next
      }
    }
  }
  
  pr_all <- pr_exp <- NA
  idx <- pmin(idx, length(probs))
  
  ## optional - add expansion cohort at safe dose
  if (idx > 0) {
    expansion <- rbinom(expansion, 1L, probs[idx])
    mat[idx, -(1:6)] <- expansion
    ## observed probs for all treated at dose and expansion cohort
    pr_all <- mean(mat[idx, ], na.rm = TRUE)
    if (length(expansion))
      pr_exp <- mean(expansion)
  }
  
  list(
    data = mat, mtd = probs[idx], level = idx, expansion = pr_exp,
    cohort = pr_all, enrolled = length(sort(mat))
  )
}

#' Simulate CRM design
#' 
#' Simulate a \link[=simCrm]{CRM} dose-escalation study.
#' 
#' @param probs vector of true rates of DLTs
#' @param expansion (optional) integer giving size of expansion cohort
#' @param target target AE rate
#' @param nptmax max number of patients used on study (not including the
#' expansion cohort)
#' @param ... additional arguments passed to \code{\link[UBCRM]{simCrm}}
#' 
#' @return
#' A list with following components:
#' 
#' \item{\code{$data}}{\code{UBCRM::simCrm()$data}}
#' \item{\code{$mtd}}{numeric dose level selected}
#' \item{\code{$level}}{integer mtd level (ie, row of \code{data} selected)}
#' \item{\code{$expansion}}{observed rate of DLT in expansion cohort}
#' \item{\code{$cohort}}{observed rate of DLT in all patients}
#' \item{\code{$enrolled}}{total number of patients enrolled on study}
#' 
#' @seealso
#' \code{\link{simp3p}}
#' 
#' @examples
#' pr <- 1:5 / 10
#' simcrm(pr)
#' 
#' x <- replicate(1000L, simcrm(pr)$level)
#' barplot(table(x), xlab = 'RP2D')
#' 
#' @export

simcrm <- function(probs, expansion = 0L, target = 1/3,
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
  
  list(data = crm$data, mtd = pr[crm$mtd], level = crm$mtd,
       expansion = exp, cohort = mean(coh), enrolled = size)
}
