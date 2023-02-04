### utils
# unexported:
# %||%, iprint, trim, trim2
###


`%||%` <- function(a, b) {
  if (!is.null(a))
    a else b
}

iprint <- function(x, copula = ' and ') {
  # iprint(1:2); iprint(1:4); iprint(1)
  if (length(x) == 2L)
    paste(x, collapse = copula)
  else if (length(x) > 2L)
    sprintf('%s,%s %s', toString(head(x, -1L)), copula, tail(x, 1L))
  else as.character(x)
}

trim <- function(x, dig = FALSE) {
  if (dig) {
    ## remove white space between digit and ,.;
    x <- gsub('(?<=\\d)\\s+(?=[,.;])', '', x, perl = TRUE)
  }
  ## remove horizontal white space
  gsub('(?<=\\v)\\h+|\\h+\\z|\\h+$|\\h+(?=\\v)', '', x, perl = TRUE)
}

trim2 <- function(x) {
  ## trim 2 or more horizontal white spaces to one space
  gsub('\\h{2,}', ' ', x, perl = TRUE)
}
