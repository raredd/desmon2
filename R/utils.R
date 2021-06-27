### utils
# unexported:
# trim, trim2
###


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
