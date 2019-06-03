### utils
# 
# 
# unexported:
# trim, trim2
###


trim <- function(x, dig = FALSE) {
  if (dig)
    x <- gsub('(?<=\\d)\\s+(?=[,.;])', '', x, perl = TRUE)
  gsub('(?<=\\v)\\h+|\\h+\\z|\\h+$|\\h+(?=\\v)', '', x, perl = TRUE)
}

trim2 <- function(x) {
  ## trim 2 or more horizontal white space to one
  gsub('\\h{2,}', ' ', x, perl = TRUE)
}
