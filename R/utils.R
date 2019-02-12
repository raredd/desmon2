### utils
# 
# 
# unexported:
# trim
###


trim <- function(x, dig = FALSE) {
  if (dig)
    x <- gsub('(?<=\\d)\\s+(?=[,.;])', '', x, perl = TRUE)
  gsub('(?<=\\v)\\h+|\\h+\\z|\\h+$|\\h+(?=\\v)', '', x, perl = TRUE)
}
