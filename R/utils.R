### utils
# 
# 
# unexported:
# trim
###


trim <- function(x) {
  gsub('(?<=\\v)\\h+|\\h+\\z|\\h+$|\\h+(?=\\v)', '', x, perl = TRUE)
}
