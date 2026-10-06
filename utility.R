

############################################################
# utility.R
#
# other functions used across multiple files
############################################################


get_first_non_missing <- function(x) {
  idx <- which(!is.na(x))[1]
  if (length(idx) == 0) return(0)
  x[idx]
}


