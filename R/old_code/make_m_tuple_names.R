# TODO: Don't export this function but have it available internally
#' Make m-tuple names
#'
#' This helper function constructs m-tuple names in the correct (alphabetical) order for a given value of m
#' @param m The size of the m-tuple. Must be an int.
#' @keywords helper
#' @export
#' @examples
#' make_m_tuple_names(1L)
#' make_m_tuple_names(2L)
#' make_m_tuple_names(3L)
#' @return A character vector
make_m_tuple_names <- function(m){
  if (!is.integer(m) || m < 1){
    stop("'m' must be an int. 'm' must be greater than 0.")
  }
  sort(do.call(paste0, expand.grid(lapply(seq_len(m), function(x){c('M', 'U')}))))
}