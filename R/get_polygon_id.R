#' Get polygon ids of a data set
#'
#' @param point_in_polygon The output of st_contains on a mesh and data set
#'
#' @return A data frame mapping each data point to the id of the polygon in which it is contained
#' @export
#'
#' @examples
get_polygon_id = function(point_in_polygon) {
  (data.frame(value = unlist(point_in_polygon),
              idx = rep(seq_along(point_in_polygon),
                        lapply(point_in_polygon, length))) %>%
     dplyr::arrange(value))$idx
}
