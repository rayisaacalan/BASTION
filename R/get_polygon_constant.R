#' Get the piecewise constant value of each polygon
#'
#' @param membership A vector with length equal to the number of polygons, the output of igraph's components function
#' @param mus A vector with length equal to the number of distinct clusters, containing the numeric piecewise constant for each cluster
#'
#' @return A vector with length equal to the number of polygons, with each value being the constant of the cluster the polygon belongs to
#' @export
#'
#' @examples
get_polygon_constant = function(membership, mus) {
  value_map = data.frame(cluster = seq(1, length(mus)),
                         mu = mus)
  value_map$mu[match(membership, value_map$cluster)]
}
