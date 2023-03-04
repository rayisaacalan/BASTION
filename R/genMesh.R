#' Generate a random mesh inside an SFC boundary
#'
#' @param sfc_bnd An SFC object, the boundary inside which the mesh should be generated
#' @param coords Optional set of interior coordinate points for the mesh to span. By default, n_ref points are uniformly sampled using st_sample
#' @param n_ref Number of interior mesh nodes to randomly generate. Larger numbers result in a finer mesh
#' @param ...
#'
#' @return An SF MULTIPOLYGON object
#' @export
#'
#' @examples
genMesh <- function(sfc_bnd,
                    coords = NULL,
                    n_ref = 100,
                    n_bnd = NULL,
                    ...) {
  # note the first and last boundary nodes are the same
  # coords are the reference points

  if (missing(coords)) {
    coords = sf::st_coordinates(sf::st_sample(sfc_bnd, size = n_ref, ...))[, 1:2]
  }
  bnd = sf::st_coordinates(sfc_bnd)[, 1:2]

  if(is.null(n_bnd)) {
    n_bnd = length(bnd[, 1]) - 1
  }

  coords_all = rbind(coords, bnd[1:n_bnd,])

  # get boundary segments
  segments = cbind((n_ref + 1):(n_ref + n_bnd), c((n_ref + 2):(n_ref + n_bnd), n_ref +
                                                    1))

  mesh = fdaPDE::create.mesh.2D(coords_all, segments = segments)  %>%
    fdaPDE::refine.mesh.2D(minimum_angle = 20,
                           verbosity = 0)
  mesh$n_int = n_ref  # number of interior nodes
  # edges that connect boundary nodes
  mesh$bnd_edges = apply(
    mesh$edges,
    1,
    FUN = function(x)
      any(x > n_ref)
  )

  meshlist = lapply(1:nrow(mesh$triangles),
                    function(x) {
                      tv <- mesh$triangles[x, , drop = TRUE]
                      sf::st_polygon(list(mesh$nodes[tv[c(1, 3, 2, 1)],
                                                     1:2,
                                                     drop = FALSE]))
                    })
  sf_mesh = sf::st_sfc(meshlist) %>% sf::st_cast('MULTIPOLYGON')
  return(sf_mesh)
}
