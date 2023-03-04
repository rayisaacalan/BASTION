#' Remove inter-cluster edges from a graph
#'
#' \code{clusterGraph} takes in a graph and the membership vector mapping each vertex to
#' a certain cluster and returns the input graph with any edges that connect vertices belonging
#' to distinct clusters removed (see \code{\link{edgeBetweenClust}}).
#'
#' @inheritParams edgeBetweenClust
#'
#' @return An object of class 'graph' from the \code{\link[igraph]{igraph}} package, the input graph with inter-cluster edges removed
#' @export
#'
#' @references
#' Luo, Z.T. (*), Sang, H. and Mallick, B.K. (2021), BAST: Bayesian Additive Regression Spanning Trees
#' for Complex Constrained Domain
#'
#' Luo, Z.T. (*), Sang, H. and Mallick, B.K. (2021), A Bayesian Contiguous Partitioning Method for
#' Learning Clustered Latent Variables, Journal of Machine Learning Research, 22, 1-52.
#'
#' @examples
#' set.seed(1)
#' coords = data.frame(lon = rnorm(50), lat = rnorm(50))
#' g = constructGraph(coords, 4)
#' clust_membership = constructClusters(g, 5, minclust = 3)$membership
#' clust_g = clusterGraph(g, clust_membership)
#' plot(clust_g, layout = as.matrix(coords), vertex.color = clust_membership, edge.arrow.mode = 0)
clusterGraph = function(graph, membership) {
  # Get vector of cluster betweenness
  between_clusters = edgeBetweenClust(graph, membership)
  # Get ids of edges which are between clusters
  edge_ids_between = which(between_clusters)
  # Delete those edges
  clustered_graph = igraph::delete.edges(graph, edge_ids_between)
  # Return clustered graph
  return(clustered_graph)
}
