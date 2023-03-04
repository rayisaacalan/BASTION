#' Get Betweenness of Edges in Clustered Graph
#'
#' \code{edgeBetweenClust} is a utility function which takes in a graph and the membership vector
#' mapping each vertex to a certain cluster and returns a vector of boolean values. A value of 'TRUE'
#' means that the corresponding edge in that edge list position connects vertices which belong to
#' distinct clusters. A value of 'FALSE' means that the edge connects vertices which belong to the same
#' cluster.
#'
#' @param graph An object of class 'graph' from the \code{\link[igraph]{igraph}} package.
#' @param membership A vector of integers of length N (see \code{\link{constructClusters}})
#'
#' @return A vector of boolean values of length \code{ecount(graph)} where each TRUE value represents an edge that connects distinct clusters
#' @export
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
#' g = constructGraph(coords, 5)
#' clust_membership = constructClusters(g, 6, minclust = 5)$membership
#' betweenness = edgeBetweenClust(g, clust_membership)
#' plot(g,
#'      layout = as.matrix(coords),
#'      edge.color = as.numeric(betweenness) + 1,
#'      edge.arrow.mode = 0)
edgeBetweenClust = function(graph, membership) {
  # Get a matrix of all edges in graph by numeric id (column 1 is source, column 2 is destination)
  edge_matrix = igraph::as_edgelist(graph, names = F)
  # Get the vector of cluster membership for all the edge sources
  membership_source = membership[edge_matrix[, 1]]
  # Get the vector of cluster membership for all the edge destinations
  membership_dest = membership[edge_matrix[, 2]]
  # Initialize all edges as being within the same cluster (between clusters is false)
  edge_status = rep(FALSE, igraph::ecount(graph))
  # If the cluster membership between the edge source & destination is not the same, between clusters is true
  edge_status[membership_source != membership_dest] = TRUE
  # Return betweenness of edges
  return(edge_status)
}
