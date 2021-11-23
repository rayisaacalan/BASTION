#' Construct a Graph Object from Coordinate Data using K-Nearest Neighbors
#'
#' @param coordinates Data frame whose first 2 columns are numeric coordinate data
#' @param k Integer, the parameter for the K-Nearest Neighbors algorithm
#'
#' @return An object of class 'graph' from the igraph package
#' @export
#'
#' @examples
constructGraph = function(coordinates, k) {
# Create a k nearest neighbor graph using first two columns of input
  # Double check that input is of appropriate dimension
  if(ncol(coordinates) < 2) {
    stop("Function requires at least 2 columns of coordinate data")
  }
  # Extract only the relevant data
  coord_data = coordinates[ , 1:2]
  if((class(coord_data[, 1]) != "numeric") || (class(coord_data[, 2]) != "numeric")) {
    stop("First 2 columns must be of class numeric")
  }
  # Construct k nearest neighbor graph
  graph = cccd::nng(coord_data, k = k)
# Measure euclidean distance for each edge (TODO: allow for output of custom distance function)
  # Get a 2 column matrix where each row is an edge from first column vertex to second column vertex
  edge_matrix = ends(graph, E(graph))
  # Calculate Euclidean distance for each edge
  edge_distances = sqrt(rowSums(coord_data[edge_matrix[ , 1], ] - coord_data[edge_matrix[ , 2], ]^2))
  # Assign the distance for each edge to the graph's edge weights
  E(graph)$weight = edge_distances
  # Return the constructed graph
  return(graph)
}

#' Construct a cluster membership list for a graph
#'
#' @param graph An object of class 'graph' from the igraph package
#' @param nclust An integer, the number of different clusters to assign points to. Must be at most N (the number of vertices in graph)
#'
#' @return A list containing two elements:
#' 'graph': The input graph with inter-cluster edges inactive
#' 'membership': A vector of integers of length N with nclust unique integers which map each vertex to a cluster
#' @export
#'
#' @examples
constructClusters = function(graph, nclust) {

}

#' Perform a cluster birth operation (split an existing cluster)
#'
#' @param graph An object of class 'graph' from the igraph package
#' @param membership A vector of integers of length N with k unique integers (k < N) which map each vertex to a cluster
#' @param clust Integer, the cluster to split. Must be between 0 and k - 1
#'
#' @return A list containing two elements:
#' 'graph': The input graph with 1 fewer active edge
#' 'membership': A vector of integers of length N with k + 1 unique integers which map each vertex to a cluster
#' @export
#'
#' @examples
graphBirth = function(graph, membership, clust) {

}

#' Perform a cluster death operation (merge an existing cluster)
#'
#' @param graph An object of class 'graph' from the igraph package
#' @param membership A vector of integers of length N with k unique integers (1 < k <= N) which map each vertex to a cluster
#' @param clust Integer, the cluster to split. Must be between 0 and k - 1
#'
#' @return A list containing two elements:
#' 'graph': The input graph with 1 additional active edge
#' 'membership': A vector of integers of length N with k - 1 unique integers which map each vertex to a cluster
#' @export
#'
#' @examples
graphDeath = function(graph, membership, clust) {

}

#' Perform a cluster death operation followed by a cluster birth operation
#'
#' @param graph An object of class 'graph' from the igraph package
#' @param membership A vector of integers of length N with k unique integers (1 < k <= N) which map each vertex to a cluster
#'
#' @return A list containing two elements:
#' 'graph': The input graph with the different set of active edges
#' 'membership': A vector of integers of length N with k unique integers which map each vertex to a cluster, likely different than input membership
#' @export
#'
#' @examples
graphChange = function(graph, membership) {

}

#' Generate new spanning trees for each cluster
#'
#' @param graph An object of class 'graph' from the igraph package
#' @param membership  A vector of integers of length N with k unique integers (1 < k <= N) which map each vertex to a cluster
#'
#' @return A list containing two elements:
#' 'graph': The input graph with each cluster having a new minimum spanning tree based on resampled edge weights
#' 'membership': A vector of integers of length N with k unique integers which map each vertex to a cluster
#' @export
#'
#' @examples
graphHyper = function(graph, membership) {

}
