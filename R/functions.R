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
