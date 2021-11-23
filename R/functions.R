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
#' @return A vector of integers of length N with nclust unique integers which map each vertex to a cluster
#' @export
#'
#' @examples
constructClusters = function(graph, nclust) {

}

#' Perform a cluster birth operation (split an existing cluster)
#'
#' @param graph An object of class 'graph' from the igraph package
#' @param membership A vector of integers of length N with k unique integers (k < N) which map each vertex to a cluster
#' @param clust Integer, the cluster to split. Must be between 0 and k
#'
#' @return A vector of integers of length N with k + 1 unique integers which map each vertex to a cluster
#' @export
#'
#' @examples
graphBirth = function(graph, membership, clust) {

}
}
