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
#' @return A vector of integers of length N with nclust unique integers
#' @export
#'
#' @examples
constructClusters = function(graph, nclust) {

}

}
