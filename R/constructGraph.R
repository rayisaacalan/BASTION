#' Construct a Graph Object from Coordinate Data using K-Nearest Neighbors
#'
#' \code{constructGraph} takes numeric data (currently it requires 2 dimensional data such as spatial data,
#' this functionality will expand in a future update) and uses the K-Nearest Neighbors algorithm to
#' construct a weighted graph object from \code{\link[igraph]{igraph}}. The K-Nearest Neighbors implementation is from
#' \code{\link[cccd]{cccd}}.
#'
#' Note: this function is not guaranteed to return a connected graph, and connectivity
#' is highly recommended for graphs to be used in other functions. If the output graph is not connected,
#' it is recommended to add more data points, increase k, or construct a graph manually.
#'
#' @param coordinates Data frame whose first 2 columns are numeric coordinate data
#' @param k Integer, the parameter for the K-Nearest Neighbors algorithm
#'
#' @return An object of class 'graph' from the \code{\link[igraph]{igraph}} package
#' @export
#'
#' @references
#' Csardi G, Nepusz T: The igraph software package for complex network research,
#' InterJournal, Complex Systems 1695. 2006. https://igraph.org
#'
#' David J. Marchette (2015). cccd: Class Cover Catch Digraphs. R package version 1.5.
#' https://CRAN.R-project.org/package=cccd
#'
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
#' plot(g, layout = as.matrix(coords), edge.arrow.mode = 0)
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
  edge_matrix = igraph::ends(graph, igraph::E(graph))
  # Calculate Euclidean distance for each edge
  edge_distances = sqrt(rowSums((coord_data[edge_matrix[ , 1], ] - coord_data[edge_matrix[ , 2], ])^2))
  # Assign the distance for each edge to the graph's edge weights
  igraph::E(graph)$weight = edge_distances
  # Return the constructed graph
  return(graph)
}

