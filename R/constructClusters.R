#' Construct a cluster membership list for a graph
#'
#' \code{constructClusters} takes in a weighted and connected graph (see \code{\link{constructGraph}}), a number of clusters
#' to assign the vertices of the graph to, and optionally the minimum number of vertices each cluster should contain.
#' If \code{minclust} is not specified, by default it will try to ensure no cluster is more than 10 times larger than any other cluster.
#' The method by which it assigns vertices to clusters is by using Prim's algorithm to find the minimum spanning tree
#' of the graph and then randomly making \code{nclust - 1} cuts to result in \code{nclust} disconnected trees. The vertices connected
#' by each tree are assigned to the same cluster, and the collection of these trees is returned as \code{spanning_forest}.
#'
#' Note that it is possible that after 100 different uniform random edge samples, no spanning forest meets the criteria
#' of the parameters (namely, not every cluster contains \code{minclust} elements). If this is the case a warning message
#' will be generated, but the algorithm will still proceed to assign cluster memberships as best it can.
#'
#' @param graph An object of class 'graph' from the \code{\link[igraph]{igraph}} package. The graph must have weights for each edge
#' @param nclust An integer, the number of different clusters to assign points to. Must be at most N (the number of vertices in graph)
#' @param minclust An integer, the smallest allowable cluster size. By default, one-tenth the ratio of vertices to nclust.
#'
#' @return A list containing three elements:
#' \item{clustered_graph}{The input graph with inter-cluster edges removed}
#' \item{spanning_forest}{The input graph with inter-cluster edges removed, and every cluster induced subgraph is a spanning tree}
#' \item{membership}{A vector of integers of length N with nclust unique integers which map each vertex to a cluster}
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
#' clust_out = constructClusters(g, 5, minclust = 3)
#' plot(clust_out$spanning_forest,
#'      layout = as.matrix(coords),
#'      vertex.color = clust_out$membership,
#'      edge.arrow.mode = 0)
constructClusters = function(graph, nclust, minclust = NULL) {
  # First, get the number of vertices
  N = igraph::vcount(graph)
  # Check that nclust is valid
  if(nclust > N) {
    stop("nclust cannot be larger than the number of vertices in graph")
  }
  if(nclust < 1) {
    stop("nclust must be positive and at least 1")
  }
  # Check that minclust is either not provided (in which case initialize)
  if(is.null(minclust)) {
    # Each cluster must contain a minimum of ~10% of the points by default
    minclust = floor((N/nclust) * 0.1)
  } else if(minclust < 1) {
    # Can't have 0 size clusters
    stop("minimum cluster size cannot be less than 1")
  } else if((minclust * nclust) > N) {
    # Can't have every cluster so large that there are insufficient vertices
    stop("minimum cluster size is too large; decrease it or nclust or add more vertices")
  }
  # Check that the input graph is connected (otherwise MST cuts will result in wrong number of clusters)
  if(!igraph::is.connected(graph)) {
    warning("Input graph is not connected; clusters may not be contiguous\n")
  }
  # Allocate the membership vector
  membership = rep(0, N)
  # Construct minimum spanning tree on the graph
  minspantree = igraph::mst(graph)
  # Do the following at least once and until all the clusters meet the minimum size requirement:
  repeat {
    # Count number of iterations in case a valid tree is not found in a 'reasonable' number of iterations
    iterations = 0
    # Choose nclust - 1 edges to delete from minspantree
    # Note: by definition a minimum spanning tree will have (N - 1) edges
    # However, if the input graph is not connected, then a 'cut' has already been implicitly done
    if(igraph::is.connected(graph)) {
      deleted_edge_ids = sample(1:(N - 1), nclust - 1, replace = FALSE)
    } else {
      disconnections = igraph::components(minspantree)$no
      deleted_edge_ids = sample(1:(N - disconnections), nclust - disconnections, replace = FALSE)
    }
    # Create a new graph with nclust spanning trees (a forest)
    spanforest = igraph::delete.edges(minspantree, deleted_edge_ids)
    # Update the membership of each vertex
    membership = igraph::components(spanforest)$membership
    # If all the clusters meet the size requirement, we are done
    if(min(table(membership)) >= minclust) {
      break
    }
    if(iterations > 100) {
      warning("Failed to find clusters meeting parameters in 100 iterations; clusters might be smaller than expected\n")
      break
    }
    iterations = iterations + 1
  }
  clustered_graph = clusterGraph(graph, membership)
  # Return list with clustered graph, spanning forest, and vertex membership
  return(list(clustered_graph = clustered_graph, spanning_forest = spanforest, membership = membership))
}
