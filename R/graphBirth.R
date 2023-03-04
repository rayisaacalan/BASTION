#' Perform a cluster birth operation (split an existing cluster)
#'
#' \code{graphBirth} takes in a spanning forest graph with k disconnected components corresponding to
#' the vector of cluster assignments 'membership' and splits an existing cluster. Optionally, the integer corresponding
#' to one of the clusters can be specified as the one to split. If it is not specified, a cluster is uniformly selected at random.
#' In the specified or selected cluster, an edge is uniformly randomly selected and cut, turning the one cluster into two.
#'
#' @param graph An object of class 'graph' from the \code{\link[igraph]{igraph}} package
#' @param membership A vector of integers of length N with k unique integers (k < N) which map each vertex to a cluster
#' @param clust (Optional) Integer, the cluster to split. Must be between 1 and k. By default, a cluster will be uniformly randomly selected
#'
#' @importFrom igraph %--%
#'
#' @return A list containing two elements:
#' \item{graph}{The input graph with 1 fewer active edge}
#' \item{membership}{A vector of integers of length N with k + 1 unique integers which map each vertex to a cluster}
#' \item{new_clust_ids}{Vertex keys of the vertices belonging to the new cluster}
#' \item{old_clust_ids}{Vertex keys of the vertices from the cluster which was split excluding the ones from \code{new_clust_ids}}
#' @export
#'
#' @seealso
#' \code{\link{constructClusters}}, \code{\link{graphDeath}}, \code{\link{graphChange}}, \code{\link{graphHyper}}
#'
#' @references
#' Luo, Z.T. (*), Sang, H. and Mallick, B.K. (2021), BAST: Bayesian Additive Regression Spanning Trees
#' for Complex Constrained Domain
#'
#' Luo, Z.T. (*), Sang, H. and Mallick, B.K. (2021), A Bayesian Contiguous Partitioning Method for
#' Learning Clustered Latent Variables, Journal of Machine Learning Research, 22, 1-52.
#'
#' @examples
#' set.seed(1234)
#' coords = data.frame(lon = rnorm(50), lat = rnorm(50))
#' g = constructGraph(coords, 5)
#' clust_out = constructClusters(g, 6, minclust = 5)
#' plot(clust_out$spanning_forest,
#'      layout = as.matrix(coords),
#'      vertex.color = clust_out$membership,
#'      edge.arrow.mode = 0)
#' g_7_clusters = graphBirth(clust_out$spanning_forest,
#'                           clust_out$membership,
#'                           4)
#' plot(g_7_clusters$graph,
#'      layout = as.matrix(coords),
#'      vertex.color = g_7_clusters$membership,
#'      edge.arrow.mode = 0)
graphBirth = function(graph, membership, clust = NULL) {
  # Test that membership is valid
  N = igraph::vcount(graph)
  if(length(membership) != N) {
    stop("membership is of incorrect length, should be same length as vcount of graph")
  }
  k = length(unique(membership))
  # Test that k is valid
  if(k >= N) {
    stop("Not enough vertices in graph to generate a new cluster")
  }
  # If clust is NULL, uniformly pick a random cluster
  if(is.null(clust)) {
    # Exclude clusters which only have one vertex
    valid_clusts = (1:k)[which(igraph::components(graph)$csize > 1)]
    if(length(valid_clusts) == 1) {
      clust = valid_clusts
    } else {
      clust = sample(valid_clusts, 1)
    }
  }
  # Test that clust is valid
  if(clust > k || clust < 1) {
    stop("Invalid cluster to split, should be an integer between 1 and k")
  }
  # Get the vertex (keys) which belong to clust from membership
  ids_by_cluster = split(1:N, membership)
  clust_ids = ids_by_cluster[[clust]]
  # Test that there is at least 2 vertices in the specified cluster
  if(length(clust_ids) < 2) {
    stop("Specified cluster has only 1 vertex; cannot split")
  }
  # Get the sequence of edges in the specified cluster
  clust_edges = igraph::E(graph)[clust_ids %--% clust_ids]
  # Check for validity of cluster (should be a tree; trees have n - 1 edges for n vertices)
  if((length(clust_ids)) - 1 != length(clust_edges)) {
    warning("Possibly unintended behavior: Input graph is not a spanning forest; clust induced subgraph is not a tree\n")
  }
  # Randomly choose an edge to remove
  cut_edge = sample(clust_edges, 1)
  # Remove the edge
  cut_graph = igraph::delete.edges(graph, cut_edge)
  # Update the membership of each cluster
  cut_membership = igraph::components(cut_graph)$membership
  # Get the vertex (keys) which belong to the 'new' and 'old' clust from membership
  # We'll define the 'new' clust as being the one with a greater cluster index
  new_clust_ids = which(cut_membership == max(cut_membership[clust_ids]))
  old_clust_ids = which(cut_membership == min(cut_membership[clust_ids]))
  # Return the updated graph and new membership
  return(list(graph = cut_graph, membership = cut_membership, new_clust_ids = new_clust_ids, old_clust_ids = old_clust_ids))
}
