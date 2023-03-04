#' Perform a cluster death operation followed by a cluster birth operation
#'
#' \code{graphChange} executes \code{\link{graphDeath}} on the inputs and then executes \code{\link{graphBirth}}
#' on the output of \code{\link{graphDeath}}. The net effect is that the output graph has the same number of clusters as
#' the input graph, but different vertex membership. This is to encourage a better mixing time of the Markov chain when this
#' function is used to implement an MCMC method.
#'
#' @param graph An object of class 'graph' from the \code{\link[igraph]{igraph}} package
#' @param membership A vector of integers of length N with k unique integers (1 < k <= N) which map each vertex to a cluster
#' @param full_graph An object of class 'graph' from the \code{\link[igraph]{igraph}} package; 'graph' should be a subgraph of full_graph
#'
#' @return A list containing two elements:
#' \item{graph}{The input graph with the different set of active edges}
#' \item{intermediate_graph}{The graph from performing the death operation before then performing the birth operation}
#' \item{membership}{A vector of integers of length N with k unique integers which map each vertex to a cluster, likely different than input membership}
#' \item{new_dclust_ids}{Vertex keys of the vertices belonging to the newly unified cluster before the birth operation}
#' \item{old_dclust_ids}{Vertex keys of the vertices belonging only to the cluster being merged which has a higher number before the birth operation}
#' \item{new_bclust_ids}{Vertex keys of the vertices belonging to the new cluster after the birth operation}
#' \item{old_bclust_ids}{Vertex keys of the vertices from the cluster which was born excluding the ones from \code{new_bclust_ids}}
#' @export
#'
#' @seealso
#' \code{\link{constructClusters}}, \code{\link{graphDeath}}, \code{\link{graphBirth}}, \code{\link{graphHyper}}
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
#' coords = data.frame(lon = rnorm(80), lat = rnorm(80))
#' g = constructGraph(coords, 5)
#' clust_out = constructClusters(g, 6, minclust = 8)
#' plot(clust_out$spanning_forest,
#'      layout = as.matrix(coords),
#'      vertex.color = clust_out$membership,
#'      edge.arrow.mode = 0)
#' g_different = graphChange(clust_out$spanning_forest,
#'                           clust_out$membership, g)
#' plot(g_different$graph,
#'      layout = as.matrix(coords),
#'      vertex.color = g_different$membership,
#'      edge.arrow.mode = 0)
graphChange = function(graph, membership, full_graph) {
  # First, perform a cluster death operation
  dead_graph = graphDeath(graph, membership, full_graph)
  # Then, perform a birth operation
  reborn_graph = graphBirth(dead_graph$graph, dead_graph$membership)
  # Return the result
  return(list(graph = reborn_graph$graph, intermediate_graph = dead_graph$graph, membership = reborn_graph$membership,
              new_dclust_ids = dead_graph$new_clust_ids, old_dclust_ids = dead_graph$old_clust_ids,
              new_bclust_ids = reborn_graph$new_clust_ids, old_bclust_ids = reborn_graph$old_clust_ids))
}
