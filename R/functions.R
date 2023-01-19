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
    warning("Input graph is not connected; clusters may not be contiguous")
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
      warning("Failed to find clusters meeting parameters in 100 iterations; clusters might be smaller than expected")
      break
    }
    iterations = iterations + 1
  }
  clustered_graph = clusterGraph(graph, membership)
  # Return list with clustered graph, spanning forest, and vertex membership
  return(list(clustered_graph = clustered_graph, spanning_forest = spanforest, membership = membership))
}




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
    warning("Possibly unintended behavior: Input graph is not a spanning forest; clust induced subgraph is not a tree")
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

#' Perform a cluster death operation (merge an existing cluster)
#'
#' \code{graphDeath} takes in a spanning forest graph with k disconnected components corresponding to
#' the vector of cluster assignments 'membership' and merges two connected clusters together. It requires the original
#' graph (as generated by \code{\link{constructGraph}}) to determine what edges in the modified graph are available to
#' be returned in order to connect disconnected components. From the list of edges which connect vertices belonging
#' to different clusters (see \code{\link{edgeBetweenClust}}), one is uniformly randomly selected and returned to
#' the output graph, reducing the number of clusters by one.
#'
#' @param graph An object of class 'graph' from the \code{\link[igraph]{igraph}} package
#' @param membership A vector of integers of length N with k unique integers (1 < k <= N) which map each vertex to a cluster
#' @param full_graph An object of class 'graph' from the \code{\link[igraph]{igraph}} package; 'graph' should be a subgraph of full_graph
#'
#' @return A list containing two elements:
#' \item{graph}{The input graph with 1 additional active edge}
#' \item{membership}{A vector of integers of length N with k - 1 unique integers which map each vertex to a cluster}
#' \item{new_clust_ids}{Vertex keys of the vertices belonging to the newly unified cluster}
#' \item{old_clust_ids}{Vertex keys of the vertices belonging only to the cluster being merged which has a higher number}
#' @export
#'
#' @seealso
#' \code{\link{constructClusters}}, \code{\link{graphBirth}}, \code{\link{graphChange}}, \code{\link{graphHyper}}
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
#' g = constructGraph(coords, 6)
#' clust_out = constructClusters(g, 8, minclust = 3)
#' plot(clust_out$spanning_forest,
#'      layout = as.matrix(coords),
#'      vertex.color = clust_out$membership,
#'      edge.arrow.mode = 0)
#' g_7_clusters = graphDeath(clust_out$spanning_forest,
#'                           clust_out$membership,
#'                           g)
#' plot(g_7_clusters$graph,
#'      layout = as.matrix(coords),
#'      vertex.color = g_7_clusters$membership,
#'      edge.arrow.mode = 0)
graphDeath = function(graph, membership, full_graph) {
  # Test that membership is valid
  N = igraph::vcount(graph)
  if(length(membership) != N) {
    stop("membership is of incorrect length, should be same length as vcount of graph")
  }
  # Test that graph has the same number of vertices as full_graph
  if(N != igraph::vcount(full_graph)) {
    stop("graph does not have the same number of vertices as full_graph")
  }
  k = length(unique(membership))
  # Test that k is valid
  if(k <= 1) {
    stop("Not enough clusters; need at least 2 clusters to be able to merge")
  }
  # Get a vector of which of the edges which are between the existing clusters
  betweenness = edgeBetweenClust(full_graph, membership)
  # Use betweenness to make a vector of the edges which are between clusters
  between_edges = igraph::E(full_graph)[betweenness]
  # Select one to return to the graph
  returning_edge = sample(between_edges, 1)
  # Get vertex IDs which the returning edge connects
  connecting_vertices = igraph::ends(full_graph, returning_edge)
  # Get vertices belonging to the higher cluster number (the 'old' cluster)
  old_clust_ids = which(membership == max(membership[connecting_vertices]))
  # Get vertices belonging to the newly merged cluster
  new_clust_ids = which(membership %in% membership[connecting_vertices])
  # Add the edge back to graph
  merged_graph = igraph::add_edges(graph, connecting_vertices)
  # Update the membership of each cluster
  merged_membership = igraph::components(merged_graph)$membership
  # Return the updated graph and new membership
  return(list(graph = merged_graph, membership = merged_membership, new_clust_ids = new_clust_ids, old_clust_ids = old_clust_ids))
}

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
  return(list(graph = reborn_graph$graph, membership = reborn_graph$membership,
              new_dclust_ids = dead_graph$new_clust_ids, old_dclust_ids = dead_graph$old_clust_ids,
              new_bclust_ids = reborn_graph$new_clust_ids, old_bclust_ids = reborn_graph$old_clust_ids))
}

#' Generate new spanning trees for each cluster
#'
#' \code{graphHyper} resamples the spanning trees which together partition the \code{full_graph} (as generated by
#' \code{\link{constructGraph}}) into the clusters specified by \code{membership}. It does this by first resampling every
#' edge weight from the full graph, with edges connecting vertices currently assigned to different clusters being drawn from
#' a uniform (0.5, 1) distribution (weighted more heavily), and edges connecting vertices currently assigned to the same
#' cluster being drawn from a uniform (0, 0.5) distribution (weighted less heavily). Then, the minimum spanning tree across the
#' entire graph is recalculated using these new weights. By removing the inter-cluster edges from this new MST (see
#' \code{\link{clusterGraph}}), the resultant graph has the same vertex membership, but the individual trees which make up
#' the spanning forest likely have a different set of edges.
#'
#' @param full_graph An object of class 'graph' from the \code{\link[igraph]{igraph}} package
#' @param membership  A vector of integers of length N with k unique integers (1 < k <= N) which map each vertex to a cluster
#'
#' @return A list containing two elements:
#' \item{graph}{The input graph with each cluster having a new minimum spanning tree based on resampled edge weights}
#' \item{membership}{A vector of integers of length N with k unique integers which map each vertex to a cluster}
#' @export
#'
#' @seealso
#' \code{\link{constructClusters}}, \code{\link{graphDeath}}, \code{\link{graphChange}}, \code{\link{graphBirth}}
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
#' coords = data.frame(lon = rnorm(100), lat = rnorm(100))
#' g = constructGraph(coords, 5)
#' clust_out = constructClusters(g, 8, minclust = 6)
#' plot(clust_out$spanning_forest,
#'      layout = as.matrix(coords),
#'      vertex.color = clust_out$membership,
#'      edge.arrow.mode = 0)
#' g_resample = graphHyper(g, clust_out$membership)
#' plot(g_resample$graph,
#'      layout = as.matrix(coords),
#'      vertex.color = g_resample$membership,
#'      edge.arrow.mode = 0)
graphHyper = function(full_graph, membership) {
  # Test that membership is valid
  N = igraph::vcount(full_graph)
  if(length(membership) != N) {
    stop("membership is of incorrect length, should be same length as vcount of graph")
  }
  # Test that full_graph is connected
  if(!igraph::is.connected(full_graph)) {
    warning("full_graph is not connected; clusters could be incorrect")
  }
  # Get the vector of edge betweenness of the full graph
  betweenness = edgeBetweenClust(full_graph, membership)
  # Get the vector of edges which are in between clusters
  between_edges = igraph::E(full_graph)[betweenness]
  # Get the vector of edges which are within clusters
  within_edges = igraph::E(full_graph)[!betweenness]
  #Initialize vector of weights
  weight = rep(0, igraph::ecount(full_graph))
  # Sample between cluster edge weights and assign them (uniform from 0.5 to 1)
  weight[betweenness] = stats::runif(length(between_edges), 0.5, 1)
  # Sample within cluster edge weights and assign them (uniform from 0 to 0.5)
  weight[!betweenness] = stats::runif(length(within_edges), 0, 0.5)
  # Recalculate the MST
  minspantree = igraph::mst(full_graph, weights = weight)
  # Remove inter-cluster edges to recover original partition
  resampled_graph = clusterGraph(minspantree, membership)
  # Return updated graph and same membership
  return(list(graph = resampled_graph, membership = membership))
}
