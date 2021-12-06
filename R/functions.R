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
#' @param graph An object of class 'graph' from the igraph package.
#' @param membership A vector of integers of length N (see \code{\link{constructClusters}})
#'
#' @return A vector of boolean values of length \code{ecount(graph)} where each TRUE value represents an edge that connects distinct clusters
#' @export
#'
#' @examples
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
#' @inheritParams edgeBetweenClust
#'
#' @return An object of class 'graph' from the igraph package, the input graph with inter-cluster edges removed
#' @export
#'
#' @examples
clusterGraph = function(graph, membership) {
  # Get vector of cluster betweenness
  between_clusters = edgeBetweenClust(graph, membership)
  # Get ids of edges which are between clusters
  edge_ids_between = which(between_clusters)
  # Delete those edges
  clustered_graph = igraph::delete.edges(graph, edge_ids_between)
  # Return clustered graph
  return(clusterGraph)
}


#' Construct a cluster membership list for a graph
#'
#' @param graph An object of class 'graph' from the igraph package. The graph must have weights for each edge
#' @param nclust An integer, the number of different clusters to assign points to. Must be at most N (the number of vertices in graph)
#' @param minclust An integer, the smallest allowable cluster size. By default, one-tenth the ratio of vertices to nclust.
#'
#' @return A list containing two elements:
#' 'clustered_graph': The input graph with inter-cluster edges removed
#' 'membership': A vector of integers of length N with nclust unique integers which map each vertex to a cluster
#' @export
#'
#' @examples
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
  # Allocate the membership vector
  membership = rep(0, N)
  # Construct minimum spanning tree on the graph
  minspantree = igraph::mst(graph)
  # Do the following at least once and until all the clusters meet the minimum size requirement:
  repeat {
    # Choose nclust - 1 edges to delete from minspantree
    # Note: by definition a minimum spanning tree will have (N - 1) edges
    deleted_edge_ids = sample(1:(N - 1), nclust - 1, replace = FALSE)
    # Create a new graph with nclust spanning trees (a forest)
    spanforest = igraph::delete.edges(minspantree, deleted_edge_ids)
    # Update the membership of each vertex
    membership = igraph::components(spanforest)$membership
    # If all the clusters meet the size requirement, we are done
    if(min(table(membership)) >= minclust) {
      break
    }
  }
  clustered_graph = clusterGraph(graph, membership)
  # Return list with clustered graph and vertex membership
  return(list(clustered_graph = clustered_graph, membership = membership))
}




#' Perform a cluster birth operation (split an existing cluster)
#'
#' @param graph An object of class 'graph' from the igraph package
#' @param membership A vector of integers of length N with k unique integers (k < N) which map each vertex to a cluster
#' @param clust Integer, the cluster to split. Must be between 1 and k
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
#' @param clust Integer, the cluster to split. Must be between 1 and k
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
