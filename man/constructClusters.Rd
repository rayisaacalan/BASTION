% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{constructClusters}
\alias{constructClusters}
\title{Construct a cluster membership list for a graph}
\usage{
constructClusters(graph, nclust, minclust = NULL)
}
\arguments{
\item{graph}{An object of class 'graph' from the \code{\link[igraph]{igraph}} package. The graph must have weights for each edge}

\item{nclust}{An integer, the number of different clusters to assign points to. Must be at most N (the number of vertices in graph)}

\item{minclust}{An integer, the smallest allowable cluster size. By default, one-tenth the ratio of vertices to nclust.}
}
\value{
A list containing three elements:
\item{clustered_graph}{The input graph with inter-cluster edges removed}
\item{spanning_forest}{The input graph with inter-cluster edges removed, and every cluster induced subgraph is a spanning tree}
\item{membership}{A vector of integers of length N with nclust unique integers which map each vertex to a cluster}
}
\description{
\code{constructClusters} takes in a weighted and connected graph (see \code{\link{constructGraph}}), a number of clusters
to assign the vertices of the graph to, and optionally the minimum number of vertices each cluster should contain.
If \code{minclust} is not specified, by default it will try to ensure no cluster is more than 10 times larger than any other cluster.
The method by which it assigns vertices to clusters is by using Prim's algorithm to find the minimum spanning tree
of the graph and then randomly making \code{nclust - 1} cuts to result in \code{nclust} disconnected trees. The vertices connected
by each tree are assigned to the same cluster, and the collection of these trees is returned as \code{spanning_forest}.
}
\details{
Note that it is possible that after 100 different uniform random edge samples, no spanning forest meets the criteria
of the parameters (namely, not every cluster contains \code{minclust} elements). If this is the case a warning message
will be generated, but the algorithm will still proceed to assign cluster memberships as best it can.
}
\examples{
set.seed(1)
coords = data.frame(lon = rnorm(50), lat = rnorm(50))
g = constructGraph(coords, 4)
clust_out = constructClusters(g, 5, minclust = 3)
plot(clust_out$spanning_forest,
     layout = as.matrix(coords),
     vertex.color = clust_out$membership,
     edge.arrow.mode = 0)
}
\references{
Luo, Z.T. (*), Sang, H. and Mallick, B.K. (2021), BAST: Bayesian Additive Regression Spanning Trees
for Complex Constrained Domain

Luo, Z.T. (*), Sang, H. and Mallick, B.K. (2021), A Bayesian Contiguous Partitioning Method for
Learning Clustered Latent Variables, Journal of Machine Learning Research, 22, 1-52.
}
