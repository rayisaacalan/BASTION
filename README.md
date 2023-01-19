# BASTION
Bayesian Additive Spanning Trees Input/Output for Non-parametric regression

The Bayesian Additive Regression Spanning Trees model is a novel (still in development) ensemble model
for non-parametric regression intended to be used on data (especially spatial) that lies on a complex or
constrained domain, or a space with irregular shape embedded in Euclidean space. Existing ensemble
models for non-parametric regression such as Bayesian Additive Regression Trees (BART) or XGBoost are
very popular and effective, but these models rely on binary decision tree partition models as their weak
learners. These binary decision tree partitions do not respect domain constraints as they only make splits
parallel to Euclidean axes. At the core of the BAST model is a novel weak learner; a random spanning tree
manifold partition model. The model is based upon 4 possible moves/graph operations which will eventually
comprise a tailored backfitting Markov chain Monte Carlo algorithm.

This package implements those 4 possible graph operations as well as some utility functions to make it easier to
utilize them.

To install this package, run the following code:

```
install.packages("devtools")
devtools::install_github("rayisaacalan/BASTION", dependencies = TRUE)
```

To install this package and view its vignette, run the following code:

```
install.packages("devtools")
devtools::install_github("rayisaacalan/BASTION", build_vignettes = TRUE, dependencies = TRUE, force = TRUE)
library(BASTION)
vignette("MCMC", package = "BASTION")
```
To get started with a simple clustered graph, try the following:

```
set.seed(1)
coords = data.frame(lon = rnorm(50), lat = rnorm(50))
g = constructGraph(coords, 4)
clust_out = constructClusters(g, 5, minclust = 3)
plot(clust_out$spanning_forest,
     layout = as.matrix(coords),
     vertex.color = clust_out$membership,
     edge.arrow.mode = 0)
```


