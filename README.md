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

To install this package, run the following code:

```
install.packages("devtools")
library(devtools)
install_github("rayisaacalan/BASTION")
```

The tasks for this semester are as follows:
- ~~Package skeleton~~
  * ~~README~~
  * ~~Description~~
- C++ graph class / functions
  * ~~Graph class (data structure, adjacency list)~~
  * ~~Graph methods (subgraph, edge status, number of vertices/edges, add edge)~~
  * Primitive graph functions (Prim's algorithm, ~~vertex connectedness~~)
- C++ graph operations
  * Birth move
  * Death move
  * Change move
  * Hyper move
- R functions
  * ~~Construct Graph (K nearest neighbors)~~
  * Construct Clusters
- R wrappers
  * Compatibility / input validity checks
  * C++ function calls
- Documentation
  * ~~Roxygen skeleton (basic i/o)~~
  * Detailed descriptions
  * Examples
