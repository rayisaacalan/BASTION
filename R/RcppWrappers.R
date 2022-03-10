

#' BASTION Regression using Rcpp
#'
#' @param Y A numeric vector of responses with length equal to the number of vertices in 'graph'
#' @param graph An object of class 'graph' from the \code{\link[igraph]{igraph}} package. The graph must have weights for each edge
#' @param init_vals A named list with the following components:
#' \item{sigmasq_y}{A numeric value}
#' \item{mu}{A list of length M where each element is a vector of length 1}
#' @param hyperpars A named list with the following components:
#' \item{k_max}{An integer, the maximum number of clusters for each weak learner}
#' \item{M}{An integer, the number of weak learners}
#' \item{lambda_k}{A numeric value}
#' \item{lambda_s}{A numeric value}
#' \item{nu}{A numeric value}
#' \item{sigmasq_mu}{A numeric value}
#' @param MCMC An integer, the number of iterations for the MCMC
#' @param BURNIN An integer, the number of initial iterations for the MCMC before recording output
#' @param THIN An integer, specifies the thinning interval
#' @param seed Optional, value to feed to R's set.seed function
#'
#' @return A named list with the following components:
#' \item{mu_out}{List of fitted values for mu for each weak learner at each output interval}
#' \item{sigmasq_y_out}{List of sigmasq_y at each output interval}
#' \item{cluster_out}{List of cluster membership of each vertex for each weak learner at each output interval}
#' \item{log_post_out}{List of the log of the posterior probability at each output interval}
#' @export
#'
#' @examples
#' test_coords = data.frame(x = c(0,1,1,0), y = c(1,0,1,0)) # observed locations
#' test_graph = constructGraph(test_coords, 3)
#' # Not run: plot this graph using ggraph package
#' # ggraph(test_graph, layout = test_coords) + geom_node_point() + geom_edge_link0()
#' n = igraph::vcount(test_graph)
#' z = c(1,2,3,4) # response
#'
#' M = 3      # number of weak learners
#' k_max = 2   # maximum number of clusters per weak learner
#' mu = list() # initial values of mu
#' cluster = matrix(1, nrow = n, ncol = M)  # initial cluster memberships
#' for(m in 1:M) {
#'   mu[[m]] = c(0)
#' }
#'
#' init_val = list()
#' init_val[['mu']] = mu
#' init_val[['sigmasq_y']] = 1
#'
#'
#' # find lambda_s
#' nu = 3; q = 0.9
#' quant = qchisq(1-q, nu)
#' lambda_s = quant * var(z) / nu
#'
#' hyperpar = c()
#' hyperpar['sigmasq_mu'] = (0.5/(2*sqrt(M)))^2
#' hyperpar['lambda_s'] = lambda_s
#' hyperpar['nu'] = nu
#' hyperpar['lambda_k'] = 4
#' hyperpar['M'] = M
#' hyperpar['k_max'] = k_max
#'
#' # MCMC parameters
#' # number of posterior samples = (MCMC - BURNIN) / THIN
#' MCMC = 100    # MCMC iterations
#' BURNIN = 50  # burnin period length
#' THIN = 5        # thinning intervals
#'
#' BASTIONfit_C(z, test_graph, init_val, hyperpar, MCMC, BURNIN, THIN, seed = 12345)
#'
#'
BASTIONfit_C = function(Y, graph, init_vals, hyperpars, MCMC, BURNIN, THIN, seed = NULL) {
  if(!is.null(seed)) {
    set.seed(seed)
  }
  if(!igraph::is.connected(graph)) {
    stop("Input graph must be connected")
  }
  if(length(Y) != igraph::vcount(graph)) {
    stop("Dimension of Y and graph are mismatched")
  }
  hyperpars_C = list(k_max = as.integer(hyperpars["k_max"]),
                     n_learners = as.integer(hyperpars["M"]),
                     lambda_k = as.double(hyperpars["lambda_k"]),
                     lambda_s = as.double(hyperpars["lambda_s"]),
                     nu = as.double(hyperpars["nu"]),
                     sigmasq_mu = as.double(hyperpars["sigmasq_mu"]))
  init_vals_C = list(sigmasq_y = as.double(init_vals["sigmasq_y"]),
                     mu = init_vals[["mu"]])
  if(length(init_vals_C[["mu"]]) != hyperpars_C["n_learners"]) {
    stop("Incorrect dimension of init_val of mu; should be of length M")
  }
  if(MCMC <= BURNIN) {
    stop("MCMC must be greater than BURNIN")
  }
  if((MCMC - BURNIN)/THIN < 1) {
    stop("(MCMC-BURNIN)/THIN must be at least 1")
  }
  if(MCMC < 1 || BURNIN < 1 || THIN < 1) {
    stop("MCMC, BURNIN, THIN must all be positive integers")
  }
  Edges = igraph::as_edgelist(graph, names = FALSE) - 1
  if(min(Edges) != 0) {
    stop("Graph's vertices are not 1-indexed")
  }
  Weights = as.numeric(igraph::get.edge.attribute(graph, "weight"))
  if(any(is.na(Weights))) {
    stop("All edge weights must be well defined (not NA)")
  }
  if(length(Weights) != nrow(Edges)) {
    stop("Not all edges have a weight")
  }
  result = BASTIONfit_cpp(Edges, Weights, Y, MCMC, BURNIN, THIN, init_vals_C, hyperpars_C)
  result$cluster_out = aperm(simplify2array(result$cluster_out), c(3,1,2)) + 1
  mode(result$cluster_out) = 'integer'
  return(result)
}
