library(BASTION)
library(ggplot2)
library(ggraph)
library(geometry)
library(igraph)
set.seed(12345)
n <- 50 # even number
x.1 = runif(n/2, min = -10, max = 10)
y.1 = runif(n/2, min = -10, max = 10)
x.2 = x.1 + 20 
y.2 = y.1 + 20 
data_full = data.frame(x = c(x.1,x.2), y = c(y.1, y.2))

Y <- c(rep(5, n/2) , rep(-5, n/2)) + rnorm(n)
d_g <-delaunayn(data_full)
#create edges
e <- c()
w <- c()
for(i in 1:dim(d_g)[1]){
  e <- rbind(e, c(d_g[i,1],d_g[i,2]))
  e <- rbind(e, c(d_g[i,2],d_g[i,3]))
  e <- rbind(e, c(d_g[i,3],d_g[i,1]))
  w <- c(w, sqrt(sum(diff(as.matrix(data_full[c(d_g[i,1],d_g[i,2]),]))^2)))
  w <- c(w, sqrt(sum(diff(as.matrix(data_full[c(d_g[i,2],d_g[i,3]),]))^2)))
  w <- c(w, sqrt(sum(diff(as.matrix(data_full[c(d_g[i,3],d_g[i,1]),]))^2)))
}
w <- w[!duplicated(e)]
e <- e[!duplicated(e),]

#edge_distances = sqrt(rowSums((coord_data[edge_matrix[ , 1], ] - coord_data[edge_matrix[ , 2], ])^2))
# Assign the distance for each edge to the graph's edge weights
#igraph::E(graph)$weight = edge_distances

data_graph <- graph_from_data_frame(data.frame(v1=e[,1], v2 = e[,2]), directed=FALSE)
data_graph <- set_edge_attr(data_graph,"weight", value = w)

M = 5       # number of weak learners
k_max = 3   # maximum number of clusters per weak learner
mu = list() # initial values of mu (piecewise constant value to fit to cluster)
n = length(Y)
cluster = matrix(1, nrow = n, ncol = M)  # initial cluster memberships
for(m in 1:M) {
  mu[[m]] = c(0)
}

init_val = list()
init_val[['mu']] = mu
init_val[['sigmasq_y']] = 1

# standardize our response
Y_std = Y / sd(Y)

# find lambda_s
nu = 3 
q = 0.9
quant = qchisq(1-q, nu)
lambda_s = quant * var(Y_std) / nu

hyperpar = c()
hyperpar['sigmasq_mu'] = (0.5/(2*sqrt(M)))^2
hyperpar['lambda_s'] = lambda_s
hyperpar['nu'] = nu
hyperpar['lambda_k'] = 4
hyperpar['M'] = M
hyperpar['k_max'] = k_max

# MCMC parameters
# number of posterior samples = (MCMC - BURNIN) / THIN
MCMC = 100   # MCMC iterations
BURNIN = 50  # burnin period length
THIN = 5        # thinning intervals

BAST_model = BASTIONfit_C(Y_std, data_graph, init_val, hyperpar, MCMC, BURNIN, THIN, seed = 1234)