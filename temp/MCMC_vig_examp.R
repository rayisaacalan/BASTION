library(BASTION)

set.seed(12345)
x = runif(1000, min = -10, max = 10)
y = runif(1000, min = -10, max = 10)
data_full = data.frame(x, y)
data_hole = data_full[x^2 + y^2 >= 16, ]

response = 1.1^(data_hole$x) +
  (data_hole$y)^2 -
  350*(data_hole$x > data_hole$y) +
  200*(-data_hole$x > data_hole$y) +
  5*data_hole$x*data_hole$y

hole_graph = constructGraph(data_hole, 5)

# function to standardize Y
standardize <- function(x) {
  xmean = mean(x)
  x = x - xmean
  xscale = 2 * max(abs(x))
  x = x / xscale
  param = c('mean' = xmean, 'scale' = xscale)
  return(list(x = x, std_par = param))
}

# function to unstandardize Y
unstandardize <- function(x, std_par, nomean = F, s2 = F) {
  if(s2) {
    x = x * std_par['scale'] ^ 2
  } else {
    x = x * std_par['scale']
  }
  if(!nomean) x = x + std_par['mean']
  return(x)
}

M = 30      # number of weak learners
k_max = 6   # maximum number of clusters per weak learner
mu = list() # initial values of mu (piecewise constant value to fit to cluster)
n = length(response)
cluster = matrix(1, nrow = n, ncol = M)  # initial cluster memberships
for(m in 1:M) {
  mu[[m]] = c(0)
}

init_val = list()
init_val[['mu']] = mu
init_val[['sigmasq_y']] = 1

# standardize our response
std_res = standardize(response)
Y_std = std_res$x
std_par = std_res$std_par

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
MCMC = 1000    # MCMC iterations
BURNIN = 500  # burnin period length
THIN = 5        # thinning intervals


BAST_model = BASTIONfit_C_new(Y_std, hole_graph, init_val, hyperpar, MCMC, BURNIN, THIN, seed = 1234)
