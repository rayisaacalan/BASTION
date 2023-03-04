// No longer using RcppClock
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins("cpp11")]]

#include "RSTlearner.hpp"

using namespace Rcpp;


NumericVector combine(const List& list);

// [[Rcpp::export]]
Rcpp::List BASTIONfit_cpp(const Rcpp::IntegerMatrix &edges,
                      const Rcpp::NumericVector &weights,
                      const Rcpp::NumericVector &Y,
                      const int MCMC_iter,
                      const int BURNIN,
                      const int THIN,
                      const List init_values,
                      const List hyperpars) {
//Rprintf("Entered function body\n");
//Rcpp::Clock clock;
// Fetch hyper parameters
  //clock.tick("init_vals");
  // Fetch the initial values and hyper parameters
    const int k_max = hyperpars["k_max"];
    const int n_learners = hyperpars["n_learners"];
    const double lambda_k = hyperpars["lambda_k"];
    const double lambda_s = hyperpars["lambda_s"];
    const double nu = hyperpars["nu"];
    const double sigmasq_mu = hyperpars["sigmasq_mu"];
// Fetch initial values
    double sigmasq_y = init_values["sigmasq_y"];
    List mu = init_values["mu"];
//Rprintf("Fetched initial values\n");
  //clock.tock("init_vals");
  //clock.tick("graph_build");
// First, construct the graph from the edge list
  int n_edges = edges.rows();
  int n_verts = Y.size();
  // Create a vector of vertex pairs
  std::vector<std::pair<int, int>> edge_vec;
  for(int i = 0; i < n_edges; ++i) {
    edge_vec.emplace_back(std::make_pair(edges(i, 0), edges(i, 1)));
  }
//Rprintf("Created edge vector\n");
  // We now have our full graph with initial weights
  Graph g(edge_vec.begin(), edge_vec.end(), as<std::vector<double>>(weights).begin(), n_verts);
//Rprintf("Constructed graph\n");
  //clock.tock("graph_build");
  //clock.tick("learners_build");
  // Initialize a matrix to store each fitted value throughout the iterations
  // Each column represents a weak learner and each row represents a vertex
  NumericMatrix fitted_mus(n_verts, n_learners);
  // Vector of k's for evaluating log posterior
  NumericVector K(n_learners);
  // Initialize each weak learner ([TODO] further consideration on whether to use
  // list vs vector to store each weak learner is warranted)
  std::vector<RSTlearner> WeakLearners;
  for (int learner = 0; learner < n_learners; ++learner) {
    WeakLearners.emplace_back(&g, learner);
    K(learner) = (WeakLearners[learner].k);
  }
//Rprintf("Initialized Weak Learners\n");
  // Initialize some values to be used in the MCMC
  // Storage for output values
  int out_length = (((MCMC_iter) - BURNIN) % THIN);
  List mu_out;
  NumericVector sigmasq_y_out(out_length);
  List cluster_out;
  NumericVector log_post_out(out_length);
  // There should be at most k_max <= n_verts clusters
  int max_clusts = (k_max <= n_verts)? k_max : n_verts;
  // 0 is a birth, 1 is a death, 2 is a change, 3 is a hyper
  NumericVector moves = {0, 1, 2, 3};
  //clock.tock("learners_build");
//Rprintf("Initialized MCMC Values\n");
  for (int iter = 0; iter < MCMC_iter; ++iter) {
    //Rprintf("Entering iteration %i\n", iter);
    NumericMatrix MembershipByLearner(n_verts, n_learners);
      for (int learner = 0; learner < n_learners; ++learner) {
        //Rprintf(" Entering Weak Learner %i,", learner);
        int k_m = WeakLearners[learner].k;
        //Rprintf(" Currently %i clusters\n", k_m);
        NumericVector LearnerResponse = Y - (rowSums(fitted_mus) - fitted_mus( _ , learner));
        NumericVector moveProbability;
        if(k_m == 1) {
          moveProbability = NumericVector::create(0.9, 0, 0, 0.1);
        } else if(k_m == max_clusts) {
          moveProbability = NumericVector::create(0, 0.6, 0.3, 0.1);
        } else {
          moveProbability = NumericVector::create(0.3, 0.3, 0.3, 0.1);
        }
        int move = sample(moves, 1, false, moveProbability)(0);
        switch(move) {
          case 0: {
            //Rprintf("   Performing a Birth step\n");
            //clock.tick("birth");
            WeakLearners[learner].birth();
          // Calculate acceptance probability acc_prob
            // Log-prior ratio
            double log_A = log(lambda_k) - log(k_m + 1);
              // If k_m is almost at k_max, rd_new is smaller
            double rd_new = (k_m == (std::min(k_max, n_verts) - 1)) ? 0.6 : 0.3;
            // Log-proposal ratio
            double log_P = log(rd_new) - log(moveProbability(0));
            // Log-likelihood ratio
            double sigma_ratio = sigmasq_y / sigmasq_mu;
            int csize_old = WeakLearners[learner].old_clust_ids.size();
            int csize_new = WeakLearners[learner].new_clust_ids.size();
            NumericVector vid_old = wrap(WeakLearners[learner].old_clust_ids);
            NumericVector vid_new = wrap(WeakLearners[learner].new_clust_ids);
            NumericVector Response_old = LearnerResponse[vid_old];
            NumericVector Response_new = LearnerResponse[vid_new];
            double sum_e_old = sum(Response_old);
            double sum_e_new = sum(Response_new);
            double logdetdiff = -0.5 * (
              log(csize_old + sigma_ratio) +
              log(csize_new + sigma_ratio) -
              log(csize_old + csize_new + sigma_ratio) -
              log(sigma_ratio)
            );
            double quaddiff = (0.5/sigmasq_y) * (
              ((sum_e_old*sum_e_old)/(csize_old+sigma_ratio)) +
              ((sum_e_new*sum_e_new)/(csize_new+sigma_ratio)) -
              ((pow(sum_e_new + sum_e_old, 2))/(csize_old+csize_new+sigma_ratio))
            );
            double log_L = logdetdiff + quaddiff;
            // Calculate acceptance probability
            double log_acc_prob = fmin(0.0, log_A + log_P + log_L);
            double acc_prob = exp(log_acc_prob);
            if(R::runif(0, 1) < acc_prob) {
              WeakLearners[learner].accept();
            } else {
              WeakLearners[learner].reject();
            }
            //clock.tock("birth");
            break;
          }
          case 1: {
            //Rprintf("   Performing a Death step\n");
            //clock.tick("death");
            WeakLearners[learner].death();
          // Calculate acceptance probability acc_prob
            // Log-prior ratio
            double log_A = log(k_m) - log(lambda_k);
            // If k_m is almost at k_max, rd_new is smaller
            double rb_new = (k_m == 2) ? 0.9 : 0.3;
            // Log-proposal ratio
            double log_P = log(moveProbability(1)) - log(rb_new);
            // Log-likelihood ratio
            double sigma_ratio = sigmasq_y / sigmasq_mu;
            int csize_old = WeakLearners[learner].old_clust_ids.size();
            int csize_new = WeakLearners[learner].new_clust_ids.size();
            NumericVector vid_old = wrap(WeakLearners[learner].old_clust_ids);
            NumericVector vid_new = wrap(WeakLearners[learner].new_clust_ids);
            NumericVector Response_old = LearnerResponse[vid_old];
            NumericVector Response_new = LearnerResponse[vid_new];
            double sum_e_old = sum(Response_old);
            double sum_e_new = sum(Response_new);
            double logdetdiff = -0.5 * (
                log(csize_new + sigma_ratio) -
                log(csize_old + sigma_ratio) -
                log(csize_new - csize_old + sigma_ratio) +
                log(sigma_ratio)
            );
            double quaddiff = (0.5/sigmasq_y) * (
                ((sum_e_new*sum_e_new)/(csize_new+sigma_ratio)) -
                ((sum_e_old*sum_e_old)/(csize_old+sigma_ratio)) -
                ((pow(sum_e_new - sum_e_old, 2))/(csize_new-csize_old+sigma_ratio))
            );
            double log_L = logdetdiff + quaddiff;
            // Calculate acceptance probability
            double log_acc_prob = fmin(0.0, log_A + log_P + log_L);
            double acc_prob = exp(log_acc_prob);
            if(R::runif(0, 1) < acc_prob) {
              WeakLearners[learner].accept();
            } else {
              WeakLearners[learner].reject();
            }
            //clock.tock("death");
            break;
          }
          case 2: {
            //Rprintf("   Performing a Change step\n");
            //clock.tick("change");
            // Store the current graph state (since it will be going forwards
            // by two steps)
            WeakLearners[learner].storeGraphState();
            // First, compute log-likelihood ratio for the death move
            WeakLearners[learner].death();

              // Log-likelihood ratio
              double sigma_ratio = sigmasq_y / sigmasq_mu;
              int csize_old = WeakLearners[learner].old_clust_ids.size();
              int csize_new = WeakLearners[learner].new_clust_ids.size();
              NumericVector vid_old = wrap(WeakLearners[learner].old_clust_ids);
              NumericVector vid_new = wrap(WeakLearners[learner].new_clust_ids);
              NumericVector Response_old = LearnerResponse[vid_old];
              NumericVector Response_new = LearnerResponse[vid_new];
              double sum_e_old = sum(Response_old);
              double sum_e_new = sum(Response_new);
              double logdetdiff = -0.5 * (
                log(csize_new + sigma_ratio) -
                  log(csize_old + sigma_ratio) -
                  log(csize_new - csize_old + sigma_ratio) +
                  log(sigma_ratio)
              );
              double quaddiff = (0.5/sigmasq_y) * (
                ((sum_e_new*sum_e_new)/(csize_new+sigma_ratio)) -
                  ((sum_e_old*sum_e_old)/(csize_old+sigma_ratio)) -
                  ((pow(sum_e_new - sum_e_old, 2))/(csize_new-csize_old+sigma_ratio))
              );
              double log_L_death = logdetdiff + quaddiff;

            // Now, compute the log-likelihood ratio for the birth move
            // We assume to have accepted the death before doing the birth
            WeakLearners[learner].accept();
            WeakLearners[learner].birth();

              // Log-likelihood ratio
              int csize_old_b = WeakLearners[learner].old_clust_ids.size();
              int csize_new_b = WeakLearners[learner].new_clust_ids.size();
              NumericVector vid_old_b = wrap(WeakLearners[learner].old_clust_ids);
              NumericVector vid_new_b = wrap(WeakLearners[learner].new_clust_ids);
              NumericVector Response_old_b = LearnerResponse[vid_old_b];
              NumericVector Response_new_b = LearnerResponse[vid_new_b];
              double sum_e_old_b = sum(Response_old_b);
              double sum_e_new_b = sum(Response_new_b);
              double logdetdiff_b = -0.5 * (
                log(csize_old_b + sigma_ratio) +
                  log(csize_new_b + sigma_ratio) -
                  log(csize_old_b + csize_new_b + sigma_ratio) -
                  log(sigma_ratio)
              );
              double quaddiff_b = (0.5/sigmasq_y) * (
                ((sum_e_old_b*sum_e_old_b)/(csize_old_b+sigma_ratio)) +
                  ((sum_e_new_b*sum_e_new_b)/(csize_new_b+sigma_ratio)) -
                  ((pow(sum_e_new_b + sum_e_old_b, 2))/(csize_old_b+csize_new_b+sigma_ratio))
              );
              double log_L_birth = logdetdiff_b + quaddiff_b;

            // Add the log likelihood ratios for the two moves
            double log_L = log_L_death + log_L_birth;
            double log_acc_prob = fmin(0.0, log_L);
            double acc_prob = exp(log_acc_prob);
            if(R::runif(0, 1) < acc_prob) {
              WeakLearners[learner].accept();
            } else {
              WeakLearners[learner].restoreGraphState();
              WeakLearners[learner].reject();
            }
            //clock.tock("change");
            break;
          }
          case 3: {
            //clock.tick("hyper");
            //Rprintf("   Performing a Hyper step\n");
            WeakLearners[learner].hyper();
            WeakLearners[learner].accept();
            //clock.tock("hyper");
            break;
          }
        }
        //clock.tick("learner_update");
        //Rprintf("   Updating values\n");
        // Update the learner column of fitted_mus w/ new fitted mus
        k_m = WeakLearners[learner].k;
        NumericVector membership = wrap(WeakLearners[learner].membership);
        if((Rcpp::max(membership) + 1) > k_m) {
          WeakLearners[learner].printGraph();
          stop("Something has gone wrong with connected clusters!\n");
        }
        NumericVector csize(WeakLearners[learner].clusts.size());
        NumericVector ClustResponse(WeakLearners[learner].clusts.size());
        for(int i = 0; i < WeakLearners[learner].clusts.size(); ++i) {
          csize(i) = (WeakLearners[learner].clusts[i].size());
          NumericVector ClustIndices = wrap(WeakLearners[learner].clusts[i]);
          NumericVector ResponseByClust = LearnerResponse[ClustIndices];
          ClustResponse(i) = (sum(ResponseByClust));
        }
        NumericVector Qinv_diag = 1/((csize/sigmasq_y) + 1/(sigmasq_mu));
        NumericVector b = (Qinv_diag * ClustResponse) / sigmasq_y;
        NumericVector NewMus(k_m);
        for(int i = 0; i < k_m; ++i) {
          NewMus(i) = (R::rnorm(b[i], sqrt(Qinv_diag)[i]));
        }
        mu[learner] = NewMus;
        std::vector<double> new_response;
        for(int i = 0; i < n_verts; ++i) {
          new_response.emplace_back(NewMus[WeakLearners[learner].membership[i]]);
        }
        NumericVector NewResponse = wrap(new_response);
        fitted_mus( _ , learner) = NewResponse;
        MembershipByLearner( _ , learner) = membership;
        K[learner] = WeakLearners[learner].k;
        //clock.tock("learner_update");
      }
  // Update sigma squared of y
    // Simplification from R logic:
    // Y_hat = fitted_mus[, n_learners - 1] + Y - (Y - (rowsums(fitted_mus) - fitted_mus[, n_learners - 1]))
    // = fitted_mus[, n_learners - 1] + Y - Y + rowsums(fitted_mus) - fitted_mus[, n_learners - 1]
    // = rowsums(fitted_mus)
    //Rprintf(" Updating Iteration Values\n");
    //clock.tick("iter_update");
    NumericVector Y_hat = rowSums(fitted_mus);
    double scale = 1/(0.5*(nu*lambda_s + sum((Y - Y_hat)*(Y - Y_hat))));
    double shape = (double)(n_verts+nu)/2.0;
    sigmasq_y = 1/(R::rgamma(shape, scale));
    // Save the result
    if(((iter+1) > BURNIN) && ((((iter+1) - BURNIN) % THIN) == 0)) {
      mu_out.push_back(clone(mu));
      sigmasq_y_out.push_back(sigmasq_y);
      cluster_out.push_back(MembershipByLearner);
    // Evaluate log posterior
      // First evaluate log-likelihood
      double log_prior = -((1+(nu/2.0))*log(sigmasq_y))
                         -((nu*lambda_s)/(2.0*sigmasq_y))
                         +sum(-1*(lchoose(n_verts-1,K-1))
                              + K*log(lambda_k)
                              - lfactorial(K))
                         -(sum(pow(combine(mu), 2))/(2.0*sigmasq_mu));
      // Then evaluate log-likelihood
      double log_like = (-n_verts/2.0)*log(sigmasq_y)
                        -(sum((Y - Y_hat)*(Y - Y_hat))/(2.0*sigmasq_y));
      // Add the two
      log_post_out.push_back(log_prior+log_like);
    }
    //clock.tock("iter_update");
  }
//Rprintf("Finished MCMC Loop\n");
  //clock.stop("BASTIONfit_C_timings");
  // Return the result
  List OutputList = List::create(
    Named("mu_out") = mu_out,
    Named("sigmasq_y_out") = sigmasq_y_out,
    Named("cluster_out") = cluster_out,
    Named("log_post_out") = log_post_out
  );
  //Rprintf("Finished computations");
  return(OutputList);
}


// Utility function (similar to unlist)
// Taken from:
// https://stackoverflow.com/a/30201460/14190296

NumericVector combine(const List& list) {
  std::size_t n = list.size();
  // Figure out the length of the output vector
  std::size_t total_length = 0;
  for (std::size_t i = 0; i < n; ++i)
    total_length += Rf_length(list[i]);
  // Allocate the vector
  NumericVector output = no_init(total_length);
  // Loop and fill
  std::size_t index = 0;
  for (std::size_t i = 0; i < n; ++i)
  {
    NumericVector el = list[i];
    std::copy(el.begin(), el.end(), output.begin() + index);

    // Update the index
    index += el.size();
  }
  return output;
}




