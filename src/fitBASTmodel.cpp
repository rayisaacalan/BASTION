#include "BAST_Sampler.hpp"
#include "SamplerEngine.hpp"

#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::List fitBASTmodel(const Rcpp::IntegerMatrix &edges,
                        const Rcpp::NumericVector &weights,
                        const Rcpp::NumericVector &Y,
                        const int MCMC_iter,
                        const int BURNIN,
                        const int THIN,
                        const double init_sigmasq_y,
                        const Rcpp::List init_mu_values,
                        const Rcpp::List hyperpars) {
try{
  // First, construct the vector of learners
  DEBUG_MSG("Entered fitBASTmodel");
  // Build the graph (constant across learners for this model)
  int n_edges = edges.rows();
  int n_verts = Y.size();
  std::vector<std::pair<int, int>> edge_vec;
  for(int i = 0; i < n_edges; ++i) {
    edge_vec.emplace_back(std::make_pair(edges(i, 0), edges(i, 1)));
  }
  Graph g(edge_vec.begin(), edge_vec.end(),
          Rcpp::as<std::vector<double>>(weights).begin(), n_verts);
  DEBUG_MSG("Constructed graph");
  // Construct PRNG seeds for each learner's generator
  const unsigned int n_learners = Rcpp::as<int>(hyperpars[1]);
  std::vector<unsigned int> inp_seeds(n_learners);
  std::iota(inp_seeds.begin(), inp_seeds.end(), Rcpp::as<unsigned int>(hyperpars[6]));
  std::seed_seq seed_seq_inp(inp_seeds.begin(), inp_seeds.end());
  std::vector<unsigned int> output_seeds(n_learners);
  seed_seq_inp.generate(output_seeds.begin(), output_seeds.end());
  DEBUG_MSG("Constructed PRNG seeds");
  // Initialize the learner vector
  // Note - doesn't reserve memory in advance because no default constructor
  std::vector<BAST_Sampler> BAST_learners;
  // Build each learner in place
  for(int i = 0; i < n_learners; ++i) {
    Rcpp::NumericVector mu_init = init_mu_values[i];
    BAST_learners.emplace_back(&g, Y, hyperpars, mu_init, init_sigmasq_y, output_seeds[i], i);
  }
  DEBUG_MSG("Constructed learners");

  // Next, initialize the sampling engine

  SamplerEngine<BAST_Sampler, BAST_output_T>
    BASTengine(BAST_learners,
               MCMC_iter, BURNIN, THIN,
               false, // No blocking is done in BAST model
               1); // BAST model is not parallel

  DEBUG_MSG("Constructed BASTengine");
  // Run the sampling engine
  DEBUG_MSG("Running BASTengine");
  BASTengine.runEngine();

  DEBUG_MSG("Finished running BASTengine");

  // Gather the outputs
  int n_output_samps = BASTengine.outputs[0].size();

  Rcpp::List mu_samps(n_output_samps);
  Rcpp::NumericVector sigmasq_y_samps(n_output_samps);
  Rcpp::List membership_samps(n_output_samps);
  Rcpp::NumericVector log_post_samps(n_output_samps);

  // BASTengine.outputs is a vector<vector<BAST_output_T>
  // Outer vector is learners; inner vector is samples
  // BAST_output_T contains vector<double> mu_out, vector<int> membership_out
  // as well as double sigmasq_y_out, double log_post_out which are constant
  // across learner wrt sample
  for(int sample = 0; sample < n_output_samps; ++sample) {
    Rcpp::List mu_samp(n_learners);
    Rcpp::List membership_samp(n_learners);
    for(int learner = 0; learner < n_learners; ++learner) {
      // Get values that are constant wrt learner
      if(learner == 0) {
        sigmasq_y_samps[sample] = BASTengine.outputs[learner][sample].sigmasq_y_out;
        log_post_samps[sample] = BASTengine.outputs[learner][sample].log_post_out;
      }
      // Get values that differ across learners

      mu_samp[learner] = (Rcpp::wrap(BASTengine.outputs[learner][sample].mu_out));
      membership_samp[learner] = (Rcpp::wrap(BASTengine.outputs[learner][sample].membership_out));
    }
    mu_samps[sample] = (Rcpp::clone(mu_samp));
    membership_samps[sample] = (Rcpp::clone(membership_samp));
  }

  DEBUG_MSG("Finished gathering outputs");

  // Return

  Rcpp::List OutputList = Rcpp::List::create(
    Rcpp::Named("mu_out") = mu_samps,
    Rcpp::Named("sigmasq_y_out") = sigmasq_y_samps,
    Rcpp::Named("cluster_out") = membership_samps,
    Rcpp::Named("log_post_out") = log_post_samps
  );
  return(OutputList);

} catch(std::exception &ex) {
  forward_exception_to_r(ex);
} catch(...) {
  Rcpp::Rcerr << "C++ exception caught; unknown reason" << std::endl;
}
  return Rcpp::List::create();
}

