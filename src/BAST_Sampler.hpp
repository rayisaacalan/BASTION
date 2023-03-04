#ifndef BAST_SAMPLER_HPP
#define BAST_SAMPLER_HPP

// This definition enables debug statements to Rcout
// #define BAST_DEBUG

#ifdef BAST_DEBUG
#define DEBUG_MSG(str) do { Rcpp::Rcout << str << std::endl; } while( false )
#else
#define DEBUG_MSG(str) do { } while ( false )
#endif

#include <cmath>
// Stop making a fool out of me; why don't you come on over valarray?
#include <valarray>
#include "RSTlearner.hpp"
#include "Sampler.hpp"
#include "util_funcs.hpp"

struct BAST_output_T {
  std::vector<double> mu_out;
  double sigmasq_y_out;
  std::vector<int> membership_out;
  double log_post_out;
};
typedef struct BAST_output_T BAST_output_T;


class BAST_Sampler : public RSTlearner, public Sampler<BAST_output_T, BAST_Sampler> {
public:
// Inherited data members
  using RSTlearner::generator;
  using Sampler<BAST_output_T, BAST_Sampler>::output;
// Hyperparameters
  int max_clusts;
  // Note that the class only cares about itself, but sometimes n_learners is used in sampling
  int learner_id;
  int n_learners;

  double lambda_k;
  double lambda_s;
  double nu;
  double sigmasq_mu;
  std::valarray<double> Y;
// Values which change across iterations
  // This maps cluster values to clusters
  std::valarray<double> mu;
  // Inherited from RSTlearner
  // std::vector<int> membership;
  // This maps individual vertices to their fitted value
  std::valarray<double> learner_response;
// Values that have knowledge of other weak learners
  // Current total fitted value for each vertex
  // std::valarray<double> Y_hat;
  // Number of clusters in each learner
  // std::valarray<double> K_learners;
  // Data variance (inverse gamma)
  double sigmasq_y;
  // For tracking model convergence to stationary dist
  double log_post;

  // Define moves as int values
  enum Move {birth_m = 0, death_m = 1, change_m = 2, hyper_m = 3};

  BAST_Sampler(const Graph* graph,
               const Rcpp::NumericVector &response,
               const Rcpp::List &hyperpars,
               const Rcpp::NumericVector &init_mu_values,
               const double init_sigmasq_y,
               const unsigned int seed,
               const int id) :
    max_clusts(Rcpp::as<int>(hyperpars[0]) <= boost::num_vertices(*graph) ?
                 Rcpp::as<int>(hyperpars[0]) :
                 boost::num_vertices(*graph)),
    n_learners(Rcpp::as<int>(hyperpars[1])),
    lambda_k(Rcpp::as<double>(hyperpars[2])),
    lambda_s(Rcpp::as<double>(hyperpars[3])),
    nu(Rcpp::as<double>(hyperpars[4])),
    sigmasq_mu(Rcpp::as<double>(hyperpars[5])),
    learner_id(id),
    RSTlearner(graph, seed),
    Sampler<BAST_output_T, BAST_Sampler>(BAST_output_T {}),
    Y((Rcpp::as<std::vector<double>>(response)).data(), response.size())
  {
    mu = std::valarray<double>((Rcpp::as<std::vector<double>>(init_mu_values)).data(), k);
    sigmasq_y = init_sigmasq_y;
    log_post = 0;
    //Y_hat = std::valarray<double>(0.0, Y.size());
    //K_learners = std::valarray<double>(1.0, n_learners);
    learner_response = std::valarray<double>(0.0, n);
    DEBUG_MSG("Constructed a BAST_Sampler");
  }


  void sampleMH() override {
    DEBUG_MSG("Entering sampleMH()");
    double unif_val = r_std_unif();
    double log_acc_prob = 0;
    double log_prior_ratio = 0;
    double log_prop_ratio = 0;
    double log_like_ratio = 0;
    int current_move = suggestMove();
    switch(current_move) {
      case birth_m: {
          DEBUG_MSG("   Selected a birth() move");
          birth();
          log_prior_ratio = std::log((double) lambda_k) -
                            std::log((double) prop_k);
          log_prop_ratio = std::log(moveProbabilities(prop_k)[death_m]) -
                           std::log(moveProbabilities(k)[birth_m]);
          log_like_ratio = logLikeRatioCalc(birth_m);
        break;
      }
      case death_m: {
          DEBUG_MSG("   Selected a death() move");
          death();
          log_prior_ratio = std::log((double) k) -
                            std::log((double) lambda_k);
          log_prop_ratio = std::log(moveProbabilities(k)[death_m]) -
                           std::log(moveProbabilities(prop_k)[birth_m]);
          log_like_ratio = logLikeRatioCalc(death_m);
        break;
      }
      case change_m: {
        DEBUG_MSG("   Selected a change() move");
        storeGraphState();
        death();
        log_like_ratio = logLikeRatioCalc(death_m);
        accept();
        birth();
        log_like_ratio += logLikeRatioCalc(birth_m);
        break;
      }
      case hyper_m: {
        DEBUG_MSG("   Selected a hyper() move");
        hyper();
        break;
      }
    }
    log_acc_prob = std::min(log_prior_ratio + log_prop_ratio + log_like_ratio, 0.0);
    if(unif_val < std::exp(log_acc_prob)) {
      DEBUG_MSG("   Accepting a move");
      accept();
    } else {
      if(current_move == change_m)
        restoreGraphState();
      DEBUG_MSG("   Rejecting a move");
      reject();
    }
    DEBUG_MSG("Finished sampleMH()");
  }

  void sampleOutputs() override {
    // Using std::valarray here; this may or may not be better than
    // operator overloading on std::vectors or using RcppEigen but for now
    // it is convenient enough
    DEBUG_MSG("Entering sampleOutputs()");
    // As a reminder, these are values that depend only on the current learner
    // First, update mu's
    std::valarray<double> clust_size(k);
    std::valarray<double> clust_response(k);
    for(int i = 0; i < k; ++i) {
      clust_size[i] = (double) clusts[i].size();
      //std::vector<double> response_by_clust = slice(learner_response, clusts[i]);
      clust_response[i] = sum(slice(learner_response, clusts[i]));
    }
    std::valarray<double> Qinv_diag = std::pow(((clust_size / sigmasq_y) + (1.0 / sigmasq_mu)), -1.0);
    std::valarray<double> b = (Qinv_diag * clust_response) / sigmasq_y;
    mu = r_std_norm(b, std::sqrt(Qinv_diag));
    output.mu_out = std::vector<double>(std::begin(mu), std::end(mu));

    // Next, update membership
    output.membership_out = membership;

    // Finally, update learner_response based on new mu values
    for(int i = 0; i < k; ++i) {
      std::valarray<std::size_t> clust_idx(clust_size[i]);
      std::copy(clusts[i].begin(), clusts[i].end(), std::begin(clust_idx));
      learner_response[clust_idx] = mu[i];
    }
    DEBUG_MSG("Finished sampleOutputs()");
  }

  void sampleBlocked(std::vector<BAST_Sampler> &learners) override {
    // Calculate Y_hat
    std::valarray<double> Y_hat(0.0, Y.size());
    for(int i = 0; i < n_learners; ++i) {
      // This will include the learner calling this method
      Y_hat += learners[i].learner_response;
    }
    // Residual backfitting: learner_response = Y - (Y_hat - learner_response)
    learner_response = Y - (Y_hat - learner_response);
  }

  void sampleConditioned(std::vector<BAST_Sampler> &learners) override {
    // These are values that depend on all the learners
    // This is only expected to be called by one learner; the rest
    // get updated using the input vector of samplers
    // The reason this is a method is so that it can use a consistent
    // PRNG generator (assuming the same learner calls this method every time)
    DEBUG_MSG("Entering sampleConditioned()");
    // First, calculate values that depend on the other learners
    std::valarray<double> Y_hat(0.0, Y.size());
    std::valarray<double> K_learners(1.0, n_learners);
    double sigmasq_y_new, log_post_new, sum_all_mu_squared = 0;

    for(int i = 0; i < n_learners; ++i) {
      // This will include the learner calling this method
      Y_hat += learners[i].learner_response;
      K_learners[i] = learners[i].k;
      sum_all_mu_squared += sum((std::valarray<double>) (std::pow(learners[i].mu, 2)));
    }

    // Calculate new sigmasq_y
    double shape = (double) (n + nu)/2.0;
    std::valarray<double> resid_sq = std::pow((Y - Y_hat), 2);
    double scale = std::pow(0.5*(nu*lambda_s + sum(resid_sq)), -1.0);
    sigmasq_y_new = std::pow(r_std_gamma(shape, scale), -1.0);

    // Calculate new log_post

    log_post_new = logPosteriorCalc(K_learners, Y_hat, sum_all_mu_squared);

    // Finally, pass the updated values to the other learners
    for(int i = 0; i < n_learners; ++i) {
      learners[i].sigmasq_y = sigmasq_y_new;
      learners[i].output.sigmasq_y_out = sigmasq_y_new;
      learners[i].log_post = log_post_new;
      learners[i].output.log_post_out = log_post_new;
    }
    DEBUG_MSG("Finished sampleConditioned()");
  }

  // Prior probability on moves as function of current number of clusters
  std::vector<double> moveProbabilities(int n_clusters) {
    if(n_clusters == 1) {
      return std::vector<double>{0.9, 0.0, 0.0, 0.1};
    } else if(n_clusters == max_clusts) {
      return std::vector<double>{0.0, 0.6, 0.3, 0.1};
    } else {
      return std::vector<double>{0.3, 0.3, 0.3, 0.1};
    }
  }
  // Generate move from discrete distribution on {0,1,2,3} for birth, death, change, hyper respectively
  int suggestMove() {
    std::vector<double> probs = moveProbabilities(k);
    return (std::discrete_distribution<int>(probs.begin(),
                                            probs.end())(generator));
  }

  // Sample vector of normal random variables with valarray of parameters
  std::valarray<double> r_std_norm(std::valarray<double> mean,
                                 std::valarray<double> sd) {
    if(mean.size() != sd.size()) {
      Rcpp::stop("Unable to draw random normals; mean and sd different length \n");
    }
    std::valarray<double> ran_vars(mean.size());
    for(int i = 0; i < mean.size(); ++i) {
      ran_vars[i] = (std::normal_distribution<double>(mean[i], sd[i]))(generator);
    }

    return ran_vars;
  }

  // Sample single gamma random variable; shape scale parameterization
  double r_std_gamma(double shape, double scale) {
    return (std::gamma_distribution<double>(shape, scale))(generator);
  }

  double logLikeRatioCalc(Move move) {
    double sigma_ratio = sigmasq_y / sigmasq_mu;
    int csize_old = old_clust_ids.size();
    int csize_new = new_clust_ids.size();
    double sum_resp_old = sum(slice(learner_response, old_clust_ids));
    double sum_resp_new = sum(slice(learner_response, new_clust_ids));
    double log_determinant_ratio = 0;
    double log_quad_ratio = 0;
    switch(move) {
      case birth_m: {
        log_determinant_ratio = -0.5 * (
          // Numerator
          std::log(csize_old + sigma_ratio) +
          std::log(csize_new + sigma_ratio) -
          // Denominator
          std::log(csize_old + csize_new + sigma_ratio) -
          std::log(sigma_ratio)
        );
        log_quad_ratio = (0.5/sigmasq_y) * (
          ((sum_resp_old*sum_resp_old)/(csize_old+sigma_ratio)) +
          ((sum_resp_new*sum_resp_new)/(csize_new+sigma_ratio)) -
          ((std::pow(sum_resp_new + sum_resp_old, 2))/(csize_old+csize_new+sigma_ratio))
        );
        break;
      }
      case death_m: {
        log_determinant_ratio = -0.5 * (
          // Numerator
          log(csize_new + sigma_ratio) +
          log(sigma_ratio) -
          // Denominator
          log(csize_old + sigma_ratio) -
          log(csize_new - csize_old + sigma_ratio)
        );
        log_quad_ratio = (0.5/sigmasq_y) * (
          ((sum_resp_new*sum_resp_new)/(csize_new+sigma_ratio)) -
          ((sum_resp_old*sum_resp_old)/(csize_old+sigma_ratio)) -
          ((std::pow(sum_resp_new - sum_resp_old, 2))/(csize_new-csize_old+sigma_ratio))
        );
        break;
      }
      default: {
        Rcpp::stop("logLikeRatioCalc() can only be called for a birth or death move\n");
      }
    }
    return(log_determinant_ratio + log_quad_ratio);
  }

  double logPosteriorCalc(std::valarray<double> K_learners,
                          std::valarray<double> Y_hat,
                          double sum_all_mu_squared) {
    std::valarray<double> log_bin_coeff = lchoose((n - 1), std::valarray<double>(K_learners - 1.0));
    double log_post_p = - ((1+(nu/2.0))*std::log(sigmasq_y))
                       - ((nu*lambda_s)/(2.0*sigmasq_y))
                       + sum((std::valarray<double>) (-1*(log_bin_coeff)
                              + (K_learners*std::log(lambda_k))
                              - lfactorial(K_learners)))
                        -(sum_all_mu_squared/(2.0*sigmasq_mu));
    double log_post_l = (-n/2.0)*log(sigmasq_y)
                      -(sum((std::valarray<double>) ((Y - Y_hat)*(Y - Y_hat)))/(2.0*sigmasq_y));
    return(log_post_p + log_post_l);
  }

};




#endif /* BAST_SAMPLER_HPP */


























