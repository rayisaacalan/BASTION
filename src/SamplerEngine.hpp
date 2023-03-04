#ifndef SAMPLERENGINE_HPP
#define SAMPLERENGINE_HPP

#include <vector>
#include "Sampler.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif
//[[Rcpp::plugins(openmp)]]

//[[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

// This definition enables debug statements to Rcout
// #define SAMPENG_DEBUG

#ifdef SAMPENG_DEBUG
#define DEBUG_MSG(str) do { Rcpp::Rcout << str << std::endl; } while( false )
#else
#define DEBUG_MSG(str) do { } while ( false )
#endif

template<typename T, typename output_T, typename = std::enable_if_t<std::is_base_of<Sampler<output_T, T>, T>::value>>
class SamplerEngine {
public:
  // Engine parameters
  const int n_learners;
  const int iterations;
  const int burn_in;
  const int thin;
  const bool blocked;
  const int threads;
  // Progress and interrupt monitor
  Progress prog_monitor;
  // Outer vector is for each learner; inner vector is for each output sample
  std::vector<std::vector<output_T>> outputs;
  std::vector<T> learners;

  SamplerEngine(std::vector<T> &inp_learners,
                const int MCMC,
                const int BURNIN,
                const int THINNINGINT,
                const bool BLOCK,
                const int THREAD) :
    n_learners(inp_learners.size()),
    iterations(MCMC),
    burn_in(BURNIN),
    thin(THINNINGINT),
    blocked(BLOCK),
    threads(THREAD),
    prog_monitor(MCMC, true)
    // May want to add another input instead of fixing at "true"
    // in case someone doesn't want progress bar displayed
  {
    DEBUG_MSG("Beginning to build a SamplerEngine");
    learners = inp_learners;
    for(int i = 0; i < n_learners; ++i) {
      outputs.push_back(std::vector<output_T>());
    }
    #ifdef _OPENMP
    if(threads > 0)
      omp_set_num_threads(std::min(threads, n_learners));
    else
      omp_set_num_threads(1);
    #endif
    DEBUG_MSG("Built a SamplerEngine with " << learners.size() << " learners");
  }

  void runEngine() {
    // Outermost loop: each MCMC iteration
    DEBUG_MSG("Starting to run SamplerEngine");
    for(int iter = 0; iter < iterations; ++iter) {
      // Check for interrupts from R
      if(Progress::check_abort())
          Rcpp::stop("Execution interrupted by R\n");
      prog_monitor.increment();
      // If using a blocked sampler, that step occurs here
      if(blocked) {
        DEBUG_MSG("****Beginning blocked sampling****");
        for(int learner = 0; learner < n_learners; ++learner) {
          // Note that this can *not* be parallelized
          learners[learner].sampleBlocked(learners);
        }
        DEBUG_MSG("****Finished blocked sampling****");
        // Repeat steps THIN times
        DEBUG_MSG("****Beginning MH/Output sampling****");
        #ifdef _OPENMP
        #pragma omp parallel
        #pragma omp for
        #endif
        for(int learner = 0; learner < n_learners; ++learner) {
          // Now do the MH updates and posterior samples (possibly in parallel)
          for(int repetition = 0; repetition < thin; ++repetition) {
            learners[learner].sampleMH();
            learners[learner].sampleOutputs();
          }
        }
        DEBUG_MSG("****Finished MH/Output sampling****");
        // If necessary, sample other constant parameters conditioned on
        // all learners here
        DEBUG_MSG("****Beginning conditioned sampling****");
        learners[0].sampleConditioned(learners);
        DEBUG_MSG("****Finished conditioned sampling****");
      } else {
        DEBUG_MSG("****Beginning MH/Output sampling****");
        for(int learner = 0; learner < n_learners; ++learner) {
          learners[learner].sampleBlocked(learners);
          learners[learner].sampleMH();
          learners[learner].sampleOutputs();
        }
        DEBUG_MSG("****Finished MH/Output sampling****");
        // If not blocked; just do one sampling (from first learner)
        // and update values for remaining learners
        DEBUG_MSG("****Beginning conditioned sampling****");
        learners[0].sampleConditioned(learners);
        DEBUG_MSG("****Finished conditioned sampling****");
      }

      // If relevant, output the samples
      // Either blocked sampler in which case output all samples after burn in
      if((blocked && (iter + 1 > burn_in)) ||
         // Or nonblocked sampler in which case output every thin iteration after burn in
         (!blocked && (iter + 1 > burn_in) && (((iter + 1 - burn_in) % thin) == 0))) {
        DEBUG_MSG("Writing a sample");
        #ifdef _OPENMP
        #pragma omp parallel
        #pragma omp for
        #endif
        for(int learner = 0; learner < n_learners; ++learner) {
          outputs[learner].push_back(learners[learner].output);
        }
        DEBUG_MSG("Wrote a sample");
      }
    }
  }
};


#endif /* SAMPLERENGINE_HPP */

