#ifndef SAMPLER_HPP
#define SAMPLER_HPP

#include <random>

template <class output_T, class learner_T>
class Sampler {
public:
  output_T output;

  virtual void sampleMH() = 0;
  virtual void sampleOutputs() = 0;
  virtual void sampleBlocked(std::vector<learner_T> &learners) {};
  virtual void sampleConditioned(std::vector<learner_T> &learners) {};

  Sampler(output_T init_output) :
    output(init_output) {}
};



#endif /* SAMPLER_HPP */


