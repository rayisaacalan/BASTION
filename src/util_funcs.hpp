#ifndef UTIL_FUNCS_CPP
#define UTIL_FUNCS_CPP

#include <vector>
#include <valarray>
#include <cmath>
// Need lgamma from boost since std::lgamma is not thread safe
#include <boost/math/special_functions/gamma.hpp>

// From https://stackoverflow.com/a/67426730 ; similar to R style vector indexing
// std::valarray does this more eloquently by default
template<typename T>
std::vector<T> slice(const std::vector<T>& v, const std::vector<int>& id) {

  size_t id_size = id.size ();
  std::vector<T> tmp (id_size);
  T *tmp_data = tmp.data ();

  const int *id_data = id.data ();
  const T* v_data = v.data ();

  for (size_t i = 0; i < id_size; ++i) {
    tmp_data [i] = v_data [id_data [i]];
  }

  return tmp;
}

// For ease of use, allow valarray slicing by int vector rather than size_t valarray
template<typename T>
std::valarray<T> slice(const std::valarray<T>& v, const std::vector<int>& id) {

  std::valarray<std::size_t> slicing_ids(id.size());
  std::copy(id.begin(), id.end(), std::begin(slicing_ids));
  std::valarray<T> tmp = v[slicing_ids];

  return tmp;
}

// Generalization of std::reduce such that it takes only one argument
template <typename T>
double sum(const std::valarray<T> &input) {
  auto start_c = std::begin(input);
  auto end_c = std::end(input);
  return(std::reduce(start_c, end_c));
}

template <typename T>
double sum(const std::vector<T> &input) {
  auto start_c = std::begin(input);
  auto end_c = std::end(input);
  return(std::reduce(start_c, end_c));
}

// Log of binomial coefficient (n choose k)
template<class T>
T lchoose(T n, T k) {
  return (boost::math::lgamma(n + 1) - boost::math::lgamma(k + 1) - boost::math::lgamma(n - k + 1));
}

// Log of factorial function (n!)
template<class T>
T lfactorial(T n) {
  return boost::math::lgamma(n + 1);
}

// lfactorial but applied elementwise to a valarray
template<class T>
std::valarray<T> lfactorial(const std::valarray<T> &n_vals) {
  std::valarray<T> tmp = n_vals;
  for(auto &i : tmp) {
    i = lfactorial(i);
  }
  return(tmp);
}

// lgamma but applied elementwise to a valarray
template<class T>
std::valarray<T> lgamma(const std::valarray<T> &n_vals) {
  std::valarray<T> tmp = n_vals;
  for(auto &i : tmp) {
    i = boost::math::lgamma(i);
  }
  return(tmp);
}

// lchoose but for a valarray of k (ex k is valarray of double and n is int)
template<class T>
std::valarray<T> lchoose(int n, const std::valarray<T> &k) {
  return (boost::math::lgamma(n + 1) -
          lgamma((std::valarray<T>)(k + 1)) -
          lgamma((std::valarray<T>)((-1*k) + n + 1)));
}


#endif /* UTIL_FUNCS_CPP */

