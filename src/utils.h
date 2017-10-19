#ifndef UTILS_H
#define UTILS_H
#include <random>

// R-seeded random engine, to pass to C++ std::shuffle

std::default_random_engine setup_random_engine();

// Check subset vector.

Rcpp::IntegerVector check_subset_vector(SEXP, size_t);

// Overloaded functions to check for NA'ness.

bool isNA(int x);
bool isNA(double x);

#endif
