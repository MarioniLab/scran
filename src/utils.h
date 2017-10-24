#ifndef UTILS_H
#define UTILS_H
#include <random>

// R-seeded random engine, to pass to C++ std::shuffle

class R_random_engine : public std::default_random_engine {
public:
    R_random_engine(bool=true);
    void reseed();
};

// Check subset vector.

Rcpp::IntegerVector check_subset_vector(SEXP, size_t);

// Overloaded functions to check for NA'ness.

bool isNA(int x);
bool isNA(double x);

#endif
