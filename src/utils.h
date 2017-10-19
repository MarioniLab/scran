#ifndef UTILS_H
#define UTILS_H

// R random number generator, to pass to C++ random_shuffle.

struct R_RNG {
    typedef size_t result_type;
    static constexpr result_type largest=result_type(-1); // Getting the maximum value.
    static constexpr result_type min() { return 0; }
    static constexpr result_type max() { return largest-1; } 
    result_type operator()() const { return unif_rand()*largest; }
};

// Check subset vector.

Rcpp::IntegerVector check_subset_vector(SEXP, size_t);

// Overloaded functions to check for NA'ness.

bool isNA(int x);
bool isNA(double x);

#endif
