#ifndef UTILS_H
#define UTILS_H

// R random number generator, to pass to C++ random_shuffle.

struct R_RNG {
    typedef size_t result_type;
    static constexpr result_type min() { return 0; }
    static constexpr result_type max() { return 999999; } // size_t probably holds 1 million okay. Setting this to size_t(-1)-1 doesn't work well.
    result_type operator()() const { 
        return std::min(R_RNG::max(), result_type(unif_rand()*double(R_RNG::max()+1))); // Enforce [min(), max()] output in case of failed truncation due to precision.
    }
};

// Check subset vector.

Rcpp::IntegerVector check_subset_vector(SEXP, size_t);

// Overloaded functions to check for NA'ness.

bool isNA(int x);
bool isNA(double x);

#endif
