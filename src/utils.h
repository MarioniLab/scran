#ifndef UTILS_H
#define UTILS_H

// R random number generator, to pass to C++ random_shuffle.

struct R_RNG {
    typedef size_t result_type;

    static constexpr result_type min() { return 0; }

    /* size_t holds 1 million on any platform where R runs, so this shouldn't be an issue. 
     * Setting this to size_t(-1)-1 seems more portable, but it doesn't work for some reason.
     * Possibly becaues it can't be held precisely upon coercion into a double()? */
    static constexpr result_type max() { return 999999; } 

    result_type operator()() const { 
        /* Enforce [min(), max()] output in case of failed truncation.
         * Theoretically possible due to numerical imprecision with doubles. */
        return std::min(R_RNG::max(), result_type(unif_rand()*double(R_RNG::max()+1))); 
    }
};

// Check subset vector.

Rcpp::IntegerVector check_subset_vector(SEXP, size_t);

// Overloaded functions to check for NA'ness.

bool isNA(int x);
bool isNA(double x);

#endif
