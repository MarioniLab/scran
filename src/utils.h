#ifndef UTILS_H
#define UTILS_H

// Shuffle that mimics C++ shuffle, but using the R random number generator.

template <typename Pt>
void Rx_shuffle (Pt start, Pt end) {
    size_t n=(end-start), i;
    for (i=n-1; i>0; --i) {
        std::swap(start[i], start[int(unif_rand()*(i+1))]); // Equal chance of getting anything from [0, i], as per random_shuffle in C++.
    }
    return;
}

// Check subset vector.

Rcpp::IntegerVector check_subset_vector(SEXP, size_t);

// Overloaded functions to check for NA'ness.

bool isNA(int x);
bool isNA(double x);

#endif
