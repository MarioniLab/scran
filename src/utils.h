#ifndef UTILS_H
#define UTILS_H

// Shuffle that mimics C++ shuffle, but using the R random number generator.

template <typename Pt>
void Rx_shuffle (Pt start, Pt end) {
    for (size_t i=(end-start)-1; i>0; --i) {
        std::swap(
                *(start+i), 
                *(start + int(R::unif_rand()*(i+1))) // Equal chance of getting anything from [0, i], as per random_shuffle in C++.
                ); 
    }
    return;
}

// Check subset vector.

Rcpp::IntegerVector check_subset_vector(SEXP, size_t);

// Overloaded functions to check for NA'ness.

bool isNA(int x);
bool isNA(double x);

// Checking for scalar inputs.

int check_integer_scalar(Rcpp::RObject, const char*);

double check_numeric_scalar(Rcpp::RObject, const char*);

bool check_logical_scalar(Rcpp::RObject, const char*);


#endif
