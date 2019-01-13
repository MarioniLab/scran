#include "utils.h"

#include <stdexcept>

Rcpp::IntegerVector check_subset_vector(SEXP subvec, size_t len) {
    Rcpp::IntegerVector sout(subvec);
    for (const auto s : sout) {
        if (isNA(s) || s < 0 || static_cast<size_t>(s) >= len) {
            throw std::runtime_error("subset indices out of range");
        }
    }
    return sout;
}

// Special function to check for NA'ness.

bool isNA(int x) {
    return x==NA_INTEGER;
}

bool isNA(double x) {
    return ISNAN(x);
}


// Checking for scalar inputs.

template<typename T, class V>
T check_scalar(Rcpp::RObject incoming, const char* arg, const char* val) {
    V vec(incoming);
    if (vec.size()!=1) {
        std::stringstream err;
        err << arg << " should be " << val;
        throw std::runtime_error(err.str());
    }
    return vec[0];
}

int check_integer_scalar(Rcpp::RObject incoming, const char* arg) {
    return check_scalar<int, Rcpp::IntegerVector>(incoming, arg, "an integer scalar");
}

double check_numeric_scalar(Rcpp::RObject incoming, const char* arg) {
    return check_scalar<double, Rcpp::NumericVector>(incoming, arg, "a numeric scalar");
}

bool check_logical_scalar(Rcpp::RObject incoming, const char* arg) {
    return check_scalar<bool, Rcpp::LogicalVector>(incoming, arg, "a logical scalar");
}
