#ifndef UTILS_H
#define UTILS_H

#include "Rcpp.h"

// Overloaded functions to check for NA'ness.

inline bool isNA(int x) {
    return x==NA_INTEGER;
}

inline bool isNA(double x) {
    return ISNAN(x);
}

#endif
