#ifndef RUN_DORMQR_H
#define RUN_DORMQR_H

#include "scran.h"

// Class to run Qy or QtY multiplications.

class run_dormqr {
public:
    run_dormqr(SEXP, SEXP, const char);
    void run();
    void run(double*);
    void solve();
    void solve(double*);

    int get_nobs() const;
    int get_ncoefs() const;

    // Where input and output goes.
    std::vector<double> rhs;
private:
    Rcpp::NumericMatrix QR;
    Rcpp::NumericVector AUX;
    const double* qrptr, * qxptr;
    const int nobs, ncoef, ncol;
    const char side, trans;
    int info, lwork;
    std::vector<double> work;
};

#endif
