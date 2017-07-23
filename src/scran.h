#ifndef SCRAN_H
#define SCRAN_H

#include "beachmat/integer_matrix.h"
#include "beachmat/numeric_matrix.h"
#include "Rcpp.h"

#include <stdexcept>
#include <algorithm>
#include <deque>
#include <vector>
#include <cmath>

extern "C" {

SEXP forge_system (SEXP, SEXP, SEXP, SEXP);

SEXP shuffle_scores (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP get_null_rho (SEXP, SEXP);

SEXP get_null_rho_design (SEXP, SEXP, SEXP);

SEXP compute_rho(SEXP, SEXP, SEXP, SEXP);

SEXP combine_corP(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP auto_shuffle(SEXP, SEXP);

SEXP get_scaled_ranks(SEXP, SEXP, SEXP);

SEXP subset_and_divide(SEXP, SEXP, SEXP); 

SEXP fit_linear_model(SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP rank_subset(SEXP, SEXP, SEXP, SEXP);

SEXP get_residuals(SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP compute_CV2(SEXP, SEXP, SEXP, SEXP);

SEXP build_snn(SEXP);

SEXP overlap_exprs(SEXP, SEXP, SEXP, SEXP);

}

#include "utils.h"

#endif 
