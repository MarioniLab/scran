#ifndef SCRAN_H
#define SCRAN_H

#include "R.h"
#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include "R_ext/Visibility.h"

#include <stdexcept>
#include <algorithm>

extern "C" {

SEXP forge_system (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP shuffle_scores (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP get_null_rho (SEXP, SEXP);

SEXP compute_rho(SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP auto_shuffle(SEXP, SEXP);

}

#include "utils.h"

#endif 
