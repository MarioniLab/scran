#include "Rcpp.h"
#include "beachmat/numeric_matrix.h"
#include "beachmat/integer_matrix.h"

#include "run_dormqr.h"

template<class M, class TRANSFORMER>
Rcpp::List compute_residual_stats(Rcpp::RObject qr, Rcpp::RObject qraux, SEXP inmat, TRANSFORMER trans) {
    auto emat=beachmat::create_matrix<M>(inmat);
    const size_t ncells=emat->get_ncol();
    const size_t ngenes=emat->get_nrow();

    residual_stats_calculator RSC(qr, qraux);
   
    // Setting up the output objects.
    Rcpp::NumericMatrix outvar(1, ngenes);
    Rcpp::NumericMatrix outmean(1, ngenes);
    Rcpp::NumericVector incoming(ncells);

    for (size_t counter=0; counter<ngenes; ++counter) {
        emat->get_row(counter, incoming.begin());
        trans(incoming.begin(), incoming.end());
        auto curvarrow=outvar.column(counter);
        auto curmeanrow=outmean.column(counter);
        RSC.compute(incoming.begin(), curmeanrow.begin(),  curvarrow.begin());
    }

    return(Rcpp::List::create(outmean, outvar));
}

/************************************************
 * Compute statistics for log-transformed counts.
 ***********************************************/

struct lognorm {
    lognorm(Rcpp::NumericVector sizefactors, double pseudo) : sf(sizefactors), ps(pseudo) {}
    void operator()(Rcpp::NumericVector::iterator start, Rcpp::NumericVector::iterator end) {
        auto sfIt=sf.begin();
        while (start!=end) {
            *start=std::log(*start/(*sfIt) + ps)/M_LN2;
            ++start;
            ++sfIt;
        }
    }
private:
    Rcpp::NumericVector sf;
    double ps;
};

// [[Rcpp::export(rng=false)]]
Rcpp::List compute_residual_stats_lognorm(Rcpp::RObject qr, Rcpp::RObject qraux, SEXP inmat,
    Rcpp::NumericVector sf, double pseudo)
{
    int rtype=beachmat::find_sexp_type(inmat);
    lognorm LN(sf, pseudo);
    if (rtype==INTSXP) {
        return compute_residual_stats<beachmat::integer_matrix>(qr, qraux, inmat, LN);
    } else {
        return compute_residual_stats<beachmat::numeric_matrix>(qr, qraux, inmat, LN);
    }
}

/***********************************************
 * Compute statistics for expression as provided.
 ***********************************************/

struct none {
    none() {}
    void operator()(Rcpp::NumericVector::iterator start, Rcpp::NumericVector::iterator end) {}
};

// [[Rcpp::export(rng=false)]]
Rcpp::List compute_residual_stats_none(Rcpp::RObject qr, Rcpp::RObject qraux, SEXP inmat) {
    int rtype=beachmat::find_sexp_type(inmat);
    none N;
    if (rtype==INTSXP) {
        return compute_residual_stats<beachmat::integer_matrix>(qr, qraux, inmat, N);
    } else {
        return compute_residual_stats<beachmat::numeric_matrix>(qr, qraux, inmat, N);
    }
}
