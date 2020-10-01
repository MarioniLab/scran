#include "Rcpp.h"
#include "beachmat3/beachmat.h"
#include "scuttle/linear_model_fit.h"

template<class TRANSFORMER>
Rcpp::List compute_residual_stats(Rcpp::NumericMatrix qr, Rcpp::NumericVector qraux, Rcpp::RObject inmat, TRANSFORMER trans) {
    auto emat = beachmat::read_lin_block(inmat);
    const size_t ncells=emat->get_ncol();
    const size_t ngenes=emat->get_nrow();

    scuttle::linear_model_fit fitter(qr, qraux);
    const size_t ncoefs=fitter.get_ncoefs();

    // Setting up the output objects.
    Rcpp::NumericMatrix outvar(1, ngenes);
    Rcpp::NumericMatrix outmean(1, ngenes);
    Rcpp::NumericVector incoming(ncells);

    for (size_t counter=0; counter<ngenes; ++counter) {
        auto iIt = incoming.begin();
        auto ptr = emat->get_row(counter, iIt);
        trans(ptr, ptr + ncells, iIt);

        auto curvarrow=outvar.column(counter);
        auto curvar=curvarrow.begin();
        auto curmeanrow=outmean.column(counter);
        auto curmean=curmeanrow.begin();

        auto iEnd = incoming.end();
        (*curmean)=std::accumulate(iIt, iEnd, 0.0)/ncells;
        fitter.multiply(iIt);

        double& v=(*curvar);
        iIt+=ncoefs;
        while (iIt != iEnd) { // only using the residual effects.
            v += (*iIt) * (*iIt);
            ++iIt;
        }
        v /= ncells - ncoefs;
    }

    return Rcpp::List::create(outmean, outvar);
}

/************************************************
 * Compute statistics for log-transformed counts.
 ***********************************************/

struct lognorm {
    lognorm(Rcpp::NumericVector sizefactors, double pseudo) : sf(sizefactors), ps(pseudo) {}

    template<class IN, class OUT>
    void operator()(IN start, IN end, OUT out) {
        auto sfIt = sf.begin();
        while (start != end) {
            *out = std::log(*start/(*sfIt) + ps)/M_LN2;
            ++start;
            ++sfIt;
            ++out;
        }
    }
private:
    Rcpp::NumericVector sf;
    double ps;
};

// [[Rcpp::export(rng=false)]]
Rcpp::List compute_residual_stats_lognorm(Rcpp::NumericMatrix qr, Rcpp::NumericVector qraux, Rcpp::RObject inmat,
    Rcpp::NumericVector sf, double pseudo)
{
    lognorm LN(sf, pseudo);
    return compute_residual_stats(qr, qraux, inmat, LN);
}

/***********************************************
 * Compute statistics for expression as provided.
 ***********************************************/

struct none {
    none() {}

    template<class IN, class OUT>
    void operator()(IN start, IN end, OUT out) {
        if (out!=start) {
            std::copy(start, end, out);
        }
    }
};

// [[Rcpp::export(rng=false)]]
Rcpp::List compute_residual_stats_none(Rcpp::NumericMatrix qr, Rcpp::NumericVector qraux, Rcpp::RObject inmat) {
    none N;
    return compute_residual_stats(qr, qraux, inmat, N);
}
