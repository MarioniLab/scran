#include "Rcpp.h"
#include "beachmat/numeric_matrix.h"
#include "beachmat/integer_matrix.h"

#include "utils.h"

#include <vector>
#include <algorithm>
#include <cmath>

/********************************************
 * Template to compute variance statistics for 
 * a given transformation of the counts.
 *********************************************/

template<class MAT, class TRANSFORMER>
Rcpp::List compute_blocked_stats(Rcpp::RObject mat, Rcpp::IntegerVector block, int nblocks, TRANSFORMER trans) {
    auto emat=beachmat::create_matrix<MAT>(mat);
    const size_t ncells=emat->get_ncol();
    const size_t ngenes=emat->get_nrow();

    // Setting up the output objects.
    Rcpp::NumericMatrix outvar(ngenes, nblocks);
    Rcpp::NumericMatrix outmean(ngenes, nblocks);
    Rcpp::NumericVector incoming(ngenes);
    std::vector<int> count(nblocks);

    // Using Welford's algorithm to compute the variance in column-major style,
    // which should be much more cache-friendly for large matrices.
    for (size_t counter=0; counter<ncells; ++counter) {
        emat->get_col(counter, incoming.begin());
        trans(counter, incoming.begin(), incoming.end());

        auto curblock=block[counter];
        if (isNA(curblock)) {
            continue;
        }

        auto M=outmean.column(curblock);
        auto S=outvar.column(curblock);
        auto& sofar=count[curblock];
        
        // Running mean.
        auto mIt=M.begin();
        auto sIt=S.begin();
        auto iIt=incoming.begin();
        ++sofar;

        for (size_t i=0; i<ngenes; ++i, ++iIt, ++mIt, ++sIt) {
            const double delta=*iIt - *mIt;
            *mIt += delta/sofar;
            *sIt += delta*(*iIt - *mIt);
        }
    }

    for (int b=0; b<nblocks; ++b) {
        auto S=outvar.column(b);
        int N=count[b]-1;
        if (N) {
            for (auto& s : S) {
                s/=N;
            }
        } else {
            std::fill(S.begin(), S.end(), R_NaReal);
        }
    }

    return(Rcpp::List::create(outmean, outvar));
}

/************************************************
 * Compute statistics for log-transformed counts.
 ***********************************************/

struct lognorm {
    lognorm(Rcpp::NumericVector sizefactors, double p) : sf(sizefactors), pseudo(p) {}
    void operator()(int i, Rcpp::NumericVector::iterator start, Rcpp::NumericVector::iterator end) {
        double target=sf[i];
        while (start!=end) {
            *start=std::log(*start/target + pseudo)/M_LN2;
            ++start;
        }
    }
private:
    Rcpp::NumericVector sf;
    double pseudo;
};


// [[Rcpp::export(rng=false)]]
Rcpp::List compute_blocked_stats_lognorm(Rcpp::RObject mat, Rcpp::IntegerVector block, 
    int nblocks, Rcpp::NumericVector sf, double pseudo) 
{
    int rtype=beachmat::find_sexp_type(mat);
    lognorm LN(sf, pseudo);
    if (rtype==INTSXP) {
        return compute_blocked_stats<beachmat::integer_matrix>(mat, block, nblocks, LN);
    } else {
        return compute_blocked_stats<beachmat::numeric_matrix>(mat, block, nblocks, LN);
    }
}

/*******************************************
 * Compute statistics for normalized counts.
 *******************************************/

struct norm {
    norm(Rcpp::NumericVector sizefactors) : sf(sizefactors) {}
    void operator()(int i, Rcpp::NumericVector::iterator start, Rcpp::NumericVector::iterator end) {
        double target=sf[i];
        while (start!=end) {
            *start/=target;
            ++start;
        }
    }
private:
    Rcpp::NumericVector sf;
};

// [[Rcpp::export(rng=false)]]
Rcpp::List compute_blocked_stats_norm(Rcpp::RObject mat, Rcpp::IntegerVector block, 
    int nblocks, Rcpp::NumericVector sf)
{
    int rtype=beachmat::find_sexp_type(mat);
    norm N(sf);
    if (rtype==INTSXP) {
        return compute_blocked_stats<beachmat::integer_matrix>(mat, block, nblocks, N);
    } else {
        return compute_blocked_stats<beachmat::numeric_matrix>(mat, block, nblocks, N);
    }
}

/***********************************************
 * Compute statistics for expression as provided.
 ***********************************************/

struct none {
    none() {}
    void operator()(int i, Rcpp::NumericVector::iterator start, Rcpp::NumericVector::iterator end) {}
};

// [[Rcpp::export(rng=false)]]
Rcpp::List compute_blocked_stats_none(Rcpp::RObject mat, Rcpp::IntegerVector block, int nblocks) {
    int rtype=beachmat::find_sexp_type(mat);
    none N;
    if (rtype==INTSXP) {
        return compute_blocked_stats<beachmat::integer_matrix>(mat, block, nblocks, N);
    } else {
        return compute_blocked_stats<beachmat::numeric_matrix>(mat, block, nblocks, N);
    }
}
