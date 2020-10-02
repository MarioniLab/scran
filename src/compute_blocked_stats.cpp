#include "Rcpp.h"
#include "beachmat3/beachmat.h"
#include "utils.h"

#include <vector>
#include <algorithm>
#include <cmath>

/********************************************
 * Template to compute variance statistics for 
 * a given transformation of the counts.
 *********************************************/

template<class TRANSFORMER>
Rcpp::List compute_blocked_stats(Rcpp::RObject mat, Rcpp::IntegerVector block, int nblocks, TRANSFORMER trans) {
    auto emat = beachmat::read_lin_block(mat);
    const size_t ncells=emat->get_ncol();
    const size_t ngenes=emat->get_nrow();

    // Setting up the output objects.
    Rcpp::NumericMatrix outvar(ngenes, nblocks);
    Rcpp::NumericMatrix outmean(ngenes, nblocks);
    std::vector<double> incoming(ngenes);
    std::vector<int> count(nblocks);

    // Using Welford's algorithm to compute the variance in column-major style,
    // which should be much more cache-friendly for large matrices.
    for (size_t counter=0; counter<ncells; ++counter) {
        auto ptr = emat->get_col(counter, incoming.data());
        trans(counter, ptr, ptr + ngenes, incoming.data());

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
        auto M=outmean.column(b);
        if (count[b] <= 0) {
            std::fill(M.begin(), M.end(), R_NaReal);
        }

        auto S=outvar.column(b);
        if (count[b] > 1) {
            const double rdf=count[b]-1;
            for (auto& s : S) {
                s/=rdf;
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

    template<class IN, class OUT>
    void operator()(int i, IN start, IN end, OUT out) {
        double target=sf[i];
        while (start!=end) {
            *out=std::log(*start/target + pseudo)/M_LN2;
            ++start;
            ++out;
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
    lognorm LN(sf, pseudo);
    return compute_blocked_stats(mat, block, nblocks, LN);
}

/*******************************************
 * Compute statistics for normalized counts.
 *******************************************/

struct norm {
    norm(Rcpp::NumericVector sizefactors) : sf(sizefactors) {}

    template<class IN, class OUT>
    void operator()(int i, IN start, IN end, OUT out) {
        double target=sf[i];
        while (start!=end) {
            (*out) = (*start)/target;
            ++start;
            ++out;
        }
    }
private:
    Rcpp::NumericVector sf;
};

// [[Rcpp::export(rng=false)]]
Rcpp::List compute_blocked_stats_norm(Rcpp::RObject mat, Rcpp::IntegerVector block, 
    int nblocks, Rcpp::NumericVector sf)
{
    norm N(sf);
    return compute_blocked_stats(mat, block, nblocks, N);
}

/***********************************************
 * Compute statistics for expression as provided.
 ***********************************************/

struct none {
    none() {}
    template<class IN, class OUT>
    void operator()(int i, IN start, IN end, OUT out) {
        if (out!=start) {
            std::copy(start, end, out);
        }
    }
};

// [[Rcpp::export(rng=false)]]
Rcpp::List compute_blocked_stats_none(Rcpp::RObject mat, Rcpp::IntegerVector block, int nblocks) {
    none N;
    return compute_blocked_stats(mat, block, nblocks, N);
}
