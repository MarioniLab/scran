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
    std::vector<int> count(nblocks);

    if (!emat->is_sparse()) { 
        std::vector<double> incoming(ngenes);

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
    } else {
        Rcpp::NumericMatrix nnzero(ngenes, nblocks);
        std::vector<double> work_x(ngenes);
        std::vector<int> work_i(ngenes);
        auto smat = beachmat::promote_to_sparse(emat);

        // Running through Welford's algorithm for only nonzero counts.
        for (size_t counter=0; counter<ncells; ++counter) {
            auto idx = smat->get_col(counter, work_x.data(), work_i.data());
            trans(counter, idx.x, idx.x + idx.n, work_x.data());

            auto curblock=block[counter];
            if (isNA(curblock)) {
                continue;
            }
            ++count[curblock];

            auto M=outmean.column(curblock);
            auto S=outvar.column(curblock);
            auto NZ=nnzero.column(curblock);
            
            auto mIt=M.begin();
            auto sIt=S.begin();
            auto nzIt=NZ.begin();

            for (size_t i=0; i<idx.n; ++i) {
                auto& curM = *(mIt + idx.i[i]);
                auto& curS = *(sIt + idx.i[i]);
                const auto& curval = work_x[i];
                auto& curNZ = *(nzIt + idx.i[i]);
                ++curNZ;
                
                const double delta = curval - curM;
                curM += delta / curNZ;
                curS += delta * (curval - curM);
            }
        }

        // Filling in the zeros afterwards. This is done by realizing that s^2
        // = sum(y^2) - n * mean(y)^2 and that the sum of the squares does not
        // change with more zeroes; this allows us to solve for a new s^2 by
        // simply adding the difference of the scaled squared means.
        for (int b=0; b<nblocks; ++b) {
            auto M=outmean.column(b);
            auto S=outvar.column(b);
            auto NZ=nnzero.column(b);
            
            auto mIt=M.begin();
            auto sIt=S.begin();
            auto nzIt=NZ.begin();
            const double blocktotal = count[b];

            for (size_t g = 0; g<ngenes; ++g, ++mIt, ++sIt, ++nzIt) {
                const double curNZ = *nzIt;
                const double ratio = curNZ / blocktotal;
                auto& curM = *mIt;
                *sIt += curM * curM * ratio * (blocktotal - curNZ);
                curM *= ratio;
            }
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
