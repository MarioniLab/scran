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

    Rcpp::NumericMatrix outvar(ngenes, nblocks);
    Rcpp::NumericMatrix outmean(ngenes, nblocks);
    Rcpp::NumericMatrix nnzero(ngenes, nblocks);
    std::vector<int> count(nblocks);

    std::vector<double> work_x(ngenes);
    std::vector<int> work_i;
    std::unique_ptr<beachmat::lin_sparse_matrix> smat;

    const bool is_sparse = emat->is_sparse(), preserve_zero = trans.preserve_zero();
    if (is_sparse && preserve_zero) {
        work_i.resize(ngenes);
        smat = beachmat::promote_to_sparse(emat);
    }

    // Using Welford's algorithm to compute the variance in column-major style,
    // which should be much more cache-friendly for large matrices. We run through
    // only non-zero counts when circumstances allow for it.
    for (size_t counter=0; counter<ncells; ++counter) {
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

        if (is_sparse && preserve_zero) {
            auto idx = smat->get_col(counter, work_x.data(), work_i.data());
            trans(counter, idx.x, idx.x + idx.n, work_x.data());

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
        } else {
            auto ptr = emat->get_col(counter, work_x.data());
            trans(counter, ptr, ptr + ngenes, work_x.data());

            auto wIt = work_x.begin();
            for (size_t i=0; i<ngenes; ++i, ++wIt, ++mIt, ++sIt, ++nzIt) {
                if (!preserve_zero || *wIt) {
                    ++(*nzIt);
                    const double delta=*wIt - *mIt;
                    *mIt += delta/(*nzIt);
                    *sIt += delta*(*wIt - *mIt);
                }
            }
        } 
    }

    for (int b=0; b<nblocks; ++b) {
        auto M=outmean.column(b);
        if (count[b] <= 0) {
            std::fill(M.begin(), M.end(), R_NaReal);
        }

        auto S=outvar.column(b);
        if (count[b] <= 1) {
            std::fill(S.begin(), S.end(), R_NaReal);
            continue; 
        }

        if (preserve_zero) {
            // Filling in the zeros afterwards. This is done by realizing that s^2
            // = sum(y^2) - n * mean(y)^2 and that the sum of the squares does not
            // change with more zeroes; this allows us to solve for a new s^2 by
            // simply adding the difference of the scaled squared means.
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

        const double rdf = count[b] - 1;
        for (auto& s : S) {
            s/=rdf;
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

    bool preserve_zero()  const { return pseudo == 1; }
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

    bool preserve_zero () const { return true; }
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

    bool preserve_zero () const { return true; }
};

// [[Rcpp::export(rng=false)]]
Rcpp::List compute_blocked_stats_none(Rcpp::RObject mat, Rcpp::IntegerVector block, int nblocks) {
    none N;
    return compute_blocked_stats(mat, block, nblocks, N);
}
