#include "scran.h"

template <class V>
double get_proportion (const V& expr, const int& minpairs, const Rcpp::IntegerVector& marker1, const Rcpp::IntegerVector& marker2) {
    int was_first=0, was_total=0;
    auto m2It=marker2.begin();
    auto eIt=expr.begin();
    for (auto m1It=marker1.begin(); m1It!=marker1.end(); ++m1It, ++m2It) {
        const auto& first=*(eIt+*m1It);
        const auto& second=*(eIt+*m2It);
        if (first > second) { ++was_first; }
        if (first != second) { ++was_total; }      
    }
    if (was_total < minpairs) { return NA_REAL; }
    return double(was_first)/was_total;
}

template <class V, class M>
SEXP shuffle_scores_internal (M mat_ptr, 
        Rcpp::IntegerVector mycells,
        Rcpp::IntegerVector marker1, Rcpp::IntegerVector marker2, Rcpp::IntegerVector used, 
        Rcpp::IntegerVector iter, Rcpp::IntegerVector miniter, Rcpp::IntegerVector minpair) {
   
    const size_t ncells=mycells.size();
    const size_t ngenes=mat_ptr->get_nrow();

    const size_t npairs=marker1.size();
    if (npairs!=marker2.size()) { throw std::runtime_error("vectors of markers must be of the same length"); }
    const size_t nused=used.size();

    if (iter.size()!=1) { throw std::runtime_error("number of iterations must be an integer scalar"); }
    const int nit=iter[0];
    if (miniter.size()!=1) { throw std::runtime_error("minimum number of iterations must be an integer scalar"); }
    const int minit=miniter[0];
    if (minpair.size()!=1) { throw std::runtime_error("minimum number of pairs must be an integer scalar"); }
    const int minp=minpair[0];

    // Checking marker sanity.    
    auto m2It=marker2.begin();
    for (auto m1It=marker1.begin(); m1It!=marker1.end(); ++m1It, ++m2It) {
        const int& m1m=(*m1It);
        if (m1m >= nused || m1m < 0) { throw std::runtime_error("first marker indices are out of range"); }
        const int& m2m=(*m2It);
        if (m2m >= nused || m2m < 0) { throw std::runtime_error("second marker indices are out of range"); }
    }

    // Checking gene index sanity.
    for (auto uIt=used.begin(); uIt!=used.end(); ++uIt) { 
        const int& usedex=(*uIt);
        if (usedex >= ngenes || usedex < 0) { throw std::runtime_error("used gene indices are out of range"); }
    }

    Rcpp::NumericVector output(ncells, NA_REAL);
    V all_exprs(ngenes), current_exprs(nused);
    Rcpp::RNGScope rng; // Place after initialization of all Rcpp objects.

    auto oIt=output.begin();
    for (auto cIt=mycells.begin(); cIt!=mycells.end(); ++cIt, ++oIt) { 

        // Extracting only the expression values that are used in at least one pair.
        auto inIt=mat_ptr->get_const_col(*cIt - 1, all_exprs.begin());
        auto curIt=current_exprs.begin();
        for (auto uIt=used.begin(); uIt!=used.end(); ++uIt, ++curIt) {
            (*curIt)=*(inIt + *uIt);
        }
        
        const double curscore=get_proportion(current_exprs, minp, marker1, marker2);
        if (ISNA(curscore)) { 
            continue;
        }

        // Iterations of shuffling to obtain a null distribution for the score.
        int below=0, total=0;
        for (int it=0; it < nit; ++it) {
            Rx_shuffle(current_exprs.begin(), current_exprs.end());
            const double newscore=get_proportion(current_exprs, minp, marker1, marker2);
            if (!ISNA(newscore)) { 
                if (newscore < curscore) { ++below; }
                ++total;
            }
        }
       
        if (total >= minit) { 
            (*oIt)=double(below)/total;
        }
    }

    return output;
}

SEXP shuffle_scores(SEXP mycells, SEXP exprs, SEXP marker1, SEXP marker2, SEXP indices, SEXP iter, SEXP miniter, SEXP minpair) {
    BEGIN_RCPP
    int rtype=beachmat::find_sexp_type(exprs);
    if (rtype==INTSXP) {
        auto mat=beachmat::create_integer_matrix(exprs);
        return shuffle_scores_internal<Rcpp::IntegerVector>(mat.get(), mycells, marker1, marker2, indices, iter, miniter, minpair);
    } else {
        auto mat=beachmat::create_numeric_matrix(exprs);
        return shuffle_scores_internal<Rcpp::NumericVector>(mat.get(), mycells, marker1, marker2, indices, iter, miniter, minpair);
    }
    END_RCPP
}

/* We could just assign ties random directions; then we'd only have to shuffle
 * once for all cells, and then we could use the same null distribution across
 * multiple cells, without worrying about whether or not one cell has more ties
 * than the other. The problem is that there's no protection from spuriously
 * high scores due to random breaking of ties; normally (for correlations),
 * we'd provide protection by controlling the type I error rate, but we're not
 * generating p-values here so it's harder to do.
 */

SEXP auto_shuffle(SEXP incoming, SEXP nits) {
    BEGIN_RCPP

    const int niters=Rcpp::IntegerVector(nits)[0];
    const Rcpp::NumericVector invec(incoming);
    const size_t N=invec.size();
    Rcpp::NumericMatrix outmat(N, niters);
    Rcpp::RNGScope rng; // Place after initialization of all Rcpp objects.

    Rcpp::NumericVector::const_iterator source=invec.begin();
    Rcpp::NumericVector::iterator oIt=outmat.begin();
    
    for (int i=0; i<niters; ++i) {
        std::copy(source, source+N, oIt);
        Rx_shuffle(oIt, oIt+N);
        source=oIt;
        oIt+=N;
    }

    return outmat;
    END_RCPP
}
