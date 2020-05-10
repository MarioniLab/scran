#include "Rcpp.h"

#include "beachmat/integer_matrix.h"
#include "beachmat/numeric_matrix.h"
#include "beachmat/utils/const_column.h"
#include "boost/range/algorithm.hpp"
#include "utils.h"
#include "rand_custom.h"

#include <stdexcept>
#include <algorithm>

template <class V>
double get_proportion (const V& expr, const int minpairs, const Rcpp::IntegerVector& marker1, const Rcpp::IntegerVector& marker2, const double threshold=NA_REAL) {
    int was_first=0, was_total=0;
    const bool short_cut = !ISNA(threshold);

    auto eIt=expr.begin();
    const size_t npairs=marker1.size();
    auto m1It=marker1.begin(), m2It=marker2.begin();

    size_t m=0;
    while (m < npairs) {
        const size_t limit=std::min(npairs, (short_cut ? m+100 : npairs)); // Checking every hundred pairs to avoid redundant calculations.
        for (; m<limit; ++m, ++m1It, ++m2It) {
            const auto& first=*(eIt+*m1It);
            const auto& second=*(eIt+*m2It);
            if (first != second) {
                if (first > second) { ++was_first; }
                ++was_total;
            }
        }

        // Returning if all we need to know is whether the score is greater than or less than 'threshold'.
        if (short_cut && was_total >= minpairs) {
            const size_t leftovers=npairs - m - 1;
            const double max_total=was_total + leftovers;
            const double n_thresh=max_total * threshold;

            // +1 to avoid incorrect early termination due to numerical imprecision upon equality.
            if (static_cast<double>(was_first + leftovers + 1) < n_thresh) {
                return -1;
            } else if (was_first && static_cast<double>(was_first - 1) > n_thresh) { // -1 for the same reason (need 'was_first' check to avoid underflow).
                return 1;
            }
        }
    }
    if (was_total < minpairs) { return NA_REAL; }
    
    const double output=static_cast<double>(was_first)/was_total;
    if (short_cut) {
        return output < threshold ? -1 : 1;
    }
    return output;
}

template <class M>
Rcpp::NumericVector cyclone_scores_internal (Rcpp::RObject input, Rcpp::IntegerVector mycells,
    Rcpp::IntegerVector marker1, Rcpp::IntegerVector marker2, Rcpp::IntegerVector used, 
    int niters, int miniters, int minpairs, Rcpp::List seeds, Rcpp::IntegerVector streams) 
{
    auto mat_ptr=beachmat::create_matrix<M>(input);
    const size_t ncells=mycells.size();
    const size_t ngenes=mat_ptr->get_nrow();
    const size_t nused=used.size();

    const size_t npairs=marker1.size();
    if (npairs!=static_cast<size_t>(marker2.size())) { 
        throw std::runtime_error("vectors of markers must be of the same length"); 
    }

    // Checking PCG setup.
    check_pcg_vectors(seeds, streams, mat_ptr->get_ncol(), "cells");

    // Checking marker sanity.    
    auto m2It=marker2.begin();
    for (auto m1It=marker1.begin(); m1It!=marker1.end(); ++m1It, ++m2It) {
        const int& m1m=(*m1It);
        if (m1m < 0 || static_cast<size_t>(m1m) >= nused) { 
            throw std::runtime_error("first marker indices are out of range"); 
        }

        const int& m2m=(*m2It);
        if (m2m < 0 || static_cast<size_t>(m2m) >= nused) {
            throw std::runtime_error("second marker indices are out of range"); 
        }
    }

    // Checking gene index sanity.
    for (auto uIt=used.begin(); uIt!=used.end(); ++uIt) { 
        const int& usedex=(*uIt);
        if (usedex < 0 || static_cast<size_t>(usedex) >= ngenes) { 
            throw std::runtime_error("used gene indices are out of range"); 
        }
    }

    Rcpp::NumericVector output(ncells, NA_REAL);
    typename M::vector current_exprs(nused);
    beachmat::const_column<M> col_holder(mat_ptr.get(), false); // need indexed access.

    auto oIt=output.begin();
    for (auto cIt=mycells.begin(); cIt!=mycells.end(); ++cIt, ++oIt) { 
        const size_t curcell=*cIt - 1;

        // Extracting only the expression values that are used in at least one pair.
        col_holder.fill(curcell);
        auto allIt=col_holder.get_values();
        auto curIt=current_exprs.begin();
        for (auto uIt=used.begin(); uIt!=used.end(); ++uIt, ++curIt) {
            (*curIt)=*(allIt + *uIt);
        }

        const double curscore=get_proportion(current_exprs, minpairs, marker1, marker2);
        if (ISNA(curscore)) { 
            continue;
        }

        // Iterations of shuffling to obtain a null distribution for the score.
        int below=0, total=0;
        auto generator=create_pcg32(seeds[curcell], streams[curcell]);
        for (int it=0; it < niters; ++it) {
            boost::range::random_shuffle(current_exprs, generator);
            const double newscore=get_proportion(current_exprs, minpairs, marker1, marker2, curscore);
            if (!ISNA(newscore)) { 
                if (newscore < 0) { ++below; }
                ++total;
            }
        }
       
        if (total >= miniters) { 
            (*oIt)=double(below)/total;
        }
    }

    return output;
}

// [[Rcpp::export(rng=false)]]
Rcpp::NumericVector cyclone_scores (Rcpp::RObject input, Rcpp::IntegerVector mycells,
    Rcpp::IntegerVector marker1, Rcpp::IntegerVector marker2, Rcpp::IntegerVector used, 
    int niters, int miniters, int minpairs, Rcpp::List seeds, Rcpp::IntegerVector streams) 
{
    int rtype=beachmat::find_sexp_type(exprs);
    if (rtype==INTSXP) {
        return cyclone_scores_internal<beachmat::integer_matrix>(exprs, mycells, marker1, marker2, indices, 
            niters, miniters, minpairs, seeds, streams);
    } else {
        return cyclone_scores_internal<beachmat::numeric_matrix>(exprs, mycells, marker1, marker2, indices, 
            niters, miniters, minpairs, seeds, streams);
    }
}

/* We could just assign ties random directions; then we'd only have to shuffle
 * once for all cells, and then we could use the same null distribution across
 * multiple cells, without worrying about whether or not one cell has more ties
 * than the other. The problem is that there's no protection from spuriously
 * high scores due to random breaking of ties; normally (for correlations),
 * we'd provide protection by controlling the type I error rate, but we're not
 * generating p-values here so it's harder to do.
 */
