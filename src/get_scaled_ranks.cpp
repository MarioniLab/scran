#include "scran.h"

#include "beachmat/numeric_matrix.h"
#include "beachmat/integer_matrix.h"
#include "beachmat/utils/const_column.h"
#include "utils.h"

#include <vector>
#include <utility>
#include <algorithm>
#include <stdexcept>

template <class M> 
SEXP average_ranks_internal(SEXP input, SEXP subset, SEXP transpose, SEXP as_sparse) { 
    auto mat=beachmat::create_matrix<M>(input);

    // Checking the subset values.
    const size_t ncells=mat->get_ncol();
    const size_t ngenes=mat->get_nrow();
    auto SS=check_subset_vector(subset, ngenes);
    const size_t slen=SS.size();

    // Checking if we should transpose or not.
    const bool do_transpose=check_logical_scalar(transpose, "transpose specification");
    int out_nr=(do_transpose ? ncells : slen);
    int out_nc=(do_transpose ? slen : ncells);

    // Creating the output matrix.
    const bool sparsify=check_logical_scalar(as_sparse, "sparse specification");
    beachmat::output_param OPARAM;
    if (sparsify) {
        OPARAM=beachmat::output_param("dgCMatrix", "Matrix");
    }
    auto omat=beachmat::create_numeric_output(out_nr, out_nc, OPARAM);

    // Various other bits and pieces.
    std::vector<std::pair<typename M::type, size_t> > collected;
    collected.reserve(slen);
    std::vector<size_t> zeroes;
    zeroes.reserve(slen);

    beachmat::const_column<M> col_holder(mat.get(), false); // no sparse, need row-level indexing.
    Rcpp::NumericVector outgoing(slen);
    const double mean_rank=static_cast<double>(slen-1)/2;

    for (size_t c=0; c<ncells; ++c) {
        col_holder.fill(c);
        auto vals=col_holder.get_values();
        collected.clear();
        zeroes.clear();

        // Sorting all subsetted values (zeroes are handled separately for greater efficiency).
        auto sIt=SS.begin();
        for (size_t s=0; s<slen; ++s, ++sIt) {
            const typename M::type curval=*(vals + *sIt);
            if (isNA(curval)) { 
                throw std::runtime_error("missing values not supported in quickCluster");
            } else if (curval==0) {
                zeroes.push_back(s);
            } else {
                collected.push_back(std::make_pair(curval, s));
            }
        }
        std::sort(collected.begin(), collected.end());

        // Computing tied ranks for negative values, then positive values.
        double accumulated_rank=0;
        size_t cur_rank=0;
        double zero_rank=0;
        auto cIt=collected.begin();

        for (int positive=0; positive<2; ++positive) {
            while (cIt!=collected.end() && (positive==1 || cIt->first < 0)) {
                auto copy=cIt;
                ++copy;
                double accumulated_rank=cur_rank;
                ++cur_rank;

                while (copy!=collected.end() && copy->first==cIt->first) {
                    accumulated_rank+=cur_rank;
                    ++cur_rank;
                    ++copy;
                }

                double mean_rank=accumulated_rank/(copy-cIt);
                while (cIt!=copy) {
                    outgoing[cIt->second]=mean_rank;
                    ++cIt;
                }
            }

            // Special handling for zeroes.
            if (positive==0 && !zeroes.empty()) {
                zero_rank=static_cast<double>(zeroes.size()-1)/2 + cur_rank;
                for (auto z : zeroes) { outgoing[z]=zero_rank; }
                cur_rank+=zeroes.size();
            }
        }

        // Mean-adjusting (unless we want to leave it sparse) and converting to cosine values.
        double sum_squares=0;
        if (sparsify) {
            for (auto& o : outgoing) {
                double tmp=o-mean_rank;
                sum_squares+=tmp*tmp;
                o-=zero_rank;
            }
        } else {
            for (auto& o : outgoing) {
                o-=mean_rank;
                sum_squares+=o*o;
            }
        }

        if (sum_squares==0) {
            throw std::runtime_error("rank variances of zero detected for a cell");
        }
        sum_squares = std::sqrt(sum_squares)*2;
        for (auto& o : outgoing) {
            o/=sum_squares;
        }

        if (do_transpose) { 
            omat->set_row(c, outgoing.begin());
        } else {
            omat->set_col(c, outgoing.begin());
        }
    }

    return omat->yield();
}

SEXP get_scaled_ranks(SEXP exprs, SEXP subset, SEXP transpose, SEXP as_sparse) {
    BEGIN_RCPP
    int rtype=beachmat::find_sexp_type(exprs);
    if (rtype==INTSXP) { 
        return average_ranks_internal<beachmat::integer_matrix>(exprs, subset, transpose, as_sparse);
    } else {
        return average_ranks_internal<beachmat::numeric_matrix>(exprs, subset, transpose, as_sparse);
    }
    END_RCPP
}
