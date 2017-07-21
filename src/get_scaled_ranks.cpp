#include "scran.h"

template <typename T, class V, class M> 
SEXP average_ranks_internal(const M mat, SEXP intype, SEXP subset, SEXP transpose) { 
    /// Checking the subset values.
    const size_t ncells=mat->get_ncol();
    const size_t ngenes=mat->get_nrow();
    auto SS=check_subset_vector(subset, ngenes);
    const size_t slen=SS.size();

    // Checking if we should transpose or not.
    Rcpp::LogicalVector tr(transpose);
    if (tr.size()!=1) { 
        throw std::runtime_error("transpose specification should be a logical scalar");
    }
    const bool do_transpose=tr[0];
    int out_nr=(do_transpose ? ncells : slen);
    int out_nc=(do_transpose ? slen : ncells);

    // Various other bits and pieces.
    std::deque<std::pair<T, int> > collected(slen);
    const double mean_adj=double(slen-1)/2;
    V incoming(ngenes);
    Rcpp::NumericVector outgoing(slen);
    auto omat=beachmat::create_numeric_output(out_nr, out_nc, beachmat::output_param(intype, true, false));

    for (size_t c=0; c<ncells; ++c) {
        mat->get_col(c, incoming.begin());

        // Sorting all subsetted values.
        auto sIt=SS.begin();
        for (size_t s=0; s<slen; ++s, ++sIt) {
            const T& curval=incoming[*sIt];
            if (isNA(curval)) { 
                throw std::runtime_error("missing values not supported in quickCluster");
            }
            collected[s].first=curval;
            collected[s].second=s;
        }
        std::sort(collected.begin(), collected.end());

        // Need a bit more effort to deal with tied ranks.
        double accumulated_rank=0, sum_squares=0;
        int n_same_rank=0;
        for (size_t s=0; s<slen; ++s) {
            ++n_same_rank;
            accumulated_rank+=s;

            if (s==slen-1 || collected[s].first!=collected[s+1].first) {
                accumulated_rank /= n_same_rank;
                accumulated_rank -= mean_adj; // getting to a cosine distance.
                sum_squares += accumulated_rank * accumulated_rank * n_same_rank;
                
                int s_same=s;
                while (n_same_rank) {
                    outgoing[collected[s_same].second]=accumulated_rank;
                    --n_same_rank;
                    --s_same;
                }
                accumulated_rank=0;
            }
        }
        
        // Converting to cosine values.
        if (sum_squares==0) {
            throw std::runtime_error("rank variances of zero detected for a cell");
        }
        sum_squares = std::sqrt(sum_squares)*2;
        for (auto oIt=outgoing.begin(); oIt!=outgoing.end(); ++oIt) {
            (*oIt)/=sum_squares;
        }

        if (do_transpose) { 
            omat->set_row(c, outgoing.begin());
        } else {
            omat->set_col(c, outgoing.begin());
        }
    }

    return omat->yield();
}

SEXP get_scaled_ranks(SEXP exprs, SEXP subset, SEXP transpose) {
    BEGIN_RCPP
    int rtype=beachmat::find_sexp_type(exprs);
    if (rtype==INTSXP) { 
        auto mat=beachmat::create_integer_matrix(exprs);
        return average_ranks_internal<int, Rcpp::IntegerVector>(mat.get(), exprs, subset, transpose);
    } else {
        auto mat=beachmat::create_numeric_matrix(exprs);
        return average_ranks_internal<double, Rcpp::NumericVector>(mat.get(), exprs, subset, transpose);
    }
    END_RCPP
}
