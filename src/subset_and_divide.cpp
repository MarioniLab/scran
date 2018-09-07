#include "scran.h"

/* A function to (a) subset by row, (b) subset by column, and (c) divide through by the library sizes. 
 * The output is equivalent to t(t(MAT[row_subset,col_subset])/lib.sizes) where lib.sizes itself is
 * computed as colSums(MAT[row_subset,col_subset]).
 */

template <class V, class M>
SEXP subset_and_divide_internal(M in, SEXP row_subset, SEXP col_subset, SEXP scaling) {
    // Checking subset vectors
    auto rsubout=check_subset_vector(row_subset, in->get_nrow());
    const size_t rslen=rsubout.size();
    auto csubout=check_subset_vector(col_subset, in->get_ncol());
    const size_t cslen=csubout.size();

    V incoming(in->get_nrow());
    size_t start_row=0, end_row=0;
    if (rslen) {
        // Cutting out extraction costs for unneeded start/end elements.
        start_row=*std::min_element(rsubout.begin(), rsubout.end());
        end_row=*std::max_element(rsubout.begin(), rsubout.end())+1;
    }

    const bool use_custom_scale=(scaling!=R_NilValue);
    Rcpp::NumericVector scale_values;
    if (use_custom_scale) {
        scale_values=Rcpp::NumericVector(scaling);
        if (scale_values.size()!=in->get_ncol()) {
            throw std::runtime_error("'length(scaling)' should be equal to 'ncol(x)'");
        }
    }

    // Setting up the output structures - only ever using in-memory matrices, as computeSumFactors() can't handle too many cells anyway.
    Rcpp::NumericVector outscale(cslen);
    Rcpp::NumericVector outgoing(rslen), averaged(rslen);

    beachmat::output_param oparam=(in->get_matrix_type()==beachmat::SPARSE ? beachmat::SPARSE_PARAM : beachmat::SIMPLE_PARAM);
    auto omat=beachmat::create_numeric_output(rslen, cslen, oparam);

    for (size_t cs=0; cs<cslen; ++cs) {
        const auto& curdex=csubout[cs];

        // Extracting the column, subsetting the rows.
        auto inIt=in->get_const_col(curdex, incoming.begin(), start_row, end_row);
        auto oIt=outgoing.begin();
        for (const auto& r : rsubout) {
            (*oIt)=*(inIt + r - start_row);
            ++oIt;
        }

        // Dividing by the library size. 
        double& curscale=outscale[cs];

        if (use_custom_scale) {
            curscale=scale_values[curdex];
            if (curscale < 0.00000001) {
                throw std::runtime_error("cells should have non-zero scaling factors");
            }
        } else {
            curscale=std::accumulate(outgoing.begin(), outgoing.end(), 0.0);
            if (curscale < 0.00000001) {
                throw std::runtime_error("cells should have non-zero library sizes");
            }
        }

        for (double& out : outgoing) {
            out/=curscale;
        }
        omat->set_col(cs, outgoing.begin());

        // Adding to the average.
        oIt=outgoing.begin();
        for (auto& a : averaged) {
            a+=(*oIt);
            ++oIt;
        }
    }

    // Scaling the sum to get the average and to adjust for the value of 'scaling'.
    const double adjustment=std::accumulate(outscale.begin(), outscale.end(), 0.0)/cslen/cslen;
    for (auto& a : averaged) {
        a*=adjustment;
    }

    return Rcpp::List::create(outscale, omat->yield(), averaged); 
}

SEXP subset_and_divide(SEXP matrix, SEXP row_subset, SEXP col_subset, SEXP scaling) {
    BEGIN_RCPP
    int rtype=beachmat::find_sexp_type(matrix);
    if (rtype==INTSXP) {
        auto input=beachmat::create_integer_matrix(matrix);
        return subset_and_divide_internal<Rcpp::IntegerVector>(input.get(), row_subset, col_subset, scaling);
    } else {
        auto input=beachmat::create_numeric_matrix(matrix);
        return subset_and_divide_internal<Rcpp::NumericVector>(input.get(), row_subset, col_subset, scaling);
    }
    END_RCPP
}


