#include "scran.h"

#include "beachmat/numeric_matrix.h"
#include "beachmat/integer_matrix.h"
#include "beachmat/utils/const_column.h"
#include "utils.h"

#include <stdexcept>
#include <algorithm>

/* A function to (a) subset by row, (b) subset by column, and (c) divide through by the library sizes. 
 * The output is equivalent to t(t(MAT[row_subset,col_subset])/lib.sizes) where lib.sizes itself is
 * computed as colSums(MAT[row_subset,col_subset]).
 */

template <class M>
SEXP subset_and_divide_internal(SEXP incoming, SEXP row_subset, SEXP col_subset, SEXP scaling) {
    auto in=beachmat::create_matrix<M>(incoming);
    beachmat::const_column<M> col_holder(in.get(), false); // need to index by rows, so turn off sparsity.

    // Checking subset vectors
    auto rsubout=check_subset_vector(row_subset, in->get_nrow());
    const size_t rslen=rsubout.size();
    auto csubout=check_subset_vector(col_subset, in->get_ncol());
    const size_t cslen=csubout.size();

    // Cutting out extraction costs for unneeded start/end elements.
    size_t start_row=0, end_row=0;
    if (rslen) {
        start_row=*std::min_element(rsubout.begin(), rsubout.end());
        end_row=*std::max_element(rsubout.begin(), rsubout.end())+1;
    }

    const bool use_custom_scale=(scaling!=R_NilValue);
    Rcpp::NumericVector scale_values;
    if (use_custom_scale) {
        scale_values=Rcpp::NumericVector(scaling);
        if (static_cast<size_t>(scale_values.size())!=in->get_ncol()) {
            throw std::runtime_error("'length(scaling)' should be equal to 'ncol(x)'");
        }
    }

    // Setting up the output structures - only ever using in-memory matrices, as computeSumFactors() can't handle too many cells anyway.
    Rcpp::NumericVector outscale(cslen);
    Rcpp::NumericVector outgoing(rslen), averaged(rslen);

    beachmat::output_param OPARAM;
    if (in->get_class()=="dgCMatrix" && in->get_package()=="Matrix") {
        OPARAM=beachmat::output_param("dgCMatrix", "Matrix");
    }
    auto omat=beachmat::create_numeric_output(rslen, cslen, OPARAM);

    for (size_t cs=0; cs<cslen; ++cs) {
        const auto& curdex=csubout[cs];
        col_holder.fill(curdex, start_row, end_row);

        // Extracting the column, subsetting the rows.
        auto val=col_holder.get_values();
        auto oIt=outgoing.begin();
        for (const auto& r : rsubout) {
            (*oIt)=*(val + r - start_row);
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
        return subset_and_divide_internal<beachmat::integer_matrix>(matrix, row_subset, col_subset, scaling);
    } else {
        return subset_and_divide_internal<beachmat::numeric_matrix>(matrix, row_subset, col_subset, scaling);
    }
    END_RCPP
}


