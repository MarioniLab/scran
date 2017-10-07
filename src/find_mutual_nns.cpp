#include "scran.h"

SEXP find_mutual_nns (SEXP left, SEXP right) {
    BEGIN_RCPP
    const Rcpp::IntegerMatrix L(left), R(right);
    const int nnR=R.ncol();
 
    std::vector<int> sortedR(R.size());
    std::deque<int> mutualL, mutualR;

    // Sorting the elements on the right.
    auto sIt=sortedR.begin();
    for (int r=0; r<R.nrow(); ++r) {
        const auto currow=R.row(r);
        std::copy(currow.begin(), currow.end(), sIt);
        std::sort(sIt, sIt+nnR);
        sIt+=nnR;
    }

    // Running through the elements on the left, and doing a binary search.
    for (int l=0; l<L.nrow(); ++l) {
        const auto currow=L.row(l);
        const int tocheck=l+1;

        for (const auto& curval : currow) {
            auto startIt=sortedR.begin() + nnR*(curval-1); // 1-indexed.
            auto endIt=startIt+nnR;
            auto closest=std::lower_bound(startIt, endIt, tocheck);

            if (closest!=endIt && *closest==tocheck) { 
                mutualL.push_back(tocheck);
                mutualR.push_back(curval);
            }
        }
    }

    return Rcpp::List::create(Rcpp::IntegerVector(mutualL.begin(), mutualL.end()),
                              Rcpp::IntegerVector(mutualR.begin(), mutualR.end()));
    END_RCPP
}

/* Performs the cosine normalization in a fairly efficient manner. */

template<class M>
SEXP cosine_norm_internal (M mat, SEXP original) {
    const size_t& nrow=mat->get_nrow();
    const size_t& ncol=mat->get_ncol();
    auto output=beachmat::create_numeric_output(nrow, ncol, beachmat::output_param(original));
    
    Rcpp::NumericVector incoming(nrow);
    for (size_t c=0; c<ncol; ++c) {
        mat->get_col(c, incoming.begin());

        double total=0;
        for (const auto& val : incoming) { 
            total+=val*val;
        }
        total=std::sqrt(total);
        total=std::max(total, 0.00000001); // avoid division by zero.

        for (auto& val : incoming) { 
            val/=total;
        }
        output->set_col(c, incoming.begin());
    }

    return output->yield();
}

SEXP cosine_norm(SEXP incoming) {
    BEGIN_RCPP
    int rtype=beachmat::find_sexp_type(incoming);
    if (rtype==INTSXP) {
        auto input=beachmat::create_integer_matrix(incoming);
        return cosine_norm_internal(input.get(), incoming);
    } else {
        auto input=beachmat::create_numeric_matrix(incoming);
        return cosine_norm_internal(input.get(), incoming);
    }
    END_RCPP
}

