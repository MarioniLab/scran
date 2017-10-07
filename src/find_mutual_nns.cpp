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
