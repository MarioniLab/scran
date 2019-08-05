#include "Rcpp.h"

#include "beachmat/integer_matrix.h"
#include "beachmat/numeric_matrix.h"
#include "utils.h"

#include <stdexcept>
#include <deque>
#include <vector>
#include <algorithm>

template <typename T, class V, class M>
Rcpp::List overlap_exprs_internal(const M mat, const Rcpp::List& groups, Rcpp::IntegerVector subset, const T tol) {
    /// Checking the subset values.
    auto SS=check_subset_vector(subset, mat->get_nrow());
    const size_t slen=SS.size();
    const size_t& ncells=mat->get_ncol();
   
    // Constructing groups. 
    const size_t ngroups=groups.size();
    std::deque<std::vector<int> > sources;
    std::deque<std::vector<T> > by_group;
    for (size_t g=0; g<ngroups; ++g) {
        Rcpp::RObject gdata=groups[g];
        if (gdata.sexp_type()!=INTSXP) {  
            throw std::runtime_error("'groups' should contain integer vectors"); 
        }
        Rcpp::IntegerVector curgroup(gdata);

        sources.push_back(std::vector<int>(curgroup.begin(), curgroup.end()));
        std::vector<int>& last_added=sources.back();
        for (auto laIt=last_added.begin(); laIt!=last_added.end(); ++laIt) {
            int& current=(*laIt);
            if (current<0 || current>=static_cast<int>(ncells)) { 
                throw std::runtime_error("indices in 'groups' out of range");
            }
        }

        by_group.push_back(std::vector<T>(curgroup.size()));;
    }
    
    // Setting up the output matrices to hold the overlap proportions, number of ties.
    Rcpp::List pout(ngroups), tout(ngroups);
    std::vector<std::vector<Rcpp::NumericMatrix::iterator> > pptrs(ngroups), tptrs(ngroups);

    for (size_t i=0; i<ngroups; ++i) {
        Rcpp::NumericMatrix tmpP(slen, i), tmpT(slen, i);
        auto pIt=tmpP.begin(), tIt=tmpT.begin();
        pout[i]=tmpP;
        tout[i]=tmpT;

        pptrs[i].resize(i);
        tptrs[i].resize(i);
        for (size_t j=0; j<i; ++j) {
            pptrs[i][j] = pIt;
            pIt += slen;
            tptrs[i][j] = tIt;
            tIt += slen;
        }
    }

    // Running through all genes and computing pairwise overlaps. 
    V tmp(ncells);
    for (auto sIt=SS.begin(); sIt!=SS.end(); ++sIt) { 
        mat->get_row(*sIt, tmp.begin());

        // Sorting expression values within each group.
        for (size_t i=0; i<ngroups; ++i) {
            const auto& cur_source=sources[i];
            auto& cur_group=by_group[i];

            auto cgIt=cur_group.begin();
            for (auto csIt=cur_source.begin(); csIt!=cur_source.end(); ++csIt, ++cgIt) { 
                (*cgIt)=tmp[*csIt];
            }
            std::sort(cur_group.begin(), cur_group.end());
        }

        for (size_t i1=0; i1<ngroups; ++i1) {
            const auto& group1=by_group[i1];
            const size_t ncells1=group1.size();
            if (ncells1==0) { continue; }

            // Comparing to every other group.
            for (size_t i2=0; i2<i1; ++i2) {
                const auto& group2=by_group[i2];
                const size_t ncells2=group2.size();
                if (ncells2==0) { continue; }

                double& score=*(pptrs[i1][i2]++); 
                double& tieval=*(tptrs[i1][i2]++);
                size_t c1=0, c2=0;

                while (1) {
                    // Each iteration of this loop should represent one set of tied ranks.
                    const bool ok1=c1 < ncells1;
                    const bool ok2=c2 < ncells2;
                    T curval;

                    if (!ok1 && !ok2) {
                        break;
                    } else if (ok1 && ok2) {
                        if (group1[c1] < group2[c2]) {
                            curval=group1[c1];
                        } else {
                            curval=group2[c2];
                        }
                    } else if (ok1) {
                        curval=group1[c1];
                    } else {
                        curval=group2[c2];
                    }

                    // Make each index point to first element outside of the range.
                    const T right=curval + tol;
                    size_t ties1=0;
                    if (ok1) { 
                        const size_t c1_old=c1;
                        while (c1 < ncells1 && group1[c1] <= right) { 
                            ++c1;
                        }
                        ties1=c1-c1_old;
                    }

                    size_t ties2=0;
                    const size_t c2_old=c2;
                    if (ok2) {
                        while (c2 < ncells2 && group2[c2] <= right) { 
                            ++c2;
                        }
                        ties2=c2 - c2_old;
                    }

                    score += (static_cast<double>(c2_old) + static_cast<double>(ties2)*0.5) * ties1;
                    if (ties1 + ties2 > 1) {
                        // To recapitulate the normal approximation in stats::wilcox.test(). 
                        const double nties=ties1 + ties2;
                        tieval += nties * (nties * nties - 1);
                    }
                }
            }
        } 
    }

    return Rcpp::List::create(pout, tout);
}

// [[Rcpp::export(rng=false)]]
Rcpp::List overlap_exprs(Rcpp::RObject exprs, Rcpp::IntegerVector subset, Rcpp::List bygroup, Rcpp::RObject tolerance) {
    int rtype=beachmat::find_sexp_type(exprs);
    if (rtype==INTSXP) {
        auto mat=beachmat::create_integer_matrix(exprs);
        const int tol=check_integer_scalar(tolerance, "tolerance");
        return overlap_exprs_internal<int, Rcpp::IntegerVector>(mat.get(), bygroup, subset, tol);
    } else {
        auto mat=beachmat::create_numeric_matrix(exprs);
        const double tol=check_numeric_scalar(tolerance, "tolerance");
        return overlap_exprs_internal<double, Rcpp::NumericVector>(mat.get(), bygroup, subset, tol);
    }
}
