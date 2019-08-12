#include "Rcpp.h"

#include "beachmat/integer_matrix.h"
#include "beachmat/numeric_matrix.h"
#include "utils.h"

#include <stdexcept>
#include <deque>
#include <vector>
#include <algorithm>

template<typename T, class V>
class wilcoxer {
public:    
    wilcoxer(Rcpp::List groups, int ncells) {
        const size_t ngroups=groups.size();
        for (size_t g=0; g<ngroups; ++g) {
            Rcpp::IntegerVector curgroup(groups[g]);
            sources.push_back(std::vector<int>(curgroup.begin(), curgroup.end()));

            std::vector<int>& last_added=sources.back();
            for (auto& current : last_added) {
                if (current<0 || current>=ncells) { 
                    throw std::runtime_error("indices in 'groups' out of range");
                }
            }

            by_group.push_back(std::vector<T>(curgroup.size()));
        }
        return;
    }

    void initialize(const V& vec) {
        // Sorting expression values within each group.
        for (size_t i=0; i<by_group.size(); ++i) {
            const auto& cur_source=sources[i];
            auto& cur_group=by_group[i];

            auto cgIt=cur_group.begin();
            for (auto csIt=cur_source.begin(); csIt!=cur_source.end(); ++csIt, ++cgIt) { 
                (*cgIt)=vec[*csIt];
            }
            std::sort(cur_group.begin(), cur_group.end());
        }
        return;
    }

    std::pair<double, int> contrast_groups(int left, int right, T shift) const {
        int c1=0, c2=0;
        const auto& group1=by_group[left];
        const int ncells1=group1.size();
        const auto& group2=by_group[right];
        const int ncells2=group2.size();

        std::pair<double, double> output;
        double& score=output.first;
        double& tieval=output.second;

        while (1) {
            // Each iteration of this loop should represent one set of tied ranks
            // (after adjusting 'group1' for the shift).
            const bool ok1=c1 < ncells1;
            const bool ok2=c2 < ncells2;
            T curval;

            if (!ok1 && !ok2) {
                break;
            } else if (ok1 && ok2) {
                curval=std::min(group1[c1] - shift, group2[c2]);
            } else if (ok1) {
                curval=group1[c1] - shift;
            } else {
                curval=group2[c2];
            }

            // Make each index point to first element outside of the range.
            const T right=curval;
            int ties1=0;
            if (ok1) { 
                const int c1_old=c1;
                while (c1 < ncells1 && group1[c1] - shift <= right) { 
                    ++c1;
                }
                ties1=c1-c1_old;
            }

            int ties2=0;
            const int c2_old=c2;
            if (ok2) {
                while (c2 < ncells2 && group2[c2] <= right) { 
                    ++c2;
                }
                ties2=c2 - c2_old;
            }

            score += (static_cast<double>(c2_old) + static_cast<double>(ties2)*0.5) * ties1;
                
            // To recapitulate the normal approximation in stats::wilcox.test(). 
            const double nties=ties1 + ties2;
            tieval += nties * (nties * nties - 1);
        }

        return output;
    }

    bool empty_group(int i) const {
        return (sources[i].size()==0);
    }
private:
    std::deque<std::vector<int> > sources;
    std::deque<std::vector<T> > by_group;
};

template <typename T, class V, class M>
Rcpp::List overlap_exprs_internal(const M mat, const Rcpp::List& groups, Rcpp::IntegerVector subset, const T lfc) {
    /// Checking the subset values.
    auto SS=check_subset_vector(subset, mat->get_nrow());
    const size_t slen=SS.size();
    const size_t ncells=mat->get_ncol();
   
    // Constructing groups. 
    const size_t ngroups=groups.size();
    wilcoxer<T, V> wilcox_calc(groups, ncells);
    
    // Setting up the output matrices to hold the overlap proportions, number of ties.
    Rcpp::List pout(ngroups), tout(ngroups);
    std::vector<std::vector<Rcpp::NumericMatrix::iterator> > pptrs(ngroups), tptrs(ngroups);

    for (size_t i=0; i<ngroups; ++i) {
        size_t ncols=(lfc==0 ? i : ngroups);
        Rcpp::NumericMatrix tmpP(slen, ncols), tmpT(slen, ncols);
        auto pIt=tmpP.begin(), tIt=tmpT.begin();
        pout[i]=tmpP;
        tout[i]=tmpT;

        pptrs[i].reserve(ncols);
        tptrs[i].reserve(ncols);
        for (size_t j=0; j<ncols; ++j, pIt+=slen, tIt+=slen) {
            pptrs[i].push_back(pIt);
            tptrs[i].push_back(tIt);
        }
    }

    // Running through all genes and computing pairwise overlaps. 
    V tmp(ncells);
    for (auto sIt=SS.begin(); sIt!=SS.end(); ++sIt) { 
        mat->get_row(*sIt, tmp.begin());
        wilcox_calc.initialize(tmp);

        for (size_t i1=0; i1<ngroups; ++i1) {
            if (wilcox_calc.empty_group(i1)) { continue; }
            for (size_t i2=0; i2<i1; ++i2) {
                if (wilcox_calc.empty_group(i2)) { continue; }

                auto left=wilcox_calc.contrast_groups(i1, i2, lfc);
                *(pptrs[i1][i2]++)=left.first;
                *(tptrs[i1][i2]++)=left.second;

                if (lfc!=0) {
                    auto right=wilcox_calc.contrast_groups(i1, i2, -lfc);
                    *(pptrs[i2][i1]++)=right.first;
                    *(tptrs[i2][i1]++)=right.second;
                }
            }
        } 
    }

    return Rcpp::List::create(pout, tout);
}

// [[Rcpp::export(rng=false)]]
Rcpp::List overlap_exprs(Rcpp::RObject exprs, Rcpp::IntegerVector subset, Rcpp::List bygroup, Rcpp::RObject lfc) {
    int rtype=beachmat::find_sexp_type(exprs);
    if (rtype==INTSXP) {
        auto mat=beachmat::create_integer_matrix(exprs);
        const int shift=check_integer_scalar(lfc, "lfc");
        return overlap_exprs_internal<int, Rcpp::IntegerVector>(mat.get(), bygroup, subset, shift);
    } else {
        auto mat=beachmat::create_numeric_matrix(exprs);
        const double shift=check_numeric_scalar(lfc, "lfc");
        return overlap_exprs_internal<double, Rcpp::NumericVector>(mat.get(), bygroup, subset, shift);
    }
}
