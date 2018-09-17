#include "scran.h"

template <typename T, class V, class M>
SEXP overlap_exprs_internal(const M mat, const Rcpp::List& groups, SEXP subset, const T tol) {
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
    
    // Setting up the output matrices to hold the overlap proportions.
    Rcpp::List pout(ngroups);
    std::vector<std::vector<Rcpp::NumericMatrix::iterator> > pptrs(ngroups);
    for (size_t i=0; i<ngroups; ++i) {
        Rcpp::NumericMatrix tmpx(slen, i);
        auto pIt=tmpx.begin();
        pout[i]=tmpx;

        pptrs[i].resize(i);
        for (size_t j=0; j<i; ++j) {
            pptrs[i][j] = pIt;
            pIt += slen;
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

        // Running through each group and comparing to each other group.
        for (size_t i1=0; i1<ngroups; ++i1) {
            const auto& group1=by_group[i1];
            const size_t ncells1=group1.size();
            if (ncells1==0) { continue; }

            for (size_t i2=0; i2<i1; ++i2) {
                const auto& group2=by_group[i2];
                const size_t ncells2=group2.size();
                if (ncells2==0) { continue; }

                double& score=(*(pptrs[i1][i2]++)=0); // Referencing and bumping it up.
                size_t c2_left=0; 
                size_t c2_right=0;

                for (size_t c1=0; c1<ncells1; ++c1) {
                    const T& cur1=group1[c1];
                    const T left=cur1 - tol;
                    const T right=cur1 + tol;
                    while (c2_left < ncells2 && group2[c2_left] <= left) { ++c2_left; } // c2_left points to first element in range.
                    while (c2_right < ncells2 && group2[c2_right] < right) { ++c2_right; } // c2_right points to first element out of range.
                    score += static_cast<double>(c2_left) + static_cast<double>(c2_right - c2_left)*0.5;
                }
            }
        } 
    }

    return pout;
}

SEXP overlap_exprs(SEXP exprs, SEXP subset, SEXP bygroup, SEXP tolerance) {
    BEGIN_RCPP
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
    END_RCPP
}
