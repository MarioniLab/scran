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
            --current; // get to 0-based indexing.
            if (current<0 || current>=int(ncells)) { 
                throw std::runtime_error("indices in 'groups' out of range");
            }
        }

        by_group.push_back(std::vector<T>(curgroup.size()));;
    }
    
    // Setting up the output matrices to hold the overlap proportions.
    Rcpp::List pout(ngroups);
    std::deque<Rcpp::NumericVector::iterator> pptrs(ngroups*ngroups);
    int counter=0;
        
    for (size_t i=0; i<ngroups; ++i) {
        Rcpp::NumericMatrix tmpx(slen, ngroups-1);
        auto pIt=tmpx.begin();
        pout[i]=tmpx;

        for (size_t j=0; j<ngroups; ++j, ++counter) {
            if (i!=j) { 
                pptrs[counter] = pIt;
                pIt += slen;
            }
        }
    }

    // Calclulating the number of cells for each pair of groups.
    Rcpp::List nout(ngroups);
    std::deque<double> ncellpairs(ngroups*ngroups);
    counter=0;

    for (size_t i=0; i<ngroups; ++i) {
        const size_t isize=sources[i].size();
        Rcpp::NumericVector tmpx(ngroups-1);
        auto nIt=tmpx.begin();
        nout[i]=tmpx;

        for (int j=0; j<ngroups; ++j, ++counter) {
            const size_t jsize=sources[j].size();
            if (i!=j) { 
                if (isize && jsize) { 
                    (*nIt) = (ncellpairs[counter] = jsize + isize);
                }
                ++nIt;
            }
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

                int counter=i1*ngroups + i2;
                double& score=(*(pptrs[counter]++)=0); // Referencing and bumping it up.
                size_t c2_left=0; 
                size_t c2_right=0;

                for (size_t c1=0; c1<ncells1; ++c1) {
                    const T& cur1=group1[c1];
                    const T left=cur1 - tol;
                    const T right=cur1 + tol;
                    while (c2_left < ncells2 && group2[c2_left] <= left) { ++c2_left; } // c2_left points to first element in range.
                    while (c2_right < ncells2 && group2[c2_right] < right) { ++c2_right; } // c2_right points to first element out of range.
                    score += double(c2_left) + double(c2_right - c2_left)*0.5;
                }
                score/=double(ncells1)*double(ncells2);

                // Accounting for the total number of cells.
                const double& total_cells=ncellpairs[counter];
                counter=i2*ngroups + i1; // Adding the symmetric value.
                *(pptrs[counter]++)=(1-score) * total_cells;
                score *= total_cells;
            }
        } 
    }

    Rcpp::List output(2);
    output[0]=pout;
    output[1]=nout;
    return output;
}

SEXP overlap_exprs(SEXP exprs, SEXP subset, SEXP bygroup, SEXP tolerance) {
    BEGIN_RCPP
    int rtype=beachmat::find_sexp_type(exprs);
    if (rtype==INTSXP) {
        auto mat=beachmat::create_integer_matrix(exprs);
        Rcpp::IntegerVector tol(tolerance);
        if (tol.size()!=1) { 
            throw std::runtime_error("tolerance should be an integer scalar");
        }
        return overlap_exprs_internal<int, Rcpp::IntegerVector>(mat.get(), bygroup, subset, tol[0]);
    } else {
        auto mat=beachmat::create_numeric_matrix(exprs);
        Rcpp::NumericVector tol(tolerance);
        if (tol.size()!=1) { 
            throw std::runtime_error("tolerance should be a double-precision scalar");
        }
        return overlap_exprs_internal<double, Rcpp::NumericVector>(mat.get(), bygroup, subset, tol[0]);
    }
    END_RCPP
}
