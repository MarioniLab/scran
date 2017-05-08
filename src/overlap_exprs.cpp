#include "scran.h"

template <typename T>
SEXP overlap_exprs_internal(const T* ptr, const matrix_info& MAT, SEXP groups, SEXP subset, const T tol) {
    /// Checking the subset values.
    subset_values SS=check_subset_vector(subset, MAT.nrow);
    const int slen=SS.first;
    const int* sptr=SS.second;
    const size_t& ngenes=MAT.nrow;
    const size_t& ncells=MAT.ncol;
   
    // Constructing groups. 
    if (!isNewList(groups)) { 
        throw std::runtime_error("'groups' should be a list");
    }
    const int ngroups=LENGTH(groups);
    std::deque<const int*> sources(ngroups);
    std::deque<int> groupsize(ngroups);
    std::deque<std::deque<T> > by_group(ngroups);
    for (int i=0; i<ngroups; ++i) {
        SEXP gdata=VECTOR_ELT(groups, i);
        if (!isInteger(gdata)) {  throw std::runtime_error("'groups' should contain integer vectors"); }
        const int* cursource=(sources[i]=INTEGER(gdata));
        const int& cursize=(groupsize[i]=LENGTH(gdata));
        for (int j=0; j<cursize; ++j) { 
            if (cursource[j]<1 || cursource[j]>ncells) { 
                throw std::runtime_error("indices in 'groups' out of range");
            }
        }
        by_group[i].resize(cursize);
    }
    
    // Setting up the output structures.
    SEXP output=PROTECT(allocVector(VECSXP, 2));
    try {
        SET_VECTOR_ELT(output, 0, allocVector(VECSXP, ngroups));
        SEXP pout=VECTOR_ELT(output, 0);
        std::deque<double*> pptrs(ngroups*ngroups, NULL);
        int counter=0;
        
        for (int i=0; i<ngroups; ++i) {
            SET_VECTOR_ELT(pout, i, allocMatrix(REALSXP, slen, ngroups-1));
            double * pptr=REAL(VECTOR_ELT(pout, i));
            for (int j=0; j<ngroups; ++j, ++counter) {
                if (i!=j) { 
                    pptrs[counter] = pptr;
                    std::fill(pptr, pptr+slen, 0);
                    pptr += slen;
                }
            }
        }

        // Calclulating the number of cells for each pair of groups.
        SET_VECTOR_ELT(output, 1, allocVector(VECSXP, ngroups));
        SEXP nout=VECTOR_ELT(output, 1);
        std::deque<double*> nptrs(ngroups*ngroups, NULL);
        counter=0;

        double isize, jsize;
        for (int i=0; i<ngroups; ++i) {
            isize=groupsize[i];
            SET_VECTOR_ELT(nout, i, allocVector(REALSXP, ngroups-1));
            double * nptr=REAL(VECTOR_ELT(nout, i));

            for (int j=0; j<ngroups; ++j, ++counter) {
                jsize=groupsize[j];
                if (i!=j) { 
                    nptrs[counter] = nptr;
                    (*nptr)= (isize && jsize) ? jsize + isize : 0;
                    ++nptr;
                }
            }
        }

        // Setting up temporaries.
        int i1, i2;
        const T* cur_ptr;
        T left, right;
        int c1, c2_left, c2_right;
        double total_cells;

        // Running through all genes and computing pairwise overlaps.
        for (int s=0; s<slen; ++s) {
            const int& g=sptr[s];

            // Sorting expression values within each group.
            cur_ptr = ptr + g - ngenes;
            for (i1=0; i1<ngroups; ++i1) {
                const int* cur_source=sources[i1];
                std::deque<T>& cur_group=by_group[i1];
                const int& cur_size=groupsize[i1];
                for (i2=0; i2<cur_size; ++i2) {
                    cur_group[i2]=cur_ptr[ngenes * cur_source[i2]];
                }
                std::sort(cur_group.begin(), cur_group.end());
            }

            // Running through each group and comparing to each other group.
            for (i1=0; i1<ngroups; ++i1) {
                const std::deque<T>& group1=by_group[i1];
                const int& ncells1=groupsize[i1];
                if (ncells1==0) { continue; }

                for (i2=0; i2<i1; ++i2) {
                    const std::deque<T>& group2=by_group[i2];
                    const int& ncells2=groupsize[i2];
                    if (ncells2==0) { continue; }

                    counter=i1*ngroups + i2;
                    double& score=(pptrs[counter][s]=0);
                    c2_left=0; 
                    c2_right=0;

                    for (c1=0; c1<ncells1; ++c1) {
                        const T& cur1=group1[c1];
                        left=cur1 - tol;
                        right=cur1 + tol;
                        while (c2_left < ncells2 && group2[c2_left] <= left) { ++c2_left; } // c2_left points to first element in range.
                        while (c2_right < ncells2 && group2[c2_right] < right) { ++c2_right; } // c2_right points to first element out of range.
                        score += double(c2_left) + double(c2_right - c2_left)*0.5;
                    }
                    score/=double(ncells1)*double(ncells2);

                    // Accounting for the total number of cells.
                    const double& total_cells=nptrs[counter][0];
                    counter=i2*ngroups + i1; // Adding the symmetric value.
                    pptrs[counter][s]=(1-score) * total_cells;
                    score *= total_cells;
                }
            } 
        }
    } catch (std::exception& e) {
        UNPROTECT(1);
        throw;
    }

    UNPROTECT(1);
    return output;
}

SEXP overlap_exprs(SEXP exprs, SEXP subset, SEXP bygroup, SEXP tolerance) try {
    matrix_info MAT=check_matrix(exprs);
    if (!isReal(tolerance) || LENGTH(tolerance)!=1) { 
        throw std::runtime_error("tolerance should be a double-precision scalar");
    }
    if (MAT.is_integer){ 
        return overlap_exprs_internal<int>(MAT.iptr, MAT, bygroup, subset, asInteger(tolerance));
    } else {
        return overlap_exprs_internal<double>(MAT.dptr, MAT, bygroup, subset, asReal(tolerance));
    }
} catch (std::exception& e) {
    return mkString(e.what());
}
