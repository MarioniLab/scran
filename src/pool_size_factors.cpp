#include "scran.h"

#include "beachmat/numeric_matrix.h"
#include "beachmat/utils/const_column.h"

#include <stdexcept>
#include <algorithm>
#include <vector>

/*** A function to estimate the pooled size factors and construct the linear equations. ***/

SEXP pool_size_factors (SEXP exprs, SEXP ref, SEXP ordering, SEXP poolsizes) {
    BEGIN_RCPP
    auto emat=beachmat::create_numeric_matrix(exprs);
    const size_t ngenes=emat->get_nrow();
    const size_t ncells=emat->get_ncol();
    if (ncells==0) { 
        throw std::runtime_error("at least one cell required for normalization"); 
    }
   
    // Checking the input sizes.
    Rcpp::IntegerVector pool_sizes(poolsizes);
    const size_t nsizes=pool_sizes.size();
    if (nsizes==0) {
        return Rcpp::List::create(Rcpp::IntegerVector(0), Rcpp::IntegerVector(0), Rcpp::NumericVector(0));
    }

    int last_size=-1, total_size=0;
    for (auto s : pool_sizes) { 
        if (s < 1 || static_cast<size_t>(s) > ncells) { throw std::runtime_error("each element of sizes should be within [1, number of cells]"); }
        if (s < last_size) { throw std::runtime_error("sizes should be sorted"); }
        total_size+=s;
        last_size=s;
    }

    // Checking pseudo cell.
    Rcpp::NumericVector pseudo_cell(ref);
    if (ngenes!=pseudo_cell.size()) { throw std::runtime_error("length of pseudo-cell vector is not the same as the number of cells"); }

    // Checking ordering.
    Rcpp::IntegerVector order(ordering);
    if (order.size() < ncells*2-1)  { throw std::runtime_error("ordering vector is too short for number of cells"); }
    for (auto o : order) { 
        if (o < 0 || static_cast<size_t>(o) >= ncells) { 
            throw std::runtime_error("elements of ordering vector are out of range");
        }
    }

    /* Setting up the storage space.
     * The first vector is unfilled as it gets dropped and refilled in the first iteration anyway.
     * We also trigger generation of indices so that each element doesn't generate its own indices.
     */
    typedef beachmat::numeric_matrix M;
    std::deque<beachmat::const_column<M> > collected(1, beachmat::const_column<M>(emat.get()));
    collected.back().get_indices(); 

    auto orIt_tail=order.begin();
    for (int s=1; s<last_size; ++s, ++orIt_tail) {
        collected.push_back(collected.front());
        collected.back().fill(*orIt_tail);
    }

    // Setting up the output vectors and other bits and pieces.
    Rcpp::IntegerVector row_num(total_size*ncells), col_num(total_size*ncells);
    Rcpp::NumericVector pool_factor(nsizes*ncells);

    std::vector<double> combined(ngenes), ratios(ngenes);
    auto rowIt=row_num.begin(), colIt=col_num.begin();
    auto orIt=order.begin();

    if (ngenes==0) {
        throw std::runtime_error("insufficient features for median calculations");
    }
    const bool is_even=bool(ngenes%2==0);
    const int halfway=int(ngenes/2);

    // Running through the sliding windows.
    for (size_t win=0; win<ncells; ++win, ++orIt, ++orIt_tail) {
        std::fill(combined.begin(), combined.end(), 0);

        /* Rotating; effectively moves the first element of 'collected' to the end.
         * The is the same as shifting the column that we've moved past to the end,
         * and then overwriting it with the next column.
         */
        std::rotate(collected.begin(), collected.begin()+1, collected.end());
        collected.back().fill(*orIt_tail);

        int index=0;
        int rownum=win; // Setting the row so that all pools with the same SIZE form consecutive equations.
        for (auto psIt=pool_sizes.begin(); psIt!=pool_sizes.end(); ++psIt, rownum+=ncells) { 
            const int& SIZE=(*psIt);
            std::fill(rowIt, rowIt+SIZE, rownum);
            rowIt+=SIZE;
            std::copy(orIt, orIt+SIZE, colIt);
            colIt+=SIZE;

            for (; index<SIZE; ++index) {
                auto& current=collected[index];
                auto val=current.get_values();

                if (current.is_sparse()) {
                    auto n=current.get_n();
                    auto idx=current.get_indices();
                    for (auto i=0; i<n; ++i, ++val, ++idx) {
                        combined[*idx]+=*val;
                    }
                } else { 
                    for (auto cIt=combined.begin(); cIt!=combined.end(); ++cIt, ++val) {
                        (*cIt)+=(*val);
                    }
                }
            }
           
            // Computing the ratio against the reference.
            auto rIt=ratios.begin(), cIt=combined.begin();
            for (auto pcIt=pseudo_cell.begin(); pcIt!=pseudo_cell.end(); ++pcIt, ++rIt, ++cIt) {
                (*rIt)=(*cIt)/(*pcIt);
            }

            // Computing the median (faster than partial sort).
            std::nth_element(ratios.begin(), ratios.begin()+halfway, ratios.end());
            if (is_even) {
                double medtmp=ratios[halfway];
                std::nth_element(ratios.begin(), ratios.begin()+halfway-1, ratios.end());
                pool_factor[rownum]=(medtmp+ratios[halfway-1])/2;
            } else {
                pool_factor[rownum]=ratios[halfway];
            }       
        }

/*
        std::partial_sort(combined, combined+halfway+1, combined+ngenes);
        if (is_even) {
            ofptr[cell]=(combined[halfway]+combined[halfway-1])/2;
        } else {
            ofptr[cell]=combined[halfway];
        } 
*/
    }    

    return Rcpp::List::create(row_num, col_num, pool_factor);
    END_RCPP
}
