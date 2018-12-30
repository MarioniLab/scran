#include "scran.h"

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
        if (s < 1 || s > ncells) { throw std::runtime_error("each element of sizes should be within [1, number of cells]"); }
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
        if (o < 0 || o > ncells) { 
            throw std::runtime_error("elements of ordering vector are out of range");
        }
    }

    // Setting up the storage space and the vector of iterators for each cell's profile.
    Rcpp::NumericVector all_collected(last_size*ngenes);
    std::deque<Rcpp::NumericVector::const_iterator> collected;
    auto acIt=all_collected.begin();
    collected.push_back(acIt); // unfilled first vector, which gets dropped and refilled in the first iteration anyway.
    acIt+=ngenes;

    auto orIt_tail=order.begin();
    Rcpp::NumericVector tmp(emat->get_nrow());
    for (int s=1; s<last_size; ++s, ++orIt_tail, acIt+=ngenes) {
        auto colIt=emat->get_const_col(*orIt_tail, acIt);
        collected.push_back(colIt); 
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
    for (size_t win=0; win<ncells; ++win, ++orIt) {
        std::fill(combined.begin(), combined.end(), 0);
        
        // Dropping the column that we've moved past, adding the next column.
        if (acIt==all_collected.end()) { 
            acIt=all_collected.begin();
        }
        collected.pop_front();
        auto curIt=emat->get_const_col(*orIt_tail, acIt);
        collected.push_back(curIt);
        ++orIt_tail;
        acIt+=ngenes;

        int index=0;
        int rownum=win; // Setting the row so that all pools with the same SIZE form consecutive equations.
        for (auto psIt=pool_sizes.begin(); psIt!=pool_sizes.end(); ++psIt, rownum+=ncells) { 
            const int& SIZE=(*psIt);
            std::fill(rowIt, rowIt+SIZE, rownum);
            rowIt+=SIZE;
            std::copy(orIt, orIt+SIZE, colIt);
            colIt+=SIZE;

            for (; index<SIZE; ++index) {
                auto ceIt=collected[index];
                for (auto cIt=combined.begin(); cIt!=combined.end(); ++cIt, ++ceIt) {
                    (*cIt)+=(*ceIt);
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
