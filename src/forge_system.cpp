#include "scran.h"

/* A function to (a) subset by row, (b) subset by column, and (c) divide through by the library sizes. 
 * The output is equivalent to t(t(MAT[row_subset,col_subset])/lib.sizes) where lib.sizes itself is
 * computed as colSums(MAT[row_subset,col_subset]).
 */

template <class V, class M>
SEXP subset_and_divide_internal(const M in, SEXP inmat, SEXP row_subset, SEXP col_subset) {
    // Checking subset vectors
    auto rsubout=check_subset_vector(row_subset, in->get_nrow());
    const size_t rslen=rsubout.size();
    auto csubout=check_subset_vector(col_subset, in->get_ncol());
    const size_t cslen=csubout.size();

    V incoming(in->get_nrow());
    size_t start_row=0, end_row=0;
    if (rslen) {
        // Cutting out extraction costs for unneeded start/end elements.
        start_row=*std::min_element(rsubout.begin(), rsubout.end());
        end_row=*std::max_element(rsubout.begin(), rsubout.end())+1;
    }

    // Setting up the output structures.
    Rcpp::NumericVector libsizes(cslen);
    Rcpp::NumericVector outgoing(rslen), averaged(rslen);

    beachmat::output_param oparam(inmat, false, true);
    oparam.set_chunk_dim(rslen, 1); // pure-column chunks for random access, if HDF5.
    auto omat=beachmat::create_numeric_output(rslen, cslen, oparam);

    auto lbIt=libsizes.begin();
    size_t cs=0;
    for (const auto& c : csubout) {

        // Extracting the column, subsetting the rows.
        auto inIt=in->get_const_col(c, incoming.begin(), start_row, end_row);
        auto oIt=outgoing.begin();
        for (const auto& r : rsubout) {
            (*oIt)=*(inIt + r - start_row);
            ++oIt;
        }
           
        // Dividing by the library size. 
        const double& curlib=((*lbIt)=std::accumulate(outgoing.begin(), outgoing.end(), 0.0));
        if (curlib < 0.00000001) {
            throw std::runtime_error("cells should have non-zero library sizes");
        }
        ++lbIt;

        for (double& out : outgoing) {
            out/=curlib;
        }
        omat->set_col(cs, outgoing.begin());
        ++cs;

        // Adding to the average.
        oIt=outgoing.begin();
        for (auto& a : averaged) {
            a+=(*oIt);
            ++oIt;
        }
    }

    // Dividing the average.
    for (auto& a : averaged) {
        a/=cslen;
    }

    return Rcpp::List::create(libsizes, omat->yield(), averaged); 
}

SEXP subset_and_divide(SEXP matrix, SEXP row_subset, SEXP col_subset) {
    BEGIN_RCPP
    int rtype=beachmat::find_sexp_type(matrix);
    if (rtype==INTSXP) {
        auto input=beachmat::create_integer_matrix(matrix);
        return subset_and_divide_internal<Rcpp::IntegerVector>(input.get(), matrix, row_subset, col_subset);
    } else {
        auto input=beachmat::create_numeric_matrix(matrix);
        return subset_and_divide_internal<Rcpp::NumericVector>(input.get(), matrix, row_subset, col_subset);
    }
    END_RCPP
}

/*** A function to estimate the pooled size factors and construct the linear equations. ***/

SEXP forge_system (SEXP exprs, SEXP ref, SEXP ordering, SEXP poolsizes) {
    BEGIN_RCPP
    auto emat=beachmat::create_numeric_matrix(exprs);
    const size_t ngenes=emat->get_nrow();
    const size_t ncells=emat->get_ncol();
    if (ncells==0) { throw std::runtime_error("at least one cell required for normalization"); }
   
    // Checking the input sizes.
    Rcpp::IntegerVector pool_sizes(poolsizes);
    const size_t nsizes=pool_sizes.size();
    if (nsizes==0) {
        throw std::runtime_error("sizes should be a non-empty integer vector"); 
    }
    int last_size=-1, total_SIZE=0;
    for (const auto& SIZE : pool_sizes) { 
        if (SIZE < 1 || SIZE > ncells) { throw std::runtime_error("each element of sizes should be within [1, number of cells]"); }
        if (SIZE < last_size) { throw std::runtime_error("sizes should be sorted"); }
        total_SIZE+=SIZE;
        last_size=SIZE;
    }

    // Checking pseudo cell.
    Rcpp::NumericVector pseudo_cell(ref);
    if (ngenes!=pseudo_cell.size()) { throw std::runtime_error("length of pseudo-cell vector is not the same as the number of cells"); }

    // Checking ordering.
    Rcpp::IntegerVector order(ordering);
    if (order.size() < ncells*2-1)  { throw std::runtime_error("ordering vector is too short for number of cells"); }
    for (const auto& o : order) { 
        if (o < 0 || o > ncells) { 
            throw std::runtime_error("elements of ordering vector are out of range");
        }
    }

    // Filling up the cell vector.
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

    // Setting up the output vectors.
    Rcpp::IntegerVector row_num(total_SIZE*ncells), col_num(total_SIZE*ncells);
    Rcpp::NumericVector pool_factor(nsizes*ncells);

    // Various other bits and pieces.
    std::vector<double> combined(ngenes), ratios(ngenes);
    const bool is_even=bool(ngenes%2==0);
    const int halfway=int(ngenes/2);
    auto rowIt=row_num.begin(), colIt=col_num.begin();
    auto orIt=order.begin();

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


