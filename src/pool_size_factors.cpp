#include "Rcpp.h"
#include "beachmat3/beachmat.h"

#include <stdexcept>
#include <algorithm>
#include <vector>

/*** A function to estimate the pooled size factors and construct the linear equations. ***/

// [[Rcpp::export(rng=false)]]
Rcpp::List pool_size_factors (Rcpp::RObject exprs, Rcpp::NumericVector pseudo_cell, 
    Rcpp::IntegerVector order, Rcpp::IntegerVector pool_sizes) 
{
    auto emat = beachmat::read_lin_block(exprs);
    const size_t ngenes = emat->get_nrow();
    const size_t ncells = emat->get_ncol();

    /******************************
     * Applying checks on inputs. *
     ******************************/

    if (ncells==0) { 
        throw std::runtime_error("at least one cell required for normalization"); 
    }
    if (ngenes==0) {
        throw std::runtime_error("insufficient features for median calculations");
    }

    // Checking the input sizes.
    const size_t nsizes=pool_sizes.size();
    if (nsizes==0) {
        return Rcpp::List::create(Rcpp::IntegerVector(0), Rcpp::IntegerVector(0), Rcpp::NumericVector(0));
    }

    int last_size=-1, total_size=0;
    for (auto s : pool_sizes) { 
        if (s < 1 || static_cast<size_t>(s) > ncells) { 
            throw std::runtime_error("each element of sizes should be within [1, number of cells]"); 
        }
        if (s < last_size) { 
            throw std::runtime_error("sizes should be sorted"); 
        }
        total_size+=s;
        last_size=s;
    }

    // Checking pseudo cell.
    if (ngenes!=pseudo_cell.size()) { 
        throw std::runtime_error("length of pseudo-cell vector is not the same as the number of cells"); 
    }

    // Checking ordering.
    if (order.size() < ncells*2-1)  { 
        throw std::runtime_error("ordering vector is too short for number of cells"); 
    }
    for (auto o : order) { 
        if (o < 0 || static_cast<size_t>(o) >= ncells) { 
            throw std::runtime_error("elements of ordering vector are out of range");
        }
    }

    /*********************************
     * Setting up the storage space. *
     *********************************/

    // Deciding whether or 
    const bool is_sparse = (emat->is_sparse());
    std::unique_ptr<beachmat::lin_sparse_matrix> smat;
    if (is_sparse) { 
        smat = beachmat::promote_to_sparse(emat);
    }

    std::vector<const double*> xptrs(last_size);
    std::vector<std::vector<double> > work_x(last_size, std::vector<double>(ngenes));

    std::vector<const int*> iptrs;
    std::vector<std::vector<int> > work_i;
    std::vector<size_t> nnzero;
    if (is_sparse) {
        work_i.resize(last_size, std::vector<int>(ngenes));
        iptrs.resize(last_size);
        nnzero.resize(last_size);
    }
     
    // The first vector is unfilled as it gets dropped and refilled in the first iteration anyway.
    // We also trigger generation of indices so that each element doesn't generate its own indices.
    auto orIt_tail=order.begin();
    for (int s=1; s<last_size; ++s, ++orIt_tail) {
        if (is_sparse) {
            auto idxs = smat->get_col(*orIt_tail, work_x[s].data(), work_i[s].data());
            xptrs[s] = idxs.x;
            iptrs[s] = idxs.i;
            nnzero[s] = idxs.n;
        } else {
            xptrs[s] = emat->get_col(*orIt_tail, work_x[s].data());
        }
    }

    // Setting up the output vectors and other bits and pieces.
    Rcpp::IntegerVector row_num(total_size*ncells), col_num(total_size*ncells);
    Rcpp::NumericVector pool_factor(nsizes*ncells);

    std::vector<double> combined(ngenes), ratios(ngenes);

    /******************************
     * Performing the iterations. *
     ******************************/

    auto rowIt=row_num.begin(), colIt=col_num.begin();
    auto orIt=order.begin();
    const bool is_even=bool(ngenes%2==0);
    const int halfway=int(ngenes/2);

    // Running through the sliding windows.
    for (size_t win=0; win<ncells; ++win, ++orIt, ++orIt_tail) {
        std::fill(combined.begin(), combined.end(), 0);

        /* Rotating; effectively moves the first element of 'collected' to the end.
         * The is the same as shifting the column that we've moved past to the end,
         * and then overwriting it with the next column.
         */
        std::rotate(xptrs.begin(), xptrs.begin()+1, xptrs.end());
        std::rotate(work_x.begin(), work_x.begin()+1, work_x.end());

        if (!is_sparse) {
            xptrs.back() = emat->get_col(*orIt_tail, work_x.back().data());
        } else {
            std::rotate(iptrs.begin(), iptrs.begin()+1, iptrs.end());
            std::rotate(work_i.begin(), work_i.begin()+1, work_i.end());
            std::rotate(nnzero.begin(), nnzero.begin()+1, nnzero.end());

            auto idxs = smat->get_col(*orIt_tail, work_x.back().data(), work_i.back().data());
            xptrs.back() = idxs.x;
            iptrs.back() = idxs.i;
            nnzero.back() = idxs.n;
        }

        int index=0;
        int rownum=win; // Setting the row so that all pools with the same SIZE form consecutive equations.
        for (auto psIt=pool_sizes.begin(); psIt!=pool_sizes.end(); ++psIt, rownum+=ncells) { 
            const int& SIZE=(*psIt);
            std::fill(rowIt, rowIt+SIZE, rownum);
            rowIt+=SIZE;
            std::copy(orIt, orIt+SIZE, colIt);
            colIt+=SIZE;

            for (; index<SIZE; ++index) {
                auto val = xptrs[index];

                if (is_sparse) {
                    auto n = nnzero[index];
                    auto rows = iptrs[index];
                    for (auto i=0; i<n; ++i, ++val, ++rows) {
                        combined[*rows]+=*val;
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
}
