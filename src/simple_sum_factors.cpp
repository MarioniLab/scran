#include "scran.h"
#include <cmath>

template<class V, class M>
void get_tricube_average(size_t i, M mat, 
    const Rcpp::IntegerMatrix& indices, const Rcpp::NumericMatrix& distances, 
    double ndist, double* out, typename V::iterator work) 
{
    auto curI=indices.column(i);
    auto curD=distances.column(i);
    
    // Computing the bandwidth.
    const size_t nneighbors=indices.nrow();
    const size_t midpoint=std::ceil(nneighbors/2)-1;
    const double bandwidth = std::max(0.00000001, curD[midpoint] * ndist);

    const size_t ngenes=mat->get_nrow();
    std::fill(out, out+ngenes, 0);
    double total_weight=0;

    // Computing the weighted average.
    for (size_t nn=0; nn<nneighbors; ++nn) {
        const double reldist=curD[nn] / bandwidth;
        if (reldist > 1) { break; }
        
        const double weight0 = 1 - reldist * reldist * reldist;
        const double weight = weight0 * weight0 * weight0;
        total_weight += weight;

        auto it=mat->get_const_col(curI[nn], work);
        auto out_copy=out;
        for (size_t g=0; g<ngenes; ++g, ++it, ++out_copy) {
            *out_copy += (*it) * weight; 
        }
    }

    // Normalizing by the total weight.
    for (size_t g=0; g<ngenes; ++g, ++out) {
        *out /= total_weight;
    }
    return;
}

template<class V, class M>
SEXP simple_sum_factors_internal(M mat, Rcpp::IntegerMatrix indices, Rcpp::NumericMatrix distances, Rcpp::IntegerVector ref_cell, SEXP min_mean, SEXP ndist) 
{
    const size_t ngenes=mat->get_nrow();
    const size_t ncells=mat->get_ncol();

    if (indices.ncol()!=ncells || distances.ncol()!=ncells) {
        throw std::runtime_error("number of columns in nearest neighbor results is not equal to the number of cells");
    }
    if (indices.nrow()!=distances.nrow()) {
        throw std::runtime_error("number of rows in nearest neighbor results are not equal");
    }

    const double mean_threshold=(min_mean!=R_NilValue ? check_numeric_scalar(min_mean, "mean threshold") : 0);
    const double num_dist=check_numeric_scalar(ndist, "bandwidth scaling");
    const size_t ref_dex=check_integer_scalar(ref_cell, "reference cell");

    // Computing the reference profile.
    Rcpp::NumericVector reference(ngenes);
    V work(ngenes);
    get_tricube_average<V>(ref_dex, mat, indices, distances, num_dist, reference.begin(), work.begin());
    const double reflibsize=std::accumulate(reference.begin(), reference.end(), 0.0);

    std::vector<double> local(ngenes), ratios(ngenes);
    Rcpp::NumericVector outfactors(ncells);

    for (size_t i=0; i<ncells; ++i) {
        get_tricube_average<V>(i, mat, indices, distances, num_dist, local.data(), work.begin());
        const double locallibsize=std::accumulate(local.begin(), local.end(), 0.0);

        // Computing ratios, with or without filtering.
        auto rIt=reference.begin();
        size_t postfilter=0;
        const double MULT1=(1 + reflibsize/locallibsize), MULT2=(1 + locallibsize/reflibsize);

        for (size_t g=0; g<ngenes; ++g, ++rIt) {
            if (local[g]==0 && (*rIt)==0) { 
                continue; 
            }

            const double mean_exprs= (local[g] * MULT1 + (*rIt) * MULT2) / 4; // same as (local/locallib + ref/reflib)/2 * (locallib + reflib)/2;
            if (mean_exprs >= mean_threshold) {
                ratios[postfilter]=local[g]/(*rIt);
                ++postfilter;
            }
        }

        if (postfilter) { 
            // Computing the median (faster than partial sort).
            bool is_even=(postfilter%2==0);
            size_t halfway=postfilter/2;
    
    		std::nth_element(ratios.begin(), ratios.begin()+halfway, ratios.begin() + postfilter);
    		if (is_even) {
    			double medtmp=ratios[halfway];
    			std::nth_element(ratios.begin(), ratios.begin()+halfway-1, ratios.begin() + postfilter);
    			outfactors[i]=(medtmp+ratios[halfway-1])/2;
    		} else {
    			outfactors[i]=ratios[halfway];
    		}
        } else {
            outfactors[i]=R_NaReal;
        }

        // Actually computing the current factor.
        auto It=mat->get_const_col(i, work.begin());
        const double curlibsize=std::accumulate(It, It+ngenes, 0.0);
        outfactors[i] *= curlibsize/locallibsize;
    }

    return Rcpp::List::create(outfactors, reference);
}

SEXP simple_sum_factors(SEXP matrix, SEXP indices, SEXP distances, SEXP ref_cell, SEXP min_mean, SEXP ndist) {
    BEGIN_RCPP
    int rtype=beachmat::find_sexp_type(matrix);
    if (rtype==INTSXP) {
        auto input=beachmat::create_integer_matrix(matrix);
        return simple_sum_factors_internal<Rcpp::IntegerVector>(input.get(), indices, distances, ref_cell, min_mean, ndist);
    } else {
        auto input=beachmat::create_numeric_matrix(matrix);
        return simple_sum_factors_internal<Rcpp::NumericVector>(input.get(),  indices, distances, ref_cell, min_mean, ndist);
    }
    END_RCPP
}

