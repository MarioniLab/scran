#include "scran.h"

SEXP find_mutual_nns (SEXP left, SEXP right) {
    BEGIN_RCPP
    const Rcpp::IntegerMatrix L(left), R(right);
    const int nnR=R.ncol();
 
    std::vector<int> sortedR(R.size());
    std::deque<int> mutualL, mutualR;

    // Sorting the elements on the right.
    auto sIt=sortedR.begin();
    for (int r=0; r<R.nrow(); ++r) {
        const auto currow=R.row(r);
        std::copy(currow.begin(), currow.end(), sIt);
        std::sort(sIt, sIt+nnR);
        sIt+=nnR;
    }

    // Running through the elements on the left, and doing a binary search.
    for (int l=0; l<L.nrow(); ++l) {
        const auto currow=L.row(l);
        const int tocheck=l+1;

        for (const auto& curval : currow) {
            auto startIt=sortedR.begin() + nnR*(curval-1); // 1-indexed.
            auto endIt=startIt+nnR;
            auto closest=std::lower_bound(startIt, endIt, tocheck);

            if (closest!=endIt && *closest==tocheck) { 
                mutualL.push_back(tocheck);
                mutualR.push_back(curval);
            }
        }
    }

    return Rcpp::List::create(Rcpp::IntegerVector(mutualL.begin(), mutualL.end()),
                              Rcpp::IntegerVector(mutualR.begin(), mutualR.end()));
    END_RCPP
}

/* Performs the cosine normalization in a fairly efficient manner. */

template<class M>
SEXP cosine_norm_internal (M mat, SEXP original, SEXP return_mat) {
    const size_t& nrow=mat->get_nrow();
    const size_t& ncol=mat->get_ncol();

    // Deciding whether or not to return the matrix.
    Rcpp::LogicalVector retmat(return_mat);
    if (retmat.size()!=1) { 
        throw std::runtime_error("return matrix specification should be a logical vector");
    }
    bool mat_return=retmat[0];
    
    beachmat::numeric_output* optr=NULL;
    std::vector<std::unique_ptr<beachmat::numeric_output> > holder;
    if (mat_return) { 
        holder.push_back(beachmat::create_numeric_output(nrow, ncol, 
                         beachmat::output_param(original, false, true)));
        optr=holder.front().get();
    }
   
    // Calculating the L2 norm of each vector and applying it. 
    Rcpp::NumericVector incoming(nrow);
    Rcpp::NumericVector l2norm(ncol);
    for (size_t c=0; c<ncol; ++c) {
        mat->get_col(c, incoming.begin());

        double& total=l2norm[c];
        for (const auto& val : incoming) { 
            total+=val*val;
        }
        total=std::sqrt(total);
        total=std::max(total, 0.00000001); // avoid division by zero.

        if (mat_return) { 
            for (auto& val : incoming) { 
                val/=total;
            }
            optr->set_col(c, incoming.begin());
        }
    }

    // Figuring out what to return.
    if (mat_return) { 
        return Rcpp::List::create(optr->yield(), l2norm);
    } else {
        return Rcpp::List::create(R_NilValue, l2norm);
    }
}

SEXP cosine_norm(SEXP incoming, SEXP getmat) {
    BEGIN_RCPP
    int rtype=beachmat::find_sexp_type(incoming);
    if (rtype==INTSXP) {
        auto input=beachmat::create_integer_matrix(incoming);
        return cosine_norm_internal(input.get(), incoming, getmat);
    } else {
        auto input=beachmat::create_numeric_matrix(incoming);
        return cosine_norm_internal(input.get(), incoming, getmat);
    }
    END_RCPP
}

/* Perform smoothing with the Gaussian kernel */

SEXP smooth_gaussian_kernel(SEXP vect, SEXP index, SEXP data, SEXP sigma) {
    BEGIN_RCPP
    const Rcpp::NumericMatrix correction_vectors(vect);
    const int npairs=correction_vectors.nrow();
    const int ngenes=correction_vectors.ncol();
    const Rcpp::IntegerVector _index(index);
    if (npairs!=_index.size()) { 
        throw std::runtime_error("number of rows in 'vect' should be equal to length of 'index'");
    }
    
    // Constructing the average vector for each MNN cell.
    std::deque<std::vector<double> > averages;
    std::deque<int> number;
    std::set<int> mnncell;

    int row=0;
    for (const auto& i : _index) {
        auto currow=correction_vectors.row(row);
        ++row;
        
        if (i >= averages.size() || averages[i].empty()) { 
            if (i >= averages.size()) { 
                averages.resize(i+1);
                number.resize(i+1);
            }
            averages[i]=std::vector<double>(currow.begin(), currow.end());
            number[i]=1;
            mnncell.insert(i);
        } else {
            auto& target=averages[i];
            auto cIt=currow.begin();
            for (auto& t : target){
                t += *cIt;
                ++cIt;
            }
            ++(number[i]);
        }
    }

    for (const auto& i : mnncell) {
        auto& target=averages[i];
        const auto& num=number[i];
        for (auto& t : target) {
            t/=num;
        }
    }

    // Setting up input constructs (including the expression matrix on which the distances are computed).
    auto mat=beachmat::create_numeric_matrix(data);
    const int ncells=mat->get_ncol();
    const int ngenes_for_dist=mat->get_nrow();

    Rcpp::NumericVector _sigma(sigma);
    if (_sigma.size()!=1) {
        throw std::runtime_error("sigma should be a double-precision scalar");
    }
    const double s2=_sigma[0];

    // Setting up output constructs.
    Rcpp::NumericMatrix output(ngenes, ncells); // yes, this is 'ngenes' not 'ngenes_for_dist'.
    std::vector<double> distances2(ncells), totalprob(ncells);
    Rcpp::NumericVector mnn_incoming(ngenes_for_dist), other_incoming(ngenes_for_dist);

    // Using distances between cells and MNN-involved cells to smooth the correction vector per cell.
    for (const auto& mnn : mnncell) {
        auto mnn_iIt=mat->get_const_col(mnn, mnn_incoming.begin());

        for (int other=0; other<ncells; ++other) {
            double& curdist2=(distances2[other]=0);
            auto other_iIt=mat->get_const_col(other, other_incoming.begin());
            auto iIt_copy=mnn_iIt;

            for (int g=0; g<ngenes_for_dist; ++g) {
                const double tmp=(*iIt_copy  - *other_iIt);
                curdist2+=tmp*tmp;
                ++other_iIt;
                ++iIt_copy;
            }
        }

        // Compute probabilities based on the squared distances (sum to get the relative MNN density with a Gaussian kernel).
        for (auto& d2 : distances2) { 
            d2=std::exp(-d2/s2); 
        }
        double density=0;
        for (const auto& other_mnn : mnncell) {
            density+=distances2[other_mnn];
        }

        // Filling in the smoothed values, after dividing by density to avoid being dominated by high-density regions.
        const auto& correction = averages[mnn];
        auto oIt=output.begin();
        for (int other=0; other<ncells; ++other) {
            const double mult=distances2[other]/density;
            totalprob[other]+=mult;

            for (const auto& corval : correction) { 
                (*oIt)+=corval*mult;
                ++oIt;
            }
        }
    }

    // Adjusting the smoothed values.
    for (int other=0; other<ncells; ++other) {
        auto curcol=output.column(other);
        const double total=totalprob[other];

        for (auto& val : curcol) {
            val/=total;
        }
    }

    return output;
    END_RCPP
}
