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
    bool mat_return=check_logical_scalar(return_mat, "return matrix specification");
    
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
    const double s2=check_numeric_scalar(sigma, "sigma");

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

        // Compute log-probabilities using a Gaussian kernel based on the squared distances.
        // We keep things logged to avoid float underflow, and just ignore the constant bit at the front.
        for (auto& d2 : distances2) { 
            d2/=-s2;
        }
        
        // We sum the probabilities to get the relative MNN density. 
        // This requires some care as the probabilities are still logged at this point.
        double density=NA_REAL;
        for (const auto& other_mnn : mnncell) {
            if (ISNA(density)) {
                density=distances2[other_mnn];
            } else {
                const double& to_add=distances2[other_mnn];
                const double larger=std::max(to_add, density), diff=std::abs(to_add - density);
                density=larger + log1pexp(-diff);
            }
        }

        // Each correction vector is weighted by the Gaussian probability (to account for distance)
        // and density (to avoid being dominated by high-density regions).
        // Summation (and then division, see below) yields smoothed correction vectors.
        const auto& correction = averages[mnn];
        auto oIt=output.begin();
        for (int other=0; other<ncells; ++other) {
            const double mult=std::exp(distances2[other] - density);
            totalprob[other]+=mult;

            for (const auto& corval : correction) { 
                (*oIt)+=corval*mult;
                ++oIt;
            }
        }
    }

    // Dividing by the total probability.
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

/* Perform variance adjustment with weighted distributions */

double sq_distance_to_line(Rcpp::NumericMatrix::Column::const_iterator ref, 
                           std::vector<double>::const_iterator grad,
                           Rcpp::NumericMatrix::Column::const_iterator point,
                           std::vector<double>& working) {
    for (auto& w : working) { // Calculating the vector difference from "point" to "ref".
        w=*ref - *point;
        ++point;
        ++ref;
    }

    // Calculating the vector difference from "point" to the line, and taking its norm.
    const double scale=std::inner_product(working.begin(), working.end(), grad, 0.0);
    double dist=0;
    for (auto& w : working) {
        w -= scale * (*grad);
        ++grad;
        dist+=w*w;
    }
   
    return dist;
}

SEXP adjust_shift_variance(SEXP data1, SEXP data2, SEXP vect, SEXP sigma) {
    BEGIN_RCPP
    const Rcpp::NumericMatrix _data1(data1), _data2(data2), _vect(vect);
    const size_t ngenes=_data1.nrow();
    if (ngenes!=_data2.nrow() || ngenes!=_vect.ncol()) { 
        throw std::runtime_error("number of genes do not match up between matrices");
    }
    const size_t ncells1=_data1.ncol(), ncells2=_data2.ncol();
    if (ncells2!=_vect.nrow()) {
        throw std::runtime_error("number of cells do not match up between matrices");
    }        
    const double s2=check_numeric_scalar(sigma, "sigma");

    std::vector<double> grad(ngenes), working(ngenes);
    std::vector<std::pair<double, double> > distance1(ncells1); 
    Rcpp::NumericVector output(ncells2);

    // Iterating through all cells.
    for (size_t cell=0; cell<ncells2; ++cell) {
        const auto curcell=_data2.column(cell);

        // Calculating the l2 norm and adjusting to a unit vector.
        double l2norm=0;
        {
            const auto curvect=_vect.row(cell);
            std::copy(curvect.begin(), curvect.end(), grad.begin());
            for (const auto& g : grad) {
                l2norm+=g*g;            
            }
            l2norm=std::sqrt(l2norm);
            for (auto& g : grad) { 
                g/=l2norm;
            }
        }
                
        const double curproj=std::inner_product(grad.begin(), grad.end(), curcell.begin(), 0.0);

        // Getting the cumulative probability of each cell in its own batch.
        double prob2=0, totalprob2=0;    
        for (size_t same=0; same<ncells2; ++same) {
            if (same==cell) { 
                prob2+=1;
                totalprob2+=1;
            } else {
                const auto samecell=_data2.column(same);
                const double sameproj=std::inner_product(grad.begin(), grad.end(), samecell.begin(), 0.0); // Projection
                const double samedist=sq_distance_to_line(curcell.begin(), grad.begin(), samecell.begin(), working); // Distance.
                
                const double sameprob=std::exp(-samedist/s2);
                if (sameproj <= curproj) {
                    prob2+=sameprob;
                }
                totalprob2+=sameprob;
            }
        }
        prob2/=totalprob2;

        // Filling up the coordinates and weights for the reference batch.
        double totalprob1=0;
        for (size_t other=0; other<ncells1; ++other) {
            const auto othercell=_data1.column(other);
            distance1[other].first=std::inner_product(grad.begin(), grad.end(), othercell.begin(), 0.0); // Projection
            const double otherdist=sq_distance_to_line(curcell.begin(), grad.begin(), othercell.begin(), working); // Distance.
            totalprob1+=(distance1[other].second=std::exp(-otherdist/s2));
        }
        std::sort(distance1.begin(), distance1.end());

        // Choosing the quantile in the projected reference coordinates that matches the cumulative probability in its own batch.
        const double target=prob2*totalprob1;
        double cumulative=0, ref_quan=R_NaReal;
        if (ncells1) { 
            ref_quan=distance1.back().first;
        }

        for (const auto& val : distance1) {
            cumulative+=val.second;
            if (cumulative >= target) { 
                ref_quan=val.first;
                break;
            }
        }

        // Distance between quantiles represents the scaling of the original vector.        
        output[cell]=(ref_quan - curproj)/l2norm;
    }
    
    return(output);
    END_RCPP
}

