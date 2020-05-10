#include "Rcpp.h"

#include "rand_custom.h"
#include "scuttle/linear_model_fit.h"
#include "boost/range/algorithm.hpp"
#include "utils.h"

#include <stdexcept>
#include <vector>
#include <algorithm>

static double rho_mult (double Ncells) {
    return 6/(Ncells*(Ncells*Ncells-1));
}

/*** Null distribution estimation without a design matrix. ***/

// [[Rcpp::export(rng=false)]]
Rcpp::NumericVector get_null_rho (int Ncells, int Niters, Rcpp::List Seeds, Rcpp::IntegerVector Streams) {
    if (Ncells <= 1) { throw std::runtime_error("number of cells should be greater than 2"); }
    if (Niters < 0) { throw std::runtime_error("number of iterations should be non-negative"); }
    check_pcg_vectors(Seeds, Streams, Niters, "iterations");

    // Filling rank vector.
    std::vector<int> rankings(Ncells);
    Rcpp::NumericVector output(Niters);
    const double mult=rho_mult(Ncells);

    for (int it=0; it<Niters; ++it) {
        std::iota(rankings.begin(), rankings.end(), 0);

        auto generator=create_pcg32(Seeds[it], Streams[it]);
        boost::range::random_shuffle(rankings, generator);

        double tmp=0;
        for (int cell=0; cell<Ncells; ++cell) {
            const double tmpdiff=rankings[cell]-cell;
            tmp+=tmpdiff*tmpdiff;
        }
        tmp*=mult;
        output[it]=1-tmp;            
    }

    return output;
}

/*** Null distribution estimation with a design matrix. ***/

// [[Rcpp::export(rng=false)]]
Rcpp::NumericVector get_null_rho_design (Rcpp::NumericMatrix qr, Rcpp::NumericVector qraux, 
    int Niters, Rcpp::List Seeds, Rcpp::IntegerVector Streams) 
{
    if (Niters <= 0) { throw std::runtime_error("number of iterations should be positive"); }
    check_pcg_vectors(Seeds, Streams, Niters, "iterations");

    // Setting up to multiply by the Q matrix.
    scuttle::QR_multiplier multQ(qr, qraux, 'N');
    const int Nobs=multQ.get_nobs();
    const int Ncoef=multQ.get_ncoefs();

    Rcpp::NumericVector output(Niters);
    std::vector<double> holding_val(Nobs);
    std::vector<int> holding_idx(Nobs);
    std::vector<int> rank1(Nobs), rank2(Nobs);
    const double mult=rho_mult(Nobs);

    /* Simulating residuals, using the Q-matrix to do it. We set the main effects to zero 
     * (hence, starting from "Ncoefs") and simulate normals for the residual effects.
     * We then use this to reconstruct the residuals themselves - twice - and then compute 
     * correlations between the two reconstructions.
     */
    for (int it=0; it<Niters; ++it) {
        auto generator=create_pcg32(Seeds[it], Streams[it]);
        boost::random::normal_distribution<double> cpp_rnorm;

        for (int mode=0; mode<2; ++mode) {
            // Simulating the residuals.
            std::fill(holding_val.begin(), holding_val.begin()+Ncoef, 0);
            for (int col=Ncoef; col<Nobs; ++col) {
                holding_val[col]=cpp_rnorm(generator);
            }
            multQ.multiply(holding_val.data());

            // Computing ranks.
            std::iota(holding_idx.begin(), holding_idx.end(), 0);
            rsort_with_index(holding_val.data(), holding_idx.data(), Nobs); 

            auto hdIt=holding_idx.begin();
            auto& rank=(mode ? rank1 : rank2);
            for (int row=0; row<Nobs; ++row, ++hdIt) {
                rank[*hdIt]=row;
            }
        }

        // Computing the squared difference in the ranks.
        double& rho=(output[it]=0);
        for (int row=0; row<Nobs; ++row) {
            const double tmpdiff=rank1[row]-rank2[row];
            rho+=tmpdiff*tmpdiff;
        }
        rho*=mult;
        rho=1-rho;
    }

    return output;    
}
