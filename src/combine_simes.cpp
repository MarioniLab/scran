#include "Rcpp.h"

#include "utils.h"

#include <vector>
#include <stdexcept>
#include <cmath>
#include <algorithm>

// [[Rcpp::export(rng=false)]]
Rcpp::NumericVector combine_simes(Rcpp::List Pvals, bool logp) {
    const size_t ncon=Pvals.size();

    std::vector<Rcpp::NumericVector> individual(ncon);
    size_t ngenes=0;
    for (size_t c=0; c<ncon; ++c) {
        auto& current=(individual[c]=Pvals[c]);
        if (c==0) {
            ngenes=current.size();
        } else if (ngenes!=static_cast<size_t>(current.size())) {
            throw std::runtime_error("p-value vectors must be of the same length");           
        }
    }
    
    // Should we process these values as if they were log-transformed?
    Rcpp::NumericVector output(ngenes, (logp ? 0 : 1));
    std::vector<double> collected(ncon);

    for (size_t g=0; g<ngenes; ++g) {
        size_t nonna=0;
        for (size_t c=0; c<ncon; ++c) {
            const auto& current=individual[c][g];
            if (!ISNA(current)) { 
                collected[nonna]=current;
                ++nonna;
            }
        }

        if (nonna==0) {
            output[g]=R_NaReal;
            continue;
        }
        std::sort(collected.begin(), collected.begin()+nonna);

        double& minval=output[g];
        for (size_t i=0; i<nonna; ++i) {
            auto P=collected[i];

            if (logp) {
                P-=std::log(i+1);
            } else {
                P/=i+1;
            }
            
            if (P<minval) { 
                minval=P;
            }
        }

        if (logp) { 
            minval+=std::log(nonna);
        } else {
            minval*=nonna;
        }
    }

    return output;
}

// [[Rcpp::export(rng=false)]]
Rcpp::NumericVector combine_holm_middle(Rcpp::List Pvals, bool logp, double prop) {
    const size_t ncon=Pvals.size();

    std::vector<Rcpp::NumericVector> individual(ncon);
    size_t ngenes=0;
    for (size_t c=0; c<ncon; ++c) {
        auto& current=(individual[c]=Pvals[c]);
        if (c==0) {
            ngenes=current.size();
        } else if (ngenes!=static_cast<size_t>(current.size())) {
            throw std::runtime_error("p-value vectors must be of the same length");           
        }
    }
    
    // Should we process these values as if they were log-transformed?
    Rcpp::NumericVector output(ngenes, (logp ? 0 : 1));
    std::vector<double> collected(ncon);

    for (size_t g=0; g<ngenes; ++g) {
        size_t nonna=0;
        for (size_t c=0; c<ncon; ++c) {
            const auto& current=individual[c][g];
            if (!ISNA(current)) { 
                collected[nonna]=current;
                ++nonna;
            }
        }

        if (nonna==0) {
            output[g]=R_NaReal;
            continue;
        }
        std::sort(collected.begin(), collected.begin()+nonna);

        double last=(logp ? R_NegInf : 0);
        double mult=nonna;
        for (size_t i=0; i<nonna; ++i, --mult) {
            auto& P=collected[i];
            if (logp) {
                P+=std::log(mult);
            } else {
                P*=mult;
            }
            if (P < last) {
                P = last;
            } else {
                last = P;
            }
        }

        double& chosen=output[g];
        const size_t jump=std::floor(nonna * prop); // -1 for zero indexing, cancels out the +1.
        std::nth_element(collected.begin(), collected.begin() + jump, collected.begin() + nonna);
        chosen=collected[jump];
        if (logp) {
            if (chosen > 0) { chosen=0; }
        } else {
            if (chosen > 1) { chosen=1; }
        }
    }

    return output;
}
