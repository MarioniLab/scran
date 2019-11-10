#include "Rcpp.h"

#include "utils.h"

#include <vector>
#include <stdexcept>
#include <cmath>
#include <algorithm>

template<class V>
size_t instantiate_list(Rcpp::List input, std::vector<V>& output, const std::string msg) {
    size_t n=0;
    for (size_t c=0; c<input.size(); ++c) {
        auto& current=(output[c]=input[c]);
        if (c==0) {
            n=current.size();
        } else if (n!=static_cast<size_t>(current.size())) {
            throw std::runtime_error(msg + " vectors must be of the same length");           
        }
    }
    return n;
}

// [[Rcpp::export(rng=false)]]
Rcpp::NumericVector combine_simes(Rcpp::List Pvals, bool logp) {
    const size_t ncon=Pvals.size();
    std::vector<Rcpp::NumericVector> individual(ncon);
    const size_t ngenes=instantiate_list(Pvals, individual, "p-value");
    
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
    const size_t ngenes=instantiate_list(Pvals, individual, "p-value");

    std::vector<double> collected(ncon);
    Rcpp::NumericVector output(ngenes, R_NaReal);

    for (size_t g=0; g<ngenes; ++g) {
        size_t nonna=0;
        for (size_t c=0; c<ncon; ++c) {
            const auto& current=individual[c][g];
            if (!ISNA(current)) { 
                collected[nonna]=current;
                ++nonna;
            }
        }

        if (nonna!=0) {
            // Computing the Holm series of p-values:
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

            // Picking a middle p-value.
            size_t jump=std::ceil(nonna * prop);
            if (jump) { --jump; } // -1 for the zero-indexing. 
            std::nth_element(collected.begin(), collected.begin() + jump, collected.begin() + nonna);

            double& chosen=(output[g]=collected[jump]);
            if (logp) {
                if (chosen > 0) { chosen=0; }
            } else {
                if (chosen > 1) { chosen=1; }
            }
        }
    }

    return output;
}

size_t define_jump (size_t ntests, double prop) {
    /* We want the p-value for the (prop*ntests)-th rejection, which rejects
     * the null that more than (1-prop) of the  nulls are true. The '-1' is
     * for the zero-indexing but obviously does not apply if prop=0.
     */
    size_t jump=std::ceil(ntests * prop);
    if (jump) { --jump; } 
    return jump;
}

// [[Rcpp::export(rng=false)]]
Rcpp::IntegerVector compute_Top_statistic_from_ranks(Rcpp::List Ranks, double prop) {
    const size_t ncon=Ranks.size();
    std::vector<Rcpp::IntegerVector> individual(ncon);
    const size_t ngenes=instantiate_list(Ranks, individual, "rank");

    std::vector<int> collected(ncon);
    Rcpp::IntegerVector output(ngenes, NA_INTEGER);

    for (size_t g=0; g<ngenes; ++g) {
        size_t nonna=0;
        for (size_t c=0; c<ncon; ++c) {
            const auto& current=individual[c][g];
            if (current!=NA_INTEGER) {
                collected[nonna]=current;
                ++nonna;
            }
        }

        if (nonna!=0) {
            const size_t jump=define_jump(nonna, prop);
            std::nth_element(collected.begin(), collected.begin() + jump, collected.begin() + nonna);
            output[g]=collected[jump];
        }
    }

    return output;
}

// [[Rcpp::export(rng=false)]]
Rcpp::NumericVector choose_middle_effect_size(Rcpp::List Pvals, Rcpp::List Effects, double prop) {
    const size_t ncon=Effects.size();
    std::vector<Rcpp::NumericVector> effects(ncon), pvals(ncon);
    const size_t neffects=instantiate_list(Effects, effects, "effect");
    const size_t npvals=instantiate_list(Pvals, pvals, "p-value");
    if (neffects!=npvals) {
        throw std::runtime_error("p-value and effect vectors should have the same length");
    }

    std::vector<std::pair<double, double> > collected(ncon);
    Rcpp::NumericVector output(neffects, R_NaReal);

    for (size_t g=0; g<neffects; ++g) {
        size_t nonna=0;

        for (size_t c=0; c<ncon; ++c) {
            const auto& curp=pvals[c][g];
            const auto& cure=effects[c][g];

            if (!ISNA(curp) && !ISNA(cure)) {
                collected[nonna]=std::make_pair(curp, cure);
                ++nonna;
            }
        }

        if (nonna!=0) {
            const size_t jump=define_jump(nonna, prop);
            std::nth_element(collected.begin(), collected.begin() + jump, collected.begin() + nonna);
            output[g]=collected[jump].second;
        }
    }

    return output;
}
