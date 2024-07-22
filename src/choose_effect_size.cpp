#include "Rcpp.h"

#include "utils.h"

#include <vector>
#include <stdexcept>
#include <string>
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
