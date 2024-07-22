#include "Rcpp.h"

#include "beachmat3/beachmat.h"
#include "utils.h"

#include <stdexcept>
#include <deque>
#include <vector>
#include <algorithm>

class wilcoxer {
public:    
    wilcoxer(Rcpp::List groups, int ncells) {
        const size_t ngroups=groups.size();
        nelements.resize(ngroups);
        nzeroes.resize(ngroups);

        for (size_t g=0; g<ngroups; ++g) {
            Rcpp::IntegerVector curgroup(groups[g]);
            sources.push_back(std::vector<int>(curgroup.begin(), curgroup.end()));

            std::vector<int>& last_added=sources.back();
            for (auto& current : last_added) {
                if (current<0 || current>=ncells) { 
                    throw std::runtime_error("indices in 'groups' out of range");
                }
            }

            by_group.push_back(std::vector<double>(curgroup.size()));
        }
        return;
    }

    void initialize(const double* vec) {
        // Sorting expression values within each group,
        // special-casing the behavior at zero.
        for (size_t i=0; i<by_group.size(); ++i) {
            auto& cur_group = by_group[i];
            auto cgIt = cur_group.begin();
            const auto& cur_source=sources[i];

            for (auto csIt=cur_source.begin(); csIt!=cur_source.end(); ++csIt) { 
                auto val = vec[*csIt];
                if (val) {
                    *(cgIt++) = val;
                }
            }

            size_t filled = cgIt - cur_group.begin();
            nzeroes[i] = cur_group.size() - filled;
            if (nzeroes[i]) {
                ++filled;
                *cgIt = 0; // adding a placeholder zero, see below.
            }

            nelements[i] = filled;
            std::sort(cur_group.begin(), cur_group.begin() + filled);
        }
        return;
    }

    std::pair<double, int> contrast_groups(int left, int right, double shift) const {
        const auto& group1 = by_group[left];
        const int nelements1 = nelements[left];
        const auto& group2 = by_group[right];
        const int nelements2 = nelements[right];

        std::pair<double, double> output;
        double& score=output.first;
        double& tieval=output.second;

        int c1=0, c2=0, below = 0;
        while (1) {
            // Each iteration of this loop should represent one set of tied ranks
            // (after adjusting 'group1' for the shift).
            const bool ok1=c1 < nelements1;
            const bool ok2=c2 < nelements2;
            double curval;

            if (!ok1 && !ok2) {
                break;
            } else if (ok1 && ok2) {
                curval=std::min(group1[c1] - shift, group2[c2]);
            } else if (ok1) {
                curval=group1[c1] - shift;
            } else {
                curval=group2[c2];
            }

            // Make each index point to first element outside of the range.
            int ties1=0;
            if (ok1) { 
                if (group1[c1]) {
                    const int c1_old=c1;
                    while (c1 < nelements1 && group1[c1] - shift == curval) {
                        ++c1;
                    } 
                    ties1 = c1-c1_old;
                } else if (-shift == curval) { // handling the placeholder zero.
                    ties1 = nzeroes[left];
                    ++c1;
                }
            }

            int ties2 = 0;
            if (ok2) {
                if (group2[c2]) {
                    const int c2_old=c2;
                    while (c2 < nelements2 && group2[c2] == curval) { 
                        ++c2;
                    }
                    ties2=c2 - c2_old;
                } else if (0 == curval) { // handling the placeholder zero.
                    ties2 = nzeroes[right];
                    ++c2;
                }
            }

            score += (static_cast<double>(below) + static_cast<double>(ties2)*0.5) * ties1;
            below += ties2;
                
            // To recapitulate the normal approximation in stats::wilcox.test(). 
            const double nties=ties1 + ties2;
            tieval += nties * (nties * nties - 1);
        }

        return output;
    }

    bool empty_group(int i) const {
        return (sources[i].size()==0);
    }
private:
    std::deque<std::vector<int> > sources;
    std::deque<std::vector<double> > by_group;
    std::deque<int> nelements, nzeroes;
};

// [[Rcpp::export(rng=false)]]
Rcpp::List overlap_exprs(Rcpp::RObject exprs, Rcpp::List groups, double lfc) {
    auto mat = beachmat::read_lin_block(exprs);
    const size_t ngenes=mat->get_nrow();
    const size_t ncells=mat->get_ncol();
   
    // Constructing groups. 
    const size_t ngroups=groups.size();
    wilcoxer wilcox_calc(groups, ncells);

    // Setting up the output matrices to hold the overlap proportions, number of ties.
    Rcpp::List pout(ngroups), tout(ngroups);
    std::vector<std::vector<Rcpp::NumericMatrix::iterator> > pptrs(ngroups), tptrs(ngroups);

    for (size_t i=0; i<ngroups; ++i) {
        size_t ncols=(lfc==0 ? i : ngroups);
        Rcpp::NumericMatrix tmpP(ngenes, ncols), tmpT(ngenes, ncols);
        auto pIt=tmpP.begin(), tIt=tmpT.begin();
        pout[i]=tmpP;
        tout[i]=tmpT;

        pptrs[i].reserve(ncols);
        tptrs[i].reserve(ncols);
        for (size_t j=0; j<ncols; ++j, pIt+=ngenes, tIt+=ngenes) {
            pptrs[i].push_back(pIt);
            tptrs[i].push_back(tIt);
        }
    }

    // Running through all genes and computing pairwise overlaps. 
    std::vector<double> tmp(ncells);
    for (size_t g=0; g<ngenes; ++g) {
        auto ptr = mat->get_row(g, tmp.data());
        wilcox_calc.initialize(ptr);

        for (size_t i1=0; i1<ngroups; ++i1) {
            if (wilcox_calc.empty_group(i1)) { continue; }
            for (size_t i2=0; i2<i1; ++i2) {
                if (wilcox_calc.empty_group(i2)) { continue; }

                auto left=wilcox_calc.contrast_groups(i1, i2, lfc);
                *(pptrs[i1][i2]++)=left.first;
                *(tptrs[i1][i2]++)=left.second;

                if (lfc!=0) {
                    auto right=wilcox_calc.contrast_groups(i1, i2, -lfc);
                    *(pptrs[i2][i1]++)=right.first;
                    *(tptrs[i2][i1]++)=right.second;
                }
            }
        } 
    }

    return Rcpp::List::create(pout, tout);
}

// [[Rcpp::export(rng=false)]]
Rcpp::NumericMatrix overlap_exprs_paired(Rcpp::RObject exprs, Rcpp::IntegerVector left, Rcpp::IntegerVector right, Rcpp::List groups, double lfc) {
    auto mat = beachmat::read_lin_block(exprs);
    const size_t ngenes=mat->get_nrow();
    const size_t ncells=mat->get_ncol();

    // Constructing groups. 
    const size_t ngroups=groups.size();
    wilcoxer wilcox_calc(groups, ncells);

    // Setting up the output matrices to hold the overlap proportions, number of ties.
    Rcpp::NumericMatrix overlaps(left.size(), ngenes);
    auto oIt = overlaps.begin();

    std::vector<double> tmp(ncells);
    for (size_t g=0; g<ngenes; ++g) {
        auto ptr = mat->get_row(g, tmp.data());
        wilcox_calc.initialize(ptr);

        for (size_t i=0; i<left.size(); ++i, ++oIt) {
            *oIt = wilcox_calc.contrast_groups(left[i] - 1, right[i] - 1, lfc).first;
        }
    }

    return overlaps;
}
