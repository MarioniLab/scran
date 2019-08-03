#include "Rcpp.h"
#include "beachmat/numeric_matrix.h"
#include "beachmat/integer_matrix.h"
#include "utils.h"
#include <vector>

/********************************
 * Central class to compute all variance statistics.
 * This should probably use Welford's algorithm to
 * encourage column-major access in the individual functions.
 ********************************/

class blocked_stats_calculator {
public:
    blocked_stats_calculator(Rcpp::List bygroup, size_t ncells) : groups(bygroup.size()) {
        for (size_t i=0; i<groups.size(); ++i) { 
            groups[i]=check_subset_vector(bygroup[i], ncells);
        }
    }

    void compute(Rcpp::NumericVector::iterator values, double* curmeanrow, double* curvarrow) const {
        for (const auto& curgroup : groups) {
            double& curmean=*curmeanrow;
            ++curmeanrow;
            double& curvar=*curvarrow;
            ++curvarrow;

            // Calculating the mean.          
            if (!curgroup.size()) {
                curmean=R_NaReal;
                curvar=R_NaReal;
                continue; 
            }

            curmean=0;
            for (const auto& index : curgroup) {
                curmean+=*(values + index);
            }
            curmean/=curgroup.size();

            // Computing the variance.
            if (curgroup.size()==1) {
                curvar=R_NaReal;
                continue;
            }

            curvar=0;
            for (const auto& index : curgroup) {
                const double tmp=*(values + index) - curmean;
                curvar += tmp * tmp;
            }
            curvar/=curgroup.size()-1;
        }
        return;
    }

    size_t ngroups() const { return groups.size(); }
private:  
    std::vector<Rcpp::IntegerVector> groups;
};

/********************************
 * Compute statistics for transformed counts.
 ********************************/

template<class M, class TRANSFORMER>
Rcpp::List compute_blocked_stats(Rcpp::List bygroup, SEXP inmat, TRANSFORMER trans) {
    auto emat=beachmat::create_matrix<M>(inmat);
    const size_t ncells=emat->get_ncol();
    const size_t ngenes=emat->get_nrow();

    blocked_stats_calculator BSC(bygroup, ncells);
    auto ngroups=BSC.ngroups();
   
    // Setting up the output objects.
    Rcpp::NumericMatrix outvar(ngroups, ngenes);
    Rcpp::NumericMatrix outmean(ngroups, ngenes);
    Rcpp::NumericVector incoming(ncells);

    for (size_t counter=0; counter<ngenes; ++counter) {
        emat->get_row(counter, incoming.begin());
        trans(incoming.begin(), incoming.end());
        auto curvarrow=outvar.column(counter);
        auto curmeanrow=outmean.column(counter);
        BSC.compute(incoming.begin(), curmeanrow.begin(),  curvarrow.begin());
    }

    return(Rcpp::List::create(outmean, outvar));
}

/********************************
 * Compute statistics for log-transformed counts.
 ********************************/

struct lognorm {
    lognorm(Rcpp::NumericVector sizefactors, double pseudo) : sf(sizefactors), ps(pseudo) {}
    void operator()(Rcpp::NumericVector::iterator start, Rcpp::NumericVector::iterator end) {
        auto sfIt=sf.begin();
        while (start!=end) {
            *start=std::log(*start/(*sfIt) + ps)/M_LN2;
            ++start;
            ++sfIt;
        }
    }
private:
    Rcpp::NumericVector sf;
    double ps;
};

// [[Rcpp::export(rng=false)]]
Rcpp::List compute_blocked_stats_lognorm(Rcpp::List bygroup, SEXP inmat, Rcpp::NumericVector sf, double pseudo) 
{
    int rtype=beachmat::find_sexp_type(inmat);
    lognorm LN(sf, pseudo);
    if (rtype==INTSXP) {
        return compute_blocked_stats<beachmat::integer_matrix>(bygroup, inmat, LN);
    } else {
        return compute_blocked_stats<beachmat::numeric_matrix>(bygroup, inmat, LN);
    }
}

/********************************
 * Compute statistics for normalized counts.
 ********************************/

struct norm {
    norm(Rcpp::NumericVector sizefactors) : sf(sizefactors) {}
    void operator()(Rcpp::NumericVector::iterator start, Rcpp::NumericVector::iterator end) {
        auto sfIt=sf.begin();
        while (start!=end) {
            *start/=*sfIt;
            ++start;
            ++sfIt;
        }
    }
private:
    Rcpp::NumericVector sf;
};

// [[Rcpp::export(rng=false)]]
Rcpp::List compute_blocked_stats_norm(Rcpp::List bygroup, SEXP inmat, Rcpp::NumericVector sf)
{
    int rtype=beachmat::find_sexp_type(inmat);
    norm N(sf);
    if (rtype==INTSXP) {
        return compute_blocked_stats<beachmat::integer_matrix>(bygroup, inmat, N);
    } else {
        return compute_blocked_stats<beachmat::numeric_matrix>(bygroup, inmat, N);
    }
}

/********************************
 * Compute statistics for expression values without further transformation.
 ********************************/

struct none {
    none() {}
    void operator()(Rcpp::NumericVector::iterator start, Rcpp::NumericVector::iterator end) {}
};

// [[Rcpp::export(rng=false)]]
Rcpp::List compute_blocked_stats_none(Rcpp::List bygroup, SEXP inmat) {
    int rtype=beachmat::find_sexp_type(inmat);
    none N;
    if (rtype==INTSXP) {
        return compute_blocked_stats<beachmat::integer_matrix>(bygroup, inmat, N);
    } else {
        return compute_blocked_stats<beachmat::numeric_matrix>(bygroup, inmat, N);
    }
}
