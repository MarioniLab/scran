#include "scran.h"

#include "utils.h"

#include <memory>
#include <stdexcept>
#include <cmath>

struct FUN {
    FUN(double lim) : limit(lim) {}
    virtual ~FUN() {};
    virtual double q(double, bool)=0;
    virtual double d(int, double)=0;
protected:
    double limit;    
};

struct pois : public FUN {
    pois(double lim) : FUN(lim) {}
    ~pois() {}
    double q(double m, bool lower) {
        return R::qpois(limit, m, lower, 0);
    }
    double d(int y, double m) {
        return R::dpois(y, m, 0);
    }
};

struct nbinom : public FUN {
    nbinom(double lim, double disp): FUN(lim), size(1/disp) {}
    ~nbinom() {}
    double q(double m, bool lower) {
        return R::qnbinom_mu(limit, size, m, lower, 0);
    }
    double d(int y, double m) {
        return R::dnbinom_mu(y, size, m, 0);
    }
protected:
    double size;
};

std::unique_ptr<FUN> choose_dist (SEXP tol, SEXP disp) {
    const double lim=check_numeric_scalar(tol, "tolerance");    
    const double phi=check_numeric_scalar(disp, "dispersion");    
    if (lim <= 0) {
        throw std::runtime_error("tolerance must be a positive double-precision value");
    }
    if (phi < 0) { 
        throw std::runtime_error("dispersion must be a non-negative double-precision value");
    }
    if (phi==0) {
        return std::unique_ptr<FUN>(new pois(lim));
    } else {
        return std::unique_ptr<FUN>(new nbinom(lim, phi));
    }
}

/********************************************/

double get_pseudo (SEXP offset) {
    const double pseudo=check_numeric_scalar(offset, "offset");    
    if (pseudo <= 0) { 
        throw std::runtime_error("offset must be a positive double-precision value");
    }
    return pseudo;
}

double get_mean (double m, double s, double p, FUN * ptr) {
    const double mu=m*s;
    const int left=ptr->q(mu, true);
    const int right=ptr->q(mu, false);
    double out_mean=0;
    double total_prob=0;
    
    for (auto i=left; i<=right; ++i) {
        const double cur_prob = ptr->d(i, mu);
        out_mean += cur_prob * std::log(i/s+p);
        total_prob += cur_prob;
    }
            
    out_mean/=total_prob;
    return out_mean;
}

double get_var (double m, double s, double p, double grand_mean, FUN * ptr) {
    const double mu=m*s;
    const int left=ptr->q(mu, true);
    const int right=ptr->q(mu, false);
    double out_var=0, total_prob=0;
    
    for (auto i=left; i<=right; ++i) {
        const double cur_prob = ptr->d(i, mu);
        const double tmp=std::log(i/s+p) - grand_mean;
        out_var+=cur_prob * tmp * tmp;
        total_prob += cur_prob;
    }
    
    out_var/=total_prob;
    return out_var;
}

/*******************************************
 * Split into three functions.
 *  - calc_log_count_stats computes the per-gene mean and variance
 *  - calc_log_expected computes the expected value for each log-count 
 *  - calc_log_sqdiff computes the expected squared difference from a per-gene constant.
 */

SEXP calc_log_count_stats (SEXP means, SEXP sf, SEXP tol, SEXP disp, SEXP offset) {
    BEGIN_RCPP
    Rcpp::NumericVector Means(means), Sizes(sf);
    const auto pseudo=get_pseudo(offset);
    auto ptr=choose_dist(tol, disp);

    // Setting up the output.
    const size_t Nmeans = Means.size();
    Rcpp::NumericVector outputm(Nmeans), outputv(Nmeans);

    for (size_t i=0; i < Nmeans; ++i) { 
        const auto& m=Means[i];
        
        // Computing the grand mean of the expectation of the log-counts for each gene.
        double& grand_mean=outputm[i];
        for (const auto& s : Sizes) {
            grand_mean += get_mean(m, s, pseudo, ptr.get());
        }
        grand_mean/=Sizes.size();

        // Computing the variance for each gene, as the expected (X - E[E[X]])^2 where E[E[X]] is the grand mean.
        double& grand_var=outputv[i];
        for (const auto& s : Sizes) {
            grand_var += get_var(m, s, pseudo, grand_mean, ptr.get());
        }
        grand_var/=Sizes.size();

        // Adjusting for the base of the log.
        grand_mean/=M_LN2;
        grand_var/=M_LN2 * M_LN2;
    }    

    return Rcpp::List::create(outputm, outputv);
    END_RCPP
}

SEXP calc_log_expected (SEXP means, SEXP sf, SEXP tol, SEXP disp, SEXP offset) {
    BEGIN_RCPP
    Rcpp::NumericVector Means(means), Sizes(sf);
    const auto pseudo=get_pseudo(offset);
    auto ptr=choose_dist(tol, disp);

    const size_t Nmeans = Means.size(), Nsizes = Sizes.size();
    Rcpp::List outputm(Nmeans);

    for (size_t i=0; i < Nmeans; ++i) { 
        const auto& m=Means[i];

        Rcpp::NumericVector curm(Nsizes);
        for (size_t j=0; j < Nsizes; ++j) {
            curm[j] = get_mean(m, Sizes[j], pseudo, ptr.get()) / M_LN2;
        }

        outputm[i]=curm;
    }
    
    return outputm;
    END_RCPP
}

SEXP calc_log_sqdiff (SEXP means, SEXP sf, SEXP tol, SEXP disp, SEXP offset, SEXP constant) {
    BEGIN_RCPP
    Rcpp::NumericVector Means(means), Sizes(sf), Constants(constant);
    if (Constants.size()!=Means.size()) {
        throw std::runtime_error("'constant' and 'means' should be of the same length");
    }
    const auto pseudo=get_pseudo(offset);
    auto ptr=choose_dist(tol, disp);
        
    const size_t Nmeans = Means.size(), Nsizes = Sizes.size();
    Rcpp::List outputv(Nmeans);

    for (size_t i=0; i < Nmeans; ++i) { 
        const auto& m=Means[i];
        const double curconst=Constants[i] * M_LN2; // get back to natural log.

        Rcpp::NumericVector curv(Nsizes);
        for (size_t j=0; j < Nsizes; ++j) {
            curv[j] = get_var(m, Sizes[j], pseudo, curconst, ptr.get()) / M_LN2 / M_LN2;
        }
        outputv[i]=curv;
    }
    
    return outputv;
    END_RCPP
}


