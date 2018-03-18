#include "scran.h"

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

SEXP calc_log_count_stats (SEXP means, SEXP sf, SEXP tol, SEXP disp, SEXP offset) {
    BEGIN_RCPP
    Rcpp::NumericVector Means(means), Sizes(sf);

    const double lim=check_numeric_scalar(tol, "tolerance");    
    const double phi=check_numeric_scalar(disp, "dispersion");    
    const double pseudo=check_numeric_scalar(offset, "offset");    

    if (lim <= 0) {
        throw std::runtime_error("tolerance must be a positive double-precision value");
    }
    if (phi < 0) { 
        throw std::runtime_error("dispersion must be a non-negative double-precision value");
    }
    if (pseudo <= 0) { 
        throw std::runtime_error("offset must be a positive double-precision value");
    }

    // Defining quantile function and PDFs.
    FUN* ptr=NULL;
    pois P(lim);
    nbinom N(lim, phi);
    if (phi==0) {
        ptr=&P;                        
    } else {
        ptr=&N;
    }

    // Setting up the output.
    Rcpp::NumericVector outputm(Means.size()), outputv(Means.size());
    auto mIt=outputm.begin(), vIt=outputv.begin();

    for (const auto& m : Means) {
        // Computing the grand mean of the expectation of the log-counts for each gene.
        double& grand_mean=(*mIt);
        ++mIt;

        for (const auto& s : Sizes) {
            const double mu=m*s;
            const int left=ptr->q(mu, true);
            const int right=ptr->q(mu, false);
            double out_mean=0;
            double total_prob=0;

            for (auto i=left; i<=right; ++i) {
                const double cur_prob = ptr->d(i, mu);
                out_mean += cur_prob *std::log(i/s+pseudo);
                total_prob += cur_prob;
            }

            out_mean/=total_prob;
            grand_mean += out_mean;
        }
        grand_mean/=Sizes.size();

        // Computing the variance for each gene, as the expected (X - E[E[X]])^2 where E[E[X]] is the grand mean.
        double& grand_var=(*vIt);
        ++vIt;

        for (const auto& s : Sizes) {
            const double mu=m*s;
            const int left=ptr->q(mu, true);
            const int right=ptr->q(mu, false);
            double out_var=0, total_prob=0;

            for (auto i=left; i<=right; ++i) {
                const double cur_prob = ptr->d(i, mu);
                const double tmp=std::log(i/s+pseudo) - grand_mean;
                out_var+=cur_prob * tmp * tmp;
                total_prob += cur_prob;
            }

            out_var/=total_prob;
            grand_var += out_var;
        }
        grand_var/=Sizes.size();

        // Adjusting for the base of the log.
        grand_mean/=M_LN2;
        grand_var/=M_LN2 * M_LN2;
    }    

    return Rcpp::List::create(outputm, outputv);
    END_RCPP
}
