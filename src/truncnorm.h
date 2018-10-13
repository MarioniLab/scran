#ifndef TRUNCNORM_H
#define TRUNCNORM_H

struct truncnorm {
    truncnorm(size_t, double, double, double, double);
    double sum, sqsum, left, right;
    size_t n;
    
    void compute(double, double);
    double loglik, dmu, dsig, ddmusig, ddmu2, ddsig2;
};

std::tuple<double, double, bool> fit_truncnorm(truncnorm&, double=0.00000001, size_t=100);

std::tuple<double, double, bool> fit_truncnorm(size_t, double, double, double, double, double=0.00000001, size_t=100);

#endif
