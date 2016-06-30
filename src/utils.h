#ifndef UTILS_H
#define UTILS_H

struct Rx_random_seed {
    Rx_random_seed() {
        GetRNGstate();
    }
    ~Rx_random_seed() {
        PutRNGstate();
    }
};

template <typename Pt>
void Rx_shuffle (Pt start, Pt end) {
    size_t n=(end-start), i;
    for (i=n-1; i>0; --i) {
        std::swap(start[i], start[int(unif_rand()*(i+1))]); // Equal chance of getting anything from [0, i], as per random_shuffle in C++.
    }
    return;
}

// Check matrix inputs.

struct matrix_info {
    matrix_info(int, int, bool);
    const size_t nrow, ncol;
    const bool is_integer;
    const int* iptr;
    const double* dptr;
};

matrix_info check_matrix(SEXP matrix);

// Check subset vector.

typedef std::pair<const int, const int*> subset_values;
subset_values check_subset_vector(SEXP, int);

// Special function to check for NA'ness.

bool isNA(int);
bool isNA(double);

// Class to run Qy or QtY multiplications.

struct run_dormqr {
    run_dormqr(const int, const int, const double*, const double*, const char);
    void run();
    void run(double*);
    const double* qr, *qrx;
    const int nobs, ncoef, ncol;
    const char side, trans;           
    int info, lwork;
    double* work;
    double* rhs;
};

#endif
