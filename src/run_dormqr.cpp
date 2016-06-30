#include "scran.h"

run_dormqr::run_dormqr(const int no, const int nc, const double* qrptr, const double* qrxptr, const char tr) :
        qr(qrptr), qrx(qrxptr), 
        nobs(no), ncoef(nc), ncol(1), 
        side('L'), trans(tr), 
        info(0), lwork(-1) {

    // Setting up the right hand side vector.
    rhs=(double*)R_alloc(nobs, sizeof(double));

    // Workspace query.            
    double tmpwork=0;
    F77_CALL(dormqr)(&side, &trans, &nobs, &ncol, &ncoef,
            qr, &nobs, qrx, rhs, &nobs,
            &tmpwork, &lwork, &info); 
    if (info) { 
        throw std::runtime_error("workspace query failed for 'dormqr'");
    }
    lwork=int(tmpwork+0.5);
    work=(double*)R_alloc(lwork, sizeof(double));
    return;
} 
       
void run_dormqr::run(double* stuff) {
    F77_CALL(dormqr)(&side, &trans, &nobs, &ncol, &ncoef,
            qr, &nobs, qrx, stuff, &nobs,
            work, &lwork, &info); 
    if (info) { 
        throw std::runtime_error("residual calculations failed for 'dormqr'");
    }
    return;
}

void run_dormqr::run() {
    run(rhs);
    return;
}


