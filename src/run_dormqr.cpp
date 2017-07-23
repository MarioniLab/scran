#include "run_dormqr.h"
#include "R_ext/BLAS.h"
#include "R_ext/Lapack.h"

run_dormqr::run_dormqr(SEXP qr, SEXP qraux, const char tr) :
        QR(qr), AUX(qraux), qrptr(NULL), qxptr(NULL), 
        nobs(QR.nrow()), ncoef(QR.ncol()), ncol(1), 
        side('L'), trans(tr), info(0), lwork(-1) {

    if (AUX.size()!=ncoef) { 
        throw std::runtime_error("QR auxiliary vector should be of length 'ncol(Q)'"); 
    }
    if (QR.size()) { 
        qrptr=&(QR[0]);
    }
    if (AUX.size()) { 
        qxptr=&(AUX[0]);
    }
    rhs.resize(nobs);
    
    // Workspace query.            
    double tmpwork=0;
    F77_CALL(dormqr)(&side, &trans, &nobs, &ncol, &ncoef,
            qrptr, &nobs, qxptr, rhs.data(), &nobs,
            &tmpwork, &lwork, &info); 
    if (info) { 
        throw std::runtime_error("workspace query failed for 'dormqr'");
    }
    lwork=int(tmpwork+0.5);
    work.resize(lwork);
    return;
} 
       
void run_dormqr::run(double* stuff) {
    F77_CALL(dormqr)(&side, &trans, &nobs, &ncol, &ncoef,
            qrptr, &nobs, qxptr, stuff, &nobs,
            work.data(), &lwork, &info); 
    if (info) { 
        throw std::runtime_error("residual calculations failed for 'dormqr'");
    }
    return;
}

void run_dormqr::run() {
    run(rhs.data());
    return;
}

void run_dormqr::solve(double* stuff) {
    const char uplo='U', xtrans='N', diag='N';
    F77_CALL(dtrtrs)(&uplo, &xtrans, &diag,
            &ncoef, &ncol, qrptr, &nobs, 
            stuff, &nobs, &info);
    if (info) { 
        throw std::runtime_error("coefficient calculations failed for 'dtrtrs'");
    }
    return;
}

void run_dormqr::solve() {
    solve(rhs.data());
}

int run_dormqr::get_nobs() const {
    return nobs;
}

int run_dormqr::get_ncoefs() const { 
    return ncoef;
}

