#include "scran.h"

/*** Combining correlated p-values for each gene into a single combined p-value. ***/

SEXP combine_rho (SEXP ng, SEXP g1, SEXP g2, SEXP rho, SEXP pval, SEXP limited, SEXP order) {
    BEGIN_RCPP

    // Checking inputs.
    const int Ngenes=check_integer_scalar(ng, "number of genes");
    if (Ngenes < 0) { throw std::runtime_error("number of genes should be non-zero"); }

    Rcpp::IntegerVector first(g1), second(g2);
    const size_t Npairs=first.size();
    if (Npairs!=second.size()) { 
        throw std::runtime_error("gene index vectors must be of the same length"); 
    }

    Rcpp::NumericVector Rho(rho);
    if (Rho.size()!=Npairs) { 
        throw std::runtime_error("'rho' must be a double precision vector of length equal to the number of pairs");
    }

    Rcpp::NumericVector Pval(pval);
    if (Pval.size()!=Npairs) {
        throw std::runtime_error("'pval' must be a double precision vector of length equal to the number of pairs");
    }

    Rcpp::LogicalVector Limited(limited);
    if (Limited.size()!=Npairs) { 
        throw std::runtime_error("'limited' must be a logical vector of length equal to the number of pairs");
    }

    Rcpp::IntegerVector Order(order);
    if (Order.size()!=Npairs) { 
        throw std::runtime_error("'order' must be an integer vector of length equal to the number of pairs");
    }
   
    // Going through and computing the combined p-value for each gene. 
    Rcpp::NumericVector pout(Ngenes), rout(Ngenes);
    Rcpp::LogicalVector lout(Ngenes);
    std::vector<int> sofar(Ngenes);

    for (auto oIt=Order.begin(); oIt!=Order.end(); ++oIt) {
        const int& curp=*oIt;
        if (curp < 0 || curp >= int(Npairs)) { 
            throw std::runtime_error("order indices out of range");
        }
        const double& currho=Rho[curp];
        const double& curpval=Pval[curp];
        const int& curlimit=Limited[curp];

        for (int i=0; i<2; ++i) {
            const int& gx=(i==0 ? first[curp] : second[curp]);
            if (gx < 0 || gx >= Ngenes) {
                throw std::runtime_error("supplied gene index is out of range");
            }

            // Checking if this is smaller than what is there, or if nothing is there yet.
            int& already_there=sofar[gx];
            ++already_there;
            const double temp_combined=curpval/already_there;
            double& combined_pval=pout[gx];
            if (already_there==1 || temp_combined < combined_pval) {
                combined_pval=temp_combined;
                lout[gx]=curlimit; // is the combined p-value computed from a limited p-value?
            }
            double& max_rho=rout[gx];
            if (already_there==1 || std::abs(max_rho) < std::abs(currho)) {
                max_rho=currho;
            }
        }
    }
    
    // Multiplying by the total number of tests for each gene.       
    auto sfIt=sofar.begin();
    for (auto poIt=pout.begin(); poIt!=pout.end(); ++poIt, ++sfIt) { 
        (*poIt)*=(*sfIt); 
    }

    return Rcpp::List::create(pout, rout, lout);
    END_RCPP
}

