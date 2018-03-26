#include "scran.h"
#include "run_dormqr.h"

double rho_mult (double Ncells) {
    return 6/(Ncells*(Ncells*Ncells-1));
}

/*** Null distribution estimation without a design matrix. ***/

SEXP get_null_rho (SEXP cells, SEXP iters) {
    BEGIN_RCPP

    // Pulling out input values.
    const int Ncells=check_integer_scalar(cells, "number of cells");
    if (Ncells <= 1) { throw std::runtime_error("number of cells should be greater than 2"); }

    const int Niters=check_integer_scalar(iters, "number of iterations");
    if (Niters <= 0) { throw std::runtime_error("number of iterations should be positive"); }

    // Filling rank vector.
    std::vector<int> rankings(Ncells);
    std::iota(rankings.begin(), rankings.end(), 0);

    Rcpp::NumericVector output(Niters);
    Rcpp::RNGScope rng;
    const double mult=rho_mult(Ncells);

    for (int it=0; it<Niters; ++it) {
        Rx_shuffle(rankings.begin(), rankings.end());
        double tmp=0;
        for (int cell=0; cell<Ncells; ++cell) {
            const double tmpdiff=rankings[cell]-cell;
            tmp+=tmpdiff*tmpdiff;
        }
        tmp*=mult;
        output[it]=1-tmp;            
    }

    return output;
    END_RCPP
}

/*** Null distribution estimation with a design matrix. ***/

SEXP get_null_rho_design(SEXP qr, SEXP qraux, SEXP iters) {
    BEGIN_RCPP
    const int Niters=check_integer_scalar(iters, "number of iterations");
    if (Niters <= 0) { throw std::runtime_error("number of iterations should be positive"); }

    // Setting up to multiply by the Q matrix.
    run_dormqr multQ(qr, qraux, 'N');
    const int Nobs=multQ.get_nobs();
    const int Ncoef=multQ.get_ncoefs();

    Rcpp::NumericVector output(Niters);
    Rcpp::RNGScope rng;

    std::deque<std::pair<double, int> > collected1(Nobs), collected2(Nobs);
    std::deque<int> rank1(Nobs), rank2(Nobs);
    const double mult=rho_mult(Nobs);

    /* Simulating residuals, using the Q-matrix to do it. We set the main effects to zero 
     * (hence, starting from "Ncoefs") and simulate normals for the residual effects.
     * We then use this to reconstruct the residuals themselves - twice - and then compute 
     * correlations between the two reconstructions.
     */
    for (int it=0; it<Niters; ++it) {
        for (int mode=0; mode<2; ++mode) {
            std::fill(multQ.rhs.begin(), multQ.rhs.begin()+Ncoef, 0);
            for (int col=Ncoef; col<Nobs; ++col) {
                multQ.rhs[col]=norm_rand();
            }

            // Computing the residuals.
            multQ.run();

            // Sorting.
            std::deque<std::pair<double, int> >& current=(mode ? collected1 : collected2);
            for (int row=0; row<Nobs; ++row) {
                current[row].first=multQ.rhs[row];
                current[row].second=row;
            }
            std::sort(current.begin(), current.end());
            std::deque<int>& rank=(mode ? rank1 : rank2);
            for (int row=0; row<Nobs; ++row) {
                rank[current[row].second]=row;
            }
        }

        // Computing the squared difference in the ranks.
        double& rho=(output[it]=0);
        for (int row=0; row<Nobs; ++row) {
            const double tmpdiff=rank1[row]-rank2[row];
            rho+=tmpdiff*tmpdiff;
        }
        rho*=mult;
        rho=1-rho;
    }

    return output;    
    END_RCPP
}

/*** This function computes error-tolerant ranks for a subset of genes in a subset of cells. ***/

template <typename T, class V, class M>
SEXP rank_subset_internal(const M mat, SEXP intype, SEXP subset_row, SEXP subset_col, const T tol) {
    // Checking subset vectors.
    auto rsubout=check_subset_vector(subset_row, mat->get_nrow());
    const size_t rslen=rsubout.size();
    auto csubout=check_subset_vector(subset_col, mat->get_ncol());
    const size_t cslen=csubout.size();
    
    // Setting up the output matrix.
    const size_t ncells=mat->get_ncol();
    beachmat::output_param oparam(intype, true, false);
    oparam.set_chunk_dim(cslen, 1); // Column chunks (tranposed, so each column is a gene now).
    auto omat=beachmat::create_integer_output(cslen, rslen, oparam);
    if (!cslen) { 
        return omat->yield();
    }

    std::vector<int> indices(cslen);
    V incoming(ncells), subsetted(cslen);
    Rcpp::NumericVector breaker(cslen);
    Rcpp::IntegerVector ranks(cslen);
    Rcpp::RNGScope rng;

    auto rsIt=rsubout.begin();
    for (size_t rs=0; rs<rslen; ++rs, ++rsIt) {
        mat->get_row(*rsIt, incoming.begin());
        std::iota(indices.begin(), indices.end(), 0);
        auto sIt=subsetted.begin();
        for (auto csIt=csubout.begin(); csIt!=csubout.end(); ++csIt, ++sIt) {
            (*sIt)=incoming[*csIt];
        }

        // First stage sorting and equalization of effective ties.
        R_orderVector1(indices.data(), cslen, SEXP(subsetted), FALSE, FALSE);
        T last_unique=subsetted[indices.front()]; // Should be okay, we've removed cases where cslen=0.
        for (auto iIt=indices.begin(); iIt!=indices.end(); ++iIt) { 
            T& val=subsetted[*iIt];
            if (val - last_unique <= tol) {
                val=last_unique;
            } else {
                last_unique=val;
            }
        }

        // Second stage sorting with broken ties. This is done in two steps, as equalization needs to be done first.
        std::iota(indices.begin(), indices.end(), 0);
        for (auto bIt=breaker.begin(); bIt!=breaker.end(); ++bIt) { 
            (*bIt)=unif_rand();
        }
        R_orderVector(indices.data(), cslen, Rf_lang2(SEXP(subsetted), SEXP(breaker)), FALSE, FALSE);

        // Filling the output matrix.
        auto iIt=indices.begin();
        for (int cs=0; cs<cslen; ++cs, ++iIt){ 
            ranks[*iIt]=cs+1;
        }
        omat->set_col(rs, ranks.begin());
    }

    return omat->yield();
}

SEXP rank_subset(SEXP exprs, SEXP subset_row, SEXP subset_col, SEXP tol) {
    BEGIN_RCPP
    int rtype=beachmat::find_sexp_type(exprs);

    if (rtype==INTSXP) { 
        auto mat=beachmat::create_integer_matrix(exprs);
        const int tolerance=check_integer_scalar(tol, "tolerance");
        return rank_subset_internal<int, Rcpp::IntegerVector>(mat.get(), exprs, subset_row, subset_col, tolerance);

    } else {
        auto mat=beachmat::create_numeric_matrix(exprs);
        const double tolerance=check_numeric_scalar(tol, "tolerance");
        return rank_subset_internal<double, Rcpp::NumericVector>(mat.get(), exprs, subset_row, subset_col, tolerance);
    }
    END_RCPP
}

/*** Estimating correlations (without expanding into a matrix to do so via 'cor'). ***/

struct cache_manager {
public:
    cache_manager(Rcpp::RObject r, Rcpp::RObject cs) : rankings(beachmat::create_integer_matrix(r)), 
            cache_size(check_integer_scalar(cs, "block size")), ncells(rankings->get_nrow()),
            cache1(ncells, cache_size), cache2(ncells, cache_size), 
            current_cache1(0), current_cache2(0) {
        if (ncells <= 1) { 
            throw std::runtime_error("number of cells should be greater than 2"); 
        }
    }

    size_t get_ncells() const {
        return ncells;
    }
   
    Rcpp::IntegerVector::iterator get1(size_t index) {
	    const size_t actual=index % cache_size;
        const size_t cache_num=index / cache_size;
        if (cache_num!=current_cache1) {
            current_cache1=cache_num;
            std::fill(cache1.filled.begin(), cache1.filled.end(), 0);
        }

        auto& loc=cache1.location[actual];
	    auto& loaded=cache1.filled[actual];
        if (!loaded) {
            loc=rankings->get_const_col(index, cache1.store.begin() + actual * ncells);
            loaded=1;
        }
	    return loc;
    }

    Rcpp::IntegerVector::iterator get2(size_t index) {
	    const size_t actual=index % cache_size;
        const size_t cache_num=index / cache_size;

        // It is important to clear the cache2.filled before skipping to cache1, as we will not get another opportunity;
        // this will result in an incorrect assumption that cache2 is loaded, the next time we need to use it.
        if (cache_num!=current_cache2) {
            current_cache2=cache_num;
            std::fill(cache2.filled.begin(), cache2.filled.end(), 0);
        }

        // Checking if the first cache has the index, in which case we can avoid loading it in the second cache.
        // Note that this assumes that get2() is always called _after_ get1(), otherwise cache1 might be set by get2() 
        // but cleared and overwritten upon get1(), which would result in incorrect values.
        rank_cache& curcache=(cache_num==current_cache1 ? cache1 : cache2);

        auto& loc=curcache.location[actual];
        auto& loaded=curcache.filled[actual];
        if (!loaded) {
            loc=rankings->get_const_col(index, curcache.store.begin() + actual*ncells);
            loaded=1;
        }
	    return loc;
    }

    std::vector<size_t> reorder(Rcpp::IntegerVector g1, Rcpp::IntegerVector g2) {
        const size_t npairs=g1.size();
        if (npairs!=g2.size()) {
            throw std::runtime_error("gene index vectors should be of the same length");
        }

        std::vector<std::pair<size_t, size_t> > pairings;
        pairings.reserve(npairs);
        const size_t Ngenes=rankings->get_nrow();

        // Running through and checking the indices for validity.
        auto it1=g1.begin(), it2=g2.begin();
        for (size_t i=0; i<npairs; ++i, ++it1, ++it2) {
            const int g1x=*it1;
            const int g2x=*it2;
            if (g1x < 0 || g1x >= Ngenes) {
                throw std::runtime_error("first gene index is out of range");
            }
            if (g2x < 0 || g2x >= Ngenes) {
                throw std::runtime_error("second gene index is out of range");
            }
            pairings.push_back(std::pair<size_t, size_t>(g1x/cache_size, g2x/cache_size));
        }

        // Sorting in a manner so that all pairs of elements in the same cache pair are together.
        std::vector<size_t> output(npairs);
        std::iota(output.begin(), output.end(), 0);
        std::sort(output.begin(), output.end(), [&](const size_t left, const size_t right) {
            const size_t left1=pairings[left].first, right1=pairings[right].first;
            if (left1 < right1) { 
                return true; 
            } else if (left1==right1) {
                return pairings[left].second < pairings[right].second;
            }
            return false;
        });
        return output;
    }
private:
    struct rank_cache {
        rank_cache(const size_t ncells, const size_t ncached) : store(ncells*ncached), filled(ncached), location(ncached) {}
        Rcpp::LogicalVector filled;
        Rcpp::IntegerVector store;
        std::vector<Rcpp::IntegerVector::iterator> location;
    };

    std::unique_ptr<beachmat::integer_matrix> rankings;
    const size_t cache_size, ncells;
    rank_cache cache1, cache2;

    size_t current_cache1, current_cache2;
};

SEXP compute_rho(SEXP g1, SEXP g2, SEXP rankings, SEXP block_size) {
    BEGIN_RCPP

    Rcpp::IntegerVector gene1(g1), gene2(g2);
    cache_manager cache(rankings, block_size);
    auto order=cache.reorder(gene1, gene2);

    const size_t Ncells=cache.get_ncells();
    const double mult=rho_mult(Ncells); 
    Rcpp::NumericVector output(order.size());

    for (const auto& o : order) {
        // get1() must be called before get2(), see above.
        auto it1=cache.get1(gene1[o]);
        auto it2=cache.get2(gene2[o]);

        // Computing the correlation.
        double working=0;
        for (size_t c=0; c<Ncells; ++c, ++it1, ++it2) { 
            const double tmp=(*it1) - (*it2);
            working+=tmp*tmp;
        }

        output[o]=1 - working*mult;
    }

    return output;
    END_RCPP
}

/*** Combining correlated p-values for each gene into a single combined p-value. ***/

SEXP combine_corP (SEXP ng, SEXP g1, SEXP g2, SEXP rho, SEXP pval, SEXP limited, SEXP order) {
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

