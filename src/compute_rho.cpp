#include "scran.h"
#include "run_dormqr.h"

double rho_mult (double Ncells) {
    return 6/(Ncells*(Ncells*Ncells-1));
}

/*** Null distribution estimation without a design matrix. ***/

SEXP get_null_rho (SEXP cells, SEXP iters) {
    BEGIN_RCPP

    // Pulling out input values.
    Rcpp::IntegerVector NC(cells);
    if (NC.size()!=1) { 
        throw std::runtime_error("number of cells should be an integer scalar"); 
    }
    const int Ncells=NC[0];
    if (Ncells <= 1) { throw std::runtime_error("number of cells should be greater than 2"); }

    Rcpp::IntegerVector Iters(iters);
    if (Iters.size()!=1)  {
        throw std::runtime_error("number of iterations should be an integer scalar"); 
    }
    const int Niters=Iters[0];
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
    Rcpp::IntegerVector Iters(iters);
    if (Iters.size()!=1)  {
        throw std::runtime_error("number of iterations should be an integer scalar"); 
    }
    const int Niters=Iters[0];
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
        Rcpp::IntegerVector tolerance(tol);
        if (tolerance.size()!=1) { 
            throw std::runtime_error("tolerance should be an integer scalar");
        }
        return rank_subset_internal<int, Rcpp::IntegerVector>(mat.get(), exprs, subset_row, subset_col, tolerance[0]);

    } else {
        auto mat=beachmat::create_numeric_matrix(exprs);
        Rcpp::NumericVector tolerance(tol);
        if (tolerance.size()!=1) { 
            throw std::runtime_error("tolerance should be an double-precision scalar");
        }
        return rank_subset_internal<double, Rcpp::NumericVector>(mat.get(), exprs, subset_row, subset_col, tolerance[0]);
    }
    END_RCPP
}

/*** Estimating correlations (without expanding into a matrix to do so via 'cor'). ***/

SEXP compute_rho(SEXP g1, SEXP g2, SEXP rankings, SEXP block_size) {
    BEGIN_RCPP

    // Checking inputs.
    auto rmat=beachmat::create_integer_matrix(rankings);
    const size_t& Ncells=rmat->get_nrow();
    if (Ncells <= 1) { throw std::runtime_error("number of cells should be greater than 2"); }
    const size_t& Ngenes=rmat->get_ncol();

    Rcpp::IntegerVector first(g1), second(g2);
    const size_t Npairs=first.size();
    if (Npairs!=second.size()) { 
        throw std::runtime_error("gene index vectors must be of the same length"); 
    }
    
    Rcpp::IntegerVector bsize(block_size);
    if (bsize.size()!=1) { 
        throw std::runtime_error("block size should be an integer scalar");
    }
    const int BLOCK=bsize[0];
    
    /* Setting up the cache, to avoid repeatedly copying/reading from file with rmat->get_col().
     * We round up the number of genes to a multiple of BLOCK to make things easier when using
     * std::fill to indicate that particular genes are no longer in memory.
     */
    Rcpp::IntegerVector cache1(BLOCK*Ncells), cache2(BLOCK*Ncells);
    const int roundedNgenes=BLOCK*int(Ngenes/BLOCK + 1); 
    std::deque<bool> filled1(roundedNgenes, false), filled2(roundedNgenes, false);
    std::vector<Rcpp::IntegerVector::const_iterator> locations1(roundedNgenes), locations2(roundedNgenes); 
    Rcpp::IntegerVector::iterator b1It=cache1.begin(), b2It=cache2.begin();
    int current_block1=0, current_block2=0;

    const double mult=rho_mult(Ncells); 
    Rcpp::NumericVector output(Npairs);

    auto fIt=first.begin(), sIt=second.begin();
    for (auto oIt=output.begin(); oIt!=output.end(); ++oIt, ++fIt, ++sIt) { 
        const int& g1x=(*fIt);
        const int& g2x=(*sIt);
        if (g1x < 0 || g1x >= Ngenes) {
            throw std::runtime_error("first gene index is out of range");
        }
        if (g2x < 0 || g2x >= Ngenes) {
            throw std::runtime_error("second gene index is out of range");
        }

        // Figuring out whether we need to update the blocks stored in the cache, starting with the first block.
        const int block1=int(g1x/BLOCK);
        const bool updated1=(block1!=current_block1);
        if (updated1) {
            if (block1 < current_block1) { 
                throw std::runtime_error("pairs should be arranged in increasing block order");
            }
            auto old_filled=filled1.begin()+current_block1*BLOCK;
            std::fill(old_filled, old_filled+BLOCK, false); 
            current_block1=block1;
            b1It=cache1.begin();          
        }
        auto& start1=locations1[g1x]; 
        if (!filled1[g1x]) {
            if (b1It==cache1.end()) { 
                throw std::runtime_error("first block cache exceeded");
            }
            start1=rmat->get_const_col(g1x, b1It);
            b1It+=Ncells;
            filled1[g1x]=true;
        }

        // Repeating for the second block. If this is the same as the first block, we just update that instead.
        const int block2=int(g2x/BLOCK);
        Rcpp::IntegerVector::const_iterator start2_copy;
        if (block2==block1) { 
            if (updated1) { 
                /* No need to discard the existing cache1, this would have already been done above. 
                 * We do, however, have to update current_block2 and discard cache2 (we can't do this
                 * later as old_start would no longer be correct with an updated current_block2).
                 */
                auto old_start=filled2.begin()+current_block2*BLOCK;
                std::fill(old_start, old_start+BLOCK, false);
                current_block2=block2;
                b2It=cache2.begin();
            }
            auto& restart1=locations1[g2x]; 
            if (!filled1[g2x]) {
                if (b1It==cache1.end()) { 
                    throw std::runtime_error("first block cache exceeded");
                }
                restart1=rmat->get_const_col(g2x, b1It);
                b1It+=Ncells;
                filled1[g2x]=true;
            }
            start2_copy=restart1;
        } else {
            if (block2!=current_block2) {
                if (!updated1 && block2 < current_block2) { 
                    throw std::runtime_error("pairs should be arranged in increasing block order");
                }
                auto old_start=filled2.begin()+current_block2*BLOCK;
                std::fill(old_start, old_start+BLOCK, false);
                current_block2=block2;
                b2It=cache2.begin();          
            } 
            auto& start2=locations2[g2x]; 
            if (!filled2[g2x]) { 
                if (b2It==cache2.end()) { 
                    throw std::runtime_error("second block cache exceeded");
                }
                start2=rmat->get_const_col(g2x, b2It);
                b2It+=Ncells;
                filled2[g2x]=true;
            }
            start2_copy=start2;
        }

        // Computing the correlation.
        double working=0;
        auto start1_copy=start1;
        for (size_t c=0; c<Ncells; ++c, ++start1_copy, ++start2_copy) { 
            const double tmp=(*start1_copy) - (*start2_copy);
            working+=tmp*tmp;
        }
        (*oIt)=1 - working*mult;
    }

    return output;
    END_RCPP
}

/*** Combining correlated p-values for each gene into a single combined p-value. ***/

SEXP combine_corP (SEXP ng, SEXP g1, SEXP g2, SEXP rho, SEXP pval, SEXP limited, SEXP order) {
    BEGIN_RCPP

    // Checking inputs.
    Rcpp::IntegerVector NG(ng);
    if (NG.size()!=1) { 
        throw std::runtime_error("number of genes must be an integer scalar");
    }
    const int Ngenes=NG[0];
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

