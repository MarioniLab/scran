#include "scran.h"

static double rho_mult (double Ncells) {
    return 6/(Ncells*(Ncells*Ncells-1));
}

/* Estimating correlations (without expanding into a matrix to do so via 'cor'). 
 * The cache manager loads blocks of genes at a time to avoid having to do multiple 
 * redundant read/writes for each gene. Some redundant reading/writing is still 
 * necessary but that is inevitable if the data cannot fit into memory.
 */

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
   
    Rcpp::IntegerVector::const_iterator get1(size_t index) {
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

    Rcpp::IntegerVector::const_iterator get2(size_t index) {
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
        const size_t Ngenes=rankings->get_ncol();

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
        Rcpp::IntegerVector store;
        Rcpp::LogicalVector filled;
        std::vector<Rcpp::IntegerVector::const_iterator> location;
    };

    std::unique_ptr<beachmat::integer_matrix> rankings;
    const size_t cache_size, ncells;
    rank_cache cache1, cache2;

    size_t current_cache1, current_cache2;
};

SEXP compute_rho_pairs(SEXP g1, SEXP g2, SEXP rankings, SEXP block_size) {
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
