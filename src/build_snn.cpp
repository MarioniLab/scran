#include "scran.h"

SEXP build_snn(SEXP neighbors) {
    BEGIN_RCPP
    auto mat=beachmat::create_integer_matrix(neighbors);
    const size_t& k=mat->get_ncol();
    const size_t& ncells=mat->get_nrow();

    // Building a host table, identifying the reverse relation from nearest neighbours to cells.
    Rcpp::IntegerVector tmp(ncells);
    std::deque<std::deque<std::pair<size_t, int> > > hosts(ncells);
    for (size_t i=1; i<=k; ++i) {
        mat->get_col(i-1, tmp.begin()); // as the first column (index=0) is the first neighbour (i=1).
        
        auto tIt=tmp.begin();
        for (size_t j=0; j<ncells; ++j, ++tIt) {
            hosts[*tIt - 1].push_back(std::make_pair(i, j)); // Getting to 0-based index, keeping 0-based ranks for now.
        }
    }

    std::deque<int> output_pairs;
    std::deque<double> output_weights;
    std::deque<size_t> current_added;
    std::deque<size_t> current_score(ncells);

    Rcpp::IntegerVector rowtmp(k);
    for (size_t j=0; j<ncells; ++j) {
        mat->get_row(j, rowtmp.begin());
        auto rtIt=rowtmp.begin();
            
        int cur_neighbor;
        for (size_t i=0; i<=k; ++i) {
            if (i==0) {
                cur_neighbor=j;
            } else {
                // Adding the actual nearest neighbors for cell 'j'.
                cur_neighbor=*rtIt - 1; 
                ++rtIt;

                size_t currank=i; // +0, as neighbour 'i' is rank 0 with respect to itself.
                size_t& existing_other=current_score[cur_neighbor];
                if (existing_other==0) { 
                    existing_other=currank;
                    current_added.push_back(cur_neighbor);
                } else if (existing_other > currank) {
                    existing_other=currank;
                }
            }

            // Adding the cells connected by shared nearest neighbors, again recording the lowest combined rank per neighbor.
            std::deque<std::pair<size_t, int> >& hosted=hosts[cur_neighbor];
            for (auto hIt=hosted.begin(); hIt!=hosted.end(); ++hIt) { 
                const int& othernode=hIt->second;
                size_t currank=hIt->first + i;

                size_t& existing_other=current_score[othernode];
                if (existing_other==0) { 
                    existing_other=currank;
                    current_added.push_back(othernode);
                } else if (existing_other > currank) {
                    existing_other=currank;
                }
            }
        }
       
        for (auto caIt=current_added.begin(); caIt!=current_added.end(); ++caIt) {
            const int& othernode=*caIt;
            size_t& otherscore=current_score[othernode];

            // Converting to edges.
            output_pairs.push_back(j + 1);
            output_pairs.push_back(othernode + 1);
            output_weights.push_back(double(k) - 0.5 * double(otherscore));

            // Resetting all those added to zero.
            otherscore=0;
        }
        current_added.clear();
    }

    Rcpp::IntegerVector pout(output_pairs.begin(), output_pairs.end());
    Rcpp::NumericVector wout(output_weights.begin(), output_weights.end());
    return Rcpp::List::create(pout, wout);
    END_RCPP
}
