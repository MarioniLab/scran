#include "scran.h"

SEXP build_snn(SEXP neighbors) try {
    matrix_info nn=check_matrix(neighbors);
    const size_t& k=nn.ncol;
    const size_t& ncells=nn.nrow;
    const int* nptr=nn.iptr;
    if (nptr==NULL) {
        throw std::runtime_error("neighbor matrix should be integer");
    }

    std::deque<std::deque<std::pair<size_t, int> > > hosts(ncells);
    size_t i, j;
    for (i=1; i<=k; ++i) {
        for (j=0; j<ncells; ++j) {
            hosts[*nptr - 1].push_back(std::make_pair(i, j)); // Getting to 0-based index, keeping 0-based ranks for now.
            ++nptr;
        }
    }

    std::deque<int> output_pairs;
    std::deque<double> output_weights;
    std::deque<int> current_added;
    std::deque<size_t> current_score(ncells);
    int cur_neighbor;
    size_t l, currank;

    nptr=nn.iptr - ncells;  
    for (j=0; j<ncells; ++j, ++nptr) {
        for (i=0; i<=k; ++i) {
            if (i==0) {
                cur_neighbor=j;
            } else {
                // Adding the actual nearest neighbors for cell 'j'.
                cur_neighbor=nptr[i*ncells] - 1; 
                currank=i;
                size_t& existing_other=current_score[cur_neighbor];
                if (existing_other==0) { 
                    existing_other=currank;
                    current_added.push_back(cur_neighbor);
                } else if (existing_other > currank) {
                    existing_other=currank;
                }
            }

            // Adding the cells connected by shared nearest neighbors.
            std::deque<std::pair<size_t, int> >& hosted=hosts[cur_neighbor];
            for (l=0; l<hosted.size(); ++l) {
                const int& othernode=hosted[l].second;
                currank=hosted[l].first + i;

                size_t& existing_other=current_score[othernode];
                if (existing_other==0) { 
                    existing_other=currank;
                    current_added.push_back(othernode);
                } else if (existing_other > currank) {
                    existing_other=currank;
                }
            }
        }
       
        for (l=0; l<current_added.size(); ++l) {
            const int& othernode=current_added[l];
            size_t& otherscore=current_score[othernode];

            // Converting to edges.
            output_pairs.push_back(j + 1);
            output_pairs.push_back(othernode + 1);
            output_weights.push_back(k - 0.5 * otherscore);

            // Resetting all those added to zero.
            otherscore=0;
        }
        current_added.clear();
    } 

    SEXP output=PROTECT(allocVector(VECSXP, 2));
    try {
        SET_VECTOR_ELT(output, 0, allocVector(INTSXP, output_pairs.size()));
        std::copy(output_pairs.begin(), output_pairs.end(), INTEGER(VECTOR_ELT(output, 0)));
        SET_VECTOR_ELT(output, 1, allocVector(REALSXP, output_weights.size()));
        std::copy(output_weights.begin(), output_weights.end(), REAL(VECTOR_ELT(output, 1)));
    } catch (std::exception& e) {
        UNPROTECT(1);
        throw;
    }
    UNPROTECT(1);   
    return output;
} catch (std::exception& e) {
    return mkString(e.what());
}
