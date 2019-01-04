#ifndef SHUFFLE_CUSTOM_H
#define SHUFFLE_CUSTOM_H
#include "boost/random.hpp"
#include <iterator>

// A possible implementation of std::shuffle, hard-coded to avoid cross-platform inconsistencies. 
// Taken from https://en.cppreference.com/w/cpp/algorithm/random_shuffle#Possible_implementation.

template <class RandomIt, class RNG>
void shuffle_custom (RandomIt first, RandomIt last, RNG&& g) {
    typedef typename std::iterator_traits<RandomIt>::difference_type diff_t;
    typedef boost::random::uniform_int_distribution<diff_t> distr_t;
    typedef typename distr_t::param_type param_t;
 
    distr_t D;
    diff_t n = last - first;
    for (diff_t i = n-1; i > 0; --i) {
        using std::swap;
        swap(first[i], first[D(g, param_t(0, i))]);
    }
    return;
}

#endif
