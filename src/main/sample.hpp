#ifndef SAMPLE_HPP
#define SAMPLE_HPP

#include <vector>
#include "rand.hpp"

template<typename T>
void sampleVals(std::vector<T> * const sample,
		std::vector<T> const * const values,
		std::vector<double> * const probs,
		int const sampleSize);
// REQUIRES: *probs add up to 1; at least sampleSize probs are non-zero;
//           values->size() == probs->size(); sampleSize <= values->size()
// MODIFIES: *sample
// EFFECTS: calls sampleVals, with *probs having all equal elements


template<typename T>
void sampleVals(std::vector<T> * const sample,
		std::vector<T> const * const values,
		int const sampleSize);
// REQUIRES: at least sampleSize probs are non-zero;
//           values->size() == probs->size(); sampleSize <= values->size()

// MODIFIES: *sample
// EFFECTS: calls sampleVals, with *probs having all equal elements


template<typename T>
void shuffle(std::vector<T> * const vec);
// REQUIRES: nothing
// MODIFIES: *vec
// EFFECTS: shuffles the elements in *vec


#endif
