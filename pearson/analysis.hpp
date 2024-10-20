/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "vector.hpp"
#include <vector>

#if !defined(ANALYSIS_HPP)
#define ANALYSIS_HPP

namespace Analysis {
std::vector<double> correlation_coefficients(std::vector<Vector> datasets);
std::vector<double> correlation_coefficients(std::vector<Vector> datasets, unsigned nthreads);
double pearson(Vector vec1, Vector vec2);
double pearson_threads(Vector vec1, Vector vec2, unsigned nthreads);
};

#endif