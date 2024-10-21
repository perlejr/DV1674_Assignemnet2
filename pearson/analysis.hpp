/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "vector.hpp"
#include <vector>

#if !defined(ANALYSIS_HPP)
#define ANALYSIS_HPP

namespace Analysis {
std::vector<double> correlation_coefficients_threads(std::vector<Vector> datasets, unsigned nthreads);
std::vector<double> correlation_coefficients(std::vector<Vector> datasets);
void* correlation_coefficients_void(void* params);
double pearson(Vector vec1, Vector vec2);
};

#endif