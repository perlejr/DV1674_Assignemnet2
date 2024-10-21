/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "analysis.hpp"
#include "dataset.hpp"
#include <iostream>
#include <cstdlib>

int main(int argc, char const* argv[])
{
    if ((argc != 4)  && (argc != 3) ) {
        std::cerr << "Usage: " << argv[0] << " [dataset] [outfile] [threads]" << std::endl;
        std::exit(1);
    }

    auto datasets { Dataset::read(argv[1]) };
    unsigned nThreads = 1;
    if (argv[3]) {
        nThreads = atoi(argv[3]);
        auto corrs { Analysis::correlation_coefficients_threads(datasets, nThreads) };
        Dataset::write(corrs, argv[2]);
    } else {
        auto corrs { Analysis::correlation_coefficients(datasets) };
        Dataset::write(corrs, argv[2]);
    }

    return 0;
}
