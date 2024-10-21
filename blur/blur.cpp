/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "matrix.hpp"
#include "ppm.hpp"
#include "filters.hpp"
#include <cstdlib>
#include <iostream>

int main(int argc, char const* argv[])
{
    if ((argc != 5)  && (argc != 4) ) {
        std::cerr << "Usage: " << argv[0] << " [radius] [infile] [outfile] [threads]" << std::endl;
        std::exit(1);
    }

    PPM::Reader reader {};
    PPM::Writer writer {};

    auto m { reader(argv[2]) };
    auto radius { static_cast<unsigned>(std::stoul(argv[1])) };
    
    if (argv[4]) {
        unsigned int nthreads = atoi(argv[4]);
        auto blurred { Filter::blur_threads(m, radius, nthreads) };
        writer(blurred, argv[3]);
    }
    else{
        auto blurred { Filter::blur(m, radius) };
        writer(blurred, argv[3]);
    }

    
    return 0;
}
