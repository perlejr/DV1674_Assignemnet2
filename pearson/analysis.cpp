/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "analysis.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <list>
#include <vector>
#include <pthread.h>


namespace Analysis {
std::vector<double> result = std::vector<double>();
struct arguments
{
    unsigned thread_nr;
    unsigned block_size;
    unsigned from_index;
    std::vector<Vector> datasets;
    std::vector<double> local_result;
};

std::vector<double> correlation_coefficients_threads(std::vector<Vector> datasets, unsigned nthreads)
{
    unsigned block = datasets.size() / nthreads;
    pthread_t threads[nthreads];   
    std::vector<arguments*> threads_args;
    for (unsigned i = 0; i < nthreads; ++i){
        arguments* new_arg = new arguments();
        new_arg->thread_nr = i;
        new_arg->local_result = std::vector<double>();
        new_arg->block_size = block;
        new_arg->datasets = datasets;
        new_arg->from_index = i * block;
        threads_args.push_back(new_arg);
        pthread_create(&threads[i], NULL,  correlation_coefficients_void, (void*) threads_args[i]);
    }
    for (int i = 0; i < nthreads; i++) {
        pthread_join(threads[i], NULL);
        result.insert(result.cend(), threads_args[i]->local_result.begin(), threads_args[i]->local_result.end());
    }
    for(auto elem : threads_args){
        delete elem;
    }
    return result;
}

void* correlation_coefficients_void(void* params)
{
    arguments* args = reinterpret_cast<arguments*>(params);

    for (unsigned sample1 = args->from_index; sample1 < (args->from_index + args->block_size); sample1++) {
        for (auto sample2 = sample1 + 1; sample2 < args->datasets.size(); sample2++) {
            auto corr { pearson(args->datasets[sample1], args->datasets[sample2]) };
            args->local_result.push_back(corr);
        }
    }
}

std::vector<double> correlation_coefficients(std::vector<Vector> datasets)
{
    for (auto sample1 { 0 }; sample1 < datasets.size() - 1; sample1++) {
        for (auto sample2 { sample1 + 1 }; sample2 < datasets.size(); sample2++) {
            auto corr { pearson(datasets[sample1], datasets[sample2]) };
            result.push_back(corr);
            
        }
    }

    return result;
}

double pearson(Vector vec1, Vector vec2)
{
    auto x_mean { vec1.mean() };
    auto y_mean { vec2.mean() };

    auto x_mm { vec1 - x_mean };
    auto y_mm { vec2 - y_mean };

    auto x_mag { x_mm.magnitude() };
    auto y_mag { y_mm.magnitude() };

    auto x_mm_over_x_mag { x_mm / x_mag };
    auto y_mm_over_y_mag { y_mm / y_mag };

    auto r { x_mm_over_x_mag.dot(y_mm_over_y_mag) };

    return std::max(std::min(r, 1.0), -1.0);
}
};
