/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "filters.hpp"
#include "matrix.hpp"
#include "ppm.hpp"
#include <cmath>
#include <list>
#include <vector>
#include <pthread.h>

namespace Filter
{
    Matrix dst = Matrix();


    struct arguments
    {
        //Matrix local_dst;
        Matrix local_m;
        int local_radius;
        unsigned thread_nr;
        unsigned local_dst_x_size;
        unsigned local_dst_y_size;
        unsigned from_index;
    };

    namespace Gauss
    {
        void get_weights(int n, double *weights_out)
        {
            for (auto i{0}; i <= n; i++)
            {
                double x{static_cast<double>(i) * max_x / n};
                weights_out[i] = exp(-x * x * pi);
            }
        }
    }

    Matrix blur_threads(Matrix m, const int radius, unsigned int nthreads)
    {
        dst = m;
        std::cout << "M x size: " << m.get_x_size() << "\n";
        std::cout << "DST x size: " << dst.get_x_size() << "\n";
        std::cout << "DST y size: " << dst.get_y_size() << "\n";
        unsigned dst_x_size_per_thread = (dst.get_x_size() / nthreads);
        unsigned dst_y_size = dst.get_y_size();
        std::vector<pthread_t> threads {nthreads};
        std::vector<arguments*> thread_args;
        for(unsigned int i = 0; i < nthreads; ++i){
            arguments* new_arg = new arguments();
            new_arg->thread_nr = i;
            //new_arg->local_dst = dst;
            new_arg->from_index = i * dst_x_size_per_thread;
            new_arg->local_dst_x_size = dst_x_size_per_thread;
            new_arg->local_dst_y_size = dst_y_size;
            new_arg->local_radius = radius;

            thread_args.push_back(new_arg);
            pthread_create(&threads[i], NULL,  blur, (void*) thread_args[i]);
        }
        for (int i = 0; i < nthreads; i++) {
            pthread_join(threads[i], NULL);
        }
        for(auto elem : thread_args){
            delete elem;
        }

        return dst;
    }


    void* blur(void* params)
    {
        arguments* args = reinterpret_cast<arguments*>(params);
            std::cout << "In thread: " << args->thread_nr <<"\n";

        Matrix scratch{PPM::max_dimension};
        
        for (auto x = args->from_index; x < (args->local_dst_x_size + args->from_index); x++)
        {
            std::cout << "Thread " << args->thread_nr << " || first iteraction: " << x <<"\n";
            for (auto y{0}; y < args->local_dst_y_size; y++)
            {
                double w[Gauss::max_radius]{};
                Gauss::get_weights(args->local_radius, w);

                // unsigned char Matrix::r(unsigned x, unsigned y) const
                // {
                //     return R[y * x_size + x];
                // }

                auto r{w[0] * dst.r(x, y)}, g{w[0] * dst.g(x, y)}, b{w[0] * dst.b(x, y)}, n{w[0]};
                
                // for (auto wi{1}; wi <= radius; wi++)    
                for (int wi{1}; wi <= args->local_radius; wi++)
                {
                    auto wc{w[wi]};
                    unsigned int x2{x - wi};
                    //if (x2 >= 0)
                    std::cout << "problem wi: " << wi << "\n";
                    std::cout << "problem x: " << x << "\n";
                    std::cout << "problem x2: " << x2 << "\n";
                    if (x2 > args->local_dst_x_size && x2 < (args->local_dst_x_size + args->from_index))
                    {
                        r += wc * dst.r(x2, y);
                        g += wc * dst.g(x2, y);
                        b += wc * dst.b(x2, y);
                        n += wc;
                    }
                    x2 = x + wi;
                    std::cout << "no problem x2: " << x2 << "\n";
                    if (x2 < (args->local_dst_x_size + args->from_index))
                    //if (x2 < dst.get_x_size())
                    {
                        r += wc * dst.r(x2, y);
                        g += wc * dst.g(x2, y);
                        b += wc * dst.b(x2, y);
                        n += wc;
                    }
                }
                scratch.r(x, y) = r / n;
                scratch.g(x, y) = g / n;
                scratch.b(x, y) = b / n;
            }
        }

        for (auto x{args->from_index}; x < (args->local_dst_x_size + args->from_index); x++)
        {
            std::cout << "Thread " << args->thread_nr << " || second iteraction: " << x <<"\n";
            for (auto y{0}; y < args->local_dst_y_size; y++)
            {
                double w[Gauss::max_radius]{};
                Gauss::get_weights(args->local_radius, w);

                auto r{w[0] * scratch.r(x, y)}, g{w[0] * scratch.g(x, y)}, b{w[0] * scratch.b(x, y)}, n{w[0]};

                for (auto wi{1}; wi <= args->local_radius; wi++)
                {
                    auto wc{w[wi]};
                    auto y2{y - wi};
                    if (y2 >= 0)
                    {
                        r += wc * scratch.r(x, y2);
                        g += wc * scratch.g(x, y2);
                        b += wc * scratch.b(x, y2);
                        n += wc;
                    }
                    y2 = y + wi;
                    if (y2 < args->local_dst_y_size)
                    {
                        r += wc * scratch.r(x, y2);
                        g += wc * scratch.g(x, y2);
                        b += wc * scratch.b(x, y2);
                        n += wc;
                    }
                }
                dst.r(x, y) = r / n;
                dst.g(x, y) = g / n;
                dst.b(x, y) = b / n;
            }
        }
    }


    Matrix blur(Matrix m, const int radius)
    {
        Matrix scratch{PPM::max_dimension};
        auto dst{m};
        unsigned dst_x_size = dst.get_x_size();
        unsigned dst_y_size = dst.get_y_size();
        for (auto x{0}; x < dst_x_size; x++)
        {
            for (auto y{0}; y < dst_y_size; y++)
            {
                double w[Gauss::max_radius]{};
                Gauss::get_weights(radius, w);

                // unsigned char Matrix::r(unsigned x, unsigned y) const
                // {
                //     return R[y * x_size + x];
                // }

                auto r{w[0] * dst.r(x, y)}, g{w[0] * dst.g(x, y)}, b{w[0] * dst.b(x, y)}, n{w[0]};

                for (auto wi{1}; wi <= radius; wi++)
                {
                    auto wc{w[wi]};
                    auto x2{x - wi};
                    if (x2 >= 0)
                    {
                        r += wc * dst.r(x2, y);
                        g += wc * dst.g(x2, y);
                        b += wc * dst.b(x2, y);
                        n += wc;
                    }
                    x2 = x + wi;
                    if (x2 < dst_x_size)
                    {
                        r += wc * dst.r(x2, y);
                        g += wc * dst.g(x2, y);
                        b += wc * dst.b(x2, y);
                        n += wc;
                    }
                }
                scratch.r(x, y) = r / n;
                scratch.g(x, y) = g / n;
                scratch.b(x, y) = b / n;
            }
        }

        for (auto x{0}; x < dst_x_size; x++)
        {
            for (auto y{0}; y < dst_y_size; y++)
            {
                double w[Gauss::max_radius]{};
                Gauss::get_weights(radius, w);

                auto r{w[0] * scratch.r(x, y)}, g{w[0] * scratch.g(x, y)}, b{w[0] * scratch.b(x, y)}, n{w[0]};

                for (auto wi{1}; wi <= radius; wi++)
                {
                    auto wc{w[wi]};
                    auto y2{y - wi};
                    if (y2 >= 0)
                    {
                        r += wc * scratch.r(x, y2);
                        g += wc * scratch.g(x, y2);
                        b += wc * scratch.b(x, y2);
                        n += wc;
                    }
                    y2 = y + wi;
                    if (y2 < dst_y_size)
                    {
                        r += wc * scratch.r(x, y2);
                        g += wc * scratch.g(x, y2);
                        b += wc * scratch.b(x, y2);
                        n += wc;
                    }
                }
                dst.r(x, y) = r / n;
                dst.g(x, y) = g / n;
                dst.b(x, y) = b / n;
            }
        }

        return dst;
    }

}
