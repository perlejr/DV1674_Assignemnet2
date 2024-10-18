/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "vector.hpp"
#include <iostream>
#include <cmath>
// #include <vector>
#include <cstring>

Vector::Vector()
    : size{0}, data{nullptr}
{
}

Vector::~Vector()
{
    if (data)
    {
        delete[] data;
    }

    size = 0;
}

Vector::Vector(unsigned size)
    : size{size}, data{new double[size]}
{
}

Vector::Vector(unsigned size, double *data)
    : size{size}, data{data}
{
}

Vector::Vector(const Vector &other)
    : Vector{other.size}
{   
    size = other.size;
    std::memcpy(data, other.data, size);
}

unsigned Vector::get_size() const
{
    return size;
}

double *Vector::get_data()
{
    return data;
}

double Vector::operator[](unsigned i) const
{
    return data[i];
}

double &Vector::operator[](unsigned i)
{
    return data[i];
}

double Vector::mean() const
{
    //  double sum{0};

    // for (auto i{0}; i < size; i++)
    // {
    //     sum += data[i];
    // }

    // return sum / static_cast<double>(size);
    double sum{0};
    unsigned t1 = 0, t2 = 0, t3 = 0, t4 = 0;
    unsigned mover = 4;
    for (auto i{0}; i < size; i += mover)
    {
        if ((i + mover) <= size) {
            t1 += data[i];
            t2 += data[i + 1];
            t3 += data[i + 2];
            t4 += data[i + 3];
        }
        else{
            for (auto z{i}; z < size; ++z){
                t1 += data[z];
            }
            break;
        }
    }

    sum = t1 + t2 + t3 + t4;

    return sum / static_cast<double>(size);
}

double Vector::magnitude() const
{
    auto dot_prod{dot(*this)};
    return std::sqrt(dot_prod);
}

Vector Vector::operator/(double div)
{
    auto result{*this};

    for (auto i{0}; i < size; i++)
    {
        result[i] /= div;
    }

    return result;
}

Vector Vector::operator-(double sub)
{
    auto result{*this};

    for (auto i{0}; i < size; i++)
    {
        result[i] -= sub;
    }

    return result;
}

double Vector::dot(Vector rhs) const
{
    double result{0};

    for (auto i{0}; i < size; i++)
    {
        result += data[i] * rhs[i];
    }

    return result;
}