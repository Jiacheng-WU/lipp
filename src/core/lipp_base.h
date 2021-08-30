#ifndef __LIPP_BASE_H__
#define __LIPP_BASE_H__

#include <limits>
#include <cmath>
#include <cstdlib>
#include <algorithm>

// Linear regression model
template <class T>
class LinearModel
{
public:
    double a = 0; // slope
    double b = 0; // intercept
    T bias = 0; // bias

    LinearModel() = default;
    LinearModel(double a, double b, T bias) : a(a), b(b), bias(bias) {}
    explicit LinearModel(const LinearModel &other) : a(other.a), b(other.b), bias(other.bias) {}

    inline int predict(T key) const
    {
        return std::floor(a * static_cast<double>(key - bias) + b);
    }

    inline double predict_double(T key) const
    {
        return a * static_cast<double>(key - bias) + b;
    }
};

#endif // __LIPP_BASE_H__
