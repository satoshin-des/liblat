#include "core.h"

#include <iostream>

template <class T>
void print(const std::vector<T> v)
{
    std::cout << "[";
    for (T w : v)
    {
        std::cout << w << ", ";
    }
    std::cout << "\b\b" << "]\n";
}

template <class U, class V>
V dot(const std::vector<U> x, const std::vector<V> y)
{
    if (x.size() != y.size())
    {
        char *err_s;
        sprintf(err_s, "An inner product of %ld-th vector and %ld-th vector cannot be defined.", x.size(), y.size());
        throw std::invalid_argument(err_s);
    }

    V S = 0;
    for (int i = 0; i < x.size(); ++i)
    {
        S += static_cast<V>(x[i]) * y[i];
    }
    return S;
}

template <class T>
bool isZero(const std::vector<T> v)
{
    for (T w : v)
    {
        if (w != 0)
        {
            return false;
        }
    }
    return true;
}

template void print<int>(const std::vector<int>);
template void print<long>(const std::vector<long>);
template int dot<int, int>(std::vector<int>, std::vector<int>);
template double dot<int, double>(std::vector<int>, std::vector<double>);
template double dot<long, double>(std::vector<long>, std::vector<double>);
template double dot<long long, double>(std::vector<long long>, std::vector<double>);
template double dot<float, double>(std::vector<float>, std::vector<double>);
template double dot<double, double>(std::vector<double>, std::vector<double>);
template long dot<long, long>(std::vector<long>, std::vector<long>);
template long long dot<long long, long long>(std::vector<long long>, std::vector<long long>);
template float dot<float, float>(std::vector<float>, std::vector<float>);
template bool isZero<int>(const std::vector<int>);
template bool isZero<long>(const std::vector<long>);
