#include "core.h"

#include <iostream>

template <class U, class V>
V dot(const std::vector<U> x, const std::vector<V> y)
{
    V S = 0;
    for (int i = 0; i < x.size(); ++i)
    {
        S += static_cast<V>(x[i]) * y[i];
    }
    return S;
}

template int dot<int, int>(std::vector<int>, std::vector<int>);
template double dot<int, double>(std::vector<int>, std::vector<double>);
template double dot<long, double>(std::vector<long>, std::vector<double>);
template double dot<long long, double>(std::vector<long long>, std::vector<double>);
template double dot<float, double>(std::vector<float>, std::vector<double>);
template double dot<double, double>(std::vector<double>, std::vector<double>);
template long dot<long, long>(std::vector<long>, std::vector<long>);
template long long dot<long long, long long>(std::vector<long long>, std::vector<long long>);
template float dot<float, float>(std::vector<float>, std::vector<float>);
