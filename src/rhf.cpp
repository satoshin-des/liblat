#include "lattice.h"

#include <iostream>
#include <cmath>

template <class T>
long double Lattice<T>::rhf()
{
    const long double hf = b1Norm() / powl(m_vol, 1.0 / m_num_rows);
    return powl(hf, 1.0 / m_num_rows);
}

template class Lattice<int>;
template class Lattice<long>;
template class Lattice<long long>;
template class Lattice<float>;
template class Lattice<double>;
