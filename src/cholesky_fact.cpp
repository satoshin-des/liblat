#include "lattice.h"

#include "core.h"

#include <iostream>

template <class T>
void Lattice<T>::choleskyFact()
{
    for (long k = 0, j, h; k < m_num_rows; ++k)
    {
        for (j = 0; j <= k; ++j)
        {
            m_r[k][j] = dot(m_basis[k], m_basis[j]);
            for (h = 0; h < j; ++h)
            {
                m_r[k][j] -= m_r[k][h] * m_mu[j][h];
            }
            m_mu[k][j] = m_r[k][j] / m_r[j][j];
        }
    }
    m_s[0] = dot(m_basis[m_num_rows - 1], m_basis[m_num_rows - 1]);
    for (long j = 0; j < m_num_rows; ++j)
    {
        m_s[j + 1] = m_s[j] - m_mu[m_num_rows - 1][j] * m_r[m_num_rows - 1][j];
    }
    m_r[m_num_rows - 1][m_num_rows - 1] = m_s[m_num_rows];
}

template class Lattice<int>;
template class Lattice<long>;
template class Lattice<long long>;
template class Lattice<float>;
template class Lattice<double>;