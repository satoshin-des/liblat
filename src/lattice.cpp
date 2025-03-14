#include "lattice.h"

#include <iostream>
#include <vector>
#include <cstdlib>
#include <sstream>

template <class T>
void Lattice<T>::setDims(const long n, const long m)
{
    m_num_rows = n;
    m_num_cols = m;
    m_basis = std::vector<std::vector<T>>(n, std::vector<T>(m, 0));
    m_b_star = std::vector<std::vector<double>>(n, std::vector<double>(m, 0));
    m_mu = std::vector<std::vector<double>>(n, std::vector<double>(n, 0));
    m_B = std::vector<double>(n, 0);
}

template <class T>
void Lattice<T>::setRandom(const long n, const long m, const T min, const T max)
{
    if ((n != 0) && (m != 0))
    {
        setDims(n, m);
    }

    m_get_rand_uni = std::uniform_real_distribution<>(min, max);

    for (int i = 0; i < n; ++i)
    {
        m_basis[i][i] = 1;
        m_basis[i][0] = static_cast<T>(m_get_rand_uni(m_mt64));
    }
}

template <class T>
void Lattice<T>::computeGSO()
{
    for (int i = 0, j, k; i < m_num_rows; ++i)
    {
        m_mu[i][i] = 1;

        for (j = 0; j < m_num_cols; ++j)
        {
            m_b_star[i][j] = static_cast<double>(m_basis[i][j]);
        }

        for (j = 0; j < i; ++j)
        {
            m_mu[i][j] = dot(m_basis[i], m_b_star[j]) / dot(m_b_star[j], m_b_star[j]);
            for (k = 0; k < m_num_cols; ++k)
            {
                m_b_star[i][k] -= m_mu[i][j] * m_b_star[j][k];
            }
        }
        m_B[i] = dot(m_b_star[i], m_b_star[i]);
    }
}

template class Lattice<int>;
template class Lattice<long>;
template class Lattice<double>;
