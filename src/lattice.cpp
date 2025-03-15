#include "lattice.h"

#include <iostream>
#include <vector>
#include <cstdlib>
#include <sstream>
#include <cmath>

#include "core.h"

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

template <class T>
void Lattice<T>::updateDeepInsGSO(const int i, const int k)
{
    int j, l;
    double t, eps;
    std::vector<double> P(m_num_rows, 0), D(m_num_rows, 0), S(m_num_rows, 0);

    P[k] = D[k] = m_B[k];
    for (j = k - 1; j >= i; --j)
    {
        P[j] = m_mu[k][j] * m_B[j];
        D[j] = D[j + 1] + m_mu[k][j] * P[j];
    }

    for (j = k; j > i; --j)
    {
        t = m_mu[k][j - 1] / D[j];
        for (l = m_num_rows - 1; l > k; --l)
        {
            S[l] += m_mu[l][j] * P[j];
            m_mu[l][j] = m_mu[l][j - 1] - t * S[l];
        }
        for (l = k; l > j; --l)
        {
            S[l] += m_mu[l - 1][j] * P[j];
            m_mu[l][j] = m_mu[l - 1][j - 1] - t * S[l];
        }
    }

    t = 1.0 / D[i];

    for (l = m_num_rows - 1; l > k; --l)
    {
        m_mu[l][i] = t * (S[l] + m_mu[l][i] * P[i]);
    }
    for (l = k; l >= i + 2; --l)
    {
        m_mu[l][i] = t * (S[l] + m_mu[l - 1][i] * P[i]);
    }

    m_mu[i + 1][i] = t * P[i];
    for (j = 0; j < i; ++j)
    {
        eps = m_mu[k][j];
        for (l = k; l > i; --l)
        {
            m_mu[l][j] = m_mu[l - 1][j];
        }
        m_mu[i][j] = eps;
    }

    for (j = k; j > i; --j)
    {
        m_B[j] = D[j] * m_B[j - 1] / D[j - 1];
    }
    m_B[i] = D[i];
}

template <class T>
void Lattice<T>::sizeReduce(const int i, const int j)
{
    if (m_mu[i][j] > 0.5 || m_mu[i][j] < -0.5)
    {
        int k;
        const long q = static_cast<long>(round(m_mu[i][j]));

        for (k = 0; k < m_num_cols; ++k)
        {
            m_basis[i][k] -= static_cast<T>(q) * m_basis[j][k];
        }
        for (k = 0; k <= j; ++k)
        {
            m_mu[i][k] -= m_mu[j][k] * static_cast<double>(q);
        }
    }
}

template class Lattice<int>;
template class Lattice<long>;
template class Lattice<long long>;
template class Lattice<float>;
template class Lattice<double>;
