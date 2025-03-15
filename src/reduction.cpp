#include "lattice.h"

#include "core.h"

template <class T>
void Lattice<T>::LLL(const double delta, const bool compute_gso)
{
    double nu, B, t;
    T tmp;

    if (compute_gso)
    {
        computeGSO();
    }

    for (int k = 1, i, j; k < m_num_rows;)
    {
        for (j = k - 1; j > -1; --j)
        {
            sizeReduce(k, j);
        }

        if (k > 0 && m_B[k] < (delta - m_mu[k][k - 1] * m_mu[k][k - 1]) * m_B[k - 1])
        {
            for (i = 0; i < m_num_cols; ++i)
            {
                tmp = m_basis[k - 1][i];
                m_basis[k - 1][i] = m_basis[k][i];
                m_basis[k][i] = tmp;
            }

            nu = m_mu[k][k - 1];
            B = m_B[k] + nu * nu * m_B[k - 1];
            m_mu[k][k - 1] = nu * m_B[k - 1] / B;
            m_B[k] *= m_B[k - 1] / B;
            m_B[k - 1] = B;

            for (i = 0; i < k - 1; ++i)
            {
                t = m_mu[k - 1][i];
                m_mu[k - 1][i] = m_mu[k][i];
                m_mu[k][i] = t;
            }
            for (i = k + 1; i < m_num_rows; ++i)
            {
                t = m_mu[i][k];
                m_mu[i][k] = m_mu[i][k - 1] - nu * t;
                m_mu[i][k - 1] = t + m_mu[k][k - 1] * m_mu[i][k];
            }

            --k;
        }
        else
        {
            ++k;
        }
    }
}

template <class T>
void Lattice<T>::deepLLL(const double delta, const bool compute_gso)
{
    double C;

    if (compute_gso)
    {
        computeGSO();
    }

    for (int k = 1, j, i, t, l; k < m_num_rows;)
    {
        for (j = k - 1; j >= 0; --j)
        {
            sizeReduce(k, j);
        }

        C = static_cast<double>(dot(m_basis[k], m_basis[k]));

        for (i = 0; i < k;)
        {
            if (C >= delta * m_B[i])
            {
                C -= m_mu[k][i] * m_mu[k][i] * m_B[i];
                ++i;
            }
            else
            {
                for (l = 0; l < m_num_cols; ++l)
                {
                    t = m_basis[k][l];
                    for (j = k; j > i; --j)
                    {
                        m_basis[j][l] = m_basis[j - 1][l];
                    }
                    m_basis[i][l] = t;
                }
                
                updateDeepInsGSO(i, k);
                k = fmax(i - 1, 0);
            }
        }
        ++k;
    }
}

template class Lattice<int>;
template class Lattice<long>;
template class Lattice<long long>;
template class Lattice<float>;
template class Lattice<double>;
