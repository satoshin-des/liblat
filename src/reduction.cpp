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

    for (long k = 1, i, j; k < m_num_rows;)
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

    for (long k = 1, j, i, t, l; k < m_num_rows;)
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
                deepInsertion(i, k);
                updateDeepInsGSO(i, k);
                
                k = fmax(i - 1, 0);
            }
        }
        ++k;
    }
}

template <class T>
void Lattice<T>::potLLL(const double delta, const bool compute_gso)
{
    double P, P_min, S;

    LLL(delta, compute_gso);

    for (long l = 0, j, i, k; l < m_num_rows;)
    {
        for (j = l - 1; j > -1; --j)
        {
            sizeReduce(l, j);
        }

        P = P_min = 1.0;
        k = 0;
        for (j = l - 1; j >= 0; --j)
        {
            S = 0;
            for (i = j; i < l; ++i)
            {
                S += m_mu[l][i] * m_mu[l][i] * m_B[i];
            }
            P *= (m_B[l] + S) / m_B[j];

            if (P < P_min)
            {
                k = j;
                P_min = P;
            }
        }

        if (delta > P_min)
        {
            deepInsertion(k, l);
            updateDeepInsGSO(k, l);
            l = k;
        }
        else
        {
            ++l;
        }
    }
}

template class Lattice<int>;
template class Lattice<long>;
template class Lattice<long long>;
template class Lattice<float>;
template class Lattice<double>;
