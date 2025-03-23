#include "lattice.h"

#include "core.h"

template <class T>
void Lattice<T>::LLL(const double delta, const bool compute_gso, long start_, long end_)
{
    try
    {
        if (delta < 0.25 || delta > 1)
        {
            throw std::out_of_range("[WARNING]The reduction parameter must be in [0.25, 1.0].");
        }
    }
    catch (const std::exception &ex)
    {
        std::cerr << ex.what() << "@ function " << __FUNCTION__ << std::endl;
    }

    if (end_ <= -2 || end_ == 0)
    {
        throw std::out_of_range("The parameter end_ must be a positive integer or -1.");
    }
    else if (end_ == -1)
    {
        if (start_ <= -1 || start_ >= m_num_rows)
        {
            throw std::out_of_range("The parameter start_ is out of index.");
        }
    }
    else
    {
        if (start_ <= -1 || start_ >= m_num_rows)
        {
            throw std::out_of_range("The parameter start_ is out of index.");
        }
        if (start_ >= end_)
        {
            throw std::invalid_argument("The parameter start_ must be less than the parameter end_.");
        }
        if (end_ > m_num_rows)
        {
            throw std::out_of_range("The parameter end_ is out of index.");
        }
    }

    long start = start_;
    long end;

    if (end_ == -1)
    {
        end = m_num_rows;
    }
    else
    {
        end = end_;
    }

    double nu, B, t;
    T tmp;

    if (compute_gso)
    {
        computeGSO();
    }

    for (long k = start + 1, i, j; k < end;)
    {
        for (j = k - 1; j > -1; --j)
        {
            sizeReduce(k, j);
        }

        if (k > start && m_B[k] < (delta - m_mu[k][k - 1] * m_mu[k][k - 1]) * m_B[k - 1])
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
            for (i = k + 1; i < end; ++i)
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
    try
    {
        if (delta < 0.25 || delta > 1)
        {
            throw std::out_of_range("[WARNING]The reduction parameter must be in [0.25, 1.0].");
        }
    }
    catch (const std::exception &ex)
    {
        std::cerr << ex.what() << "@ function " << __FUNCTION__ << std::endl;
    }

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

                k = std::max(i - 1, static_cast<long>(0));
            }
        }
        ++k;
    }
}

template <class T>
void Lattice<T>::potLLL(const double delta, const bool compute_gso)
{
    try
    {
        if (delta < 0.25 || delta > 1)
        {
            throw std::out_of_range("[WARNING]The reduction parameter must be in [0.25, 1.0].");
        }
    }
    catch (const std::exception &ex)
    {
        std::cerr << ex.what() << "@ function " << __FUNCTION__ << std::endl;
    }

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

template <class T>
void Lattice<T>::BKZ(const long beta, const double delta, const bool compute_gso)
{
    try
    {
        if (delta < 0.25 || delta > 1)
        {
            throw std::out_of_range("[WARNING]The reduction parameter must be in [0.25, 1.0].");
        }
    }
    catch (const std::exception &ex)
    {
        std::cerr << ex.what() << "@ function " << __FUNCTION__ << std::endl;
    }

    if (beta < 2 || beta > m_num_rows)
    {
        char err_s[100];
        sprintf(err_s, "The blocksize is %ld. The blocksize must be in [2, %ld].", beta, m_num_rows);
        throw std::out_of_range(err_s);
    }

    std::vector<T> v(m_num_cols);
    std::vector<long> w(m_num_rows);

    LLL(delta, true);

    for (long z = 0, j, t, num_tour = 0, i, k = 0, h, d, l; z < m_num_rows - 1;)
    {
        if (num_tour >= m_max_loop)
        {
            break;
        }

        if (k == m_num_rows - 1)
        {
            k = 0;
            ++num_tour;
        }
        ++k;
        l = std::min(k + beta - 1, m_num_rows);
        h = std::min(l + 1, m_num_rows);
        d = l - k + 1;

        w = enumShortVec_(false, k - 1, l);
        --w[0];

        if (!isZero(w))
        {
            ++w[0];
            z = 0;

            for (i = 0; i < m_num_cols; ++i)
            {
                v[i] = 0;
                for (j = 0; j < d; ++j)
                {
                    v[i] += w[j] * m_basis[j + k - 1][i];
                }
            }

            for (i = d - 1; i >= 0; --i)
            {
                if (w[i] == 1 || w[i] == -1)
                {
                    for (j = 0; j < m_num_cols; ++j)
                    {
                        m_basis[i + k - 1][j] = v[j];
                    }

                    deepInsertion(k - 1, i + k - 1);
                    break;
                }
            }

            LLL(delta, true, 0, h);
        }
        else
        {
            ++z;
            LLL(delta, false, 0, h);
        }
    }
}

template<class T>
void Lattice<T>::HKZ(const double delta, const bool compute_gso)
{
    BKZ(m_num_rows, delta, compute_gso);
}

template class Lattice<int>;
template class Lattice<long>;
template class Lattice<long long>;
template class Lattice<float>;
template class Lattice<double>;
