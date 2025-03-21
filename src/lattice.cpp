#include "lattice.h"

#include <iostream>
#include <vector>
#include <cstdlib>
#include <sstream>
#include <cmath>

#include "core.h"

template <class T>
std::vector<long> Lattice<T>::ENUM_(const double R, const long start, const long end)
{
    long n = end - start;
    long i, r[n + 1];
    long last_nonzero = 0; // index of last non-zero elements
    double temp;
    std::vector<long> weight(n, 0);
    std::vector<long> coeff_vector(n, 0);
    std::vector<double> center(n, 0);
    std::vector<std::vector<double>> sigma(n + 1, std::vector<double>(n, 0));
    std::vector<double> rho(n + 1, 0);

    coeff_vector[0] = 1;
    for (i = 0; i < n; ++i)
    {
        r[i] = i;
    }

    for (long k = 0;;)
    {
        temp = static_cast<double>(coeff_vector[k]) - center[k];
        temp *= temp;
        rho[k] = rho[k + 1] + temp * m_B[k + start]; // rho[k]=∥πₖ(shortest_vec)∥
        if (rho[k] <= R)
        {
            if (k == 0)
            {
                return coeff_vector;
            }
            else
            {
                --k;
                if (r[k + 1] >= r[k])
                {
                    r[k] = r[k + 1];
                }
                for (i = r[k]; i > k; --i)
                {
                    sigma[i][k] = sigma[i + 1][k] + m_mu[i + start][k + start] * coeff_vector[i];
                }
                center[k] = -sigma[k + 1][k];
                coeff_vector[k] = round(center[k]);
                weight[k] = 1;
            }
        }
        else
        {
            ++k;
            if (k == n)
            { // no solution
                coeff_vector = std::vector<long>(n, 0);
                return coeff_vector;
            }
            else
            {
                r[k] = k;
                if (k >= last_nonzero)
                {
                    last_nonzero = k;
                    ++coeff_vector[k];
                }
                else
                {
                    if (coeff_vector[k] > center[k])
                    {
                        coeff_vector[k] -= weight[k];
                    }
                    else
                    {
                        coeff_vector[k] += weight[k];
                    }

                    ++weight[k];
                }
            }
        }
    }
}

template <class T>
std::vector<long> Lattice<T>::enumShortVec_(const bool compute_gso, const long start, const long end)
{
    const long n = end - start;
    std::vector<long> enum_vector(n);
    std::vector<long> old_enum_vector(n);

    if (compute_gso)
    {
        computeGSO();
    }

    for (double R = m_B[start];; R *= 0.99)
    {
        for (long i = 0; i < n; ++i)
        {
            old_enum_vector[i] = enum_vector[i];
        }

        enum_vector = ENUM_(R, start, end);
        if (isZero(enum_vector))
        {
            return old_enum_vector;
        }
    }
}

template <class T>
void Lattice<T>::setMaxLoop(const long max_loop)
{
    try
    {
        if (max_loop <= 0)
        {
            throw std::out_of_range("The augment max_loop must be a positive integer.");
        }
    }
    catch (const std::exception &ex)
    {
        std::cerr << ex.what() << "@ function " << __FUNCTION__ << std::endl;
        exit(EXIT_FAILURE);
    }

    m_max_loop = max_loop;
}

template <class T>
long Lattice<T>::numRows() const
{
    return m_num_rows;
}

template <class T>
long Lattice<T>::numCols() const
{
    return m_num_cols;
}

template <class T>
void Lattice<T>::setDims(const long n, const long m)
{
    try
    {
        if (n <= 0)
        {
            throw std::invalid_argument("The number of rows of basis must be positive integer.");
        }

        if (n <= 0)
        {
            throw std::invalid_argument("The number of columns of basis must be positive integer.");
        }
    }
    catch (const std::exception &ex)
    {
        std::cerr << ex.what() << "@ function " << __FUNCTION__ << std::endl;
        exit(EXIT_FAILURE);
    }

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
    if (n <= 0)
    {
        throw std::invalid_argument("The number of rows of basis must be positive integer.");
    }

    if (n <= 0)
    {
        throw std::invalid_argument("The number of columns of basis must be positive integer.");
    }

    if (min > max)
    {
        throw std::invalid_argument("The augment min must be less than the augment max.");
    }

    if ((m_num_rows != n) && (m_num_cols != m))
    {
        setDims(n, m);
    }

    m_get_rand_uni = std::uniform_real_distribution<>(min, max);

    for (long i = 0; i < n; ++i)
    {
        m_basis[i][i] = 1;
        m_basis[i][0] = static_cast<T>(m_get_rand_uni(m_mt64));
    }
}

template <class T>
std::vector<T> Lattice<T>::mulVecBasis(const std::vector<long> v)
{
    std::vector<T> w(m_num_cols);
    long sum = 0;
    for (long j = 0, i; j < m_num_cols; ++j)
    {
        sum = 0;
        for (i = 0; i < m_num_rows; ++i)
        {
            sum += v[i] * static_cast<long>(m_basis[i][j]);
        }
        w[j] = static_cast<T>(sum);
    }
    return w;
}

template <class T>
void Lattice<T>::deepInsertion(const long k, const long l)
{
    if (k < 0 || k >= m_num_rows)
    {
        char err_s[100];
        sprintf(err_s, "%ld is out of index.", k);
        throw std::out_of_range(err_s);
    }

    if (l < 0 || l >= m_num_rows)
    {
        char err_s[100];
        sprintf(err_s, "%ld is out of index.", l);
        throw std::out_of_range(err_s);
    }

    T t;
    for (long i = 0, j; i < m_num_cols; ++i)
    {
        t = m_basis[l][i];
        for (j = l; j > k; --j)
        {
            m_basis[j][i] = m_basis[j - 1][i];
        }
        m_basis[k][i] = t;
    }
}

template <class T>
void Lattice<T>::computeGSO()
{
    for (long i = 0, j, k; i < m_num_rows; ++i)
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
void Lattice<T>::updateDeepInsGSO(const long i, const long k)
{
    long j, l;
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
void Lattice<T>::sizeReduce(const long i, const long j)
{
    if (m_mu[i][j] > 0.5 || m_mu[i][j] < -0.5)
    {
        long k;
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

template <class T>
std::vector<long> Lattice<T>::ENUM(const double R)
{
    long i, r[m_num_rows + 1];
    long last_nonzero = 0; // index of last non-zero elements
    double temp;
    std::vector<long> weight(m_num_rows, 0);
    std::vector<long> coeff_vector(m_num_rows, 0);
    std::vector<double> center(m_num_rows, 0);
    std::vector<std::vector<double>> sigma(m_num_rows + 1, std::vector<double>(m_num_rows, 0));
    std::vector<double> rho(m_num_rows + 1, 0);

    coeff_vector[0] = 1;
    for (i = 0; i < m_num_rows; ++i)
    {
        r[i] = i;
    }

    for (long k = 0;;)
    {
        temp = static_cast<double>(coeff_vector[k]) - center[k];
        temp *= temp;
        rho[k] = rho[k + 1] + temp * m_B[k]; // rho[k]=∥πₖ(shortest_vec)∥
        if (rho[k] <= R)
        {
            if (k == 0)
            {
                return coeff_vector;
            }
            else
            {
                --k;
                if (r[k + 1] >= r[k])
                {
                    r[k] = r[k + 1];
                }
                for (i = r[k]; i > k; --i)
                {
                    sigma[i][k] = sigma[i + 1][k] + m_mu[i][k] * coeff_vector[i];
                }
                center[k] = -sigma[k + 1][k];
                coeff_vector[k] = round(center[k]);
                weight[k] = 1;
            }
        }
        else
        {
            ++k;
            if (k == m_num_rows)
            { // no solution
                coeff_vector = std::vector<long>(m_num_rows, 0);
                return coeff_vector;
            }
            else
            {
                r[k] = k;
                if (k >= last_nonzero)
                {
                    last_nonzero = k;
                    ++coeff_vector[k];
                }
                else
                {
                    if (coeff_vector[k] > center[k])
                    {
                        coeff_vector[k] -= weight[k];
                    }
                    else
                    {
                        coeff_vector[k] += weight[k];
                    }

                    ++weight[k];
                }
            }
        }
    }
}

template <class T>
std::vector<T> Lattice<T>::enumShortVec(const bool compute_gso)
{
    std::vector<long> enum_vector(m_num_rows);
    std::vector<long> old_enum_vector(m_num_rows);

    if (compute_gso)
    {
        computeGSO();
    }

    for (double R = m_B[0];; R *= 0.99)
    {
        for (long i = 0; i < m_num_rows; ++i)
        {
            old_enum_vector[i] = enum_vector[i];
        }
        enum_vector = ENUM(R);
        if (isZero(enum_vector))
        {
            return mulVecBasis(old_enum_vector);
        }
    }
}

template class Lattice<int>;
template class Lattice<long>;
template class Lattice<long long>;
template class Lattice<float>;
template class Lattice<double>;
