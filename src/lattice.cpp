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
    if ((start < 0) || (start >= m_num_rows))
    {
        throw std::out_of_range("The argument start is out of index. @ function ENUM_.");
    }
    if ((end < 0) || (end > m_num_rows))
    {
        throw std::out_of_range("The argument end is out of index. @ function ENUM_.");
    }
    if (start >= end)
    {
        throw std::invalid_argument("The arguments start and end is invalid. @ function ENUM_.");
    }

    if (R <= 0)
    {
        return std::vector<long>(m_num_rows, 0);
    }

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
    if ((start < 0) || (start >= m_num_rows))
    {
        char err_s[ERR_STR_LEN];
        sprintf(err_s, "%ld is out of index. @ function enumShortVec_.", start);
        throw std::out_of_range(err_s);
    }
    if ((end < 0) || (end >= m_num_rows))
    {
        char err_s[ERR_STR_LEN];
        sprintf(err_s, "%ld is out of index. @ function enumShortVec_.", end);
        throw std::out_of_range(err_s);
    }
    if (start >= end)
    {
        throw std::invalid_argument("The arguments start and end is invalid. @ function enumShortVec_.");
    }

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
std::vector<long> Lattice<T>::dualENUM_(const double R, const long start, const long end)
{
    if ((start < 0) || (start >= m_num_rows))
    {
        throw std::out_of_range("The argument start is out of index. @ function dualENUM_.");
    }
    if ((end < 0) || (end > m_num_rows))
    {
        throw std::out_of_range("The argument end is out of index. @ function dualENUM_.");
    }
    if (start >= end)
    {
        throw std::invalid_argument("The arguments start and end is invalid. @ function dualENUM_.");
    }

    if (R <= 0)
    {
        return std::vector<long>(m_num_rows, 0);
    }

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
        rho[k] = rho[k + 1] + temp / m_B[k + start]; // rho[k]=∥πₖ(shortest_vec)∥
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
                    sigma[i][k] = sigma[i + 1][k] + m_dual_mu[i + start][k + start] * coeff_vector[i];
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
std::vector<long> Lattice<T>::dualEnumShortVec_(const bool compute_gso, const long start, const long end)
{
    if ((start < 0) || (start >= m_num_rows))
    {
        throw std::out_of_range("The argument start is out of index. @ function dualEnumShortVec_.");
    }
    if ((end < 0) || (end > m_num_rows))
    {
        throw std::out_of_range("The argument end is out of index. @ function dualEnumShortVec_.");
    }
    if (start >= end)
    {
        throw std::invalid_argument("The arguments start and end is invalid. @ function dualEnumShortVec_.");
    }

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
        enum_vector = dualENUM_(R, start, end);

        if (isZero(enum_vector))
        {
            return old_enum_vector;
        }
    }
}

template <class T>
void Lattice<T>::setMaxLoop(const long max_loop)
{
    if (max_loop <= 0)
    {
        throw std::out_of_range("The augment max_loop must be a positive integer.");
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
    if (n <= 0)
    {
        throw std::invalid_argument("The number of rows of basis must be positive integer.");
    }
    if (m <= 0)
    {
        throw std::invalid_argument("The number of columns of basis must be positive integer.");
    }

    m_num_rows = n;
    m_num_cols = m;
    m_basis = std::vector<std::vector<T>>(n, std::vector<T>(m, 0));
    m_b_star = std::vector<std::vector<double>>(n, std::vector<double>(m, 0));
    m_dual_b_star = std::vector<std::vector<double>>(n, std::vector<double>(m, 0));
    m_mu = std::vector<std::vector<double>>(n, std::vector<double>(n, 0));
    m_dual_mu = std::vector<std::vector<double>>(n, std::vector<double>(n, 0));
    m_B = std::vector<double>(n, 0);
    m_dual_B = std::vector<double>(n, 0);
}

template <class T>
void Lattice<T>::setRandom(const long n, const long m, const T min, const T max)
{
    if (n <= 0)
    {
        throw std::invalid_argument("The number of rows of basis must be positive integer.");
    }

    if (m <= 0)
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
void Lattice<T>::setBasis(const std::vector<std::vector<T>> basis_mat)
{
    if ((m_num_rows != basis_mat.size()) && (m_num_cols != basis_mat.at(0).size()))
    {
        setDims(basis_mat.size(), basis_mat[0].size());
    }

    for (long i = 0, j; i < m_num_rows; ++i)
    {
        for (j = 0; j < m_num_cols; ++j)
        {
            m_basis[i][j] = basis_mat[i][j];
        }
    }
}

template <class T>
void Lattice<T>::setGoldesteinMayerLattice(const T p, const T q)
{
    if ((q < 1) || (p < 1))
    {
        throw std::invalid_argument("The argument p and q must be greater than or equal to 1.");
    }

    if (q > p)
    {
        throw std::invalid_argument("The argument p must be greater than or equal to q.");
    }

    for (long i = 0, j; i < m_num_rows; ++i)
    {
        for (j = 0; j < m_num_cols; ++j)
        {
            m_basis[i][j] = 0;
        }
    }

    for (long i = 0; i < m_num_rows; ++i)
    {
        m_basis[i][i] = 1;
        m_basis[i][0] = q;
    }
    m_basis[0][0] = p;
}

template <class T>
void Lattice<T>::setSchnorrLattice(const long N, const double c)
{
    if (m_num_rows + 1 != m_num_cols)
    {
        throw std::invalid_argument("Schnorr lattice can be defined if and only if the size of a lattice basis matrix is (N, N + 1).  @ function setSchnorrLattice");
    }
    if (N <= 0)
    {
        throw std::invalid_argument("The augment N must be a positive integer. @ function setSchnorrLattice");
    }

    for (long i = 0, j; i < m_num_rows; ++i)
    {
        for (j = 0; j < m_num_cols; ++j)
        {
            m_basis[i][j] = 0;
        }

        m_basis[i][i] = sqrt(log(prime(i + 1)));
        m_basis[i][m_num_cols - 1] = std::pow(N, c) * log(prime(i + 1));
    }
}

template <class T>
T Lattice<T>::volume(const bool compute_gso)
{
    if (compute_gso)
    {
        computeGSO();
    }

    double vol = 1.0;
    for (long i = 0; i < m_num_rows; ++i)
    {
        vol *= m_B[i];
    }
    return sqrt(vol);
}

template <class T>
double Lattice<T>::potential(const bool compute_gso)
{
    if (compute_gso)
    {
        computeGSO();
    }

    double pot = 1.0;
    for (long i = 0; i < m_num_rows; ++i)
    {
        pot *= pow(m_B[i], m_num_rows - i);
    }
    return pot;
}

template <class T>
double Lattice<T>::logPotential(const bool compute_gso)
{
    if (compute_gso)
    {
        computeGSO();
    }

    double pot = 0.0;
    for (long i = 0; i < m_num_rows; ++i)
    {
        pot += (m_num_rows - i) * log(m_B[i]);
    }
    return pot;
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
    if ((k < 0) || (k >= m_num_rows))
    {
        char err_s[ERR_STR_LEN];
        sprintf(err_s, "%ld is out of index. @ function deepInsertion.", k);
        throw std::out_of_range(err_s);
    }

    if ((l < 0) || (l >= m_num_rows))
    {
        char err_s[ERR_STR_LEN];
        sprintf(err_s, "%ld is out of index. @ function deepInsertion.", l);
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
void Lattice<T>::dualDeepInsertion(const long k, const long l)
{
    if ((k < 0) || (k >= m_num_rows))
    {
        char err_s[ERR_STR_LEN];
        sprintf(err_s, "%ld is out of index. @ function dualDeepInsertion.", k);
        throw std::out_of_range(err_s);
    }
    if ((l < 0) || (l >= m_num_rows))
    {
        char err_s[ERR_STR_LEN];
        sprintf(err_s, "%ld is out of index. @ function dualDeepInsertion.", l);
        throw std::out_of_range(err_s);
    }

    T t;
    for (long i = 0; i < m_num_cols; ++i)
    {
        t = m_basis[k][i];
        for (int j = k; j < l; ++j)
        {
            m_basis[j][i] = m_basis[j + 1][i];
        }
        m_basis[l][i] = t;
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
void Lattice<T>::computeDualGSO()
{
    double sum;
    for (long i = 0, j, k; i < m_num_rows; ++i)
    {
        for (j = i + 1; j < m_num_cols; ++j)
        {
            sum = 0;
            for (k = 0; k < j - i; ++k)
            {
                sum += m_mu[j][i + k] * m_dual_mu[i][i + k];
            }
        }
    }
}

template <class T>
void Lattice<T>::updateDeepInsGSO(const long i, const long k, const long start, const long end)
{
    if (k <= i)
    {
        return;
    }

    if ((start < 0) || (start >= m_num_rows))
    {
        throw std::out_of_range("The argument start is out of index. @ function updateDeepInsGSO.");
    }
    if ((end < 0) || (end > m_num_rows))
    {
        throw std::out_of_range("The argument end is out of index. @ function updateDeepInsGSO.");
    }
    if (start >= end)
    {
        throw std::invalid_argument("The arguments start and end is invalid. @ function updateDeepInsGSO.");
    }

    long j, l;
    long n = end - start;
    double t, eps;
    std::vector<double> P(n, 0), D(n, 0), S(n, 0);

    P[k] = D[k] = m_B[k];
    for (j = k - 1; j >= i; --j)
    {
        P[j] = m_mu[k][j] * m_B[j];
        D[j] = D[j + 1] + m_mu[k][j] * P[j];
    }

    for (j = k; j > i; --j)
    {
        t = m_mu[k][j - 1] / D[j];
        for (l = end - 1; l > k; --l)
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

    for (l = end - 1; l > k; --l)
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
void Lattice<T>::updateDualDeepInsGSO(const long k, const long l, const std::vector<double> dual_D)
{
    long i, j, h;
    double sum;
    std::vector<std::vector<double>> xi(m_num_rows, std::vector<double>(m_num_rows, 0));

    for (i = 0; i < m_num_rows; ++i)
    {
        for (j = 0; j < m_num_rows; ++j)
        {
            xi[i][j] = m_mu[i][j];
        }
    }

    for (i = l + 1; i < m_num_rows; ++i)
    {
        sum = 0.0;
        for (h = k; h <= l; ++h)
        {
            sum += m_dual_mu[k][h] * m_mu[i][h];
        }
        xi[i][l] = sum;
    }

    for (j = k; j < l; ++j)
    {
        for (i = j + 1; i < l; ++i)
        {
            sum = 0.0;
            for (h = k; h <= j; ++h)
            {
                sum += m_dual_mu[k][h] * m_mu[i + 1][h];
            }

            xi[i][j] = m_mu[i + 1][j + 1] * dual_D[j] / dual_D[j + 1] - m_dual_mu[k][j + 1] / (dual_D[j + 1] * m_B[j + 1]) * sum;
        }
        xi[l][j] = -m_dual_mu[k][j + 1] / (dual_D[j + 1] * m_B[j + 1]);
        for (i = l + 1; i < m_num_rows; ++i)
        {
            sum = 0;
            for (h = k; h <= j; ++h)
            {
                sum += m_dual_mu[k][h] * m_mu[i][h];
            }

            xi[i][j] = m_mu[i][j + 1] * dual_D[j] / dual_D[j + 1] - m_dual_mu[k][j + 1] / (dual_D[j + 1] * m_B[j + 1]) * sum;
        }
    }

    for (j = 0; j < k; ++j)
    {
        for (i = k; i < l; ++i)
        {
            xi[i][j] = m_mu[i + 1][j];
        }
        xi[l][j] = m_mu[k][j];
    }

    for (i = 0; i < m_num_rows; ++i)
    {
        for (j = 0; j < m_num_rows; ++j)
        {
            m_mu[i][j] = xi[i][j];
        }
    }

    for (j = k; j < l; ++j)
    {
        m_B[j] = dual_D[j + 1] * m_B[j + 1] / dual_D[j];
    }
    m_B[l] = 1.0 / dual_D[l];
}

template <class T>
void Lattice<T>::sizeReduce(const long i, const long j)
{
    if ((m_mu[i][j] > 0.5) || (m_mu[i][j] < -0.5))
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
    if (R <= 0)
    {
        return std::vector<long>(m_num_rows, 0);
    }

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

template <class T>
std::vector<long> Lattice<T>::potENUM(const long start, const long n)
{
    if((start < 0) || (start >= m_num_rows))
    {
        throw std::out_of_range("The argument start is out of index. @ function potENUM.");
    }
    if((n <= 0) || (n > m_num_rows))
    {
        throw std::invalid_argument("The argument n is invalid. @ function potENUM.");
    }

    long i, r[n + 1];
    double R = log(m_B[start]), P = 0, temp;
    std::vector<long> w(n, 0), v(n, 0);
    std::vector<double> c(n, 0), D(n + 1, 0);
    std::vector<std::vector<double>> sigma(n + 1, std::vector<double>(n, 0));

    v[0] = 1;

    for (i = 0; i <= n; ++i)
    {
        r[i] = i;
    }

    for (long k = 0, last_nonzero = 0;;)
    {
        temp = static_cast<double>(v[k]) - c[k];
        temp *= temp;
        D[k] = D[k + 1] + temp * m_B[k + start];

        if ((k + 1) * log(D[k]) + P < (k + 1) * log(0.99) + R)
        {
            if (k == 0)
            {
                return v;
            }
            else
            {
                P += log(D[k]);
                --k;

                if (r[k] <= r[k + 1])
                {
                    r[k] = r[k + 1];
                }

                for (i = r[k]; i > k; --i)
                {
                    sigma[i][k] = sigma[i + 1][k] + m_mu[i + start][k + start] * v[i];
                }
                c[k] = -sigma[k + 1][k];
                v[k] = static_cast<long>(round(c[k]));
                w[k] = 1;
            }
        }
        else
        {
            ++k;
            if (k == n)
            {
                std::fill(v.begin(), v.end(), 0);
                return v;
            }
            else
            {
                r[k - 1] = k;
                if (k >= last_nonzero)
                {
                    last_nonzero = k;
                    ++v[k];
#if 1
                    if (v[last_nonzero] >= 2)
                    {
                        ++k;
                        if (k == n)
                        {
                            std::fill(v.begin(), v.end(), 0);
                            return v;
                        }
                        else
                        {
                            r[k - 1] = k;
                            last_nonzero = k;
                            v[last_nonzero] = 1;
                        }
                    }
#endif
                    P = 0;
                    R = 0;
                    for (i = 0; i <= last_nonzero; ++i)
                    {
                        R += log(m_B[i + start]);
                    }
                }
                else
                {
                    if (v[k] > c[k])
                    {
                        v[k] -= w[k];
                    }
                    else
                    {
                        v[k] += w[k];
                    }

                    ++w[k];
                    P -= log(D[k]);
                }
            }
        }
    }
}

template <class T>
void Lattice<T>::insertToDualBasis(const std::vector<long> x, const long dim)
{
    if((dim <= 0) || (dim > m_num_rows))
    {
        throw std::invalid_argument("The augment dim is invalid. @ function insertToDualBasis.");
    }
    if (x.size() != dim)
    {
        char err_s[ERR_STR_LEN];
        sprintf(err_s, "%ld-th vector cannot insert into the dual basis( its size is %ld times %ld). @function insertToDualBasis.", x.size(), dim, m_num_cols);
        throw std::invalid_argument(err_s);
    }

    long i, j, k;
    const double beta = 1.16247638743819280699046939662833968000372102125550;
    double gamma, temp;
    Lattice<T> U(dim, m_num_cols);
    Lattice<T> temp_basis(dim, dim + 1);

    temp = dot(x, x);
    temp *= pow(beta, dim - 2);
    gamma = round(temp + temp);

    for (i = 0; i < dim; ++i)
    {
        for (j = 0; j < dim; ++j)
        {
            temp_basis.m_basis[i][j] = 0;
        }
        temp_basis.m_basis[i][i] = 1;
        temp_basis.m_basis[i][dim] = gamma * x[i];
    }

    temp_basis.LLL(0.99, true, 0, dim);

    for (i = 0; i < dim; ++i)
    {
        for (j = 0; j < m_num_cols; ++j)
        {
            for (k = 0; k < dim; ++k)
            {
                U.m_basis[i][j] = static_cast<T>(temp_basis.m_basis[i][k]) * m_basis[k][j];
            }
        }
    }

    for (i = 0; i < dim; ++i)
    {
        for (j = 0; j < m_num_cols; ++j)
        {
            m_basis[i][j] = U.m_basis[i][j];
        }
    }

    computeGSO();
}

template class Lattice<int>;
template class Lattice<long>;
template class Lattice<long long>;
template class Lattice<float>;
template class Lattice<double>;
