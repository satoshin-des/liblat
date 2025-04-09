#include "lattice.h"

#include "core.h"

template <class T>
void Lattice<T>::sizeReduce(const bool compute_gso)
{
    if (compute_gso)
    {
        computeGSO();
    }

    for (long i = 1, j; i < m_num_rows; ++i)
    {
        for (j = i - 1; j >= 0; --j)
        {
            sizeReduce(i, j);
        }
    }
}

template <class T>
void Lattice<T>::LLL(const double delta, const bool compute_gso, long start_, long end_)
{
    try
    {
        if ((delta < 0.25) || (delta > 1))
        {
            throw std::out_of_range("[WARNING]The reduction parameter must be in [0.25, 1.0].");
        }
    }
    catch (const std::exception &ex)
    {
        std::cerr << ex.what() << "@ function " << __FUNCTION__ << std::endl;
    }

    if ((end_ <= -2) || (end_ == 0))
    {
        throw std::out_of_range("The parameter end_ must be a positive integer or -1. @ function LLL");
    }
    else if (end_ == -1)
    {
        if ((start_ <= -1) || (start_ >= m_num_rows))
        {
            throw std::out_of_range("The parameter start_ is out of index. @ function LLL");
        }
    }
    else
    {
        if ((start_ <= -1) || (start_ >= m_num_rows))
        {
            throw std::out_of_range("The parameter start_ is out of index. @ function LLL");
        }
        if (start_ >= end_)
        {
            throw std::invalid_argument("The parameter start_ must be less than the parameter end_. @ function LLL");
        }
        if (end_ > m_num_rows)
        {
            throw std::out_of_range("The parameter end_ is out of index. @ function LLL");
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

        if ((k > start) && (m_B[k] < (delta - m_mu[k][k - 1] * m_mu[k][k - 1]) * m_B[k - 1]))
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
void Lattice<T>::deepLLL(const double delta, const bool compute_gso, long start_, long end_)
{
    try
    {
        if ((delta < 0.25) || (delta > 1))
        {
            throw std::out_of_range("[WARNING]The reduction parameter must be in [0.25, 1.0].");
        }
    }
    catch (const std::exception &ex)
    {
        std::cerr << ex.what() << "@ function " << __FUNCTION__ << std::endl;
    }

    if ((end_ <= -2) || (end_ == 0))
    {
        throw std::out_of_range("The parameter end_ must be a positive integer or -1. @ function deepLLL");
    }
    else if (end_ == -1)
    {
        if ((start_ <= -1) || (start_ >= m_num_rows))
        {
            throw std::out_of_range("The parameter start_ is out of index. @ function deepLLL");
        }
    }
    else
    {
        if ((start_ <= -1) || (start_ >= m_num_rows))
        {
            throw std::out_of_range("The parameter start_ is out of index. @ function deepLLL");
        }
        if (start_ >= end_)
        {
            throw std::invalid_argument("The parameter start_ must be less than the parameter end_. @ function deepLLL");
        }
        if (end_ > m_num_rows)
        {
            throw std::out_of_range("The parameter end_ is out of index. @ function deepLLL");
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

    double C;

    if (compute_gso)
    {
        computeGSO();
    }

    for (long k = start + 1, j, i, t, l; k < end;)
    {
        for (j = k - 1; j >= start; --j)
        {
            sizeReduce(k, j);
        }

        C = static_cast<double>(dot(m_basis[k], m_basis[k]));

        for (i = start; i < k;)
        {
            if (C >= delta * m_B[i])
            {
                C -= m_mu[k][i] * m_mu[k][i] * m_B[i];
                ++i;
            }
            else
            {
                deepInsertion(i, k);
                updateDeepInsGSO(i, k, start, end);

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
        if ((delta < 0.25) || (delta > 1))
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
            updateDeepInsGSO(k, l, 0, m_num_rows);
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
        if ((delta < 0.25) || (delta > 1))
        {
            throw std::out_of_range("[WARNING]The reduction parameter must be in [0.25, 1.0].");
        }
    }
    catch (const std::exception &ex)
    {
        std::cerr << ex.what() << "@ function " << __FUNCTION__ << std::endl;
    }

    if ((beta < 2) || (beta > m_num_rows))
    {
        char err_s[ERR_STR_LEN];
        sprintf(err_s, "The blocksize is %ld. The blocksize must be in [2, %ld]. @ function BKZ.", beta, m_num_rows);
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
                if ((w[i] == 1) || (w[i] == -1))
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

template <class T>
void Lattice<T>::HKZ(const double delta, const bool compute_gso)
{
    BKZ(m_num_rows, delta, compute_gso);
}

template <class T>
void Lattice<T>::deepBKZ(const long beta, const double delta, const bool compute_gso)
{
    try
    {
        if ((delta < 0.25) || (delta > 1))
        {
            throw std::out_of_range("[WARNING]The reduction parameter must be in [0.25, 1.0].");
        }
    }
    catch (const std::exception &ex)
    {
        std::cerr << ex.what() << "@ function " << __FUNCTION__ << std::endl;
    }

    if ((beta < 2) || (beta > m_num_rows))
    {
        char err_s[ERR_STR_LEN];
        sprintf(err_s, "The blocksize is %ld. The blocksize must be in [2, %ld]. @ function deepBKZ.", beta, m_num_rows);
        throw std::out_of_range(err_s);
    }

    std::vector<T> v(m_num_cols);
    std::vector<long> w(m_num_rows);

    deepLLL(delta, compute_gso);

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
                if ((w[i] == 1) || (w[i] == -1))
                {
                    for (j = 0; j < m_num_cols; ++j)
                    {
                        m_basis[i + k - 1][j] = v[j];
                    }

                    deepInsertion(k - 1, i + k - 1);
                    break;
                }
            }

            deepLLL(delta, true, 0, h);
        }
        else
        {
            ++z;
            deepLLL(delta, false, 0, h);
        }
    }
}

template <class T>
void Lattice<T>::potBKZ(const long beta, const double delta, const bool compute_gso)
{
    try
    {
        if ((delta < 0.25) || (delta > 1))
        {
            throw std::out_of_range("[WARNING]The reduction parameter must be in [0.25, 1.0].");
        }
    }
    catch (const std::exception &ex)
    {
        std::cerr << ex.what() << "@ function " << __FUNCTION__ << std::endl;
    }

    if ((beta < 2) || (beta > m_num_rows))
    {
        char err_s[ERR_STR_LEN];
        sprintf(err_s, "The blocksize is %ld. The blocksize must be in [2, %ld]. @ function deepBKZ.", beta, m_num_rows);
        throw std::out_of_range(err_s);
    }

    std::vector<long> v;
    std::vector<T> w(m_num_cols, 0);

    if (compute_gso)
    {
        computeGSO();
    }

    for (long z = 0, j = 0, i, k, l, d; z < m_num_rows - 1;)
    {
        if (j == m_num_rows - 2)
        {
            j = 0;
        }
        ++j;

        k = std::min(j + beta - 1, m_num_rows - 1);

        d = k - j + 1;
        v.resize(d);

        v = potENUM(j - 1, d);
        if (!isZero(v))
        {
            z = 0;

            for (i = 0; i < m_num_cols; ++i)
            {
                w[i] = 0;
                for (l = 0; l < d; ++l)
                {
                    w[i] += static_cast<T>(v[l]) * m_basis[l + j - 1][i];
                }
            }

            for (i = d - 1; i >= 0; --i)
            {
                if ((v[i] == 1) || (v[i] == -1))
                {
                    for (l = 0; l < m_num_cols; ++l)
                    {
                        m_basis[i + j - 1][l] = w[l];
                    }

                    deepInsertion(j, i + j - 1);
                    break;
                }
            }

            potLLL(delta);
        }
        else
        {
            ++z;
        }
    }
}

template <class T>
void Lattice<T>::dualLLL(const double delta, const bool compute_gso)
{
    try
    {
        if ((delta < 0.25) || (delta > 1))
        {
            throw std::out_of_range("[WARNING]The reduction parameter must be in [0.25, 1.0].");
        }
    }
    catch (const std::exception &ex)
    {
        std::cerr << ex.what() << "@ function " << __FUNCTION__ << std::endl;
    }

    long j, i;
    double q;
    T tmp;

    if (compute_gso)
    {
        computeGSO();
    }

    m_dual_mu[m_num_rows - 1][m_num_rows - 1] = 1.0;

    for (long k = m_num_rows - 2; k >= 0;)
    {
        m_dual_mu[k][k] = 1.0;

        for (j = k + 1; j < m_num_rows; ++j)
        {
            m_dual_mu[k][j] = 0;
            for (i = k; i < j; ++i)
            {
                m_dual_mu[k][j] -= m_mu[j][i] * m_dual_mu[k][i];
            }

            if (m_dual_mu[k][j] > 0.5 || m_dual_mu[k][j] < -0.5)
            {
                q = round(m_dual_mu[k][j]);
                for (i = 0; i < m_num_cols; ++i)
                {
                    m_basis[j][i] += static_cast<T>(q) * m_basis[k][i];
                }
                for (i = j; i < m_num_rows; ++i)
                {
                    m_dual_mu[k][i] -= q * m_dual_mu[j][i];
                }
                for (i = 0; i <= k; ++i)
                {
                    m_mu[j][i] += q * m_mu[k][i];
                }
            }
        }

        if (k < m_num_rows - 1 && (delta - m_dual_mu[k][k + 1] * m_dual_mu[k][k + 1]) * m_B[k] > m_B[k + 1])
        {
            for (i = 0; i < m_num_cols; ++i)
            {
                tmp = m_basis[k + 1][i];
                m_basis[k + 1][i] = m_basis[k][i];
                m_basis[k][i] = tmp;
            }
            updateDeepInsGSO(k, k + 1, 0, m_num_rows);
            ++k;
        }
        else
        {
            --k;
        }
    }
}

template <class T>
void Lattice<T>::dualDeepLLL(const double delta, const bool compute_gso)
{
    try
    {
        if ((delta < 0.25) || (delta > 1))
        {
            throw std::out_of_range("[WARNING]The reduction parameter must be in [0.25, 1.0].");
        }
    }
    catch (const std::exception &ex)
    {
        std::cerr << ex.what() << "@ function " << __FUNCTION__ << std::endl;
    }

    long j, i, l, h;
    double q, d, D;
    std::vector<double> dual_D(m_num_rows);

    if (compute_gso)
    {
        computeGSO();
    }

    m_dual_mu[m_num_rows - 1][m_num_rows - 1] = 1.0;

    for (long k = m_num_rows - 2; k >= 0;)
    {
        m_dual_mu[k][k] = 1.0;

        for (j = k + 1; j < m_num_rows; ++j)
        {
            m_dual_mu[k][j] = 0.0;
            for (i = k; i < j; ++i)
            {
                m_dual_mu[k][j] -= m_mu[j][i] * m_dual_mu[k][i];
            }

            if ((m_dual_mu[k][j] > 0.5) || (m_dual_mu[k][j] < -0.5))
            {
                q = round(m_dual_mu[k][j]);
                for (i = 0; i < m_num_cols; ++i)
                {
                    m_basis[j][i] += static_cast<T>(q) * m_basis[k][i];
                }
                for (i = j; i < m_num_rows; ++i)
                {
                    m_dual_mu[k][i] -= q * m_dual_mu[j][i];
                }
                for (i = 0; i <= k; ++i)
                {
                    m_mu[j][i] += q * m_mu[k][i];
                }
            }
        }

        d = 0.0;
        l = m_num_rows - 1;
        for (j = k; j < m_num_rows; ++j)
        {
            d += m_dual_mu[k][j] * m_dual_mu[k][j] / m_B[j];
        }

        while (l > k)
        {
            if (m_B[l] * d < delta)
            {
                D = 1.0 / m_B[k];

                std::fill(dual_D.begin(), dual_D.end(), 0.0);
                dual_D[k] = D;
                for (h = k + 1; h < m_num_rows; ++h)
                {
                    D += m_dual_mu[k][h] * m_dual_mu[k][h] / m_B[h];
                    dual_D[h] = D;
                }

                dualDeepInsertion(k, l);
                updateDualDeepInsGSO(k, l, dual_D);

                if (l < m_num_rows - 2)
                {
                    k = l + 1;
                }
                else
                {
                    k = m_num_rows - 1;
                }
            }
            else
            {
                d -= m_dual_mu[k][l] * m_dual_mu[k][l] / m_B[l];
                --l;
            }
        }
        --k;
    }
}

template <class T>
void Lattice<T>::dualPotLLL(const double delta, const bool compute_gso)
{
    try
    {
        if ((delta < 0.25) || (delta > 1))
        {
            throw std::out_of_range("[WARNING]The reduction parameter must be in [0.25, 1.0].");
        }
    }
    catch (const std::exception &ex)
    {
        std::cerr << ex.what() << "@ function " << __FUNCTION__ << std::endl;
    }

    double P;
    double P_min;
    double s;
    double D;
    std::vector<double> dual_D(m_num_rows);

    LLL(delta, compute_gso);

    for (long k = m_num_rows - 1, j, i, l, q, h; k >= 0;)
    {
        m_dual_mu[k][k] = 1.0;

        for (j = k + 1; j < m_num_rows; ++j)
        {
            m_dual_mu[k][j] = 0.0;
            for (i = k; i < j; ++i)
            {
                m_dual_mu[k][j] -= m_mu[j][i] * m_dual_mu[k][i];
            }

            if ((m_dual_mu[k][j] > 0.5) || (m_dual_mu[k][j] < -0.5))
            {
                q = static_cast<long>(round(m_dual_mu[k][j]));
                for (i = 0; i < m_num_cols; ++i)
                {
                    m_basis[j][i] += static_cast<T>(q) * m_basis[k][i];
                }
                for (i = j; i < m_num_rows; ++i)
                {
                    m_dual_mu[k][i] -= q * m_dual_mu[j][i];
                }
                for (i = 0; i <= k; ++i)
                {
                    m_mu[j][i] += q * m_mu[k][i];
                }
            }
        }

        P = 1.0;
        P_min = 1.0;
        l = m_num_rows - 1;
        for (j = k + 1; j < m_num_rows; ++j)
        {
            s = 0.0;
            for (i = k; i <= j; ++i)
            {
                s += m_dual_mu[k][i] * m_dual_mu[k][i] / m_B[i];
            }
            P *= m_B[j];
            P *= s;

            if (P < P_min)
            {
                l = j;
                P_min = P;
            }
        }

        if (delta > P_min)
        {
            D = 1.0 / m_B[k];
            std::fill(dual_D.begin(), dual_D.end(), 0);

            dual_D[k] = D;
            for (h = k + 1; h < m_num_rows; ++h)
            {
                D += m_dual_mu[k][h] * m_dual_mu[k][h] / m_B[h];
                dual_D[h] = D;
            }

            dualDeepInsertion(k, l);
            updateDualDeepInsGSO(k, l, dual_D);

            k = l;
        }
        else
        {
            --k;
        }
    }
}

template<class T>
void Lattice<T>::dualBKZ(const long beta, const double delta, const bool compute_gso)
{
    dualLLL(delta, compute_gso);

    long h, d;
    std::vector<long> x;

    for(long flag = 1, j; flag >= 1;)
    {
        flag = 0;
        for(j = m_num_rows - 1; j >= 1; --j)
        {
            h = std::max(j - beta, static_cast<long>(0));
        }
        d = j - h + 1;

        computeDualGSO();
        
        x = dualEnumShortVec_(false, h, j + 1);
        --x[x.size() - 1];
        if (!isZero(x))
        {
            ++x[x.size() - 1];
            ++flag;

            insertToDualBasis(x, d);
            dualLLL(0.99, true);
        }
    }
}

template class Lattice<int>;
template class Lattice<long>;
template class Lattice<long long>;
template class Lattice<float>;
template class Lattice<double>;
