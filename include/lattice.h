#ifndef LATTICE_H_
#define LATTICE_H_

#include <time.h>

#include <iostream>
#include <vector>
#include <cstdlib>
#include <sstream>
#include <random>

/**
 * @brief 格子のクラス
 *
 * @tparam T int, long, long long, float, double
 */
template <class T>
class Lattice
{
private:
    bool m_gso_is_computed = false;
    long m_num_rows = 0;
    long m_num_cols = 0;
    std::vector<std::vector<T>> m_basis;
    std::vector<double> m_B;
    std::vector<std::vector<double>> m_b_star;
    std::vector<std::vector<double>> m_mu;
    std::mt19937_64 m_mt64 = std::mt19937_64((unsigned int)time(nullptr));
    std::uniform_real_distribution<> m_get_rand_uni;

public:
    friend std::ostream &operator<<(std::ostream &os, const Lattice<T> &lat)
    {
        os << "[" << std::endl;
        for (std::vector<T> bb : lat.m_basis)
        {
            os << "[";
            for (T b : bb)
            {
                os << b << ", ";
            }
            os << "\b\b]" << std::endl;
        }
        os << "]" << std::endl;
        return os;
    }

    /// @brief 格子のサイズの設定
    /// @param n 格子次元
    /// @param m 格子ベクトルのサイズ
    void setDims(const long n, const long m);

    /// @brief ランダムなSVP-challenge型格子基底の生成
    /// @param n 格子次元
    /// @param m 格子ベクトルのサイズ
    /// @param min 下限
    /// @param max 上限
    void setRandom(const long n, const long m, const T min, const T max);

    /**
     * @brief GSO情報の計算
     *
     */
    void computeGSO();

    /**
     * @brief deep insertion後のGSO情報の効率的な更新
     * 
     * @param i 
     * @param k 
     */
    void updateDeepInsGSO(const int i, const int k);

    /**
     * @brief 部分サイズ基底簡約
     *
     * @param i
     * @param j
     */
    void sizeReduce(const int i, const int j);

    /**
     * @brief LLL簡約
     * 
     * @param delta 簡約パラメタ
     * @param compute_gso LLL前にGSOを計算するか
     */
    void LLL(const double delta, const bool compute_gso = true);

    void deepLLL(const double delta, const bool compute_gso = true);
};

#endif // !LATTICE_H_