#ifndef LATTICE_H_
#define LATTICE_H_

#include <iostream>
#include <vector>
#include <cstdlib>
#include <sstream>
#include <random>
#include <time.h>

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

    /// @brief 内積を計算する関数
    /// @tparam U 何らかの数のクラス
    /// @tparam V 何らかの数のクラス
    /// @param x ベクトル
    /// @param y ベクトル
    /// @return 内積
    template<class U, class V>
    V dot(std::vector<U> x, std::vector<V> y)
    {
        V S = 0;
        for(int i = 0; i < x.size(); ++i)
        {
            S += static_cast<V>(x[i]) * y[i];
        }
        return S;
    }

public:
    friend std::ostream &operator<<(std::ostream &os, const Lattice<T> &lat)
    {
        for (long i = 0, j; i < lat.m_num_rows; ++i)
        {
            for (j = 0; j < lat.m_num_cols; ++j)
            {
                os << lat.m_basis[i][j] << " ";
            }
            os << std::endl;
        }
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

    /// @brief GSO情報の計算
    void computeGSO();
};

#endif // !LATTICE_H_