/**
 * @file lattice.h
 * @author Arata Sato
 * @brief
 * @version 0.1
 * @date 2025-03-16
 *
 * @copyright Copyright (c) 2025
 *
 */

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
    long m_num_rows = 0;                                                   // 格子の次元
    long m_num_cols = 0;                                                   // 格子ベクトルのサイズ
    std::vector<std::vector<T>> m_basis;                                   // 格子基底
    std::vector<double> m_B;                                               // GSOベクトルの二乗ノルム
    std::vector<std::vector<double>> m_b_star;                             // GSOベクトル
    std::vector<std::vector<double>> m_mu;                                 // GSO係数行列
    std::mt19937_64 m_mt64 = std::mt19937_64((unsigned int)time(nullptr)); // 乱数生成機
    std::uniform_real_distribution<> m_get_rand_uni;                       // 一様ランダムに生成

public:
    /**
     * @brief Construct a new Lattice object
     * 
     * @param n 基底行列の行数
     * @param m 基底行列の列数
     */
    Lattice(const long n = 0, const long m = 0) : m_num_rows(n), m_num_cols(m), m_basis(std::vector<std::vector<T>>(n, std::vector<T>(m, 0))), m_B(std::vector<double>(n, 0)), m_b_star(std::vector<std::vector<double>>(n, std::vector<double>(m, 0))), m_mu(std::vector<std::vector<double>>(n, std::vector<double>(n, 0))) {};

    /**
     * @brief ストリームへの出力
     * 
     * @param os ストリーム
     * @param lat 格子
     * @return std::ostream& 出力ストリーム
     */
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
     * @brief deep insertion
     *
     * @param k
     * @param l
     */
    void deepInsertion(const long k, const long l);

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
    void updateDeepInsGSO(const long i, const long k);

    /**
     * @brief 部分サイズ基底簡約
     *
     * @param i
     * @param j
     */
    void sizeReduce(const long i, const long j);

    /**
     * @brief LLL簡約
     *
     * @param delta 簡約パラメタ
     * @param compute_gso LLL前にGSOを計算するか
     * @cite A. K. Lenstra, H. W. Lenstra, L. Lovasz. Factoring polynomials with rational coefficients. 1982
     */
    void LLL(const double delta, const bool compute_gso = true);

    /**
     * @brief DeepLLL簡約
     *
     * @param delta 簡約パラメタ
     * @param compute_gso DeepLLL前にGSOを計算するか
     * @cite C. P. Schnorr, M. Euchner. Lattice basis reduction: Improved practical algorithms and solving subset sum problem.(1994)
     */
    void deepLLL(const double delta, const bool compute_gso = true);

    /**
     * @brief PotLLL簡約
     *
     * @param delta 簡約パラメタ
     * @param compute_gso PotLLL前にGSOを計算するか
     * @cite F. Fontein, M. Schneider, U. Wagner. PotLLL: A polynomial time version of LLL with deep insertions.(2014)
     */
    void potLLL(const double delta, const bool compute_gso = true);

    /**
     * @brief 二乗ノルムがR以下であるような格子ベクトルの数え上げ
     * 
     * @param R 探索半径
     * @return std::vector<long> 格子ベクトル
     */
    std::vector<long> ENUM(const double R);

    std::vector<long> enumShortVec(const bool compute_gso = true);
};

#endif // !LATTICE_H_