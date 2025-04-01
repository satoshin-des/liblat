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
    long m_max_loop = 10;                                                  // BKZの最大ループ数
    std::vector<std::vector<T>> m_basis;                                   // 格子基底
    std::vector<double> m_B;                                               // GSOベクトルの二乗ノルム
    std::vector<double> m_dual_B;                                          // 双対基底のGSOベクトルの二乗ノルム
    std::vector<std::vector<double>> m_b_star;                             // GSOベクトル
    std::vector<std::vector<double>> m_dual_b_star;                        // 双対基底のGSOベクトル
    std::vector<std::vector<double>> m_mu;                                 // GSO係数行列
    std::vector<std::vector<double>> m_dual_mu;                            // 双対基底のGSO係数行列
    std::mt19937_64 m_mt64 = std::mt19937_64((unsigned int)time(nullptr)); // 乱数生成機
    std::uniform_real_distribution<> m_get_rand_uni;                       // 一様ランダムに生成

    /**
     * @brief 局所射影ブロック格子上のノルムがR未満である格子ベクトルの数え上げ
     *
     * @param R 数え上げ半径
     * @param start
     * @param end
     * @return std::vector<long>
     */
    std::vector<long> ENUM_(const double R, const long start, const long end);

    std::vector<long> enumShortVec_(const bool compute_gso, const long start, const long end);

public:
    /**
     * @brief Construct a new Lattice object
     *
     * @param n 基底行列の行数
     * @param m 基底行列の列数
     */
    Lattice(const long n = 0, const long m = 0) : m_num_rows(n), m_num_cols(m), m_basis(std::vector<std::vector<T>>(n, std::vector<T>(m, 0))), m_B(std::vector<double>(n, 0)), m_b_star(std::vector<std::vector<double>>(n, std::vector<double>(m, 0))), m_mu(std::vector<std::vector<double>>(n, std::vector<double>(n, 0))), m_dual_B(std::vector<double>(n, 0)), m_dual_b_star(std::vector<std::vector<double>>(n, std::vector<double>(m, 0))), m_dual_mu(std::vector<std::vector<double>>(n, std::vector<double>(n, 0))) {};

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

    /**
     * @brief Set the Max Loop object
     *
     * @param max_loop BKZのループ回数の上限の設定
     */
    void setMaxLoop(const long max_loop);

    /**
     * @brief 格子基底行列の行数
     *
     * @return long 行数
     */
    long numRows() const;

    /**
     * @brief 格子基底行列の列数
     *
     * @return long 列数
     */
    long numCols() const;

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
     * @brief 基底行列の設定
     *
     * @param basis_mat 基底行列
     */
    void setBasis(const std::vector<std::vector<T>> basis_mat);

    /**
     * @brief Goldestei-Mayer格子に設定
     * 
     * @param p 
     * @param q 
     */
    void setGoldesteinMayerLattice(const T p, const T q);

    /**
     * @brief Schnorr格子に設定
     * 
     * 
     * @param N 合成数
     * @param c 精度パラメータ
     */
    void setSchnorrLattice(const long N, const double c);

    /**
     * @brief 係数ベクトルと基底の積
     *
     * @param v
     * @return std::vector<T>
     */
    std::vector<T> mulVecBasis(const std::vector<long> v);

    /**
     * @brief deep insertion
     *
     * @param k
     * @param l
     */
    void deepInsertion(const long k, const long l);

    /**
     * @brief dual deep insertion
     *
     * @param k
     * @param l
     */
    void dualDeepInsertion(const long k, const long l);

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
    void updateDeepInsGSO(const long i, const long k, const long start, const long end);

    /**
     * @brief 双対型deep insertion後のGSO情報の効率的な更新
     *
     * @param k
     * @param l
     * @param dual_D
     */
    void updateDualDeepInsGSO(const long k, const long l, const std::vector<double> dual_D);

    /**
     * @brief 部分サイズ基底簡約
     *
     * Y. Aono, M. Yasuda. 格子暗号解読のための数学的基礎.(2019)
     *
     * @param i
     * @param j
     */
    void sizeReduce(const long i, const long j);

    /**
     * @brief サイズ基底簡約
     *
     * Y. Aono, M. Yasuda. 格子暗号解読のための数学的基礎.(2019)
     *
     * @param compute_gso サイズ基底簡約前にGSOを計算するか
     */
    void sizeReduce(const bool compute_gso = true);

    /**
     * @brief LLL簡約
     *
     * A. K. Lenstra, H. W. Lenstra, L. Lovasz. Factoring polynomials with rational coefficients. 1982
     *
     * @param delta 簡約パラメタ
     * @param compute_gso LLL前にGSOを計算するか
     */
    void LLL(const double delta = 0.75, const bool compute_gso = true, const long start_ = 0, const long end_ = -1);

    /**
     * @brief DeepLLL簡約
     *
     * C. P. Schnorr, M. Euchner. Lattice basis reduction: Improved practical algorithms and solving subset sum problem.(1994)
     *
     * @param delta 簡約パラメタ
     * @param compute_gso DeepLLL前にGSOを計算するか
     */
    void deepLLL(const double delta = 0.75, const bool compute_gso = true, const long start_ = 0, const long end_ = -1);

    /**
     * @brief PotLLL簡約
     *
     * F. Fontein, M. Schneider, U. Wagner. PotLLL: A polynomial time version of LLL with deep insertions.(2014)
     *
     * @param delta 簡約パラメタ
     * @param compute_gso PotLLL前にGSOを計算するか
     */
    void potLLL(const double delta = 0.75, const bool compute_gso = true);

    /**
     * @brief 二乗ノルムがR以下であるような格子ベクトルの数え上げ
     *
     * N. Gama, P. Q. Nguyen, O. Regev. Lattice enumeration using extreme pruning.(2010)
     *
     * @param R 探索半径
     * @return std::vector<long> 格子ベクトル
     */
    std::vector<long> ENUM(const double R);

    /**
     * @brief 格子上の最短ベクトルの数え上げ
     *
     * @param compute_gso 数え上げ前にGSOを計算するかどうか
     * @return std::vector<T> 最短ベクトル
     */
    std::vector<T> enumShortVec(const bool compute_gso = true);

    /**
     * @brief BKZ簡約アルゴリズム
     *
     * C. P. Schnorr, M. Euchner. Lattice basis reduction: Improved practical algorithms and solving subset sum problem.(1994)
     *
     * @param beta ブロックサイズ
     * @param delta 簡約パラメタ
     * @param compute_gso BKZ前にGSO情報を更新するか
     */
    void BKZ(const long beta, const double delta = 0.75, const bool compute_gso = true);

    /**
     * @brief HKZ簡約
     *
     * C. P. Schnorr, M. Euchner. Lattice basis reduction: Improved practical algorithms and solving subset sum problem.(1994)
     *
     * @param delta 簡約パラメタ
     * @param compute_gso HKZ前にGSO情報を更新するか
     */
    void HKZ(const double delta = 0.75, const bool compute_gso = true);

    /**
     * @brief DeepBKZ簡約アルゴリズム
     *
     * J. Yamaguchi, M. Yasuda. Explicit formula for Gram-Schmidt vectors in LLL with deep insertions and its applications.(2017)
     *
     * @param beta ブロックサイズ
     * @param delta 簡約パラメタ
     * @param compute_gso DeepBKZ前にGSO情報を更新するか
     */
    void deepBKZ(const long beta, const double delta = 0.75, const bool compute_gso = true);

    /**
     * @brief PotENUMアルゴリズム
     *
     * A. Sato, M. Yasuda. 自己双対型PotBKZ基底簡約の提案とBKZとの比較.(2025)
     *
     * @param start
     * @param n 局所射影ブロック格子の次元
     * @return std::vector<long>
     */
    std::vector<long> potENUM(const long start, const long n);

    /**
     * @brief PotBKZ簡約アルゴリズム
     *
     * A. Sato, M. Yasuda. 自己双対型PotBKZ基底簡約の提案とBKZとの比較.(2025)
     *
     * @param beta ブロックサイズ
     * @param delta 簡約パラメタ
     * @param compute_gso
     */
    void potBKZ(const long beta, const double delta = 0.75, const bool compute_gso = true);

    /**
     * @brief 双対型DeepLLL簡約アルゴリズム
     *
     * M. Yasuda, J. Yamaguchi, M. Ooka, S. Nakamura. Development of a dual version of DeepBKZ and its application to solving the LWE challenge.(2018)
     *
     * @param delta 簡約パラメタ
     * @param compute_gso 双対型DeepLLL前にGSO情報を更新するか
     */
    void dualDeepLLL(const double delta = 0.75, const bool compute_gso = true);

    /**
     * @brief 双対型PotLLL簡約アルゴリズム
     *
     * A. Sato, M. Yasuda. 自己双対型PotBKZ基底簡約の提案とBKZとの比較.(2025)
     *
     * @param delta 簡約パラメタ
     * @param compute_gso 双対型DeepLLL前にGSO情報を更新するか
     */
    void dualPotLLL(const double delta = 0.75, const bool compute_gso = true);

    void insertToDualBasis(const std::vector<double> x);
};

#endif // !LATTICE_H_