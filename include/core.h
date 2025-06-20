#ifndef CORE_H_
#define CORE_H_

#include <vector>

#define ERR_STR_LEN 200

/**
 * @brief ベクトルの表示
 * 
 * @tparam T 
 * @param v ベクトル
 */
template<class T>
void print(const std::vector<T> v);

/**
 * @brief 素数生成
 * 
 * @param n 正整数
 * @return long n番目の素数
 */
long prime(const long n);

/// @brief 内積を計算する関数
/// @tparam U 何らかの数のクラス
/// @tparam V 何らかの数のクラス
/// @param x ベクトル
/// @param y ベクトル
/// @return 内積
template <class U, class V>
V dot(const std::vector<U> x, const std::vector<V> y);

/**
 * @brief ベクトルが零ベクトルかどうかを判定する関数
 * 
 * @tparam T 
 * @param v ベクトル
 * @return true vが零ベクトル
 * @return false vが非零ベクトル
 */
template<class T>
bool isZero(const std::vector<T> v);

#endif // !CORE_H_
