#ifndef CORE_H_
#define CORE_H_

#include <iostream>
#include <vector>

/// @brief 内積を計算する関数
/// @tparam U 何らかの数のクラス
/// @tparam V 何らかの数のクラス
/// @param x ベクトル
/// @param y ベクトル
/// @return 内積
template <class U, class V>
V dot(const std::vector<U> x, const std::vector<V> y);

#endif // !CORE_H_
