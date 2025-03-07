#ifndef quantity_hpp
#define quantity_hpp

#include "neighbors.hpp"
#include "neighbors.cpp"
#include "initialize.hpp"
#include <array>
#include <algorithm> // std::max_element

namespace quantity
{
    // 函数声明
    double get_energy(const int spins[SIZE + 1]);
    void distance_square(const int spins[SIZE + 1], double sum[4]);
    void distance_diced(const int spins[SIZE + 1], double sum[4]);
    void distance_union_jack(const int spins[SIZE + 1], double sum[4]);
    // int area(int spins[SIZE + 1]);
    // void add_to_cluster(int spins_copy[SIZE], int site[SIZE], int pos);
}

#endif // QUANTITY_H
