#ifndef initialize_hpp
#define initialize_hpp

#include "constant.hpp"

// 使用 constexpr 代替宏定义
constexpr int DATA = 21;
constexpr int NBIN = 53;

// 使用 enum class 提升类型安全性
enum class Quantity
{
    // 物理量
    MAG_PLANE,
    MAG_PLANE2,
    MAG,
    MAG2,
    ENE,
    ENE2,
    ORDER,
    ORDER2,
    DIS,
    DIS2,
    AREA,
    AREA2,
    EULER,
    EULER2,
    SINGLE_VALUE,
    SINGLE_VALUE2,
    NUCLEAR_NORM,
    NUCLEAR_NORM2,
    SHANNON,
    SHANNON2,
    MAG4,

    // 拿来计算bins中的涨落平均值
    MAG_PLANE2_bin,
    MAG2_bin,
    ENE2_bin,
    ORDER2_bin,
    DIS2_bin,
    AREA2_bin,
    EULER2_bin,
    SINGLE_VALUE2_bin,
    NUCLEAR_NORM2_bin,
    SHANNON2_bin,
    MAG4_bin,

    // 拿来计算误差
    MAG_PLANERR,
    MAG_PLANE2RR,
    MAGRR,
    MAG2RR,
    ENERR,
    ENE2RR,
    ORDERRR,
    ORDER2RR,
    DISRR,
    DIS2RR,
    AREARR,
    AREA2RR,
    EULERRR,
    EULER2RR,
    SINGLE_VALUERR,
    SINGLE_VALUE2RR,
    NUCLEAR_NORMRR,
    NUCLEAR_NORM2RR,
    SHANNONRR,
    SHANNON2RR,
    MAG4RR
};

namespace Initialization
{
    extern double m[TN][DATA], bins[TN][NBIN];
    extern int up[Parameters::Q];

    // 函数声明
    void initialize();
    void initialize_spin_up();
    void initialize_matrices(double *private_m, double *private_bins);
    void initialize_all_matrices();
    void initialize_spins(int spins[SIZE + 1], std::mt19937 &gen, std::uniform_int_distribution<int> &brandom);
}

// 不再需要重复声明外部变量，因为已经在 Initialization 命名空间中声明过

#endif // INITIALIZE_H_INCLUDED
