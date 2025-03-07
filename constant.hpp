#ifndef constant_hpp
#define constant_hpp

#include "_parameter.hpp" // 确保路径正确
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <random>
#include <fstream>

// 使用 constexpr 代替宏定义
constexpr int q = (Parameters::Q - 1);
constexpr int SIZE = (Parameters::LATTICE_SIZE * Parameters::LATTICE_SIZE);

constexpr int TN = ((Parameters::T_MAX - Parameters::T_CRIT_UP + Parameters::T_CRIT_DOWN - Parameters::T_MIN + 0.001) / Parameters::DELTA_T + (Parameters::T_CRIT_UP - Parameters::T_CRIT_DOWN + 0.001) / Parameters::DELTA_T_CRIT + 1);

namespace Constants
{
    // 全局常量声明
    extern double doubleprecision;
    extern double tmincompare;
    extern double tcrit_upCompare;
    extern double tcrit_downCompare;
    extern double eta[Parameters::Q][Parameters::Q];
    extern double center[Parameters::Q];
    extern double r_max;
    extern double temperature[TN];

    // 函数声明
    void constant();
    void c_eta();
    void c_center();
    void distance_max();
    void get_temperature();
}

// 如果函数仅在常量命名空间中使用，考虑将它们移到 Constants 命名空间内

#endif // CONSTANT_H_INCLUDED
