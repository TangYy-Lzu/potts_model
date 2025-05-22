#ifndef matrix_hpp
#define matrix_hpp

#include "constant.hpp"
#include <Eigen/Dense>
#include <Eigen/SVD> // 包含SVD相关操作的头文件

namespace matrix
{
    // 函数声明
    void Single_Value(const int spins[SIZE + 1], double sv[3]); // 计算奇异值之和
}

#endif // matrix_hpp
