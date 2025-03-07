#ifndef neighbors_hpp
#define neighbors_hpp

#include <string>
#include "constant.hpp" // 包含常数

// 声明全局变量和函数

namespace neighbors
{
    extern std::vector<std::vector<int>> neighs;
    extern int nei;

    void get_neighbors();
    void get_neighbors_square();
    void get_neighbors_triangle();
    void get_neighbors_diced();
    void get_neighbors_union_jack();
}
#endif // NEIGHBORS_H
