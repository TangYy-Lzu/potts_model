// neighbors.cpp
#include "neighbors.hpp"

namespace neighbors
{
    std::vector<std::vector<int>> neighs; // 定义在命名空间中的变量
    int nei;

    // 设置邻居关系
    void get_neighbors()
    {
        if (Parameters::LATTICE_TYPE == "square")
        {
            get_neighbors_square();
        }
        else if (Parameters::LATTICE_TYPE == "triangle")
        {
            get_neighbors_triangle();
        }
        else if (Parameters::LATTICE_TYPE == "diced")
        {
            get_neighbors_diced();
        }
        else if (Parameters::LATTICE_TYPE == "union_jack")
        {
            get_neighbors_union_jack();
        }
        else
        {
            // 使用 std::cerr 输出错误信息
            std::cerr << "Error: No such lattice type \"" << Parameters::LATTICE_TYPE << "\"." << std::endl;

            // 终止程序
            std::exit(EXIT_FAILURE);
        }
    }

    void get_neighbors_square()
    {
        // 假设 square 网格每个点有 4 个邻居
        nei = 4;
        neighs.resize(SIZE, std::vector<int>(4, -1)); // 初始化为 -1 表示没有邻居

        int UP = 0, DOWN = 1, LEFT = 2, RIGHT = 3;

        for (int i = 0; i < Parameters::LATTICE_SIZE; ++i)
        {
            for (int j = 0; j < Parameters::LATTICE_SIZE; ++j)
            {
                int index = i + j * Parameters::LATTICE_SIZE;

                // 计算周期边界条件下的邻居
                int u = (j + 1) % Parameters::LATTICE_SIZE;
                int d = (j - 1 + Parameters::LATTICE_SIZE) % Parameters::LATTICE_SIZE;
                int r = (i + 1) % Parameters::LATTICE_SIZE;
                int l = (i - 1 + Parameters::LATTICE_SIZE) % Parameters::LATTICE_SIZE;

                // 设置邻居
                neighs[index][UP] = i + u * Parameters::LATTICE_SIZE;
                neighs[index][DOWN] = i + d * Parameters::LATTICE_SIZE;
                neighs[index][LEFT] = l + j * Parameters::LATTICE_SIZE;
                neighs[index][RIGHT] = r + j * Parameters::LATTICE_SIZE;
            }
        }
    }

    void get_neighbors_triangle()
    {
        // 假设 diced 网格每个点有 6 个邻居
        nei = 6;
        neighs.resize(SIZE, std::vector<int>(6, -1)); // 初始化为 -1 表示没有邻居

        int UP = 0, DOWN = 1, LEFT = 2, RIGHT = 3, UP_LEFT = 4, DOWN_RIGHT = 5;

        for (int i = 0; i < Parameters::LATTICE_SIZE; ++i)
        {
            for (int j = 0; j < Parameters::LATTICE_SIZE; ++j)
            {
                int index = i + j * Parameters::LATTICE_SIZE;

                // 计算周期边界条件下的邻居
                int u = (j + 1) % Parameters::LATTICE_SIZE;
                int d = (j - 1 + Parameters::LATTICE_SIZE) % Parameters::LATTICE_SIZE;
                int r = (i + 1) % Parameters::LATTICE_SIZE;
                int l = (i - 1 + Parameters::LATTICE_SIZE) % Parameters::LATTICE_SIZE;

                // 设置邻居

                neighs[index][UP] = i + u * Parameters::LATTICE_SIZE;
                neighs[index][DOWN] = i + d * Parameters::LATTICE_SIZE;
                neighs[index][LEFT] = l + j * Parameters::LATTICE_SIZE;
                neighs[index][RIGHT] = r + j * Parameters::LATTICE_SIZE;
                neighs[index][UP_LEFT] = l + u * Parameters::LATTICE_SIZE;
                neighs[index][DOWN_RIGHT] = r + d * Parameters::LATTICE_SIZE;
            }
        }
    }

    void get_neighbors_diced()
    {
        // 假设 diced 网格每个点有 6 个邻居
        nei = 6;
        neighs.resize(SIZE, std::vector<int>(6, -1)); // 初始化为 -1 表示没有邻居

        int UP = 0, DOWN = 1, LEFT = 2, RIGHT = 3, UP_LEFT = 4, DOWN_RIGHT = 5;

        for (int i = 0; i < Parameters::LATTICE_SIZE; ++i)
        {
            for (int j = 0; j < Parameters::LATTICE_SIZE; ++j)
            {
                int index = i + j * Parameters::LATTICE_SIZE;

                // 计算周期边界条件下的邻居
                int u = (j + 1) % Parameters::LATTICE_SIZE;
                int d = (j - 1 + Parameters::LATTICE_SIZE) % Parameters::LATTICE_SIZE;
                int r = (i + 1) % Parameters::LATTICE_SIZE;
                int l = (i - 1 + Parameters::LATTICE_SIZE) % Parameters::LATTICE_SIZE;

                // 设置邻居
                if ((i % 3) == (j % 3))
                {
                    neighs[index][UP] = i + u * Parameters::LATTICE_SIZE;
                    neighs[index][DOWN] = i + d * Parameters::LATTICE_SIZE;
                    neighs[index][LEFT] = l + j * Parameters::LATTICE_SIZE;
                    neighs[index][RIGHT] = r + j * Parameters::LATTICE_SIZE;
                    neighs[index][UP_LEFT] = l + u * Parameters::LATTICE_SIZE;
                    neighs[index][DOWN_RIGHT] = r + d * Parameters::LATTICE_SIZE;
                }
                else if (((i + 2) % 3) == (j % 3))
                {
                    neighs[index][UP] = i + u * Parameters::LATTICE_SIZE;
                    neighs[index][DOWN] = SIZE;
                    neighs[index][LEFT] = l + j * Parameters::LATTICE_SIZE;
                    neighs[index][RIGHT] = SIZE;
                    neighs[index][UP_LEFT] = SIZE;
                    neighs[index][DOWN_RIGHT] = r + d * Parameters::LATTICE_SIZE;
                }
                else
                {
                    neighs[index][UP] = SIZE;
                    neighs[index][DOWN] = i + d * Parameters::LATTICE_SIZE;
                    neighs[index][LEFT] = SIZE;
                    neighs[index][RIGHT] = r + j * Parameters::LATTICE_SIZE;
                    neighs[index][UP_LEFT] = l + u * Parameters::LATTICE_SIZE;
                    neighs[index][DOWN_RIGHT] = SIZE;
                }
            }
        }
    }

    void get_neighbors_union_jack()
    {
        // 假设 union_jack 网格每个点有 8 个邻居
        nei = 8;
        neighs.resize(SIZE, std::vector<int>(8, -1)); // 初始化为 -1 表示没有邻居

        int UP = 0, DOWN = 1, LEFT = 2, RIGHT = 3, UP_LEFT = 4, UP_RIGHT = 5, DOWN_RIGHT = 6, DOWN_LEFT = 7;

        for (int i = 0; i < Parameters::LATTICE_SIZE; ++i)
        {
            for (int j = 0; j < Parameters::LATTICE_SIZE; ++j)
            {
                int index = i + j * Parameters::LATTICE_SIZE;

                // 计算周期边界条件下的邻居
                int u = (j + 1) % Parameters::LATTICE_SIZE;
                int d = (j - 1 + Parameters::LATTICE_SIZE) % Parameters::LATTICE_SIZE;
                int r = (i + 1) % Parameters::LATTICE_SIZE;
                int l = (i - 1 + Parameters::LATTICE_SIZE) % Parameters::LATTICE_SIZE;

                // 设置邻居
                if ((i + j) % 2 == 0)
                {
                    neighs[index][UP] = i + u * Parameters::LATTICE_SIZE;
                    neighs[index][DOWN] = i + d * Parameters::LATTICE_SIZE;
                    neighs[index][LEFT] = l + j * Parameters::LATTICE_SIZE;
                    neighs[index][RIGHT] = r + j * Parameters::LATTICE_SIZE;
                    neighs[index][UP_LEFT] = SIZE;
                    neighs[index][UP_RIGHT] = SIZE;
                    neighs[index][DOWN_RIGHT] = SIZE;
                    neighs[index][DOWN_LEFT] = SIZE;
                }
                else
                {
                    neighs[index][UP] = i + u * Parameters::LATTICE_SIZE;
                    neighs[index][DOWN] = i + d * Parameters::LATTICE_SIZE;
                    neighs[index][LEFT] = l + j * Parameters::LATTICE_SIZE;
                    neighs[index][RIGHT] = r + j * Parameters::LATTICE_SIZE;
                    neighs[index][UP_LEFT] = l + u * Parameters::LATTICE_SIZE;
                    neighs[index][UP_RIGHT] = r + u * Parameters::LATTICE_SIZE;
                    neighs[index][DOWN_RIGHT] = r + d * Parameters::LATTICE_SIZE;
                    neighs[index][DOWN_LEFT] = l + d * Parameters::LATTICE_SIZE;
                }
            }
        }
    }
}

// int main()
// {
//     neighbors::get_neighbors();
//     std::cout << neighbors::nei << std::endl;
//     // 遍历所有的行
//     for (size_t i = 0; i < neighbors::neighs.size(); ++i) // neighs.size()就是SIZE
//     {
//         std::cout << "Index " << i << ": ";
//         // 遍历每个元素
//         for (size_t j = 0; j < neighbors::neighs[i].size(); ++j) // neighs[i].size()就是nei
//         {
//             std::cout << neighbors::neighs[i][j] << " ";
//         }
//         std::cout << std::endl;
//     }

//     return 0;
// }
