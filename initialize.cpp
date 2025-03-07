#include "initialize.hpp"

namespace Initialization
{
    double m[TN][DATA], bins[TN][NBIN];
    int up[Parameters::Q];

    // 初始化所有变量
    void initialize()
    {
        initialize_spin_up();
        initialize_all_matrices();
    }

    // 初始化自旋状态
    void initialize_spin_up()
    {
        for (int i = 0; i < Parameters::Q; ++i)
        {
            Initialization::up[i] = (i + 1) % Parameters::Q;
        }
    }

    // 初始化物理量的矩阵
    void initialize_matrices(double *private_m, double *private_bins)
    {
        std::fill(private_m, private_m + DATA, 0.0);
        std::fill(private_bins, private_bins + NBIN, 0.0);
    }

    void initialize_all_matrices()
    {
        for (int i = 0; i < TN; ++i)
        {
            std::fill(m[i], m[i] + DATA, 0.0); // 将每行的元素设置为0
        }
        for (int i = 0; i < TN; ++i)
        {
            std::fill(bins[i], bins[i] + NBIN, 0.0); // 将每行的元素设置为0
        }
    }

    // 随机初始化自旋
    void initialize_spins(int spins[SIZE + 1], std::mt19937 &gen, std::uniform_int_distribution<int> &brandom)
    {
        for (int i = 0; i < SIZE; ++i)
        {
            spins[i] = brandom(gen);
        }
        spins[SIZE] = Parameters::Q;
    }
}

// int main()
// {
//     Initialization::initialize();
//     std::cout << "自旋状态全部向上加一后变成" << std::endl;
//     for (int i = 0; i < Parameters::Q; i++)
//     {
//         std::cout << Initialization::up[i];
//     }
//     std::cout << std::endl;
//     std::cout << "初始记录物理量的矩阵是" << std::endl;
//     for (int i = 0; i < TN; i++)
//     {
//         for (int j = 0; j < DATA; j++)
//         {
//             std::cout << Initialization::m[i][j];
//         }
//         std::cout << std::endl;
//     }
//     // std::cout << std::endl;
//     // std::cout << std::endl;
//     // std::cout << "初始构型空间" << std::endl;

//     // for (size_t j = 0; j < (SIZE + 1); j++)
//     // {
//     //     for (size_t k = 0; k < (SIZE + 1); k++)
//     //     {
//     //         std::cout << Initialization::space[j * (SIZE + 1) + k] << " "; // 例如，初始化为索引值
//     //     }
//     //     std::cout << std::endl;
//     // }
//     return 0;
// }