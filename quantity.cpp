#include "quantity.hpp"

namespace quantity
{
    // 计算系统能量
    double get_energy(const int spins[SIZE + 1])
    {
        double energy = 0;

        // 对于所有自旋
        for (int index = 0; index < SIZE; ++index)
        {
            for (int i = 0; i < neighbors::nei; ++i)
            {
                if (spins[index] == spins[neighbors::neighs[index][i]])
                {
                    energy += Parameters::J;
                }
            }
        }

        return 0.5 * energy / (1.0 * SIZE); // 返回平均每个格点的能量值
    }

    // 计算square格子构型空间的距离
    void distance_square(const int spins[SIZE + 1], double sum[4])
    {
        std::array<double, Parameters::Q> dist_total = {0.0};

        // 计算序参量
        for (int i = 0; i < Parameters::LATTICE_SIZE; ++i)
        {
            for (int j = 0; j < Parameters::LATTICE_SIZE; ++j)
            {
                int spin_val = spins[i + j * Parameters::LATTICE_SIZE];
                dist_total[spin_val]++;
            }
        }
        auto max_it = std::max_element(dist_total.begin(), dist_total.end()); // 序参量

        // 将最大元素赋值给 sum[3]，这个判断的目的是为了确保 max_it 有效，避免访问无效的迭代器或空容器的情况。
        if (max_it != dist_total.end())
        {
            sum[0] = *max_it / SIZE; // 将最大元素赋值给 sum[3]
        }

        // 计算离原点距离
        for (int i = 0; i < Parameters::Q; ++i)
        {
            for (int j = 0; j < Parameters::Q; ++j)
            {
                sum[3] += dist_total[i] * dist_total[j] * Constants::eta[i][j];
            }
        }
        sum[3] = std::sqrt(sum[3]) / SIZE;

        // 计算构型距离（立体的磁化强度）
        for (int i = 0; i < Parameters::Q; ++i)
        {
            dist_total[i] -= Constants::center[i];
        }
        for (int i = 0; i < Parameters::Q; ++i)
        {
            for (int j = 0; j < Parameters::Q; ++j)
            {
                sum[2] += dist_total[i] * dist_total[j] * Constants::eta[i][j];
            }
        }
        sum[2] = std::sqrt(sum[2]) / Constants::r_max;
        sum[1] = sum[2];
    }

    // 计算diced格子构型空间的距离
    void distance_diced(const int spins[SIZE + 1], double sum[4])
    {
        std::array<double, Parameters::Q> dist_total = {0.0};
        std::array<int, Parameters::Q> sum_Q = {0};

        // 计算序参量
        for (int i = 0; i < Parameters::LATTICE_SIZE; ++i)
        {
            for (int j = 0; j < Parameters::LATTICE_SIZE; ++j)
            {
                int spin_val = spins[i + j * Parameters::LATTICE_SIZE];
                dist_total[spin_val]++;

                if ((i % 3) == (j % 3)) // hi前的代码乘取得是diced格子上的三分之一格点计算的，所以少了3*3=9倍/这个结果要更小。
                {
                    sum_Q[spin_val]++;
                }
            }
        }
        auto max_it = std::max_element(sum_Q.begin(), sum_Q.end()); // 序参量

        // 将最大元素赋值给 sum[3]，这个判断的目的是为了确保 max_it 有效，避免访问无效的迭代器或空容器的情况。
        if (max_it != sum_Q.end())
        {
            sum[0] = *max_it * 3.0 / SIZE; // 将最大元素赋值给 sum[3]
        }

        // 计算磁化强度，这是在一个平面内的
        double mag;
        double mag_max = 0.0;
        double theta = 2 * M_PI / Parameters::Q;

        for (int i = 0; i < Parameters::Q; ++i) // 表示以某个方向为0夹角开始计算
        {
            mag = 0;
            for (int j = 0; j < Parameters::Q; ++j) // 表示计算对于某个方向的磁化强度分量
            {
                mag += cos(abs(i - j) * theta) * dist_total[j];
            }
            if (mag > mag_max)
            {
                mag_max = mag;
            }
        }

        sum[1] = mag_max / (1.0 * SIZE);

        // 计算离原点距离
        for (int i = 0; i < Parameters::Q; ++i)
        {
            for (int j = 0; j < Parameters::Q; ++j)
            {
                sum[3] += dist_total[i] * dist_total[j] * Constants::eta[i][j];
            }
        }
        sum[3] = std::sqrt(sum[3]) / SIZE;

        // 计算构型距离（立体的磁化强度）
        for (int i = 0; i < Parameters::Q; ++i)
        {
            dist_total[i] -= Constants::center[i];
        }
        for (int i = 0; i < Parameters::Q; ++i)
        {
            for (int j = 0; j < Parameters::Q; ++j)
            {
                sum[2] += dist_total[i] * dist_total[j] * Constants::eta[i][j];
            }
        }
        sum[2] = std::sqrt(sum[2]) / Constants::r_max;
    }

    // 计算diced格子构型空间的距离
    void distance_union_jack(const int spins[SIZE + 1], double sum[4])
    {
        std::array<double, Parameters::Q> dist_total = {0.0};
        std::array<int, Parameters::Q> sum_Q_A = {0};
        std::array<int, Parameters::Q> sum_Q_B = {0};

        // 计算序参量
        for (int i = 0; i < Parameters::LATTICE_SIZE; ++i)
        {
            for (int j = 0; j < Parameters::LATTICE_SIZE; ++j)
            {
                int spin_val = spins[i + j * Parameters::LATTICE_SIZE];
                dist_total[spin_val]++;

                if (((i % 2) != 0) && ((j % 2) == 0))
                {
                    sum_Q_A[spin_val]++;
                }
                else if (((i % 2) == 0) && ((j % 2) != 0))
                {
                    sum_Q_B[spin_val]++;
                }
            }
        }
        auto max_it_A = std::max_element(sum_Q_A.begin(), sum_Q_A.end()); // 序参量
        auto max_it_B = std::max_element(sum_Q_B.begin(), sum_Q_B.end()); // 序参量
        auto max_it = std::max(max_it_A, max_it_B);

        // 将最大元素赋值给 sum[3]，这个判断的目的是为了确保 max_it 有效，避免访问无效的迭代器或空容器的情况。
        if (max_it != sum_Q_A.end())
        {
            sum[0] = *max_it * 4.0 / SIZE; // 将最大元素赋值给 sum[3]
        }

        // 计算磁化强度，这是在一个平面内的
        double mag;
        double mag_max = 0.0;
        double theta = 2 * M_PI / Parameters::Q;

        for (int i = 0; i < Parameters::Q; ++i) // 表示以某个方向为0夹角开始计算
        {
            mag = 0;
            for (int j = 0; j < Parameters::Q; ++j) // 表示计算对于某个方向的磁化强度分量
            {
                mag += cos(abs(i - j) * theta) * dist_total[j];
            }
            if (mag > mag_max)
            {
                mag_max = mag;
            }
        }

        sum[1] = mag_max / (1.0 * SIZE);

        // 计算离原点距离
        for (int i = 0; i < Parameters::Q; ++i)
        {
            for (int j = 0; j < Parameters::Q; ++j)
            {
                sum[3] += dist_total[i] * dist_total[j] * Constants::eta[i][j];
            }
        }
        sum[3] = std::sqrt(sum[3]) / SIZE;

        // 计算构型距离（立体的磁化强度）
        for (int i = 0; i < Parameters::Q; ++i)
        {
            dist_total[i] -= Constants::center[i];
        }
        for (int i = 0; i < Parameters::Q; ++i)
        {
            for (int j = 0; j < Parameters::Q; ++j)
            {
                sum[2] += dist_total[i] * dist_total[j] * Constants::eta[i][j];
            }
        }
        sum[2] = std::sqrt(sum[2]) / Constants::r_max;
    }

    // // 计算欧拉示性数
    // int area(int spins[SIZE + 1])
    // {

    //     int i;                // Counter
    //     int n;                // 拿来记录区域数量
    //     int site[SIZE] = {0}; // 记录被涂色过的格点

    //     int spins_copy[SIZE];

    //     for (i = 0; i < SIZE; i++)
    //     {
    //         spins_copy[i] = spins[i];
    //     }

    //     n = 0;

    //     for (i = 0; i < SIZE; ++i)
    //     {
    //         if (site[i] == 0)
    //         {
    //             add_to_cluster(spins_copy, site, i);
    //             n++;
    //         }
    //     }

    //     return n;
    // }

    // // 计算欧拉示性数
    // void add_to_cluster(int spins_copy[SIZE], int site[SIZE], int pos)
    // {

    //     int i; // Counter
    //     int n; // Neighbour position近邻位置

    //     spins_copy[pos] += 99; // 用某种颜色替换格子
    //     site[pos]++;           // 记录已经涂色过的格子

    //     // For every neighbour,
    //     for (i = 0; i < neighbors::nei; i++)
    //     {
    //         if (neighbors::neighs[pos][i] != SIZE)
    //         {
    //             n = neighbors::neighs[pos][i];               // Get its position取近邻位置
    //             if ((spins_copy[n] + 99) == spins_copy[pos]) // Check if it the same
    //             {
    //                 add_to_cluster(spins_copy, site, n); // 注意这里把n的取值传递给pos了
    //             }
    //         }
    //     }
    // }

    // void spin_correlation(const int spins[SIZE + 1])
    // {
    //     int i, j;
    //     int co_length=Parameters::LATTICE_SIZE/3;
    //     int spin_corr[Parameters::LATTICE_SIZE];
    //     int co_spins[co_length][co_length];

    //     for (int i = 0; i < Parameters::LATTICE_SIZE; ++i)
    //     {
    //         for (int j = 0; j < Parameters::LATTICE_SIZE; ++j)
    //         {
    //             if ((i % 3) == (j % 3))
    //             {
    //                 co_spins[i / 3][j / 3] = spins[i][j];
    //             }
    //         }
    //     }

    //     for (i = 0; i < co_length; i++)
    //     {
    //         co_spins[i/3]=
    //     }

    // }
}