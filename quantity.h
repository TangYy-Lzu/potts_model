#ifndef quantity_h
#define quantity_h

#include "neighbors.h"
#include "write.h"

double get_energy(int spins[SIZE + 1], int neighs[SIZE][nei]);
double magnetization(int spins[SIZE + 1]);
double distance(int spins[SIZE + 1], double eta[q][q], double center[Q], double r_max, int count);

// int main()
// {
//     int spins[SIZE];
//     double energy = get_energy(spins, neighs);
//     cout << energy;
// }

// Computes the energy of the system
double get_energy(int spins[SIZE + 1], int neighs[SIZE][nei])
{
    int index, i; // Counters，计数
    int sum_neigh;

    int energy = 0; // Sum，总和

    // For every spin,对于所有自旋
    for (index = 0; index < SIZE; index++)
    {
        for (i = 0; i < nei; i++)
        {
            if (spins[index] == spins[neighs[index][i]])
            {
                energy -= J;
            }
        }
    }

    return 0.5 * energy / (1.0 * SIZE); // Return the energy，返回平均每个格点能量值。这里没有必要再乘以二了
}

// Compute the magnetization，计算磁化强度
double magnetization(int spins[SIZE + 1])
{
    if (lattice == diced)
    {
        int i;
        int j;

        double sum = 0.0;
        double sum_0 = 0.0;
        double sum_1 = 0.0;
        double sum_2 = 0.0;
        // Sum all the values of the spins，计算自旋的总值，计算三号格点的自旋为0，为1，为2的数量并且选出最大的那个。
        for (i = 0; i < L; i++) // i列
        {
            for (j = 0; j < L; j++) // j行
            {
                if ((i % 3) == (j % 3))
                {
                    if (spins[i + j * L] == 0)
                    {
                        sum_0++;
                    }
                    else if (spins[i + j * L] == 1)
                    {
                        sum_1++;
                    }
                    else
                    {
                        sum_2++;
                    }
                }
            }
        }
        sum = max(sum_0, sum_1);
        sum = max(sum, sum_2);

        return sum / (SIZE); // 因为自旋实际上是-1和1，而这里用0表示了，磁化强度定义为媒质微小体元ΔV内的全部分子磁矩矢量和与ΔV之比，实际上这里只计算磁矩就行，不考虑轨道角动量，磁矩与总自旋之间相差一个因子（参胡安p239）
    }
    else
    {
        int i;

        double sum = 0.0;
        // Sum all the values of the spins，计算自旋的总值
        for (i = 0; i < SIZE; i++)
        {
            sum += spins[i];
        }
        // And then return them
        return (2.0 * sum - SIZE) / (SIZE);
        // 因为自旋实际上是-1和1，而这里用0表示了，磁化强度定义为媒质微小体元ΔV内的全部分子磁矩矢量和与ΔV之比，实际上这里只计算磁矩就行，不考虑轨道角动量，磁矩与总自旋之间相差一个因子（参胡安p239）
    }
}

// Compute the magnetization，计算构型空间的距离
double distance(int spins[SIZE + 1], double eta[q][q], double center[Q], double r_max, int count)
{
    int i, j;             // 拿来计数
    int site;             // 拿这个计算对应构型在space中的编号
    double r = 0.0;       // 计算总距离
    double dist[Q] = {0}; // 存Q个值，但是我们只在切面看，只取前Q-1个值。最后一个值作为零点

    for (i = 0; i < SIZE; i++) // 这屋里计算出的时里原点的距离而不是就中心的距离，要减去中心坐标。
    {
        dist[spins[i]]++; // 还是计算出了三个方向的距离，只是第三个方向我们不用
    }

    site = 0;
    for (i = 0; i < q; i++)
    {
        int d_site = static_cast<int>(dist[i]); // 在第i个维度上的坐标
        site += d_site * pow((SIZE + 1), i);    // 第i个维度就要先加上(SIZE + 1)^i*d_site
    }
    space[count][site]++;

    for (i = 0; i < Q; i++) // 减去中心坐标得到了真正各个构型距离中心的距离，这里得到
    {
        dist[i] -= center[i];
    }

    for (i = 0; i < q; i++) // 计算此构型的距离
    {
        for (j = 0; j < q; j++)
        {
            r += dist[i] * dist[j] * eta[i][j];
        }
    }

    r = sqrt(r);
    r /= r_max;

    return r;
}

#endif