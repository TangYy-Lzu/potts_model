#ifndef quantity_h
#define quantity_h

#include "neighbors.h"
#include "write.h"

double get_energy(int spins[SIZE + 1], int neighs[SIZE][nei]);
void magnetization(int spins[SIZE + 1], double dist_total[Q], double sum[2]);
double distance(int spins[SIZE + 1], double eta[q][q], double center[Q], double r_max, int count, double dist_total[Q]);

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
void magnetization(int spins[SIZE + 1], double dist_total[Q], double sum[2]) // sum存序参量和磁化强度
{
    int i, j;

    // 计算序参量，sum[0]
    double sum_Q[Q] = {0.0};
    double sum_max; // 用来选出最大的那个

    for (i = 0; i < L; i++) // i列(如第0行第i个)，注意是列
    {
        for (j = 0; j < L; j++) // j行
        {
            dist_total[spins[i + j * L]]++; // 还是计算出了三个方向的距离，只是第三个方向我们不用

            if ((i % 3) == (j % 3)) // 计算序参量
            {
                sum_Q[spins[i + j * L]]++;
            }
        }
    }

    for (i = 0; i < Q; i++)
    {
        sum[0] = max(sum[0], sum_Q[i]);
    }

    sum[0] /= (1.0 * SIZE);

    // 计算磁化强度，sum[1]
    double mag;
    double mag_max = 0.0;
    double theta = 2 * M_PI / Q;

    for (i = 0; i < Q; i++) // 表示以某个方向为0夹角开始计算
    {
        mag = 0;
        for (j = 0; j < Q; j++) // 表示计算对于某个方向的磁化强度分量
        {
            mag += cos(abs(i - j) * theta) * dist_total[j];
        }
        if (mag > mag_max)
        {
            mag_max = mag;
        }
    }

    sum[1] = mag_max / (1.0 * SIZE);

    // And then return them
    return;
    // 因为自旋实际上是-1和1，而这里用0表示了，磁化强度定义为媒质微小体元ΔV内的全部分子磁矩矢量和与ΔV之比，实际上这里只计算磁矩就行，不考虑轨道角动量，磁矩与总自旋之间相差一个因子（参胡安p239）
}

// Compute the magnetization，计算构型空间的距离
double distance(int spins[SIZE + 1], double eta[q][q], double center[Q], double r_max, int count, double dist_total[Q])
{
    int i, j;       // 拿来计数
    int site;       // 拿这个计算对应构型在space中的编号
    double r = 0.0; // 计算总距离

    if (WRITE)
    {
        site = 0;
        for (i = 0; i < q; i++)
        {
            int d_site = static_cast<int>(dist_total[i]); // 在第i个维度上的坐标,将一个 dist[i] 类型转换为 int 类型
            site += d_site * pow((SIZE + 1), i);          // 第i个维度就要先加上(SIZE + 1)^i*d_site
        }
        space[count][site]++;
    }

    for (i = 0; i < Q; i++) // 减去中心坐标得到了真正各个构型距离中心的距离，这里得到
    {
        dist_total[i] -= center[i];
    }

    for (i = 0; i < q; i++) // 计算此构型的距离
    {
        for (j = 0; j < q; j++)
        {
            r += dist_total[i] * dist_total[j] * eta[i][j];
        }
    }

    r = sqrt(r);
    r /= r_max;

    return r;
}

#endif