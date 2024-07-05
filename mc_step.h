#ifndef mc_step_h
#define mc_step_h

#include "initialize.h"
#include "quantity.h"

using namespace std; // 默认用库是std，这样可能造成冲突，比如用std里的函数当变量名

void do_step(int spins[SIZE + 1], int neighs[SIZE][nei], double tstar, double &energy, mt19937 &gen, uniform_real_distribution<double> &ran_u, double m[TN][DATA], double bins[TN][NBIN][B], double old[TN][NBIN], uniform_int_distribution<int> &brandom, int up[Q], int b, int count, double eta[q][q], double center[Q], double r_max);
void flip_spin(int spins[SIZE + 1], int neighs[SIZE][nei], double &energy, mt19937 &gen, uniform_real_distribution<double> &ran_u, int index, double tstar, uniform_int_distribution<int> &brandom); // 采样翻转

void do_step(int spins[SIZE + 1], int neighs[SIZE][nei], double tstar, double &energy, mt19937 &gen, uniform_real_distribution<double> &ran_u, double m[TN][DATA], double bins[TN][NBIN][B], double old[TN][NBIN], uniform_int_distribution<int> &brandom, int up[Q], int b, int count, double eta[q][q], double center[Q], double r_max)
// double m[DATA]这个数组包含了一些磁矩和能量，进行普通操作
{

    int i, j;                                         // Counters，计数
    int MCn;                                          // 每个构型之间的马尔可夫步数
    double sum[2];                                    // To compute the sum of spins，算序参量和磁化强度
    double energysum;                                 // 总内能
    double dist;                                      // 构型空间的距离，磁化率，距离均方差，比热
    double mag_chi, energy_chi, order_chi, dist_chi;  // 各个物理量的涨落
    double old_sum[2] = {0};                          // 拿来计算序参量和磁化强度的关联
    double old_energy, old_dist;                      // 拿来计算构型之间能量和距离的自关联
    double old_mag_chi, old_energy_chi, old_dist_chi; // 拿来计算构型之间磁化率和比热和距离涨落的自关联

    // Thermalize the state，预热
    for (i = 0; i < 200 * Q; i++)
    {
        for (j = 0; j < SIZE; j++)
        {
            flip_spin(spins, neighs, energy, gen, ran_u, j, tstar, brandom);
            // 进行SIZE*q次自旋翻转(初始随机挑选1个格点进行翻转，操作N次
        }
    }

    ///----- TODO: optimize the number of steps for thermalization/measures，优化
    old_mag_chi = 0.0; // 两个物理量序参量和磁化强度已经被初始化了
    old_energy = 0.0;
    old_energy_chi = 0.0;
    old_dist = 0.0;
    old_dist_chi = 0.0;

    for (i = 0; i < N; i++) // 抽取1000个样本
    {
        // Make changes and then average变化后计算平均
        for (j = 0; j < SIZE; j++) // 必须每个自旋都给他翻一下，这样避免采样出来的构型都差不多，然后由于我们保证每次都向不同状态翻转，可以保证均匀分配
        {
            flip_spin(spins, neighs, energy, gen, ran_u, j, tstar, brandom); // 以h[(i + 4) / 2] = min(1.0, exp(-2.0 * i / tstar))概率翻转。且计算能量。
        }

        double dist_total[Q] = {0};

        // Compute quantities at time j，计算j时刻的数量
        magnetization(spins, dist_total, sum);                         // 计算序参量和磁化强度
        mag_chi = sum[1] * sum[1];                                     // 用来计算磁化强度涨落
        energy_chi = energy * energy;                                  // 能量涨落
        order_chi = sum[0] * sum[0];                                   // 计算序参量涨落
        dist = distance(spins, eta, center, r_max, count, dist_total); // 构型空间的距离
        dist_chi = dist * dist;                                        // 用来计算构型空间距离涨落

        // Add all the quantities，#define MAG 0，#define MAG2 1，#define MAG4 2，#define MAGERR 3，#define SUSERR 4，#define ENE 5，#define ENE2 6，#define ENE4 7，#define ENERR 8，#define CHERR 9
        m[count][MAG] += sum[1];                   // Magnetization，磁化强度#define MAG 0//#define MAG2 1
        m[count][MAG2] += mag_chi;                 // For the susceptibility，对于磁化率
        m[count][MAG4] += mag_chi * mag_chi;       // For the Binder cumulant and also variance of susceptibility拿来计算binder ratio
        m[count][ENE] += energy;                   // Energy，能量
        m[count][ENE2] += energy_chi;              // For specific heat
        m[count][ENE4] += energy_chi * energy_chi; // For the variance of specific heat
        m[count][ORDER] += sum[0];                 // 序参量
        m[count][ORDER2] += order_chi;             // 序参量涨落
        m[count][DIS] += dist;                     // 距离
        m[count][DIS2] += dist_chi;                // 距离涨落
        // This are used for errors,
        m[count][MAGERR] += old_sum[1] * sum[1];          // in magnetization，计算磁化强度关联
        m[count][SUSERR] += old_mag_chi * mag_chi;        // in susceptibility，计算磁化率关联
        m[count][ENERR] += old_energy * energy;           // in energy，计算能量关联
        m[count][CHERR] += old_energy_chi * energy_chi;   // in specific heat，计算比热关联
        m[count][ORDERRR] += old_energy_chi * energy_chi; // in specific heat，计算序参量关联
        m[count][DISRR] += old_dist * dist;               // in specific heat，计算距离关联
        m[count][DIS2RR] += old_dist_chi * dist_chi;      // in specific heat，计算距离涨落关联

        // p[count][distance(spins, eta)] += 1;

        // Get the value for the next iteration，取值进行下一个迭代
        old_sum[0] = sum[0];
        old_mag_chi = mag_chi;
        old_energy = energy;
        old_energy_chi = energy_chi;
        old_sum[1] = sum[1];
        old_dist = dist;
        old_dist_chi = dist_chi;

        for (j = 0; j < SIZE; j++)
        {
            spins[j] = up[spins[j]]; // 所有自旋都向上加一，防止出现由于自发对称破缺导致停留在一个构型而取不到破缺到的另一个构型
        }
    }

    bins[count][0][b] = (sum[1] - old_sum[1]) / (1.0 * N);
    bins[count][1][b] = (mag_chi - old_mag_chi) / (1.0 * N);
    bins[count][2][b] = (energy - old_energy) / (1.0 * N);
    bins[count][3][b] = (energy_chi - old_energy_chi) / (1.0 * N);
    bins[count][4][b] = (sum[0] - old_sum[0]) / (1.0 * N);
    bins[count][5][b] = (dist - old_dist) / (1.0 * N);
    bins[count][6][b] = (dist_chi - old_dist_chi) / (1.0 * N);

    old_sum[1] = sum[1];
    old_mag_chi = mag_chi;
    old_energy = energy;
    old_energy_chi = energy_chi;
    old_sum[0] = sum[0];
    old_dist = dist;
    old_dist_chi = m[count][DIS2];

    return;
}

// Flip a spin via Metropolis
void flip_spin(int spins[SIZE + 1], int neighs[SIZE][nei], double &energy, mt19937 &gen, uniform_real_distribution<double> &ran_u, int index, double tstar, uniform_int_distribution<int> &brandom) // 采样翻转
// 自旋翻转
{
    int change = 0;
    int state = brandom(gen);
    while (state == spins[index]) // 防止brandom生成的随机数就等于state
    {
        state = brandom(gen);
    }

    for (int i = 0; i < nei; i++)
    {
        if (spins[index] == spins[neighs[index][i]])
        {
            change += J;
        }
        else if (state == spins[neighs[index][i]])
        {
            change -= J;
        }
    }
    // Apply Metropolis skim，用Metropolis接受准则
    if (ran_u(gen) < min(1.0, exp(-change / tstar))) // 如果生成的随机数小于这个式子，则进行翻转，否则不翻转//Metropolis接受准则，tstar为温度
    {
        spins[index] = state;            // 自旋状态变为q
        energy += change / (1.0 * SIZE); // 计算变化后能量
        // cout << change << "  " << (2.0*change)/(1.0*SIZE) << "  " << energy << endl;//可以输出一下能量
    }
    return;
}

#endif