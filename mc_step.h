#ifndef mc_step_h
#define mc_step_h

#include "initialize.h"
#include "quantity.h"

using namespace std; // 默认用库是std，这样可能造成冲突，比如用std里的函数当变量名

void do_step(int spins[SIZE + 1], int neighs[SIZE][nei], double tstar, double &energy, mt19937 &gen, uniform_real_distribution<double> &ran_u, double m[TN][DATA], double bins[TN][6][B], double old[TN][6], uniform_int_distribution<int> &brandom, int up[Q], int b, int count, double eta[q][q], double center[Q], double r_max);
void flip_spin(int spins[SIZE + 1], int neighs[SIZE][nei], double &energy, mt19937 &gen, uniform_real_distribution<double> &ran_u, int index, double tstar, uniform_int_distribution<int> &brandom); // 采样翻转

void do_step(int spins[SIZE + 1], int neighs[SIZE][nei], double tstar, double &energy, mt19937 &gen, uniform_real_distribution<double> &ran_u, double m[TN][DATA], double bins[TN][6][B], double old[TN][6], uniform_int_distribution<int> &brandom, int up[Q], int b, int count, double eta[q][q], double center[Q], double r_max)
// double m[DATA]这个数组包含了一些磁矩和能量，进行普通操作
{

    int i, j;                                      // Counters，计数
    int MCn;                                       // 每个构型之间的马尔可夫步数
    double sum;                                    // To compute the sum of spins，算总自旋
    double energysum;                              // 总内能
    double dist;                                   // 构型空间的距离
    double chi, order, heat;                       // 磁化率，距离均方差，比热
    double old_sum, old_chi, old_heat, old_energy; // 拿来计算构型之间的自关联

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
    old_sum = 0.0;
    old_chi = 0.0;
    old_heat = 0.0;
    old_energy = 0.0;

    for (i = 0; i < N; i++) // 抽取1000个样本
    {
        // Make changes and then average变化后计算平均
        for (j = 0; j < SIZE; j++) // 必须每个自旋都给他翻一下，这样避免采样出来的构型都差不多，然后由于我们保证每次都向不同状态翻转，可以保证均匀分配
        {
            flip_spin(spins, neighs, energy, gen, ran_u, j, tstar, brandom); // 以h[(i + 4) / 2] = min(1.0, exp(-2.0 * i / tstar))概率翻转。且计算能量。
        }

        // Compute quantities at time j，计算j时刻的数量
        sum = magnetization(spins);                        // abs是求数据的绝对值
        dist = distance(spins, eta, center, r_max, count); // 构型空间的距离
        chi = sum * sum;                                   // 用来计算磁化强度涨落
        order = dist * dist;                               // 用来计算构型空间距离均方差
        heat = energy * energy;                            // 能量涨落

        // Add all the quantities，#define MAG 0，#define MAG2 1，#define MAG4 2，#define MAGERR 3，#define SUSERR 4，#define ENE 5，#define ENE2 6，#define ENE4 7，#define ENERR 8，#define CHERR 9
        m[count][MAG] += sum;          // Magnetization，磁化强度#define MAG 0//#define MAG2 1
        m[count][MAG2] += chi;         // For the susceptibility，对于磁化率
        m[count][MAG4] += chi * chi;   // For the Binder cumulant and also variance of susceptibility
        m[count][ENE] += energy;       // Energy，能量
        m[count][ENE2] += heat;        // For specific heat
        m[count][ENE4] += heat * heat; // For the variance of specific heat
        // This are used for errors,
        m[count][MAGERR] += old_sum * sum;      // in magnetization
        m[count][SUSERR] += old_chi * chi;      // in susceptibility
        m[count][ENERR] += old_energy * energy; // in energy
        m[count][CHERR] += old_heat * heat;     // in specific heat，比热

        m[count][DIS] += dist;
        m[count][DIS2] += order;

        // p[count][distance(spins, eta)] += 1;

        // Get the value for the next iteration，取值进行下一个迭代
        old_sum = sum;
        old_energy = energy;
        old_chi = chi;
        old_heat = heat;
        for (j = 0; j < SIZE; j++)
        {
            spins[j] = up[spins[j]]; // 所有自旋都向上加一，防止出现由于自发对称破缺导致停留在一个构型而取不到破缺到的另一个构型
        }
    }

    bins[count][0][b] = (m[count][MAG] - old[count][0]) / (1.0 * N);
    bins[count][1][b] = (m[count][MAG2] - old[count][1]) / (1.0 * N);
    bins[count][2][b] = (m[count][ENE] - old[count][2]) / (1.0 * N);
    bins[count][3][b] = (m[count][ENE2] - old[count][3]) / (1.0 * N);
    bins[count][4][b] = (m[count][ENE2] - old[count][4]) / (1.0 * N);
    bins[count][5][b] = (m[count][ENE2] - old[count][5]) / (1.0 * N);

    old[count][0] = m[count][MAG];
    old[count][1] = m[count][MAG2];
    old[count][2] = m[count][ENE];
    old[count][3] = m[count][ENE2];
    old[count][4] = m[count][DIS];
    old[count][5] = m[count][DIS2];

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