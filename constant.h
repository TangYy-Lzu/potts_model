#ifndef constant_H
#define constant_H

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <random>
#include <fstream>

#include "_parameter.h"

#define q (Q - 1)
#define SIZE L *L
#define square 4
#define diced 6
#define union_jack 8
#define nei lattice

using namespace std; // 默认用库是std，这样可能造成冲突，比如用std里的函数当变量名

double deltat = 0.1;
double deltat_crit = 0.01; // 相变点附近精确计算
double tmin = 0.0;
double tmax = ((TN - 1 - (tcrit_up - tcrit_down) / deltat_crit) * deltat + tcrit_up - tcrit_down + tmin); // 最大温度为
                                                                                                          // 1 + (tcrit_up - tcrit_down) / deltat_crit + (tmax - tcrit_up + tcrit_down - tmin) / deltat = TN
double doubleprecision = 0.0001;
double tmincompare = tmin + doubleprecision; // 最小温度为0
double tcrit_upCompare = tcrit_up + doubleprecision;
double tcrit_downCompare = tcrit_down + doubleprecision; // Max and min temperature, and interval where we apply Wolff，最大和最小温度，我们用Wolff算法，一种一次翻转一大堆分子的方式

double eta[q][q]; // 度规
double center[Q]; // 还是Q个方向的距离，只是我们以第Q个方向为基准。实际计算只用此空间的一个切面，前q-1个态就行，最后一个作为零点
double r_max;
double T[TN]; // 储存TN个要计算的温度

void constant();
void c_eta(double eta[q][q]);
void c_center(double center[Q]);
double distance_max(double eta[q][q], double center[Q]);
void get_temperature(double T[TN]); // 这个数组存放要计算的各个温度

// int main()
// {
//     constant();
//     cout << nei << endl;
//     cout << endl;
//     cout << "各个要输出的温度是:";
//     for (int i = 0; i < TN; i++)
//     {
//         cout << T[i] << " ";
//     }
//     cout << endl;
// cout << endl;
// cout << "度规是：" << endl;
// for (int i = 0; i < q; i++)
// {
//     for (int j = 0; j < q; j++)
//     {
//         cout << eta[i][j] << " ";
//     }
//     cout << endl;
// }
// cout << endl;
// cout << "中心位置距离是：";
// for (int i = 0; i < q; i++)
// {
//     cout << center[i] << " ";
// }
// cout << endl;
// cout << endl;
// cout << "构型空间最大距离是：" << r_max << endl;
//     return 0;
// }

void constant()
{
    get_temperature(T);
    c_eta(eta);
    c_center(center);
    r_max = distance_max(eta, center);

    return;
}

void c_eta(double eta[q][q])
{
    int i, j, k;
    // 先计算度规
    double c_eta[q][q] = {0};
    c_eta[0][0] = 1;

    if (Q == 2)
    {
        for (i = 0; i < q; i++)
        {
            for (j = 0; j < q; j++)
            {
                for (k = 0; k < q; k++)
                {
                    eta[i][j] += c_eta[i][k] * c_eta[j][k];
                }
            }
        }
        return;
    }

    c_eta[1][0] = 1.0 / 2;
    c_eta[1][1] = sqrt(3) / 2;
    if (Q == 3)
    {
        for (i = 0; i < q; i++)
        {
            for (j = 0; j < q; j++)
            {
                for (k = 0; k < q; k++)
                {
                    eta[i][j] += c_eta[i][k] * c_eta[j][k];
                }
            }
        }
        return;
    }
}

void c_center(double center[Q])
{
    int i;
    for (i = 0; i < Q; i++)
    {
        center[i] = SIZE / Q;
    }
    return;
}

double distance_max(double eta[q][q], double center[Q])
{
    int i, j;
    double r_max;
    double dist_max[Q] = {0};

    r_max = 0;              // 初始化
    dist_max[0] = SIZE;     // 表示有SIZE个格点都是第一个态的构型，这时候距离最大
    for (i = 0; i < Q; i++) // 减去中心坐标得到了真正距离中心的距离，这里得到
    {
        dist_max[i] -= center[i];
    }

    for (i = 0; i < q; i++) // 计算计算构型的最大距离
    {
        for (j = 0; j < q; j++)
        {
            r_max += dist_max[i] * dist_max[j] * eta[i][j];
        }
    }

    r_max = sqrt(r_max);
    return r_max;
}

void get_temperature(double T[TN]) // 这个数组存放要计算的各个温度
{
    double tstar; // Control parameter，控制参数
    int count = 0;
    for (tstar = tmax; tstar > tcrit_upCompare; tstar -= deltat)
    // 由于double类型的精度问题导致i储存进去要比本来的值多一点，比如tstar应该是4.9，在计算机储存时变为了4.900000004或者更小一点，导致这里用大于号得到的结果一样
    // double类型比较大小时要时刻注意
    {
        T[count] = tstar;
        count++;
    }
    for (tstar = tcrit_up; tstar > tcrit_downCompare; tstar -= deltat_crit) // 相变点附近精确计算
    {
        T[count] = tstar;
        count++;
    }
    for (tstar = tcrit_down; tstar > tmincompare; tstar -= deltat)
    {
        T[count] = tstar;
        count++;
    }

    T[TN - 1] = 0.001;

    return;
}

#endif // CONSTANT_H