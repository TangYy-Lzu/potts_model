#ifndef intialize_h
#define intialize_h

#include "constant.h"

#define DATA 17
#define NBIN 7

#define MAG 0
#define MAG2 1
#define MAG4 2
#define MAGERR 3
#define SUSERR 4
#define ENE 5
#define ENE2 6
#define ENE4 7
#define ENERR 8
#define CHERR 9
#define ORDER 10
#define ORDER2 11
#define ORDERRR 12
#define DIS 13
#define DIS2 14
#define DISRR 15
#define DIS2RR 16

double m[TN][DATA], bins[TN][NBIN][B], old[TN][NBIN]; // This array contains several moments of the magnetization and energy，这个数组包含了一些磁矩和能量，bins和old拿来计算误差
double private_m[TN][DATA], private_bins[TN][NBIN][B], private_old[TN][NBIN];
int up[Q];

// 全局变量space拿来记录空间各个位置的构型数目
vector<vector<int>> space;

void initialize();
void initialize_up(int up[Q]);
void initialize_m(double m[TN][DATA], double private_m[TN][DATA]);
void initialize_Array();

// int main()
// {
//     initialize();
//     cout << "自旋状态全部向上加一后变成" << endl;
//     for (int i = 0; i < Q; i++)
//     {
//         cout << up[i];
//     }
//     cout << endl;
//     cout << "初始记录各个温度下物理量的矩阵是" << endl;
//     for (int i = 0; i < TN; i++)
//     {
//         for (int j = 0; j < DATA; j++)
//         {
//             cout << m[i][j];
//         }
//         cout << endl;
//     }
//     cout << endl;
//     cout << "初始构型空间" << endl;

//     for (size_t i = 0; i < TN; i++)
//     {
//         for (size_t j = 0; j < (SIZE + 1); j++)
//         {
//             for (size_t k = 0; k < (SIZE + 1); k++)
//             {
//                 cout << k << ":" << space[i][j * (SIZE + 1) + k] << " "; // 例如，初始化为索引值
//             }
//             cout << endl;
//         }
//     }
//     return 0;
// }

void initialize() // 存放状态全部向上变化一和统计出来的各个TN温度下的物理量
{
    initialize_up(up);
    initialize_m(m, private_m);

    if (WRITE)
    {
        initialize_Array();
    }

    return;
}

void initialize_up(int up[Q])
{
    int i, j;
    for (i = 0; i < Q; i++) // 存放所近邻数少于nei的空近邻, 最近邻，全部向上加一自旋
    {
        up[i] = (i + 1) % Q;
    }
    return;
}

void initialize_m(double m[TN][DATA], double private__m[TN][DATA])
{
    int i, j;
    for (i = 0; i < TN; i++)
    {
        for (j = 0; j < DATA; j++)
        {
            m[i][j] = 0.0; // Init the values，181，#define DATA 9，初始化值
        }
    }
    for (i = 0; i < TN; i++)
    {
        for (j = 0; j < DATA; j++)
        {
            private_m[i][j] = 0.0; // Init the values，181，#define DATA 9，初始化值
        }
    }
    return;
}

// Initialices the grid in which we are going to do the Ising, using random values，初始化我们要做Ising的网格，使用随机值
void initialize_spins(int spins[SIZE + 1], mt19937 &gen, uniform_int_distribution<int> &brandom) // spin size为bool类型，
// 接收一个值引用给gen，gen数据类型为mt19937，接收一个值引用给brandom，brandom数据类型为uniform_int_distribution
{
    int i, j; // j有什么用？

    // Init spins with a random distribution，Init自旋以随机分布
    for (i = 0; i < SIZE; i++) // 由于是小于，所以i从0到SIZE一共SIZE个数字
    {
        spins[i] = brandom(gen); // Generate numbers，生成数字， uniform_int_distribution<int> brandom(0, 1); 获取任意随机0到q之间的整数值
    }
    spins[SIZE] = Q; // 空近邻设置为它

    return;
}

void initialize_Array()
{
    size_t arraySize = pow((SIZE + 1), q);
    // 初始化二维数组，外层大小为 TN，内层大小为 arraySize，所有元素初始化为0
    space = vector<vector<int>>(TN, vector<int>(arraySize, 0));
}

#endif