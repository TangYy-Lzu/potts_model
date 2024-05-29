#ifndef write_h
#define write_h

#include "constant.h"
#include "initialize.h"

using namespace std; // 默认用库是std，这样可能造成冲突，比如用std里的函数当变量名

#define dimentions SIZE + 1

void write(ofstream &output, double tstar, int spins[SIZE + 1]);
void bin(ofstream &outputBining, double T[TN], double bins[TN][6][B], int count);
void w_output(ofstream &file, double tstar, double m[TN][DATA], int count); // output传给file，file就等于write文件中的output了
void w_distance(ofstream &outputdistance);
void w_d_minus_one(ofstream &outputdistance, int i, int d, int sum);

void write(ofstream &output, double tstar, int spins[SIZE + 1])
{
    int i;
    // ofstream output;          // 创建流对象为output，输入输出流的名字就叫output
    // output.open("check.txt"); // 创建文件check.txt，文件打开路径为check.txt（这里直接为文件名称），不指定打开方式
    output << tstar << " ";
    for (i = 0; i < SIZE; i++)
    {
        if (i % SIZE == 0)
            output << endl;        // /每当i取16的整数倍时换行，比如第一行为i取0到15，endl是换行
        output << spins[i] << " "; // 显示16个格点的自旋
    }
    output << endl;
    output << endl;
    return;
}

void bin(ofstream &outputBining, double T[TN], double bins[TN][6][B], int count)
{
    int i, j, k;
    // ofstream output;          // 创建流对象为output，输入输出流的名字就叫output
    // output.open("check.txt"); // 创建文件check.txt，文件打开路径为check.txt（这里直接为文件名称），不指定打开方式

    for (i = 0; i < count; i++)
    {
        outputBining << T[i] << " " << endl;
        for (j = 0; j < 6; j++)
        {
            for (k = 0; k < B; k++)
            {
                outputBining << bins[i][j][k] << " "; // 显示16个格点的自旋
            }
            outputBining << endl;
        }
        outputBining << endl;
    }

    outputBining << endl;
    return;
}

void w_output(ofstream &file, double tstar, double m[TN][DATA], int count) // output传给file，file就等于write文件中的output了
{
    file << tstar << " " << 1.0 / tstar << " "; // Write T and B，写出温度tstar和β=1/kt，这里k取1

    // We here take in account residual errors, which, for low T, makes the quantities chi, ch, etc.
    // to diverge. This must be substracted. That's why we use an abs for correlation time and also
    // a check to avoid zero value of variances.

    // Then write the quantities and the corresponding errors to a file. The four things are equal,
    // but each one referred to a different quantity.

    double chi = (m[count][MAG2] - m[count][MAG] * m[count][MAG]) * SIZE;                    // Magnetic susceptibility (up to T factor)，磁场强度涨落，磁化率
    double rhom = chi != 0 ? (m[count][MAGERR] - m[count][MAG] * m[count][MAG]) / chi : 0.0; // Rho magnetization, computed if chi != 0，磁化强度的自关联函数
    double taugm = rhom != 1.0 ? rhom / (1.0 - rhom) : 0.0;                                  // Taug magnetization, computed if rhom != 0，Taug磁化
    file << m[count][MAG] << " " << sqrt(chi * abs(2.0 * taugm + 1) / (1.0 * N)) << " ";     // Write everything，把这些物理量全写出来，sqrt为开根号，abs为绝对值

    double fourth = m[count][MAG4] - m[count][MAG2] * m[count][MAG2];                                  // Susceptibility variance，磁化率涨落
    double rhos = fourth != 0.0 ? (m[count][SUSERR] - m[count][MAG2] * m[count][MAG2]) / fourth : 0.0; // Rho susceptibility，磁化率的自关联函数
    double taugs = rhos != 1.0 ? rhos / (1.0 - rhos) : 0.0;                                            // Taug susceptibility，Taug极化率
    double error_sq = sqrt(fourth * abs(2.0 * taugs + 1) / (1.0 * N));
    file << " " << chi << " " << error_sq << " ";

    double heat = (m[count][ENE2] - m[count][ENE] * m[count][ENE]) * SIZE;                      // Specific heat (up to T^2 factor)，比热
    double rhoe = heat != 0.0 ? (m[count][ENERR] - m[count][ENE] * m[count][ENE]) / heat : 0.0; // 能量的自关联函数
    double tauge = rhoe != 1.0 ? rhoe / (1.0 - rhoe) : 0.0;
    file << " " << m[count][ENE] << " " << sqrt(heat * abs(2.0 * tauge + 1) / (1.0 * N)) << " ";

    double fourth_ene = m[count][ENE4] - m[count][ENE2] * m[count][ENE2];
    double rhoc = fourth_ene != 0.0 ? (m[count][CHERR] - m[count][ENE2] * m[count][ENE2]) / fourth_ene : 0.0; // 比热的自关联函数
    double taugc = rhoc != 1.0 ? rhoc / (1.0 - rhoc) : 0.0;
    file << " " << heat << " " << sqrt(fourth_ene * abs(2.0 * taugc + 1) / (1.0 * N)) << " ";

    // Binder cumulant
    double binder = 1.0 - m[count][MAG4] / (3.0 * m[count][MAG2] * m[count][MAG2]); // Computes 4th cumulant minus one, b-1.计算binder 4th cumulant，用它来得到临界温度
    file << binder << " " << 2.0 * (1.0 - binder) * (error_sq / m[count][MAG2]);

    double order = (m[count][DIS2] - m[count][DIS] * m[count][DIS]);
    file << " " << m[count][DIS] * SIZE << " " << order << " " << endl;

    return;
}

void w_distance(ofstream &outputdistance)
{
    int i, j;
    // ofstream output;          // 创建流对象为output，输入输出流的名字就叫output
    // output.open("check.txt"); // 创建文件check.txt，文件打开路径为check.txt（这里直接为文件名称），不指定打开方式

    int d = q;
    int sum; // 统计该输出到哪个维度要加上的前置元素个数
    sum = 0;
    for (i = 0; i < TN; i++)
    {
        outputdistance << T[i] << endl;
        w_d_minus_one(outputdistance, i, d - 1, sum);
        outputdistance << endl;
    }

    return;
}

void w_d_minus_one(ofstream &outputdistance, int i, int d, int sum)
{
    int k;

    for (k = 0; k < dimentions; k++)
    {
        if (d)
        {
            w_d_minus_one(outputdistance, i, d - 1, sum);
            sum += pow((dimentions), d);
            outputdistance << endl;
        }
        else
        {
            outputdistance << space[i][sum] << " ";
            sum++;
        }
    }

    return;
}

#endif