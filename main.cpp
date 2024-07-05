#include <sstream>
#include <chrono>

#include "mc_step.h"
#include "quantity.h"
#include "write.h"

int spins[SIZE + 1]; // 储存自旋信息

using namespace std; // 默认用库是std，这样可能造成冲突，比如用std里的函数当变量名

int main(void)
{
    auto start = chrono::high_resolution_clock::now(); // 获取开始时间
    int i, j;

    initialize();

    string Lattice;
    if (nei == 4)
    {
        Lattice = "square"; // lattice已经被定义过了
    }
    else if (nei == 6)
    {
        Lattice = "diced"; // lattice已经被定义过了
    }
    else if (nei == 8)
    {
        Lattice = "union_jack"; // lattice已经被定义过了
    }

    string mag;
    if (J > 0)
    {
        mag = "fm"; // MAG已经被定义过了
    }
    else
    {
        mag = "afm";
    }

    stringstream quantityStream;
    quantityStream << "../quantity/quantity_" << Lattice << "_" << mag << "_Q" << Q << "_L" << L << ".txt";
    string quantity = quantityStream.str();

    stringstream imageStream;
    imageStream << "../spins/spins_" << Lattice << "_" << mag << "_Q" << Q << "_L" << L << ".txt";
    string image = imageStream.str();

    stringstream BiningStream;
    BiningStream << "../bining/bining_" << Lattice << "_" << mag << "_Q" << Q << "_L" << L << ".txt";
    string Bining = BiningStream.str();

    stringstream distanceStream;
    distanceStream << "../distance/distance_" << Lattice << "_" << mag << "_Q" << Q << "_L" << L << ".txt";
    string distance = distanceStream.str();

    mt19937 gen(958431198);                      // Mersenne Twister RNG，马特赛特旋转演算法，是伪随机数发生器之一，958431198是伪随机数种子
    uniform_int_distribution<int> brandom(0, q); // Get any random integer，获取任意0到q之间随机整数值，也让他决定自旋往哪个态调整
    // uniform_int_distribution<int> ran_pos(0, SIZE - 1); // Get any random integer，获取任意0到SIZE-1随机整数值
    uniform_real_distribution<double> ran_u(0.0, 1.0); // Our uniform variable generator，产生均匀分布

    ofstream output, outputspins, outputBining, outputdistance; // Output of the stream，创建流对象为output，再看87行

    if (WRITE)
    {
        ofstream outputdistance;
    }

    output.open(quantity);
    outputspins.open(image);
    outputBining.open(Bining);

    // 初始化
    constant();
    get_neighbors(neighs); // Get neighbour table，取近邻

    if (WRITE)
    {
        outputdistance.open(distance);
    }

#pragma omp parallel num_threads(12) firstprivate(spins, private_m, private_bins, private_old)
    {
#pragma omp for
        for (int b = 0; b < B; b++)
        {
            double energy;
            initialize_spins(spins, gen, brandom); // Init randomly，初始化
            energy = get_energy(spins, neighs);    // Compute initial energy，计算初始内能

            double tstar; // Control parameter，控制参数
            int count = 0;

            for (tstar = T[count]; tstar > tcrit_up; tstar -= deltat)
            // 由于double类型的精度问题导致i储存进去要比本来的值多一点，比如tstar应该是4.9，在计算机储存时变为了4.900000004或者更小一点，导致这里用大于号得到的结果一样
            // double类型比较大小时要时刻注意
            {
                do_step(spins, neighs, tstar, energy, gen, ran_u, private_m, private_bins, private_old, brandom, up, b, count, eta, center, r_max);
                if (b == 99)
                {
                    write(outputspins, tstar, spins);
                }
                count++;
            }
            for (tstar = T[count]; tstar > tcrit_downCompare; tstar -= deltat_crit) // 相变点附近精确计算
            {
                do_step(spins, neighs, tstar, energy, gen, ran_u, private_m, private_bins, private_old, brandom, up, b, count, eta, center, r_max);
                if (b == 99)
                {
                    write(outputspins, tstar, spins);
                }
                count++;
            }
            for (tstar = T[count]; tstar > tmincompare; tstar -= deltat)
            {
                do_step(spins, neighs, tstar, energy, gen, ran_u, private_m, private_bins, private_old, brandom, up, b, count, eta, center, r_max);
                if (b == 99)
                {
                    write(outputspins, tstar, spins);
                }
                count++;
            }
            tstar = T[count];
            do_step(spins, neighs, tstar, energy, gen, ran_u, private_m, private_bins, private_old, brandom, up, b, count, eta, center, r_max);
            if (b == 99)
            {
                write(outputspins, tstar, spins);
            }
            cout << b << endl;
        }

#pragma omp critical
        {
            for (int i = 0; i < TN; ++i)
            {
                for (int j = 0; j < DATA; ++j)
                {
                    m[i][j] += private_m[i][j];
                }
            }
            for (int i = 0; i < TN; ++i)
            {
                for (int j = 0; j < NBIN; ++j)
                {
                    old[i][j] += private_old[i][j];
                }
            }
            for (int i = 0; i < TN; ++i)
            {
                for (int j = 0; j < NBIN; ++j)
                {
                    for (int k = 0; k < B; ++k)
                    {

                        bins[i][j][B] += private_bins[i][j][B];
                    }
                }
            }
        }
    }

    // Finish the average完成平均
    for (i = 0; i < TN; i++)
    {
        for (j = 0; j < DATA; j++)
        {
            m[i][j] /= (1.0 * N * B); // a /=b 的意思是 a = a / b，运算“/”在C++中默认向下取整，若想设为向上取整可设为 a = ceil(a / b)，b亦可指一个表达式
        }
        w_output(output, T[i], m, i);
    }

    bin(outputBining, T, bins, i);

    if (WRITE)
    {
        w_distance(outputdistance);
    }

    output.close();
    outputspins.close();
    outputBining.close();
    if (WRITE)
    {
        outputdistance.close();
    }

    // 获取结束时间
    auto end = chrono::high_resolution_clock::now();

    // 计算持续时间
    chrono::duration<double> duration = end - start;

    cout << "程序执行时间: " << duration.count() << " 秒" << endl;

    return 0;
}