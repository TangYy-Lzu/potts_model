#include <sstream>
#include <chrono>
#include <omp.h>

#include "constant.cpp"
#include "mc_step.hpp"
#include "mc_step.cpp"
#include "quantity.hpp"
#include "quantity.cpp"
#include "matrix.hpp"
#include "matrix.cpp"
#include "write.hpp"
#include "write.cpp"

int main(void)
{
    auto start = std::chrono::high_resolution_clock::now(); // 获取开始时间

    std::string mag;
    if (Parameters::J < 0)
    {
        mag = "fm"; // MAG已经被定义过了
    }
    else
    {
        mag = "afm";
    }

    std::stringstream quantityStream;
    quantityStream << "../quantity/quantity_" << Parameters::LATTICE_TYPE << "_" << mag << "_Q" << Parameters::Q << "_L" << Parameters::LATTICE_SIZE << ".txt";
    std::string quantity = quantityStream.str();

    std::stringstream imageStream;
    imageStream << "../spins/spins_" << Parameters::LATTICE_TYPE << "_" << mag << "_Q" << Parameters::Q << "_L" << Parameters::LATTICE_SIZE << ".txt";
    std::string image = imageStream.str();

    std::stringstream BiningStream;
    BiningStream << "../bining/bining_" << Parameters::LATTICE_TYPE << "_" << mag << "_Q" << Parameters::Q << "_L" << Parameters::LATTICE_SIZE << ".txt";
    std::string Bining = BiningStream.str();

    std::stringstream correlationStream;
    correlationStream << "../correlation/correlation_" << Parameters::LATTICE_TYPE << "_" << mag << "_Q" << Parameters::Q << "_L" << Parameters::LATTICE_SIZE << ".txt";
    std::string correlation = correlationStream.str();

    std::ofstream output, outputspins, outputBining, outputdistance, outputcorrelation; // Output of the stream，创建流对象为output，再看87行

    output.open(quantity);
    outputspins.open(image);
    outputBining.open(Bining);
    outputcorrelation.open(correlation);

    Constants::constant();
    Initialization::initialize();
    neighbors::get_neighbors(); // Get neighbour table，取近邻

#pragma omp parallel num_threads(8)
    {
        std::uniform_int_distribution<int> brandom(0, q);        // Get any random integer，获取任意0到q之间随机整数值，也让他决定自旋往哪个态调整
        std::uniform_int_distribution<int> ran_pos(0, SIZE - 1); // Get any random integer，获取任意0到SIZE-1随机整数值
        std::uniform_real_distribution<double> ran_u(0.0, 1.0);  // Our uniform variable generator，产生均匀分布

        thread_local std::mt19937 gen(958431198 + omp_get_thread_num()); // 线程独立的随机数生成器

        int spins[SIZE + 1]; // 储存自旋信息

        double tstar, energy; // Control parameter，控制参数

        double private_m[DATA], private_bins[NBIN], private_corr[DATA]; // 暂时储存物理量

#pragma omp for
        for (int i = 0; i < Parameters::BROKEN_SIZE; ++i)
        {
            Initialization::initialize_spins(spins, gen, brandom); // Init randomly，初始化
            energy = quantity::get_energy(spins);                  // Compute initial energy，计算初始内能

            int count = 0;

            for (tstar = Parameters::T_MAX; tstar > Constants::tcrit_upCompare; tstar -= Parameters::DELTA_T)
            // 由于double类型的精度问题导致i储存进去要比本来的值多一点，比如tstar应该是4.9，在计算机储存时变为了4.900000004或者更小一点，导致这里用大于号得到的结果一样
            // double类型比较大小时要时刻注意
            {
                mc_step::do_step(spins, tstar, energy, gen, brandom, ran_pos, ran_u, private_m, private_bins, private_corr);
#pragma omp critical
                {
                    if (i == (Parameters::BROKEN_SIZE - 1))
                    {
                        FileUtils::write(outputspins, tstar, spins);
                    }
                    for (int j = 0; j < DATA; ++j)
                    {
                        Initialization::m[count][j] += private_m[j];
                        Initialization::bins[count][j] += private_bins[j];
                        Initialization::corr[count][j] += private_corr[j];
                    }
                    for (int j = DATA; j < NBIN; ++j)
                    {
                        Initialization::bins[count][j] += private_bins[j];
                    }
                    // for (int j = 0; j < NBIN; ++j)
                    // {
                    //     std::cout << Initialization::bins[count][j] << " " << private_bins[j] << std::endl;
                    // }
                }
                count++;
            }
            for (tstar = Parameters::T_CRIT_UP; tstar > Constants::tcrit_downCompare; tstar -= Parameters::DELTA_T_CRIT) // 相变点附近精确计算
            {
                if (Parameters::J < 0)
                {
                    mc_step::do_step_wolff(spins, tstar, energy, gen, brandom, ran_pos, ran_u, private_m, private_bins, private_corr);
                }
                else if (Parameters::J > 0)
                {
                    mc_step::do_step(spins, tstar, energy, gen, brandom, ran_pos, ran_u, private_m, private_bins, private_corr);
                }
#pragma omp critical
                {
                    if (i == (Parameters::BROKEN_SIZE - 1))
                    {
                        FileUtils::write(outputspins, tstar, spins);
                    }
                    for (int j = 0; j < DATA; ++j)
                    {
                        Initialization::m[count][j] += private_m[j];
                        Initialization::bins[count][j] += private_bins[j];
                        Initialization::corr[count][j] += private_corr[j];
                    }
                    for (int j = DATA; j < NBIN; ++j)
                    {
                        Initialization::bins[count][j] += private_bins[j];
                    }
                }
                count++;
            }
            for (tstar = Parameters::T_CRIT_DOWN; tstar > Constants::tmincompare; tstar -= Parameters::DELTA_T)
            {
                mc_step::do_step(spins, tstar, energy, gen, brandom, ran_pos, ran_u, private_m, private_bins, private_corr);
#pragma omp critical
                {
                    if (i == (Parameters::BROKEN_SIZE - 1))
                    {
                        FileUtils::write(outputspins, tstar, spins);
                    }
                    for (int j = 0; j < DATA; ++j)
                    {
                        Initialization::m[count][j] += private_m[j];
                        Initialization::bins[count][j] += private_bins[j];
                        Initialization::corr[count][j] += private_corr[j];
                    }
                    for (int j = DATA; j < NBIN; ++j)
                    {
                        Initialization::bins[count][j] += private_bins[j];
                    }
                }
                count++;
            }
            tstar = 0.00001;
            mc_step::do_step(spins, tstar, energy, gen, brandom, ran_pos, ran_u, private_m, private_bins, private_corr);
#pragma omp critical
            {
                if (i == (Parameters::BROKEN_SIZE - 1))
                {
                    FileUtils::write(outputspins, tstar, spins);
                }
                for (int j = 0; j < DATA; ++j)
                {
                    Initialization::m[count][j] += private_m[j];
                    Initialization::bins[count][j] += private_bins[j];
                    Initialization::corr[count][j] += private_corr[j];
                }
                for (int j = DATA; j < NBIN; ++j)
                {
                    Initialization::bins[count][j] += private_bins[j];
                }
            }
        }
    }

    // Finish the average 完成平均 and error bar
    for (int i = 0; i < TN; i++)
    {
        for (int j = 0; j < DATA; j++)
        {
            Initialization::m[i][j] /= 1.0 * Parameters::BROKEN_SIZE; // a /=b 的意思是 a = a / b，运算“/”在C++中默认向下取整，若想设为向上取整可设为 a = ceil(a / b)，b亦可指一个表达式
            Initialization::corr[i][j] /= 1.0 * Parameters::BROKEN_SIZE;
            Initialization::bins[i][j] /= 1.0 * Parameters::BROKEN_SIZE;
        }
        FileUtils::w_output(output, Constants::temperature[i], Initialization::m[i]);
        for (int j = DATA; j < NBIN; j++)
        {
            Initialization::bins[i][j] /= 1.0 * Parameters::BROKEN_SIZE; // a /=b 的意思是 a = a / b，运算“/”在C++中默认向下取整，若想设为向上取整可设为 a = ceil(a / b)，b亦可指一个表达式
        }
        FileUtils::w_error(outputBining, Constants::temperature[i], Initialization::m[i], Initialization::bins[i]);
        FileUtils::w_corr(outputcorrelation, Constants::temperature[i], Initialization::corr[i]);
    }

    output.close();
    outputspins.close();
    outputBining.close();
    outputcorrelation.close();

    // 获取结束时间
    auto end = std::chrono::high_resolution_clock::now();

    // 计算持续时间
    std::chrono::duration<double> duration = end - start;

    std::cout << Parameters::LATTICE_TYPE << "_J" << Parameters::J << "_Q" << Parameters::Q << "程序执行时间: " << duration.count() << " 秒" << std::endl;
    std::cout << std::endl;

    return 0;
}