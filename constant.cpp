#include "constant.hpp"

namespace Constants
{
    double doubleprecision = 0.0001;
    double tmincompare = Parameters::T_MIN + doubleprecision;
    double tcrit_upCompare = Parameters::T_CRIT_UP + doubleprecision;
    double tcrit_downCompare = Parameters::T_CRIT_DOWN + doubleprecision;
    double eta[Parameters::Q][Parameters::Q] = {0};
    double center[Parameters::Q] = {0};
    double r_max = 0.0;
    double temperature[TN];

    void constant()
    {
        c_eta();
        c_center();
        distance_max();
        get_temperature();
    }

    void c_eta()
    {
        int i, j, k;
        double c_eta[Parameters::Q][Parameters::Q] = {0};
        c_eta[0][0] = 1;

        if (Parameters::Q == 1)
        {
            for (i = 0; i < Parameters::Q; i++)
            {
                for (j = 0; j < Parameters::Q; j++)
                {
                    eta[i][j] = (i == j) ? 1.0 : 0.0;
                }
            }
            return;
        }

        c_eta[1][0] = 1.0 / 2;
        c_eta[1][1] = sqrt(3) / 2;
        if (Parameters::Q == 2)
        {
            for (i = 0; i < Parameters::Q; i++)
            {
                for (j = 0; j < Parameters::Q; j++)
                {
                    for (k = 0; k < Parameters::Q; k++)
                    {
                        eta[i][j] += c_eta[i][k] * c_eta[j][k];
                    }
                }
            }
            return;
        }

        double normalize = c_eta[1][0] * c_eta[1][0];

        for (i = 2; i < Parameters::Q - 1; i++)
        {
            for (j = 0; j < i - 1; j++)
            {
                c_eta[i][j] = c_eta[i - 1][j];
            }
            c_eta[i][i - 1] = (0.5 - normalize) / c_eta[i - 1][i - 1];
            normalize += c_eta[i][i - 1] * c_eta[i][i - 1];
            c_eta[i][i] = sqrt(1 - normalize);
        }

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
    }

    void c_center()
    {
        for (int i = 0; i < Parameters::Q; i++)
        {
            center[i] = SIZE / static_cast<double>(Parameters::Q);
        }
    }

    void distance_max()
    {
        double dist_max[Parameters::Q] = {0};
        dist_max[0] = SIZE;

        for (int i = 0; i < Parameters::Q; i++)
        {
            dist_max[i] -= center[i];
        }

        for (int i = 0; i < Parameters::Q; i++)
        {
            for (int j = 0; j < Parameters::Q; j++)
            {
                r_max += dist_max[i] * dist_max[j] * eta[i][j];
            }
        }

        r_max = sqrt(r_max);
    }

    void get_temperature() // 这个数组存放要计算的各个温度
    {
        double tstar; // Control parameter，控制参数
        int count = 0;
        for (tstar = Parameters::T_MAX; tstar > Constants::tcrit_upCompare; tstar -= Parameters::DELTA_T)
        // 由于double类型的精度问题导致i储存进去要比本来的值多一点，比如tstar应该是4.9，在计算机储存时变为了4.900000004或者更小一点，导致这里用大于号得到的结果一样
        // double类型比较大小时要时刻注意
        {
            temperature[count] = tstar;
            count++;
        }
        for (tstar = Parameters::T_CRIT_UP; tstar > Constants::tcrit_downCompare; tstar -= Parameters::DELTA_T_CRIT) // 相变点附近精确计算
        {
            temperature[count] = tstar;
            count++;
        }
        for (tstar = Parameters::T_CRIT_DOWN; tstar > Constants::tmincompare; tstar -= Parameters::DELTA_T)
        {
            temperature[count] = tstar;
            count++;
        }

        temperature[TN - 1] = 0.001;

        return;
    }
}

// int main()
// {
//     Constants::constant();
//     std::cout << std::endl;
//     std::cout << "度规是：" << std::endl;
//     for (int i = 0; i < Parameters::Q; i++)
//     {
//         for (int j = 0; j < Parameters::Q; j++)
//         {
//             std::cout << Constants::eta[i][j] << " ";
//         }
//         std::cout << std::endl;
//     }
//     std::cout << std::endl;
//     std::cout << "中心位置距离是：";
//     for (int i = 0; i < Parameters::Q; i++)
//     {
//         std::cout << Constants::center[i] << " ";
//     }
//     std::cout << std::endl;
//     std::cout << std::endl;
//     std::cout << "构型空间最大距离是：" << Constants::r_max << std::endl;
//     std::cout << std::endl;
//     std::cout << "温度个数是：";
//     std::cout << TN << std::endl;
//     std::cout << std::endl;
//     std::cout << "温度分别是：";
//     for (int i = 0; i < TN; i++)
//     {
//         std::cout << Constants::temperature[i] << ",";
//     }
//     return 0;
// }