#include "mc_step.hpp"

namespace mc_step
{
    void do_step(int spins[SIZE + 1], double tstar, double &energy, std::mt19937 &gen,
                 std::uniform_int_distribution<int> &brandom, std::uniform_int_distribution<int> ran_pos, std::uniform_real_distribution<double> &ran_u,
                 double *private_m, double *private_bins, double *private_corr)
    {
        int i, j, k;
        double mag, mag_plane, dist, order, n_area, euler_number, single_value, nuclera_norm, shannon;
        double energy_chi, mag_chi, mag_plane_chi, dist_chi, order_chi, n_area_chi, euler_number_chi, single_value_chi, nuclera_norm_chi, shannon_chi;

        Initialization::initialize_matrices(private_m, private_bins, private_corr);

        // 预热步骤
        for (i = 0; i < 200 * Parameters::Q; i++)
        {
            for (j = 0; j < SIZE; j++)
            {
                flip_spin(spins, energy, gen, tstar, brandom, ran_pos, ran_u);
            }
        }

        double old[DATA];
        for (i = 0; i < DATA; i++)
        {
            old[i] = 0;
        }

        double private_m_corr[DATA];
        double private_bins_corr[NBIN];

        for (i = 0; i < Parameters::SAMPLE_SIZE; i++)
        {
            for (j = 0; j < DATA; j++)
            {
                private_bins[j] = 0;
            }

            for (j = 0; j < Parameters::GROUP_SIZE; j++)
            {
                for (k = 0; k < SIZE; k++)
                {
                    flip_spin(spins, energy, gen, tstar, brandom, ran_pos, ran_u);
                }

                double sum[4] = {0}; // 空间距离
                double sv[3] = {0};  // 范数

                // energy = quantity::get_energy(spins);
                if (Parameters::LATTICE_TYPE == "square")
                {
                    quantity::distance_square(spins, sum);
                }
                else if (Parameters::LATTICE_TYPE == "diced")
                {
                    quantity::distance_diced(spins, sum);
                }
                else
                {
                    quantity::distance_union_jack(spins, sum);
                }
                matrix::Single_Value(spins, sv);

                order = sum[0];
                mag_plane = sum[1];
                mag = sum[2];
                dist = sum[3];
                n_area = 0; // 块数
                euler_number = 0;
                single_value = sv[0];
                nuclera_norm = sv[1];
                shannon = sv[2];

                order_chi = order * order;
                mag_plane_chi = mag_plane * mag_plane;
                mag_chi = mag * mag;
                energy_chi = energy * energy;
                dist_chi = dist * dist;
                n_area_chi = n_area * n_area;
                euler_number_chi = euler_number * euler_number;
                single_value_chi = single_value * single_value;
                nuclera_norm_chi = nuclera_norm * nuclera_norm;
                shannon_chi = shannon * shannon;

                // 存储物理量
                private_bins[static_cast<int>(Quantity::MAG_PLANE)] += mag_plane;
                private_bins[static_cast<int>(Quantity::MAG_PLANE2)] += mag_plane_chi;
                private_bins[static_cast<int>(Quantity::MAG)] += mag;
                private_bins[static_cast<int>(Quantity::MAG2)] += mag_chi;
                private_bins[static_cast<int>(Quantity::ENE)] += energy;
                private_bins[static_cast<int>(Quantity::ENE2)] += energy_chi;
                private_bins[static_cast<int>(Quantity::ORDER)] += order;
                private_bins[static_cast<int>(Quantity::ORDER2)] += order_chi;
                private_bins[static_cast<int>(Quantity::DIS)] += dist;
                private_bins[static_cast<int>(Quantity::DIS2)] += dist_chi;
                private_bins[static_cast<int>(Quantity::AREA)] += n_area;
                private_bins[static_cast<int>(Quantity::AREA2)] += n_area_chi;
                private_bins[static_cast<int>(Quantity::EULER)] += euler_number;
                private_bins[static_cast<int>(Quantity::EULER2)] += euler_number_chi;
                private_bins[static_cast<int>(Quantity::SINGLE_VALUE)] += single_value;
                private_bins[static_cast<int>(Quantity::SINGLE_VALUE2)] += single_value_chi;
                private_bins[static_cast<int>(Quantity::NUCLEAR_NORM)] += nuclera_norm;
                private_bins[static_cast<int>(Quantity::NUCLEAR_NORM2)] += nuclera_norm_chi;
                private_bins[static_cast<int>(Quantity::SHANNON)] += shannon;
                private_bins[static_cast<int>(Quantity::SHANNON2)] += shannon_chi;
                private_bins[static_cast<int>(Quantity::MAG4)] += mag_chi * mag_chi;

                // for (k = 0; k < SIZE; k++) // 各项同性ptts模型可以加上
                // {
                //     spins[j] = up[spins[j]];
                // }
            }

            for (j = 0; j < DATA; j++)
            {
                private_bins[j] /= (1.0 * Parameters::GROUP_SIZE);
            }

            update_bins(private_m, private_bins, private_corr, old);

            if (i == Parameters::SAMPLE_SIZE - 2)
            {
                for (int j = 0; j < DATA; j++)
                {
                    private_m_corr[j] = private_m[j];
                    private_bins_corr[j] = private_bins[j];
                }
                for (int j = DATA; j < NBIN; j++)
                {
                    private_bins_corr[j] = private_bins[j];
                }
            }
        }

        for (i = 0; i < DATA; i++)
        {
            private_m[i] /= (1.0 * Parameters::SAMPLE_SIZE);
            private_m_corr[i] /= (1.0 * Parameters::SAMPLE_SIZE - 1);
            private_bins[i] /= (1.0 * Parameters::SAMPLE_SIZE);
            private_bins_corr[i] /= (1.0 * Parameters::SAMPLE_SIZE - 1);
            private_corr[i] /= (1.0 * Parameters::SAMPLE_SIZE - 1);
        }
        for (i = DATA; i < NBIN; i++)
        {
            private_bins[i] /= (1.0 * Parameters::SAMPLE_SIZE); // a /=b 的意思是 a = a / b，运算“/”在C++中默认向下取整，若想设为向上取整可设为 a = ceil(a / b)，b亦可指一个表达式
            private_bins_corr[i] /= (1.0 * Parameters::SAMPLE_SIZE - 1);
        }

        count_corr(private_m_corr, private_bins_corr, private_corr);
    }

    void flip_spin(int spins[SIZE + 1], double &energy, std::mt19937 &gen, double tstar,
                   std::uniform_int_distribution<int> &brandom, std::uniform_int_distribution<int> ran_pos, std::uniform_real_distribution<double> &ran_u)
    {
        int index;
        index = ran_pos(gen);

        int change = 0;
        int state = brandom(gen);
        while (state == spins[index])
        {
            state = brandom(gen);
        }

        for (int i = 0; i < neighbors::nei; i++)
        {
            if (spins[index] == spins[neighbors::neighs[index][i]])
            {
                change -= Parameters::J;
            }
            else if (state == spins[neighbors::neighs[index][i]])
            {
                change += Parameters::J;
            }
        }

        if (ran_u(gen) < std::min(1.0, exp(-change / tstar)))
        {
            spins[index] = state;
            energy += change / (1.0 * SIZE);
        }
    }

    // Do all the things needed for a certain temperature
    void do_step_wolff(int spins[SIZE + 1], double tstar, double &energy, std::mt19937 &gen,
                       std::uniform_int_distribution<int> &brandom, std::uniform_int_distribution<int> ran_pos, std::uniform_real_distribution<double> &ran_u,
                       double *private_m, double *private_bins, double *private_corr)
    {
        int i, j, k;
        double mag, mag_plane, dist, order, n_area, euler_number, single_value, nuclera_norm, shannon;
        double energy_chi, mag_chi, mag_plane_chi, dist_chi, order_chi, n_area_chi, euler_number_chi, single_value_chi, nuclera_norm_chi, shannon_chi;

        Initialization::initialize_matrices(private_m, private_bins, private_corr);

        double pa = 1.0 - exp(Parameters::J / tstar); // TODO change添加种子的概率，这里用温度代替了β，tstar=1/kβ，J在计算能量处选为1

        // 预热步骤
        for (i = 0; i < Parameters::Q; i++)
        {
            for (j = 0; j < SIZE; j++)
            {
                area(spins, energy, gen, tstar, brandom, ran_pos, ran_u, pa);
            }
        }

        double old[DATA];
        for (i = 0; i < DATA; i++)
        {
            old[i] = 0;
        }

        double private_m_corr[DATA];
        double private_bins_corr[NBIN];

        for (i = 0; i < Parameters::SAMPLE_SIZE; i++)
        {
            for (j = 0; j < DATA; j++)
            {
                private_bins[j] = 0;
            }

            for (j = 0; j < Parameters::GROUP_SIZE; j++)
            {
                area(spins, energy, gen, tstar, brandom, ran_pos, ran_u, pa);

                double sum[4] = {0}; // 空间距离
                double sv[3] = {0};  // 范数

                // energy = quantity::get_energy(spins);
                if (Parameters::LATTICE_TYPE == "square")
                {
                    quantity::distance_square(spins, sum);
                }
                else if (Parameters::LATTICE_TYPE == "diced")
                {
                    quantity::distance_diced(spins, sum);
                }
                else
                {
                    quantity::distance_union_jack(spins, sum);
                }
                matrix::Single_Value(spins, sv);

                order = sum[0];
                mag_plane = sum[1];
                mag = sum[2];
                dist = sum[3];
                n_area = 0; // 块数
                euler_number = 0;
                single_value = sv[0];
                nuclera_norm = sv[1];
                shannon = sv[2];

                order_chi = order * order;
                mag_plane_chi = mag_plane * mag_plane;
                mag_chi = mag * mag;
                energy_chi = energy * energy;
                dist_chi = dist * dist;
                n_area_chi = n_area * n_area;
                euler_number_chi = euler_number * euler_number;
                single_value_chi = single_value * single_value;
                nuclera_norm_chi = nuclera_norm * nuclera_norm;
                shannon_chi = shannon * shannon;

                // 存储物理量
                private_bins[static_cast<int>(Quantity::MAG_PLANE)] += mag_plane;
                private_bins[static_cast<int>(Quantity::MAG_PLANE2)] += mag_plane_chi;
                private_bins[static_cast<int>(Quantity::MAG)] += mag;
                private_bins[static_cast<int>(Quantity::MAG2)] += mag_chi;
                private_bins[static_cast<int>(Quantity::ENE)] += energy;
                private_bins[static_cast<int>(Quantity::ENE2)] += energy_chi;
                private_bins[static_cast<int>(Quantity::ORDER)] += order;
                private_bins[static_cast<int>(Quantity::ORDER2)] += order_chi;
                private_bins[static_cast<int>(Quantity::DIS)] += dist;
                private_bins[static_cast<int>(Quantity::DIS2)] += dist_chi;
                private_bins[static_cast<int>(Quantity::AREA)] += n_area;
                private_bins[static_cast<int>(Quantity::AREA2)] += n_area_chi;
                private_bins[static_cast<int>(Quantity::EULER)] += euler_number;
                private_bins[static_cast<int>(Quantity::EULER2)] += euler_number_chi;
                private_bins[static_cast<int>(Quantity::SINGLE_VALUE)] += single_value;
                private_bins[static_cast<int>(Quantity::SINGLE_VALUE2)] += single_value_chi;
                private_bins[static_cast<int>(Quantity::NUCLEAR_NORM)] += nuclera_norm;
                private_bins[static_cast<int>(Quantity::NUCLEAR_NORM2)] += nuclera_norm_chi;
                private_bins[static_cast<int>(Quantity::SHANNON)] += shannon;
                private_bins[static_cast<int>(Quantity::SHANNON2)] += shannon_chi;
                private_bins[static_cast<int>(Quantity::MAG4)] += mag_chi * mag_chi;

                // for (k = 0; k < SIZE; k++) // 各项同性ptts模型可以加上
                // {
                //     spins[j] = up[spins[j]];
                // }
            }

            for (j = 0; j < DATA; j++)
            {
                private_bins[j] /= (1.0 * Parameters::GROUP_SIZE);
            }

            if (i == Parameters::SAMPLE_SIZE - 2)
            {
                for (int j = 0; j < DATA; j++)
                {
                    private_m_corr[j] = private_m[j];
                    private_bins_corr[j] = private_bins[j];
                }
                for (int j = DATA; j < NBIN; j++)
                {
                    private_bins_corr[j] = private_bins[j];
                }
            }
        }

        for (i = 0; i < DATA; i++)
        {
            private_m[i] /= (1.0 * Parameters::SAMPLE_SIZE);
            private_m_corr[i] /= (1.0 * Parameters::SAMPLE_SIZE - 1);
            private_bins[i] /= (1.0 * Parameters::SAMPLE_SIZE);
            private_bins_corr[i] /= (1.0 * Parameters::SAMPLE_SIZE - 1);
            private_corr[i] /= (1.0 * Parameters::SAMPLE_SIZE - 1);
        }
        for (i = DATA; i < NBIN; i++)
        {
            private_bins[i] /= (1.0 * Parameters::SAMPLE_SIZE); // a /=b 的意思是 a = a / b，运算“/”在C++中默认向下取整，若想设为向上取整可设为 a = ceil(a / b)，b亦可指一个表达式
            private_bins_corr[i] /= (1.0 * Parameters::SAMPLE_SIZE - 1);
        }

        count_corr(private_m_corr, private_bins_corr, private_corr);
    }

    void area(int spins[SIZE + 1], double &energy, std::mt19937 &gen, double tstar,
              std::uniform_int_distribution<int> &brandom, std::uniform_int_distribution<int> ran_pos, std::uniform_real_distribution<double> &ran_u, double pa)
    {
        int pos;
        int state_be;

        pos = ran_pos(gen);

        int state_af = brandom(gen);
        while (state_af == spins[pos])
        {
            state_af = brandom(gen);
        }

        state_be = spins[pos];

        add_to_cluster(spins, pos, state_be, state_af, energy, gen, tstar, brandom, ran_pos, ran_u, pa);
    }

    // Executes Wolff algorithm in a recursive way，用递归执行wolff算法。
    void add_to_cluster(int spins[SIZE + 1], int pos, int state_be, int state_af, double &energy, std::mt19937 &gen, double tstar,
                        std::uniform_int_distribution<int> &brandom, std::uniform_int_distribution<int> ran_pos, std::uniform_real_distribution<double> &ran_u, double pa) // 注意吧pa的值传给p了
    // 这里pos取值已经在260行给了一个随机数，添加一个cluster，并且在这个过程中翻转好了
    {
        int i; // Counter
        int n; // Neighbour position近邻位置

        int change = 0;

        // Compute sum of neighbours and change in energy，计算近邻和能量变化
        for (i = 0; i < neighbors::nei; i++)
        {
            if (state_be == spins[neighbors::neighs[pos][i]])
            {
                change -= Parameters::J;
            }
            else if (state_af == spins[neighbors::neighs[pos][i]])
            {
                change += Parameters::J;
            }
        }

        // Then modify the energy
        energy += (1.0 * change) / (1.0 * SIZE); // 计算总能量

        spins[pos] = state_af; // Flip the spin自旋翻转

        // For every neighbour,
        for (i = 0; i < neighbors::nei; i++)
        {
            n = neighbors::neighs[pos][i]; // Get its position取近邻位置
            if (state_be == spins[n])      // Check if it the same (remember we flipped)检测自旋是否相同注意我们在上面翻转了自旋。
            {
                // Add it to the cluster with certain probability.
                if (ran_u(gen) < pa) // uniform_int_distribution<int> ran_pos(0, SIZE - 1); Get any random integer，获取任意0到SIZE-1随机整数值，按照这个几率把它加入到cluster，进行翻转
                {
                    add_to_cluster(spins, n, state_be, state_af, energy, gen, tstar, brandom, ran_pos, ran_u, pa); // 注意这里把n的取值传递给pos了
                }
            }
        }

        return;
    }

    // Do all the things needed for a certain temperature
    void do_step_wolff_afm(int spins[SIZE + 1], double tstar, double &energy, std::mt19937 &gen,
                           std::uniform_int_distribution<int> &brandom, std::uniform_int_distribution<int> ran_pos, std::uniform_real_distribution<double> &ran_u,
                           double *private_m, double *private_bins, double *private_corr)
    {
        int i, j, k;
        double mag, mag_plane, dist, order, n_area, euler_number, single_value, nuclera_norm, shannon;
        double energy_chi, mag_chi, mag_plane_chi, dist_chi, order_chi, n_area_chi, euler_number_chi, single_value_chi, nuclera_norm_chi, shannon_chi;

        Initialization::initialize_matrices(private_m, private_bins, private_corr);

        double pa = 1.0 - exp(-Parameters::J / tstar); // TODO change添加种子的概率，这里用温度代替了β，tstar=1/kβ，J在计算能量处选为1

        // 预热步骤
        for (i = 0; i < Parameters::Q; i++)
        {
            for (j = 0; j < SIZE; j++)
            {
                area_afm(spins, energy, gen, tstar, brandom, ran_pos, ran_u, pa);
            }
        }

        double old[DATA];
        for (i = 0; i < DATA; i++)
        {
            old[i] = 0;
        }

        double private_m_corr[DATA];
        double private_bins_corr[NBIN];

        for (i = 0; i < Parameters::SAMPLE_SIZE; i++)
        {
            for (j = 0; j < DATA; j++)
            {
                private_bins[j] = 0;
            }

            for (j = 0; j < Parameters::GROUP_SIZE; j++)
            {
                area_afm(spins, energy, gen, tstar, brandom, ran_pos, ran_u, pa);

                double sum[4] = {0}; // 空间距离
                double sv[3] = {0};  // 范数

                // energy = quantity::get_energy(spins);
                if (Parameters::LATTICE_TYPE == "square")
                {
                    quantity::distance_square(spins, sum);
                }
                else if (Parameters::LATTICE_TYPE == "diced")
                {
                    quantity::distance_diced(spins, sum);
                }
                else
                {
                    quantity::distance_union_jack(spins, sum);
                }
                matrix::Single_Value(spins, sv);

                order = sum[0];
                mag_plane = sum[1];
                mag = sum[2];
                dist = sum[3];
                n_area = 0; // 块数
                euler_number = 0;
                single_value = sv[0];
                nuclera_norm = sv[1];
                shannon = sv[2];

                order_chi = order * order;
                mag_plane_chi = mag_plane * mag_plane;
                mag_chi = mag * mag;
                energy_chi = energy * energy;
                dist_chi = dist * dist;
                n_area_chi = n_area * n_area;
                euler_number_chi = euler_number * euler_number;
                single_value_chi = single_value * single_value;
                nuclera_norm_chi = nuclera_norm * nuclera_norm;
                shannon_chi = shannon * shannon;

                // 存储物理量
                private_bins[static_cast<int>(Quantity::MAG_PLANE)] += mag_plane;
                private_bins[static_cast<int>(Quantity::MAG_PLANE2)] += mag_plane_chi;
                private_bins[static_cast<int>(Quantity::MAG)] += mag;
                private_bins[static_cast<int>(Quantity::MAG2)] += mag_chi;
                private_bins[static_cast<int>(Quantity::ENE)] += energy;
                private_bins[static_cast<int>(Quantity::ENE2)] += energy_chi;
                private_bins[static_cast<int>(Quantity::ORDER)] += order;
                private_bins[static_cast<int>(Quantity::ORDER2)] += order_chi;
                private_bins[static_cast<int>(Quantity::DIS)] += dist;
                private_bins[static_cast<int>(Quantity::DIS2)] += dist_chi;
                private_bins[static_cast<int>(Quantity::AREA)] += n_area;
                private_bins[static_cast<int>(Quantity::AREA2)] += n_area_chi;
                private_bins[static_cast<int>(Quantity::EULER)] += euler_number;
                private_bins[static_cast<int>(Quantity::EULER2)] += euler_number_chi;
                private_bins[static_cast<int>(Quantity::SINGLE_VALUE)] += single_value;
                private_bins[static_cast<int>(Quantity::SINGLE_VALUE2)] += single_value_chi;
                private_bins[static_cast<int>(Quantity::NUCLEAR_NORM)] += nuclera_norm;
                private_bins[static_cast<int>(Quantity::NUCLEAR_NORM2)] += nuclera_norm_chi;
                private_bins[static_cast<int>(Quantity::SHANNON)] += shannon;
                private_bins[static_cast<int>(Quantity::SHANNON2)] += shannon_chi;
                private_bins[static_cast<int>(Quantity::MAG4)] += mag_chi * mag_chi;

                // for (k = 0; k < SIZE; k++) // 各项同性ptts模型可以加上
                // {
                //     spins[j] = up[spins[j]];
                // }
            }

            for (j = 0; j < DATA; j++)
            {
                private_bins[j] /= (1.0 * Parameters::GROUP_SIZE);
            }

            if (i == Parameters::SAMPLE_SIZE - 2)
            {
                for (int j = 0; j < DATA; j++)
                {
                    private_m_corr[j] = private_m[j];
                    private_bins_corr[j] = private_bins[j];
                }
                for (int j = DATA; j < NBIN; j++)
                {
                    private_bins_corr[j] = private_bins[j];
                }
            }
        }

        for (i = 0; i < DATA; i++)
        {
            private_m[i] /= (1.0 * Parameters::SAMPLE_SIZE);
            private_m_corr[i] /= (1.0 * Parameters::SAMPLE_SIZE - 1);
            private_bins[i] /= (1.0 * Parameters::SAMPLE_SIZE);
            private_bins_corr[i] /= (1.0 * Parameters::SAMPLE_SIZE - 1);
            private_corr[i] /= (1.0 * Parameters::SAMPLE_SIZE - 1);
        }
        for (i = DATA; i < NBIN; i++)
        {
            private_bins[i] /= (1.0 * Parameters::SAMPLE_SIZE); // a /=b 的意思是 a = a / b，运算“/”在C++中默认向下取整，若想设为向上取整可设为 a = ceil(a / b)，b亦可指一个表达式
            private_bins_corr[i] /= (1.0 * Parameters::SAMPLE_SIZE - 1);
        }

        count_corr(private_m_corr, private_bins_corr, private_corr);
    }

    void area_afm(int spins[SIZE + 1], double &energy, std::mt19937 &gen, double tstar,
                  std::uniform_int_distribution<int> &brandom, std::uniform_int_distribution<int> ran_pos, std::uniform_real_distribution<double> &ran_u, double pa)
    {
        int pos;

        pos = ran_pos(gen);

        add_to_cluster_afm(spins, pos, energy, gen, tstar, brandom, ran_pos, ran_u, pa);

        // for (int i = 0; i < Parameters::LATTICE_SIZE; i++)
        // {
        //     for (int i = 0; i < Parameters::LATTICE_SIZE; i++)
        //     {
        //         std::cout << spins[i] << " ";
        //     }
        //     std::cout << std::endl;
        // }
        // std::cout << spins[SIZE] << std::endl;

        for (int i = 0; i < SIZE; i++)
        {
            spins[i] = spins[i] % Parameters::Q;
        }
    }

    // Executes Wolff algorithm in a recursive way，用递归执行wolff算法。
    void add_to_cluster_afm(int spins[SIZE + 1], int pos, double &energy, std::mt19937 &gen, double tstar,
                            std::uniform_int_distribution<int> &brandom, std::uniform_int_distribution<int> ran_pos, std::uniform_real_distribution<double> &ran_u, double pa) // 注意吧pa的值传给p了
    // 这里pos取值已经在260行给了一个随机数，添加一个cluster，并且在这个过程中翻转好了
    {
        int i; // Counter
        int n; // Neighbour position近邻位置

        int change = 0;

        int state_be = spins[pos];

        int state_mark = spins[pos] + Parameters::Q + 1; // 反铁磁自旋向上翻转1，这里先标记

        int state_af = Initialization::up[spins[pos]];

        // Compute sum of neighbours and change in energy，计算近邻和能量变化
        for (i = 0; i < neighbors::nei; i++)
        {
            if (state_be == (spins[neighbors::neighs[pos][i]] % Parameters::Q))
            {
                change -= Parameters::J;
            }
            else if (state_af == (spins[neighbors::neighs[pos][i]] % Parameters::Q))
            {
                change += Parameters::J;
            }
        }

        // Then modify the energy
        energy += (1.0 * change) / (1.0 * SIZE); // 计算总能量

        spins[pos] = state_mark; // Flip the spin自旋翻转

        // for (int i = 0; i < Parameters::LATTICE_SIZE; i++)
        // {
        //     for (int j = 0; j < Parameters::LATTICE_SIZE; j++)
        //     {
        //         std::cout << spins[i * Parameters::LATTICE_SIZE + j] << " ";
        //     }
        //     std::cout << std::endl;
        // }
        // std::cout << spins[SIZE] << std::endl;

        // For every neighbour,
        for (i = 0; i < neighbors::nei; i++)
        {
            n = neighbors::neighs[pos][i]; // Get its position取近邻位置
            if (state_be == spins[n])      // Check if it the same (remember we flipped)检测自旋是否相同注意我们在上面翻转了自旋。
            {
                spins[n] = spins[n] + Parameters::Q; // 标记周围状态相同的格点
            }
        }

        // For every neighbour,
        for (i = 0; i < neighbors::nei; i++)
        {
            n = neighbors::neighs[pos][i];                                               // Get its position取近邻位置
            if ((state_be != spins[n]) && (spins[n] > -1) && (spins[n] < Parameters::Q)) // Check if it the same (remember we flipped)检测自旋是否相同注意我们在上面翻转了自旋。
            {
                // Add it to the cluster with certain probability.
                if (ran_u(gen) < pa) // uniform_int_distribution<int> ran_pos(0, SIZE - 1); Get any random integer，获取任意0到SIZE-1随机整数值，按照这个几率把它加入到cluster，进行翻转
                {
                    add_to_cluster_afm(spins, n, energy, gen, tstar, brandom, ran_pos, ran_u, pa); // 注意这里把n的取值传递给pos了
                }
            }
        }

        return;
    }

    void update_bins(double *private_m, double *private_bins, double *private_corr, double *old)
    {
        private_m[static_cast<int>(Quantity::MAG_PLANE)] += private_bins[static_cast<int>(Quantity::MAG_PLANE)];
        private_m[static_cast<int>(Quantity::MAG_PLANE2)] += private_bins[static_cast<int>(Quantity::MAG_PLANE2)];
        private_m[static_cast<int>(Quantity::MAG)] += private_bins[static_cast<int>(Quantity::MAG)];
        private_m[static_cast<int>(Quantity::MAG2)] += private_bins[static_cast<int>(Quantity::MAG2)];
        private_m[static_cast<int>(Quantity::ENE)] += private_bins[static_cast<int>(Quantity::ENE)];
        private_m[static_cast<int>(Quantity::ENE2)] += private_bins[static_cast<int>(Quantity::ENE2)];
        private_m[static_cast<int>(Quantity::ORDER)] += private_bins[static_cast<int>(Quantity::ORDER)];
        private_m[static_cast<int>(Quantity::ORDER2)] += private_bins[static_cast<int>(Quantity::ORDER2)];
        private_m[static_cast<int>(Quantity::DIS)] += private_bins[static_cast<int>(Quantity::DIS)];
        private_m[static_cast<int>(Quantity::DIS2)] += private_bins[static_cast<int>(Quantity::DIS2)];
        private_m[static_cast<int>(Quantity::AREA)] += private_bins[static_cast<int>(Quantity::AREA)];
        private_m[static_cast<int>(Quantity::AREA2)] += private_bins[static_cast<int>(Quantity::AREA2)];
        private_m[static_cast<int>(Quantity::EULER)] += private_bins[static_cast<int>(Quantity::EULER)];
        private_m[static_cast<int>(Quantity::EULER2)] += private_bins[static_cast<int>(Quantity::EULER2)];
        private_m[static_cast<int>(Quantity::SINGLE_VALUE)] += private_bins[static_cast<int>(Quantity::SINGLE_VALUE)];
        private_m[static_cast<int>(Quantity::SINGLE_VALUE2)] += private_bins[static_cast<int>(Quantity::SINGLE_VALUE2)];
        private_m[static_cast<int>(Quantity::NUCLEAR_NORM)] += private_bins[static_cast<int>(Quantity::NUCLEAR_NORM)];
        private_m[static_cast<int>(Quantity::NUCLEAR_NORM2)] += private_bins[static_cast<int>(Quantity::NUCLEAR_NORM2)];
        private_m[static_cast<int>(Quantity::SHANNON)] += private_bins[static_cast<int>(Quantity::SHANNON)];
        private_m[static_cast<int>(Quantity::SHANNON2)] += private_bins[static_cast<int>(Quantity::SHANNON2)];
        private_m[static_cast<int>(Quantity::MAG4)] += private_bins[static_cast<int>(Quantity::MAG4)];

        // 单个bin计算的磁化率，samplesize个bin的的磁化率，磁化强度平方，磁化率平方总和
        private_bins[static_cast<int>(Quantity::MAG_PLANE2)] = (private_bins[static_cast<int>(Quantity::MAG_PLANE2)] - private_bins[static_cast<int>(Quantity::MAG_PLANE)] * private_bins[static_cast<int>(Quantity::MAG_PLANE)]) * SIZE;
        private_bins[static_cast<int>(Quantity::MAG_PLANE2_bin_sus)] += private_bins[static_cast<int>(Quantity::MAG_PLANE2)];
        private_bins[static_cast<int>(Quantity::MAG_PLANERR)] += (private_bins[static_cast<int>(Quantity::MAG_PLANE)] * private_bins[static_cast<int>(Quantity::MAG_PLANE)]);
        private_bins[static_cast<int>(Quantity::MAG_PLANE2RR)] += (private_bins[static_cast<int>(Quantity::MAG_PLANE2)] * private_bins[static_cast<int>(Quantity::MAG_PLANE2)]);

        private_bins[static_cast<int>(Quantity::MAG2)] = (private_bins[static_cast<int>(Quantity::MAG2)] - private_bins[static_cast<int>(Quantity::MAG)] * private_bins[static_cast<int>(Quantity::MAG)]) * SIZE;
        private_bins[static_cast<int>(Quantity::MAG2_bin_sus)] += private_bins[static_cast<int>(Quantity::MAG2)];
        private_bins[static_cast<int>(Quantity::MAGRR)] += (private_bins[static_cast<int>(Quantity::MAG)] * private_bins[static_cast<int>(Quantity::MAG)]);
        private_bins[static_cast<int>(Quantity::MAG2RR)] += (private_bins[static_cast<int>(Quantity::MAG2)] * private_bins[static_cast<int>(Quantity::MAG2)]);

        private_bins[static_cast<int>(Quantity::ENE2)] = (private_bins[static_cast<int>(Quantity::ENE2)] - private_bins[static_cast<int>(Quantity::ENE)] * private_bins[static_cast<int>(Quantity::ENE)]) * SIZE;
        private_bins[static_cast<int>(Quantity::ENE2_bin_sus)] += private_bins[static_cast<int>(Quantity::ENE2)];
        private_bins[static_cast<int>(Quantity::ENERR)] += (private_bins[static_cast<int>(Quantity::ENE)] * private_bins[static_cast<int>(Quantity::ENE)]);
        private_bins[static_cast<int>(Quantity::ENE2RR)] += (private_bins[static_cast<int>(Quantity::ENE2)] * private_bins[static_cast<int>(Quantity::ENE2)]);

        private_bins[static_cast<int>(Quantity::ORDER2)] = (private_bins[static_cast<int>(Quantity::ORDER2)] - private_bins[static_cast<int>(Quantity::ORDER)] * private_bins[static_cast<int>(Quantity::ORDER)]) * SIZE;
        private_bins[static_cast<int>(Quantity::ORDER2_bin_sus)] += private_bins[static_cast<int>(Quantity::ORDER2)];
        private_bins[static_cast<int>(Quantity::ORDERRR)] += (private_bins[static_cast<int>(Quantity::ORDER)] * private_bins[static_cast<int>(Quantity::ORDER)]);
        private_bins[static_cast<int>(Quantity::ORDER2RR)] += (private_bins[static_cast<int>(Quantity::ORDER2)] * private_bins[static_cast<int>(Quantity::ORDER2)]);

        private_bins[static_cast<int>(Quantity::DIS2)] = (private_bins[static_cast<int>(Quantity::DIS2)] - private_bins[static_cast<int>(Quantity::DIS)] * private_bins[static_cast<int>(Quantity::DIS)]) * SIZE;
        private_bins[static_cast<int>(Quantity::DIS2_bin_sus)] += private_bins[static_cast<int>(Quantity::DIS2)];
        private_bins[static_cast<int>(Quantity::DISRR)] += (private_bins[static_cast<int>(Quantity::DIS)] * private_bins[static_cast<int>(Quantity::DIS)]);
        private_bins[static_cast<int>(Quantity::DIS2RR)] += (private_bins[static_cast<int>(Quantity::DIS2)] * private_bins[static_cast<int>(Quantity::DIS2)]);

        private_bins[static_cast<int>(Quantity::AREA2)] = (private_bins[static_cast<int>(Quantity::AREA2)] - private_bins[static_cast<int>(Quantity::AREA)] * private_bins[static_cast<int>(Quantity::AREA)]);
        private_bins[static_cast<int>(Quantity::AREA2_bin_sus)] += private_bins[static_cast<int>(Quantity::AREA2)];
        private_bins[static_cast<int>(Quantity::AREARR)] += (private_bins[static_cast<int>(Quantity::AREA)] * private_bins[static_cast<int>(Quantity::AREA)]);
        private_bins[static_cast<int>(Quantity::AREA2RR)] += (private_bins[static_cast<int>(Quantity::AREA2)] * private_bins[static_cast<int>(Quantity::AREA2)]);

        private_bins[static_cast<int>(Quantity::EULER2)] = (private_bins[static_cast<int>(Quantity::EULER2)] - private_bins[static_cast<int>(Quantity::EULER)] * private_bins[static_cast<int>(Quantity::EULER)]);
        private_bins[static_cast<int>(Quantity::EULER2_bin_sus)] += private_bins[static_cast<int>(Quantity::EULER2)];
        private_bins[static_cast<int>(Quantity::EULERRR)] += (private_bins[static_cast<int>(Quantity::EULER)] * private_bins[static_cast<int>(Quantity::EULER)]);
        private_bins[static_cast<int>(Quantity::EULER2RR)] += (private_bins[static_cast<int>(Quantity::EULER2)] * private_bins[static_cast<int>(Quantity::EULER2)]);

        private_bins[static_cast<int>(Quantity::SINGLE_VALUE2)] = (private_bins[static_cast<int>(Quantity::SINGLE_VALUE2)] - private_bins[static_cast<int>(Quantity::SINGLE_VALUE)] * private_bins[static_cast<int>(Quantity::SINGLE_VALUE)]) * SIZE;
        private_bins[static_cast<int>(Quantity::SINGLE_VALUE2_bin_sus)] += private_bins[static_cast<int>(Quantity::SINGLE_VALUE2)];
        private_bins[static_cast<int>(Quantity::SINGLE_VALUERR)] += (private_bins[static_cast<int>(Quantity::SINGLE_VALUE)] * private_bins[static_cast<int>(Quantity::SINGLE_VALUE)]);
        private_bins[static_cast<int>(Quantity::SINGLE_VALUE2RR)] += (private_bins[static_cast<int>(Quantity::SINGLE_VALUE2)] * private_bins[static_cast<int>(Quantity::SINGLE_VALUE2)]);

        private_bins[static_cast<int>(Quantity::NUCLEAR_NORM2)] = (private_bins[static_cast<int>(Quantity::NUCLEAR_NORM2)] - private_bins[static_cast<int>(Quantity::NUCLEAR_NORM)] * private_bins[static_cast<int>(Quantity::NUCLEAR_NORM)]) * SIZE;
        private_bins[static_cast<int>(Quantity::NUCLEAR_NORM2_bin_sus)] += private_bins[static_cast<int>(Quantity::NUCLEAR_NORM2)];
        private_bins[static_cast<int>(Quantity::NUCLEAR_NORMRR)] += (private_bins[static_cast<int>(Quantity::NUCLEAR_NORM)] * private_bins[static_cast<int>(Quantity::NUCLEAR_NORM)]);
        private_bins[static_cast<int>(Quantity::NUCLEAR_NORM2RR)] += (private_bins[static_cast<int>(Quantity::NUCLEAR_NORM2)] * private_bins[static_cast<int>(Quantity::NUCLEAR_NORM2)]);

        private_bins[static_cast<int>(Quantity::SHANNON2)] = (private_bins[static_cast<int>(Quantity::SHANNON2)] - private_bins[static_cast<int>(Quantity::SHANNON)] * private_bins[static_cast<int>(Quantity::SHANNON)]) * SIZE;
        private_bins[static_cast<int>(Quantity::SHANNON2_bin_sus)] += private_bins[static_cast<int>(Quantity::SHANNON2)];
        private_bins[static_cast<int>(Quantity::SHANNONRR)] += (private_bins[static_cast<int>(Quantity::SHANNON)] * private_bins[static_cast<int>(Quantity::SHANNON)]);
        private_bins[static_cast<int>(Quantity::SHANNON2RR)] += (private_bins[static_cast<int>(Quantity::SHANNON2)] * private_bins[static_cast<int>(Quantity::SHANNON2)]);

        private_bins[static_cast<int>(Quantity::MAG4RR)] += (private_bins[static_cast<int>(Quantity::MAG4)] * private_bins[static_cast<int>(Quantity::MAG4)]);

        private_corr[static_cast<int>(Quantity::MAG_PLANE)] += private_bins[static_cast<int>(Quantity::MAG_PLANE)] * old[static_cast<int>(Quantity::MAG_PLANE)];
        private_corr[static_cast<int>(Quantity::MAG_PLANE2)] += private_bins[static_cast<int>(Quantity::MAG_PLANE2)] * old[static_cast<int>(Quantity::MAG_PLANE2)];
        private_corr[static_cast<int>(Quantity::MAG)] += private_bins[static_cast<int>(Quantity::MAG)] * old[static_cast<int>(Quantity::MAG)];
        private_corr[static_cast<int>(Quantity::MAG2)] += private_bins[static_cast<int>(Quantity::MAG2)] * old[static_cast<int>(Quantity::MAG2)];
        private_corr[static_cast<int>(Quantity::ENE)] += private_bins[static_cast<int>(Quantity::ENE)] * old[static_cast<int>(Quantity::ENE)];
        private_corr[static_cast<int>(Quantity::ENE2)] += private_bins[static_cast<int>(Quantity::ENE2)] * old[static_cast<int>(Quantity::ENE2)];
        private_corr[static_cast<int>(Quantity::ORDER)] += private_bins[static_cast<int>(Quantity::ORDER)] * old[static_cast<int>(Quantity::ORDER)];
        private_corr[static_cast<int>(Quantity::ORDER2)] += private_bins[static_cast<int>(Quantity::ORDER2)] * old[static_cast<int>(Quantity::ORDER2)];
        private_corr[static_cast<int>(Quantity::DIS)] += private_bins[static_cast<int>(Quantity::DIS)] * old[static_cast<int>(Quantity::DIS)];
        private_corr[static_cast<int>(Quantity::DIS2)] += private_bins[static_cast<int>(Quantity::DIS2)] * old[static_cast<int>(Quantity::DIS2)];
        private_corr[static_cast<int>(Quantity::AREA)] += private_bins[static_cast<int>(Quantity::AREA)] * old[static_cast<int>(Quantity::AREA)];
        private_corr[static_cast<int>(Quantity::AREA2)] += private_bins[static_cast<int>(Quantity::AREA2)] * old[static_cast<int>(Quantity::AREA2)];
        private_corr[static_cast<int>(Quantity::EULER)] += private_bins[static_cast<int>(Quantity::EULER)] * old[static_cast<int>(Quantity::EULER)];
        private_corr[static_cast<int>(Quantity::EULER2)] += private_bins[static_cast<int>(Quantity::EULER2)] * old[static_cast<int>(Quantity::EULER2)];
        private_corr[static_cast<int>(Quantity::SINGLE_VALUE)] += private_bins[static_cast<int>(Quantity::SINGLE_VALUE)] * old[static_cast<int>(Quantity::SINGLE_VALUE)];
        private_corr[static_cast<int>(Quantity::SINGLE_VALUE2)] += private_bins[static_cast<int>(Quantity::SINGLE_VALUE2)] * old[static_cast<int>(Quantity::SINGLE_VALUE2)];
        private_corr[static_cast<int>(Quantity::NUCLEAR_NORM)] += private_bins[static_cast<int>(Quantity::NUCLEAR_NORM)] * old[static_cast<int>(Quantity::NUCLEAR_NORM)];
        private_corr[static_cast<int>(Quantity::NUCLEAR_NORM2)] += private_bins[static_cast<int>(Quantity::NUCLEAR_NORM2)] * old[static_cast<int>(Quantity::NUCLEAR_NORM2)];
        private_corr[static_cast<int>(Quantity::SHANNON)] += private_bins[static_cast<int>(Quantity::SHANNON)] * old[static_cast<int>(Quantity::SHANNON)];
        private_corr[static_cast<int>(Quantity::SHANNON2)] += private_bins[static_cast<int>(Quantity::SHANNON2)] * old[static_cast<int>(Quantity::SHANNON2)];
        private_corr[static_cast<int>(Quantity::MAG4)] += private_bins[static_cast<int>(Quantity::MAG4)] * old[static_cast<int>(Quantity::MAG4)];

        old[static_cast<int>(Quantity::MAG_PLANE)] = private_bins[static_cast<int>(Quantity::MAG_PLANE)];
        old[static_cast<int>(Quantity::MAG_PLANE2)] = private_bins[static_cast<int>(Quantity::MAG_PLANE2)];
        old[static_cast<int>(Quantity::MAG)] = private_bins[static_cast<int>(Quantity::MAG)];
        old[static_cast<int>(Quantity::MAG2)] = private_bins[static_cast<int>(Quantity::MAG2)];
        old[static_cast<int>(Quantity::ENE)] = private_bins[static_cast<int>(Quantity::ENE)];
        old[static_cast<int>(Quantity::ENE2)] = private_bins[static_cast<int>(Quantity::ENE2)];
        old[static_cast<int>(Quantity::ORDER)] = private_bins[static_cast<int>(Quantity::ORDER)];
        old[static_cast<int>(Quantity::ORDER2)] = private_bins[static_cast<int>(Quantity::ORDER2)];
        old[static_cast<int>(Quantity::DIS)] = private_bins[static_cast<int>(Quantity::DIS)];
        old[static_cast<int>(Quantity::DIS2)] = private_bins[static_cast<int>(Quantity::DIS2)];
        old[static_cast<int>(Quantity::AREA)] = private_bins[static_cast<int>(Quantity::AREA)];
        old[static_cast<int>(Quantity::AREA2)] = private_bins[static_cast<int>(Quantity::AREA2)];
        old[static_cast<int>(Quantity::EULER)] = private_bins[static_cast<int>(Quantity::EULER)];
        old[static_cast<int>(Quantity::EULER2)] = private_bins[static_cast<int>(Quantity::EULER2)];
        old[static_cast<int>(Quantity::SINGLE_VALUE)] = private_bins[static_cast<int>(Quantity::SINGLE_VALUE)];
        old[static_cast<int>(Quantity::SINGLE_VALUE2)] = private_bins[static_cast<int>(Quantity::SINGLE_VALUE2)];
        old[static_cast<int>(Quantity::NUCLEAR_NORM)] = private_bins[static_cast<int>(Quantity::NUCLEAR_NORM)];
        old[static_cast<int>(Quantity::NUCLEAR_NORM2)] = private_bins[static_cast<int>(Quantity::NUCLEAR_NORM2)];
        old[static_cast<int>(Quantity::SHANNON)] = private_bins[static_cast<int>(Quantity::SHANNON)];
        old[static_cast<int>(Quantity::SHANNON2)] = private_bins[static_cast<int>(Quantity::SHANNON2)];
        old[static_cast<int>(Quantity::MAG4)] = private_bins[static_cast<int>(Quantity::MAG4)];

        // private_corr[static_cast<int>(Quantity::MAG_PLANE)] += private_bins[static_cast<int>(Quantity::MAG_PLANE)] * private_bins[static_cast<int>(Quantity::MAG_PLANE)];
        // private_corr[static_cast<int>(Quantity::MAG_PLANE2)] += private_bins[static_cast<int>(Quantity::MAG_PLANE2)] * private_bins[static_cast<int>(Quantity::MAG_PLANE2)];
        // private_corr[static_cast<int>(Quantity::MAG)] += private_bins[static_cast<int>(Quantity::MAG)] * private_bins[static_cast<int>(Quantity::MAG)];
        // private_corr[static_cast<int>(Quantity::MAG2)] += private_bins[static_cast<int>(Quantity::MAG2)] * private_bins[static_cast<int>(Quantity::MAG2)];
        // private_corr[static_cast<int>(Quantity::ENE)] += private_bins[static_cast<int>(Quantity::ENE)] * private_bins[static_cast<int>(Quantity::ENE)];
        // private_corr[static_cast<int>(Quantity::ENE2)] += private_bins[static_cast<int>(Quantity::ENE2)] * private_bins[static_cast<int>(Quantity::ENE2)];
        // private_corr[static_cast<int>(Quantity::ORDER)] += private_bins[static_cast<int>(Quantity::ORDER)] * private_bins[static_cast<int>(Quantity::ORDER)];
        // private_corr[static_cast<int>(Quantity::ORDER2)] += private_bins[static_cast<int>(Quantity::ORDER2)] * private_bins[static_cast<int>(Quantity::ORDER2)];
        // private_corr[static_cast<int>(Quantity::DIS)] += private_bins[static_cast<int>(Quantity::DIS)] * private_bins[static_cast<int>(Quantity::DIS)];
        // private_corr[static_cast<int>(Quantity::DIS2)] += private_bins[static_cast<int>(Quantity::DIS2)] * private_bins[static_cast<int>(Quantity::DIS2)];
        // private_corr[static_cast<int>(Quantity::AREA)] += private_bins[static_cast<int>(Quantity::AREA)] * private_bins[static_cast<int>(Quantity::AREA)];
        // private_corr[static_cast<int>(Quantity::AREA2)] += private_bins[static_cast<int>(Quantity::AREA2)] * private_bins[static_cast<int>(Quantity::AREA2)];
        // private_corr[static_cast<int>(Quantity::EULER)] += private_bins[static_cast<int>(Quantity::EULER)] * private_bins[static_cast<int>(Quantity::EULER)];
        // private_corr[static_cast<int>(Quantity::EULER2)] += private_bins[static_cast<int>(Quantity::EULER2)] * private_bins[static_cast<int>(Quantity::EULER2)];
        // private_corr[static_cast<int>(Quantity::SINGLE_VALUE)] += private_bins[static_cast<int>(Quantity::SINGLE_VALUE)] * private_bins[static_cast<int>(Quantity::SINGLE_VALUE)];
        // private_corr[static_cast<int>(Quantity::SINGLE_VALUE2)] += private_bins[static_cast<int>(Quantity::SINGLE_VALUE2)] * private_bins[static_cast<int>(Quantity::SINGLE_VALUE2)];
        // private_corr[static_cast<int>(Quantity::NUCLEAR_NORM)] += private_bins[static_cast<int>(Quantity::NUCLEAR_NORM)] * private_bins[static_cast<int>(Quantity::NUCLEAR_NORM)];
        // private_corr[static_cast<int>(Quantity::NUCLEAR_NORM2)] += private_bins[static_cast<int>(Quantity::NUCLEAR_NORM2)] * private_bins[static_cast<int>(Quantity::NUCLEAR_NORM2)];
        // private_corr[static_cast<int>(Quantity::SHANNON)] += private_bins[static_cast<int>(Quantity::SHANNON)] * private_bins[static_cast<int>(Quantity::SHANNON)];
        // private_corr[static_cast<int>(Quantity::SHANNON2)] += private_bins[static_cast<int>(Quantity::SHANNON2)] * private_bins[static_cast<int>(Quantity::SHANNON2)];
        // private_corr[static_cast<int>(Quantity::MAG4)] += private_bins[static_cast<int>(Quantity::MAG4)] * private_bins[static_cast<int>(Quantity::MAG4)];
    }

    void count_corr(double *private_m, double *private_bins, double *private_corr)
    {
        // MAG_PLANE
        double MAG_PLANE_square = private_m[static_cast<int>(Quantity::MAG_PLANE)] * private_m[static_cast<int>(Quantity::MAG_PLANE)];
        double MAG_PLANE_error = (private_bins[static_cast<int>(Quantity::MAG_PLANERR)] - MAG_PLANE_square);
        private_corr[static_cast<int>(Quantity::MAG_PLANE)] =
            (MAG_PLANE_error != 0) ? (private_corr[static_cast<int>(Quantity::MAG_PLANE)] - MAG_PLANE_square) / MAG_PLANE_error : 0.0;

        double MAG_PLANE2_square = private_bins[static_cast<int>(Quantity::MAG_PLANE2_bin_sus)] * private_bins[static_cast<int>(Quantity::MAG_PLANE2_bin_sus)];
        double MAG_PLANE2_error = (private_bins[static_cast<int>(Quantity::MAG_PLANE2RR)] - MAG_PLANE2_square);
        private_corr[static_cast<int>(Quantity::MAG_PLANE2)] =
            (MAG_PLANE2_error != 0) ? (private_corr[static_cast<int>(Quantity::MAG_PLANE2)] - MAG_PLANE2_square) / MAG_PLANE2_error : 0.0;

        // MAG
        double MAG_square = private_m[static_cast<int>(Quantity::MAG)] * private_m[static_cast<int>(Quantity::MAG)];
        double MAG_error = (private_bins[static_cast<int>(Quantity::MAGRR)] - MAG_square);
        private_corr[static_cast<int>(Quantity::MAG)] =
            (MAG_error != 0) ? (private_corr[static_cast<int>(Quantity::MAG)] - MAG_square) / MAG_error : 0.0;

        double MAG2_square = private_bins[static_cast<int>(Quantity::MAG2_bin_sus)] * private_bins[static_cast<int>(Quantity::MAG2_bin_sus)];
        double MAG2_error = (private_bins[static_cast<int>(Quantity::MAG2RR)] - MAG2_square);
        private_corr[static_cast<int>(Quantity::MAG2)] =
            (MAG2_error != 0) ? (private_corr[static_cast<int>(Quantity::MAG2)] - MAG2_square) / MAG2_error : 0.0;

        // ENE
        double ENE_square = private_m[static_cast<int>(Quantity::ENE)] * private_m[static_cast<int>(Quantity::ENE)];
        double ENE_error = (private_bins[static_cast<int>(Quantity::ENERR)] - ENE_square);
        private_corr[static_cast<int>(Quantity::ENE)] =
            (ENE_error != 0) ? (private_corr[static_cast<int>(Quantity::ENE)] - ENE_square) / ENE_error : 0.0;

        double ENE2_square = private_bins[static_cast<int>(Quantity::ENE2_bin_sus)] * private_bins[static_cast<int>(Quantity::ENE2_bin_sus)];
        double ENE2_error = (private_bins[static_cast<int>(Quantity::ENE2RR)] - ENE2_square);
        private_corr[static_cast<int>(Quantity::ENE2)] =
            (ENE2_error != 0) ? (private_corr[static_cast<int>(Quantity::ENE2)] - ENE2_square) / ENE2_error : 0.0;

        // ORDER
        double ORDER_square = private_m[static_cast<int>(Quantity::ORDER)] * private_m[static_cast<int>(Quantity::ORDER)];
        double ORDER_error = (private_bins[static_cast<int>(Quantity::ORDERRR)] - ORDER_square);
        private_corr[static_cast<int>(Quantity::ORDER)] =
            (ORDER_error != 0) ? (private_corr[static_cast<int>(Quantity::ORDER)] - ORDER_square) / ORDER_error : 0.0;

        double ORDER2_square = private_bins[static_cast<int>(Quantity::ORDER2_bin_sus)] * private_bins[static_cast<int>(Quantity::ORDER2_bin_sus)];
        double ORDER2_error = (private_bins[static_cast<int>(Quantity::ORDER2RR)] - ORDER2_square);
        private_corr[static_cast<int>(Quantity::ORDER2)] =
            (ORDER2_error != 0) ? (private_corr[static_cast<int>(Quantity::ORDER2)] - ORDER2_square) / ORDER2_error : 0.0;

        // DIS
        double DIS_square = private_m[static_cast<int>(Quantity::DIS)] * private_m[static_cast<int>(Quantity::DIS)];
        double DIS_error = (private_bins[static_cast<int>(Quantity::DISRR)] - DIS_square);
        private_corr[static_cast<int>(Quantity::DIS)] =
            (DIS_error != 0) ? (private_corr[static_cast<int>(Quantity::DIS)] - DIS_square) / DIS_error : 0.0;

        double DIS2_square = private_bins[static_cast<int>(Quantity::DIS2_bin_sus)] * private_bins[static_cast<int>(Quantity::DIS2_bin_sus)];
        double DIS2_error = (private_bins[static_cast<int>(Quantity::DIS2RR)] - DIS2_square);
        private_corr[static_cast<int>(Quantity::DIS2)] =
            (DIS2_error != 0) ? (private_corr[static_cast<int>(Quantity::DIS2)] - DIS2_square) / DIS2_error : 0.0;

        // AREA
        double AREA_square = private_m[static_cast<int>(Quantity::AREA)] * private_m[static_cast<int>(Quantity::AREA)];
        double AREA_error = (private_bins[static_cast<int>(Quantity::AREARR)] - AREA_square);
        private_corr[static_cast<int>(Quantity::AREA)] =
            (AREA_error != 0) ? (private_corr[static_cast<int>(Quantity::AREA)] - AREA_square) / AREA_error : 0.0;

        double AREA2_square = private_bins[static_cast<int>(Quantity::AREA2_bin_sus)] * private_bins[static_cast<int>(Quantity::AREA2_bin_sus)];
        double AREA2_error = (private_bins[static_cast<int>(Quantity::AREA2RR)] - AREA2_square);
        private_corr[static_cast<int>(Quantity::AREA2)] =
            (AREA2_error != 0) ? (private_corr[static_cast<int>(Quantity::AREA2)] - AREA2_square) / AREA2_error : 0.0;

        // EULER
        double EULER_square = private_m[static_cast<int>(Quantity::EULER)] * private_m[static_cast<int>(Quantity::EULER)];
        double EULER_error = (private_bins[static_cast<int>(Quantity::EULERRR)] - EULER_square);
        private_corr[static_cast<int>(Quantity::EULER)] =
            (EULER_error != 0) ? (private_corr[static_cast<int>(Quantity::EULER)] - EULER_square) / EULER_error : 0.0;

        double EULER2_square = private_bins[static_cast<int>(Quantity::EULER2_bin_sus)] * private_bins[static_cast<int>(Quantity::EULER2_bin_sus)];
        double EULER2_error = (private_bins[static_cast<int>(Quantity::EULER2RR)] - EULER2_square);
        private_corr[static_cast<int>(Quantity::EULER2)] =
            (EULER2_error != 0) ? (private_corr[static_cast<int>(Quantity::EULER2)] - EULER2_square) / EULER2_error : 0.0;

        // SINGLE_VALUE
        double SINGLE_VALUE_square = private_m[static_cast<int>(Quantity::SINGLE_VALUE)] * private_m[static_cast<int>(Quantity::SINGLE_VALUE)];
        double SINGLE_VALUE_error = (private_bins[static_cast<int>(Quantity::SINGLE_VALUERR)] - SINGLE_VALUE_square);
        // std::cout << private_corr[static_cast<int>(Quantity::SINGLE_VALUE)] << " " << private_bins[static_cast<int>(Quantity::SINGLE_VALUERR)] << std::endl;
        private_corr[static_cast<int>(Quantity::SINGLE_VALUE)] =
            (SINGLE_VALUE_error != 0) ? (private_corr[static_cast<int>(Quantity::SINGLE_VALUE)] - SINGLE_VALUE_square) / SINGLE_VALUE_error : 0.0;

        double SINGLE_VALUE2_square = private_bins[static_cast<int>(Quantity::SINGLE_VALUE2_bin_sus)] * private_bins[static_cast<int>(Quantity::SINGLE_VALUE2_bin_sus)];
        double SINGLE_VALUE2_error = (private_bins[static_cast<int>(Quantity::SINGLE_VALUE2RR)] - SINGLE_VALUE2_square);
        private_corr[static_cast<int>(Quantity::SINGLE_VALUE2)] =
            (SINGLE_VALUE2_error != 0) ? (private_corr[static_cast<int>(Quantity::SINGLE_VALUE2)] - SINGLE_VALUE2_square) / SINGLE_VALUE2_error : 0.0;

        // NUCLEAR_NORM
        double NUCLEAR_NORM_square = private_m[static_cast<int>(Quantity::NUCLEAR_NORM)] * private_m[static_cast<int>(Quantity::NUCLEAR_NORM)];
        double NUCLEAR_NORM_error = (private_bins[static_cast<int>(Quantity::NUCLEAR_NORMRR)] - NUCLEAR_NORM_square);
        private_corr[static_cast<int>(Quantity::NUCLEAR_NORM)] =
            (NUCLEAR_NORM_error != 0) ? (private_corr[static_cast<int>(Quantity::NUCLEAR_NORM)] - NUCLEAR_NORM_square) / NUCLEAR_NORM_error : 0.0;

        double NUCLEAR_NORM2_square = private_bins[static_cast<int>(Quantity::NUCLEAR_NORM2_bin_sus)] * private_bins[static_cast<int>(Quantity::NUCLEAR_NORM2_bin_sus)];
        double NUCLEAR_NORM2_error = (private_bins[static_cast<int>(Quantity::NUCLEAR_NORM2RR)] - NUCLEAR_NORM2_square);
        private_corr[static_cast<int>(Quantity::NUCLEAR_NORM2)] =
            (NUCLEAR_NORM2_error != 0) ? (private_corr[static_cast<int>(Quantity::NUCLEAR_NORM2)] - NUCLEAR_NORM2_square) / NUCLEAR_NORM2_error : 0.0;

        // SHANNON
        double SHANNON_square = private_m[static_cast<int>(Quantity::SHANNON)] * private_m[static_cast<int>(Quantity::SHANNON)];
        double SHANNON_error = (private_bins[static_cast<int>(Quantity::SHANNONRR)] - SHANNON_square);
        private_corr[static_cast<int>(Quantity::SHANNON)] =
            (SHANNON_error != 0) ? (private_corr[static_cast<int>(Quantity::SHANNON)] - SHANNON_square) / SHANNON_error : 0.0;

        double SHANNON2_square = private_bins[static_cast<int>(Quantity::SHANNON2_bin_sus)] * private_bins[static_cast<int>(Quantity::SHANNON2_bin_sus)];
        double SHANNON2_error = (private_bins[static_cast<int>(Quantity::SHANNON2RR)] - SHANNON2_square);
        private_corr[static_cast<int>(Quantity::SHANNON2)] =
            (SHANNON2_error != 0) ? (private_corr[static_cast<int>(Quantity::SHANNON2)] - SHANNON2_square) / SHANNON2_error : 0.0;

        // MAG4
        double MAG4_square = private_m[static_cast<int>(Quantity::MAG4)] * private_m[static_cast<int>(Quantity::MAG4)];
        double MAG4_error = (private_bins[static_cast<int>(Quantity::MAG4RR)] - MAG4_square);
        private_corr[static_cast<int>(Quantity::MAG4)] =
            (MAG4_error != 0) ? (private_corr[static_cast<int>(Quantity::MAG4)] - MAG4_square) / MAG4_error : 0.0;
    }
}