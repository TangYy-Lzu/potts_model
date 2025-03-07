#include "mc_step.hpp"

namespace mc_step
{
    void do_step(int spins[SIZE + 1], double tstar, double &energy, std::mt19937 &gen,
                 std::uniform_int_distribution<int> &brandom, std::uniform_int_distribution<int> ran_pos, std::uniform_real_distribution<double> &ran_u,
                 double *private_m, double *private_bins)
    {
        int i, j, k;
        double mag, mag_plane, dist, order, n_area, euler_number, single_value, nuclera_norm, shannon;
        double energy_chi, mag_chi, mag_plane_chi, dist_chi, order_chi, n_area_chi, euler_number_chi, single_value_chi, nuclera_norm_chi, shannon_chi;

        Initialization::initialize_matrices(private_m, private_bins);

        // 预热步骤
        for (i = 0; i < 200 * Parameters::Q; i++)
        {
            for (j = 0; j < SIZE; j++)
            {
                flip_spin(spins, energy, gen, tstar, brandom, ran_pos, ran_u);
            }
        }

        for (i = 0; i < Parameters::SAMPLE_SIZE; i++)
        {
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

            for (j = 0; j < DATA; j++)
            {
                private_bins[j] /= (1.0 * Parameters::GROUP_SIZE);
            };

            update_bins(private_bins);
        }

        for (i = 0; i < DATA; i++)
        {
            private_m[i] = private_m[i] / (Parameters::SAMPLE_SIZE * Parameters::GROUP_SIZE);
        }
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
    void do_step_woff(int spins[SIZE + 1], double tstar, double &energy, std::mt19937 &gen,
                      std::uniform_int_distribution<int> &brandom, std::uniform_int_distribution<int> ran_pos, std::uniform_real_distribution<double> &ran_u,
                      double *private_m, double *private_bins)
    {
        int i, j, k;
        double mag, mag_plane, dist, order, n_area, euler_number, single_value, nuclera_norm, shannon;
        double energy_chi, mag_chi, mag_plane_chi, dist_chi, order_chi, n_area_chi, euler_number_chi, single_value_chi, nuclera_norm_chi, shannon_chi;

        Initialization::initialize_matrices(private_m, private_bins);

        double pa = 1.0 - exp(Parameters::J / tstar); // TODO change添加种子的概率，这里用温度代替了β，tstar=1/kβ，J在计算能量处选为1

        // 预热步骤
        for (i = 0; i < Parameters::Q; i++)
        {
            for (j = 0; j < SIZE; j++)
            {
                area(spins, energy, gen, tstar, brandom, ran_pos, ran_u, pa);
            }
        }

        for (i = 0; i < Parameters::SAMPLE_SIZE; i++)
        {
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

            for (j = 0; j < DATA; j++)
            {
                private_bins[j] /= (1.0 * Parameters::GROUP_SIZE);
            };

            update_bins(private_bins);
        }

        for (i = 0; i < DATA; i++)
        {
            private_m[i] = private_m[i] / (Parameters::SAMPLE_SIZE * Parameters::GROUP_SIZE);
        }
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

    void update_bins(double *private_bins)
    {
        private_bins[static_cast<int>(Quantity::MAG_PLANE2)] = (private_bins[static_cast<int>(Quantity::MAG_PLANE2)] - private_bins[static_cast<int>(Quantity::MAG_PLANE)] * private_bins[static_cast<int>(Quantity::MAG_PLANE)]) * SIZE;
        private_bins[static_cast<int>(Quantity::MAG_PLANE2_bin)] += private_bins[static_cast<int>(Quantity::MAG_PLANE2)];
        private_bins[static_cast<int>(Quantity::MAG_PLANERR)] += (private_bins[static_cast<int>(Quantity::MAG_PLANE)] * private_bins[static_cast<int>(Quantity::MAG_PLANE)]);
        private_bins[static_cast<int>(Quantity::MAG_PLANE2RR)] += (private_bins[static_cast<int>(Quantity::MAG_PLANE2)] * private_bins[static_cast<int>(Quantity::MAG_PLANE2)]);

        private_bins[static_cast<int>(Quantity::MAG2)] = (private_bins[static_cast<int>(Quantity::MAG2)] - private_bins[static_cast<int>(Quantity::MAG)] * private_bins[static_cast<int>(Quantity::MAG)]) * SIZE;
        private_bins[static_cast<int>(Quantity::MAG2_bin)] += private_bins[static_cast<int>(Quantity::MAG2)];
        private_bins[static_cast<int>(Quantity::MAGRR)] += (private_bins[static_cast<int>(Quantity::MAG)] * private_bins[static_cast<int>(Quantity::MAG)]);
        private_bins[static_cast<int>(Quantity::MAG2RR)] += (private_bins[static_cast<int>(Quantity::MAG2)] * private_bins[static_cast<int>(Quantity::MAG2)]);

        private_bins[static_cast<int>(Quantity::ENE2)] = (private_bins[static_cast<int>(Quantity::ENE2)] - private_bins[static_cast<int>(Quantity::ENE)] * private_bins[static_cast<int>(Quantity::ENE)]) * SIZE;
        private_bins[static_cast<int>(Quantity::ENE2_bin)] += private_bins[static_cast<int>(Quantity::ENE2)];
        private_bins[static_cast<int>(Quantity::ENERR)] += (private_bins[static_cast<int>(Quantity::ENE)] * private_bins[static_cast<int>(Quantity::ENE)]);
        private_bins[static_cast<int>(Quantity::ENE2RR)] += (private_bins[static_cast<int>(Quantity::ENE2)] * private_bins[static_cast<int>(Quantity::ENE2)]);

        private_bins[static_cast<int>(Quantity::ORDER2)] = (private_bins[static_cast<int>(Quantity::ORDER2)] - private_bins[static_cast<int>(Quantity::ORDER)] * private_bins[static_cast<int>(Quantity::ORDER)]) * SIZE;
        private_bins[static_cast<int>(Quantity::ORDER2_bin)] += private_bins[static_cast<int>(Quantity::ORDER2)];
        private_bins[static_cast<int>(Quantity::ORDERRR)] += (private_bins[static_cast<int>(Quantity::ORDER)] * private_bins[static_cast<int>(Quantity::ORDER)]);
        private_bins[static_cast<int>(Quantity::ORDER2RR)] += (private_bins[static_cast<int>(Quantity::ORDER2)] * private_bins[static_cast<int>(Quantity::ORDER2)]);

        private_bins[static_cast<int>(Quantity::DIS2)] = (private_bins[static_cast<int>(Quantity::DIS2)] - private_bins[static_cast<int>(Quantity::DIS)] * private_bins[static_cast<int>(Quantity::DIS)]) * SIZE;
        private_bins[static_cast<int>(Quantity::DIS2_bin)] += private_bins[static_cast<int>(Quantity::DIS2)];
        private_bins[static_cast<int>(Quantity::DISRR)] += (private_bins[static_cast<int>(Quantity::DIS)] * private_bins[static_cast<int>(Quantity::DIS)]);
        private_bins[static_cast<int>(Quantity::DIS2RR)] += (private_bins[static_cast<int>(Quantity::DIS2)] * private_bins[static_cast<int>(Quantity::DIS2)]);

        private_bins[static_cast<int>(Quantity::AREA2)] = (private_bins[static_cast<int>(Quantity::AREA2)] - private_bins[static_cast<int>(Quantity::AREA)] * private_bins[static_cast<int>(Quantity::AREA)]);
        private_bins[static_cast<int>(Quantity::AREA2_bin)] += private_bins[static_cast<int>(Quantity::AREA2)];
        private_bins[static_cast<int>(Quantity::AREARR)] += (private_bins[static_cast<int>(Quantity::AREA)] * private_bins[static_cast<int>(Quantity::AREA)]);
        private_bins[static_cast<int>(Quantity::AREA2RR)] += (private_bins[static_cast<int>(Quantity::AREA2)] * private_bins[static_cast<int>(Quantity::AREA2)]);

        private_bins[static_cast<int>(Quantity::EULER2)] = (private_bins[static_cast<int>(Quantity::EULER2)] - private_bins[static_cast<int>(Quantity::EULER)] * private_bins[static_cast<int>(Quantity::EULER)]);
        private_bins[static_cast<int>(Quantity::EULER2_bin)] += private_bins[static_cast<int>(Quantity::EULER2)];
        private_bins[static_cast<int>(Quantity::EULERRR)] += (private_bins[static_cast<int>(Quantity::EULER)] * private_bins[static_cast<int>(Quantity::EULER)]);
        private_bins[static_cast<int>(Quantity::EULER2RR)] += (private_bins[static_cast<int>(Quantity::EULER2)] * private_bins[static_cast<int>(Quantity::EULER2)]);

        private_bins[static_cast<int>(Quantity::SINGLE_VALUE2)] = (private_bins[static_cast<int>(Quantity::SINGLE_VALUE2)] - private_bins[static_cast<int>(Quantity::SINGLE_VALUE)] * private_bins[static_cast<int>(Quantity::SINGLE_VALUE)]) * SIZE;
        private_bins[static_cast<int>(Quantity::SINGLE_VALUE2_bin)] += private_bins[static_cast<int>(Quantity::SINGLE_VALUE2)];
        private_bins[static_cast<int>(Quantity::SINGLE_VALUERR)] += (private_bins[static_cast<int>(Quantity::SINGLE_VALUE)] * private_bins[static_cast<int>(Quantity::SINGLE_VALUE)]);
        private_bins[static_cast<int>(Quantity::SINGLE_VALUE2RR)] += (private_bins[static_cast<int>(Quantity::SINGLE_VALUE2)] * private_bins[static_cast<int>(Quantity::SINGLE_VALUE2)]);

        private_bins[static_cast<int>(Quantity::NUCLEAR_NORM2)] = (private_bins[static_cast<int>(Quantity::NUCLEAR_NORM2)] - private_bins[static_cast<int>(Quantity::NUCLEAR_NORM)] * private_bins[static_cast<int>(Quantity::NUCLEAR_NORM)]) * SIZE;
        private_bins[static_cast<int>(Quantity::NUCLEAR_NORM2_bin)] += private_bins[static_cast<int>(Quantity::NUCLEAR_NORM2)];
        private_bins[static_cast<int>(Quantity::NUCLEAR_NORMRR)] += (private_bins[static_cast<int>(Quantity::NUCLEAR_NORM)] * private_bins[static_cast<int>(Quantity::NUCLEAR_NORM)]);
        private_bins[static_cast<int>(Quantity::NUCLEAR_NORM2RR)] += (private_bins[static_cast<int>(Quantity::NUCLEAR_NORM2)] * private_bins[static_cast<int>(Quantity::NUCLEAR_NORM2)]);

        private_bins[static_cast<int>(Quantity::SHANNON2)] = (private_bins[static_cast<int>(Quantity::SHANNON2)] - private_bins[static_cast<int>(Quantity::SHANNON)] * private_bins[static_cast<int>(Quantity::SHANNON)]) * SIZE;
        private_bins[static_cast<int>(Quantity::SHANNON2_bin)] += private_bins[static_cast<int>(Quantity::SHANNON2)];
        private_bins[static_cast<int>(Quantity::SHANNONRR)] += (private_bins[static_cast<int>(Quantity::SHANNON)] * private_bins[static_cast<int>(Quantity::SHANNON)]);
        private_bins[static_cast<int>(Quantity::SHANNON2RR)] += (private_bins[static_cast<int>(Quantity::SHANNON2)] * private_bins[static_cast<int>(Quantity::SHANNON2)]);

        private_bins[static_cast<int>(Quantity::MAG4RR)] += (private_bins[static_cast<int>(Quantity::MAG4)] * private_bins[static_cast<int>(Quantity::MAG4)]);
    }
}