#ifndef mc_step_hpp
#define mc_step_hpp

#include "initialize.hpp"
#include "matrix.hpp"
#include "quantity.hpp"

namespace mc_step
{
    void do_step(int spins[SIZE + 1], double tstar, double &energy, std::mt19937 &gen,
                 std::uniform_int_distribution<int> &brandom, std::uniform_int_distribution<int> ran_pos, std::uniform_real_distribution<double> &ran_u,
                 double *private_m, double *private_bins, double *private_corr);

    void flip_spin(int spins[SIZE + 1], double &energy, std::mt19937 &gen, double tstar,
                   std::uniform_int_distribution<int> &brandom, std::uniform_int_distribution<int> ran_pos, std::uniform_real_distribution<double> &ran_u);

    void do_step_wolff(int spins[SIZE + 1], double tstar, double &energy, std::mt19937 &gen,
                       std::uniform_int_distribution<int> &brandom, std::uniform_int_distribution<int> ran_pos, std::uniform_real_distribution<double> &ran_u,
                       double *private_m, double *private_bins, double *private_corr);

    void area(int spins[SIZE + 1], double &energy, std::mt19937 &gen, double tstar,
              std::uniform_int_distribution<int> &brandom, std::uniform_int_distribution<int> ran_pos, std::uniform_real_distribution<double> &ran_u, double pa);

    void add_to_cluster(int spins[SIZE + 1], int pos, int state_be, int state_af, double &energy, std::mt19937 &gen, double tstar,
                        std::uniform_int_distribution<int> &brandom, std::uniform_int_distribution<int> ran_pos, std::uniform_real_distribution<double> &ran_u, double pa);

    void do_step_wolff_afm(int spins[SIZE + 1], double tstar, double &energy, std::mt19937 &gen,
                           std::uniform_int_distribution<int> &brandom, std::uniform_int_distribution<int> ran_pos, std::uniform_real_distribution<double> &ran_u,
                           double *private_m, double *private_bins, double *private_corr);

    void area_afm(int spins[SIZE + 1], double &energy, std::mt19937 &gen, double tstar,
                  std::uniform_int_distribution<int> &brandom, std::uniform_int_distribution<int> ran_pos, std::uniform_real_distribution<double> &ran_u, double pa);

    void add_to_cluster_afm(int spins[SIZE + 1], int pos, double &energy, std::mt19937 &gen, double tstar,
                            std::uniform_int_distribution<int> &brandom, std::uniform_int_distribution<int> ran_pos, std::uniform_real_distribution<double> &ran_u, double pa);

    void update_bins(double *private_m, double *private_bins, double *private_corr, double *old);

    void count_corr(double *private_m, double *private_bins, double *private_corr);
}

#endif // MC_STEP_H