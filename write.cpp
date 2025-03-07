#include "write.hpp"

namespace FileUtils
{
    int dimentions = SIZE + 1;
    void write(std::ofstream &output, double tstar, int spins[SIZE + 1])
    {
        output << tstar << " ";
        for (int i = 0; i < SIZE; i++)
        {
            if (i % Parameters::LATTICE_SIZE == 0)
                output << std::endl;
            output << spins[i] << " ";
        }
        output << std::endl
               << std::endl;
    }

    void w_output(std::ofstream &file, double tstar, double *m)
    {
        file << tstar << " " << 1.0 / tstar << " ";

        double mag_plane_chi = (m[static_cast<int>(Quantity::MAG_PLANE2)] - m[static_cast<int>(Quantity::MAG_PLANE)] * m[static_cast<int>(Quantity::MAG_PLANE)]) * SIZE;
        file << m[static_cast<int>(Quantity::MAG_PLANE)] << " " << mag_plane_chi << " ";

        double mag_chi = (m[static_cast<int>(Quantity::MAG2)] - m[static_cast<int>(Quantity::MAG)] * m[static_cast<int>(Quantity::MAG)]) * SIZE;
        file << m[static_cast<int>(Quantity::MAG)] << " " << mag_chi << " ";

        double heat = (m[static_cast<int>(Quantity::ENE2)] - m[static_cast<int>(Quantity::ENE)] * m[static_cast<int>(Quantity::ENE)]) * SIZE;
        file << m[static_cast<int>(Quantity::ENE)] << " " << heat << " ";

        double order_chi = (m[static_cast<int>(Quantity::ORDER2)] - m[static_cast<int>(Quantity::ORDER)] * m[static_cast<int>(Quantity::ORDER)]) * SIZE;
        file << m[static_cast<int>(Quantity::ORDER)] << " " << order_chi << " ";

        double dis_chi = (m[static_cast<int>(Quantity::DIS2)] - m[static_cast<int>(Quantity::DIS)] * m[static_cast<int>(Quantity::DIS)]) * SIZE;
        file << m[static_cast<int>(Quantity::DIS)] << " " << dis_chi << " ";

        double area_chi = (m[static_cast<int>(Quantity::AREA2)] - m[static_cast<int>(Quantity::AREA)] * m[static_cast<int>(Quantity::AREA)]) * SIZE;
        file << m[static_cast<int>(Quantity::AREA)] << " " << area_chi << " ";

        double euler_chi = (m[static_cast<int>(Quantity::EULER2)] - m[static_cast<int>(Quantity::EULER)] * m[static_cast<int>(Quantity::EULER)]) * SIZE;
        file << m[static_cast<int>(Quantity::EULER)] << " " << euler_chi << " ";

        double single_value_chi = (m[static_cast<int>(Quantity::SINGLE_VALUE2)] - m[static_cast<int>(Quantity::SINGLE_VALUE)] * m[static_cast<int>(Quantity::SINGLE_VALUE)]) * SIZE;
        file << m[static_cast<int>(Quantity::SINGLE_VALUE)] << " " << single_value_chi << " ";

        double nuclear_chi = (m[static_cast<int>(Quantity::NUCLEAR_NORM2)] - m[static_cast<int>(Quantity::NUCLEAR_NORM)] * m[static_cast<int>(Quantity::NUCLEAR_NORM)]) * SIZE;
        file << m[static_cast<int>(Quantity::NUCLEAR_NORM)] << " " << nuclear_chi << " ";

        double shannon_chi = (m[static_cast<int>(Quantity::SHANNON2)] - m[static_cast<int>(Quantity::SHANNON)] * m[static_cast<int>(Quantity::SHANNON)]);
        file << m[static_cast<int>(Quantity::SHANNON)] << " " << shannon_chi << " ";

        double binder = 1.0 - m[static_cast<int>(Quantity::MAG4)] / (3.0 * m[static_cast<int>(Quantity::MAG2)] * m[static_cast<int>(Quantity::MAG2)]);
        file << binder << std::endl;
    }

    void w_error(std::ofstream &file, double tstar, double *bins, double *m)
    {
        int nbin = Parameters::BRONKEN_SIZE * Parameters::SAMPLE_SIZE;
        file << tstar << " " << 1.0 / tstar << " ";

        // 方差，bins只有rr的数据可以用
        double mag_plane_error = (bins[static_cast<int>(Quantity::MAG_PLANERR)] - m[static_cast<int>(Quantity::MAG_PLANE)] * m[static_cast<int>(Quantity::MAG_PLANE)]) / (nbin - 1);
        double mag_plane_chi_error = (bins[static_cast<int>(Quantity::MAG_PLANE2RR)] - bins[static_cast<int>(Quantity::MAG_PLANE2_bin)] * bins[static_cast<int>(Quantity::MAG_PLANE2_bin)]) / (nbin - 1);
        // 标准误差
        // mag_plane_error = sqrt(mag_plane_error);
        // mag_plane_chi_error = sqrt(mag_plane_chi_error);
        file << mag_plane_error << " " << mag_plane_chi_error << " ";

        double mag_error = (bins[static_cast<int>(Quantity::MAGRR)] - m[static_cast<int>(Quantity::MAG)] * m[static_cast<int>(Quantity::MAG)]) / (nbin - 1);
        double mag_chi_error = (bins[static_cast<int>(Quantity::MAG2RR)] - bins[static_cast<int>(Quantity::MAG2_bin)] * bins[static_cast<int>(Quantity::MAG2_bin)]) / (nbin - 1);
        // mag_error = sqrt(mag_error);
        // mag_chi_error = sqrt(mag_chi_error);
        file << mag_error << " " << mag_chi_error << " ";

        double energy_erro = (bins[static_cast<int>(Quantity::ENERR)] - m[static_cast<int>(Quantity::ENE)] * m[static_cast<int>(Quantity::ENE)]) / (nbin - 1);
        double heat_erro = (bins[static_cast<int>(Quantity::ENE2RR)] - bins[static_cast<int>(Quantity::ENE2_bin)] * bins[static_cast<int>(Quantity::ENE2_bin)]) / (nbin - 1);
        // energy_erro = sqrt(energy_erro);
        // heat_erro = sqrt(heat_erro);
        file << energy_erro << " " << heat_erro << " ";

        double order_erro = (bins[static_cast<int>(Quantity::ORDERRR)] - m[static_cast<int>(Quantity::ORDER)] * m[static_cast<int>(Quantity::ORDER)]) / (nbin - 1);
        double order_chi_erro = (bins[static_cast<int>(Quantity::ORDER2RR)] - bins[static_cast<int>(Quantity::ORDER2_bin)] * bins[static_cast<int>(Quantity::ORDER2_bin)]) / (nbin - 1);
        // order_erro = sqrt(order_erro);
        // order_chi_erro = sqrt(order_chi_erro);
        file << order_erro << " " << order_chi_erro << " ";

        double dis_erro = (bins[static_cast<int>(Quantity::DISRR)] - m[static_cast<int>(Quantity::DIS)] * m[static_cast<int>(Quantity::DIS)]) / (nbin - 1);
        double dis_chi_erro = (bins[static_cast<int>(Quantity::DIS2RR)] - bins[static_cast<int>(Quantity::DIS2_bin)] * bins[static_cast<int>(Quantity::DIS2_bin)]) / (nbin - 1);
        // dis_erro = sqrt(dis_erro);
        // dis_chi_erro = sqrt(dis_chi_erro);
        file << dis_erro << " " << dis_chi_erro << " ";

        double area_erro = (bins[static_cast<int>(Quantity::AREARR)] - m[static_cast<int>(Quantity::AREA)] * m[static_cast<int>(Quantity::AREA)]) / (nbin - 1);
        double area_chi_erro = (bins[static_cast<int>(Quantity::AREA2RR)] - bins[static_cast<int>(Quantity::AREA2_bin)] * bins[static_cast<int>(Quantity::AREA2_bin)]) / (nbin - 1);
        // area_erro = sqrt(area_erro);
        // area_chi_erro = sqrt(area_chi_erro);
        file << area_erro << " " << area_chi_erro << " ";

        double euler_erro = (bins[static_cast<int>(Quantity::EULERRR)] - m[static_cast<int>(Quantity::EULER)] * m[static_cast<int>(Quantity::EULER)]) / (nbin - 1);
        double euler_chi_erro = (bins[static_cast<int>(Quantity::EULER2RR)] - bins[static_cast<int>(Quantity::EULER2_bin)] * bins[static_cast<int>(Quantity::EULER2_bin)]) / (nbin - 1);
        // euler_erro = sqrt(euler_erro);
        // euler_chi_erro = sqrt(euler_chi_erro);
        file << euler_erro << " " << euler_chi_erro << " ";

        double single_value_erro = (bins[static_cast<int>(Quantity::SINGLE_VALUERR)] - m[static_cast<int>(Quantity::SINGLE_VALUE)] * m[static_cast<int>(Quantity::SINGLE_VALUE)]) / (nbin - 1);
        double single_value_chi_erro = (bins[static_cast<int>(Quantity::SINGLE_VALUE2RR)] - bins[static_cast<int>(Quantity::SINGLE_VALUE2_bin)] * bins[static_cast<int>(Quantity::SINGLE_VALUE2_bin)]) / (nbin - 1);
        // single_value_erro = sqrt(single_value_erro);
        // single_value_chi_erro = sqrt(single_value_chi_erro);
        file << single_value_erro << " " << single_value_chi_erro << " ";

        double nuclear_norm_erro = (bins[static_cast<int>(Quantity::NUCLEAR_NORMRR)] - m[static_cast<int>(Quantity::NUCLEAR_NORM)] * m[static_cast<int>(Quantity::NUCLEAR_NORM)]) / (nbin - 1);
        double nuclear_norm_chi_erro = (bins[static_cast<int>(Quantity::NUCLEAR_NORM2RR)] - bins[static_cast<int>(Quantity::NUCLEAR_NORM2_bin)] * bins[static_cast<int>(Quantity::NUCLEAR_NORM2_bin)]) / (nbin - 1);
        // nuclear_norm_erro = sqrt(nuclear_norm_erro);
        // nuclear_norm_chi_erro = sqrt(nuclear_norm_chi_erro);
        file << nuclear_norm_erro << " " << nuclear_norm_chi_erro << " ";

        double shannon_erro = (bins[static_cast<int>(Quantity::SHANNONRR)] - m[static_cast<int>(Quantity::SHANNON)] * m[static_cast<int>(Quantity::SHANNON)]) / (nbin - 1);
        double shannon_chi_erro = (bins[static_cast<int>(Quantity::SHANNON2RR)] - bins[static_cast<int>(Quantity::SHANNON2_bin)] * bins[static_cast<int>(Quantity::SHANNON2_bin)]) / (nbin - 1);
        // shannon_erro = sqrt(shannon_erro);
        // shannon_chi_erro = sqrt(shannon_chi_erro);
        file << shannon_erro << " " << shannon_chi_erro << " ";

        double binder_erro = (bins[static_cast<int>(Quantity::MAG4RR)] - m[static_cast<int>(Quantity::MAG4)] * m[static_cast<int>(Quantity::MAG4)]) / (nbin - 1);
        binder_erro = sqrt(binder_erro);
        file << binder_erro << std::endl;
    }
}