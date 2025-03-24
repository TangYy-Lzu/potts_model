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

    void w_error(std::ofstream &outputBining, double tstar, double *m, double *bins)
    {
        int nbin = Parameters::BROKEN_SIZE * Parameters::SAMPLE_SIZE;
        outputBining << tstar << " " << 1.0 / tstar << " ";

        // 方差，bins只有rr的数据可以用
        double mag_plane_error = (bins[static_cast<int>(Quantity::MAG_PLANERR)] - m[static_cast<int>(Quantity::MAG_PLANE)] * m[static_cast<int>(Quantity::MAG_PLANE)]) / (nbin - 1);
        double mag_plane_chi_error = (bins[static_cast<int>(Quantity::MAG_PLANE2RR)] - bins[static_cast<int>(Quantity::MAG_PLANE2_bin_sus)] * bins[static_cast<int>(Quantity::MAG_PLANE2_bin_sus)]) / (nbin - 1);
        // 标准误差
        // mag_plane_error = sqrt(mag_plane_error);
        // mag_plane_chi_error = sqrt(mag_plane_chi_error);
        outputBining << mag_plane_error << " " << mag_plane_chi_error << " ";

        double mag_error = (bins[static_cast<int>(Quantity::MAGRR)] - m[static_cast<int>(Quantity::MAG)] * m[static_cast<int>(Quantity::MAG)]) / (nbin - 1);
        double mag_chi_error = (bins[static_cast<int>(Quantity::MAG2RR)] - bins[static_cast<int>(Quantity::MAG2_bin_sus)] * bins[static_cast<int>(Quantity::MAG2_bin_sus)]) / (nbin - 1);
        // mag_error = sqrt(mag_error);
        // mag_chi_error = sqrt(mag_chi_error);
        outputBining << mag_error << " " << mag_chi_error << " ";

        double energy_error = (bins[static_cast<int>(Quantity::ENERR)] - m[static_cast<int>(Quantity::ENE)] * m[static_cast<int>(Quantity::ENE)]) / (nbin - 1);
        double heat_error = (bins[static_cast<int>(Quantity::ENE2RR)] - bins[static_cast<int>(Quantity::ENE2_bin_sus)] * bins[static_cast<int>(Quantity::ENE2_bin_sus)]) / (nbin - 1);
        // energy_error = sqrt(energy_error);
        // heat_error = sqrt(heat_error);
        outputBining << energy_error << " " << heat_error << " ";

        double order_error = (bins[static_cast<int>(Quantity::ORDERRR)] - m[static_cast<int>(Quantity::ORDER)] * m[static_cast<int>(Quantity::ORDER)]) / (nbin - 1);
        double order_chi_error = (bins[static_cast<int>(Quantity::ORDER2RR)] - bins[static_cast<int>(Quantity::ORDER2_bin_sus)] * bins[static_cast<int>(Quantity::ORDER2_bin_sus)]) / (nbin - 1);
        // order_error = sqrt(order_error);
        // order_chi_error = sqrt(order_chi_error);
        outputBining << order_error << " " << order_chi_error << " ";

        double dis_error = (bins[static_cast<int>(Quantity::DISRR)] - m[static_cast<int>(Quantity::DIS)] * m[static_cast<int>(Quantity::DIS)]) / (nbin - 1);
        double dis_chi_error = (bins[static_cast<int>(Quantity::DIS2RR)] - bins[static_cast<int>(Quantity::DIS2_bin_sus)] * bins[static_cast<int>(Quantity::DIS2_bin_sus)]) / (nbin - 1);
        // dis_error = sqrt(dis_error);
        // dis_chi_error = sqrt(dis_chi_error);
        outputBining << dis_error << " " << dis_chi_error << " ";

        double area_error = (bins[static_cast<int>(Quantity::AREARR)] - m[static_cast<int>(Quantity::AREA)] * m[static_cast<int>(Quantity::AREA)]) / (nbin - 1);
        double area_chi_error = (bins[static_cast<int>(Quantity::AREA2RR)] - bins[static_cast<int>(Quantity::AREA2_bin_sus)] * bins[static_cast<int>(Quantity::AREA2_bin_sus)]) / (nbin - 1);
        // area_error = sqrt(area_error);
        // area_chi_error = sqrt(area_chi_error);
        outputBining << area_error << " " << area_chi_error << " ";

        double euler_error = (bins[static_cast<int>(Quantity::EULERRR)] - m[static_cast<int>(Quantity::EULER)] * m[static_cast<int>(Quantity::EULER)]) / (nbin - 1);
        double euler_chi_error = (bins[static_cast<int>(Quantity::EULER2RR)] - bins[static_cast<int>(Quantity::EULER2_bin_sus)] * bins[static_cast<int>(Quantity::EULER2_bin_sus)]) / (nbin - 1);
        // euler_error = sqrt(euler_error);
        // euler_chi_error = sqrt(euler_chi_error);
        outputBining << euler_error << " " << euler_chi_error << " ";

        double single_value_error = (bins[static_cast<int>(Quantity::SINGLE_VALUERR)] - m[static_cast<int>(Quantity::SINGLE_VALUE)] * m[static_cast<int>(Quantity::SINGLE_VALUE)]) / (nbin - 1);
        double single_value_chi_error = (bins[static_cast<int>(Quantity::SINGLE_VALUE2RR)] - bins[static_cast<int>(Quantity::SINGLE_VALUE2_bin_sus)] * bins[static_cast<int>(Quantity::SINGLE_VALUE2_bin_sus)]) / (nbin - 1);
        // single_value_erro = sqrt(single_value_error);
        // single_value_chi_erro = sqrt(single_value_chi_error);
        outputBining << single_value_error << " " << single_value_chi_error << " ";

        double nuclear_norm_error = (bins[static_cast<int>(Quantity::NUCLEAR_NORMRR)] - m[static_cast<int>(Quantity::NUCLEAR_NORM)] * m[static_cast<int>(Quantity::NUCLEAR_NORM)]) / (nbin - 1);
        double nuclear_norm_chi_error = (bins[static_cast<int>(Quantity::NUCLEAR_NORM2RR)] - bins[static_cast<int>(Quantity::NUCLEAR_NORM2_bin_sus)] * bins[static_cast<int>(Quantity::NUCLEAR_NORM2_bin_sus)]) / (nbin - 1);
        // nuclear_norm_error = sqrt(nuclear_norm_error);
        // nuclear_norm_chi_error = sqrt(nuclear_norm_chi_error);
        outputBining << nuclear_norm_error << " " << nuclear_norm_chi_error << " ";

        double shannon_error = (bins[static_cast<int>(Quantity::SHANNONRR)] - m[static_cast<int>(Quantity::SHANNON)] * m[static_cast<int>(Quantity::SHANNON)]) / (nbin - 1);
        double shannon_chi_error = (bins[static_cast<int>(Quantity::SHANNON2RR)] - bins[static_cast<int>(Quantity::SHANNON2_bin_sus)] * bins[static_cast<int>(Quantity::SHANNON2_bin_sus)]) / (nbin - 1);
        // shannon_error = sqrt(shannon_error);
        // shannon_chi_error = sqrt(shannon_chi_error);
        outputBining << shannon_error << " " << shannon_chi_error << " ";

        double binder_error = (bins[static_cast<int>(Quantity::MAG4RR)] - m[static_cast<int>(Quantity::MAG4)] * m[static_cast<int>(Quantity::MAG4)]) / (nbin - 1);
        // binder_error = sqrt(binder_erro);
        outputBining << binder_error << std::endl;
    }

    void w_corr(std::ofstream &outputcorrelation, double tstar, double *corr)
    {
        outputcorrelation << tstar << " " << 1.0 / tstar << " ";

        // 关联
        outputcorrelation << corr[static_cast<int>(Quantity::MAG_PLANE)] << " " << corr[static_cast<int>(Quantity::MAG_PLANE2)] << " ";
        outputcorrelation << corr[static_cast<int>(Quantity::MAG)] << " " << corr[static_cast<int>(Quantity::MAG2)] << " ";
        outputcorrelation << corr[static_cast<int>(Quantity::ENE)] << " " << corr[static_cast<int>(Quantity::ENE2)] << " ";
        outputcorrelation << corr[static_cast<int>(Quantity::ORDER)] << " " << corr[static_cast<int>(Quantity::ORDER2)] << " ";
        outputcorrelation << corr[static_cast<int>(Quantity::DIS)] << " " << corr[static_cast<int>(Quantity::DIS2)] << " ";
        outputcorrelation << corr[static_cast<int>(Quantity::AREA)] << " " << corr[static_cast<int>(Quantity::AREA2)] << " ";
        outputcorrelation << corr[static_cast<int>(Quantity::EULER)] << " " << corr[static_cast<int>(Quantity::EULER2)] << " ";
        outputcorrelation << corr[static_cast<int>(Quantity::SINGLE_VALUE)] << " " << corr[static_cast<int>(Quantity::SINGLE_VALUE2)] << " ";
        outputcorrelation << corr[static_cast<int>(Quantity::NUCLEAR_NORM)] << " " << corr[static_cast<int>(Quantity::NUCLEAR_NORM2)] << " ";
        outputcorrelation << corr[static_cast<int>(Quantity::SHANNON)] << " " << corr[static_cast<int>(Quantity::SHANNON2)] << " ";
        outputcorrelation << corr[static_cast<int>(Quantity::MAG4)] << std::endl;
    }
}