#include "matrix.hpp"

namespace matrix
{
    void Single_Value(const int spins[SIZE + 1], double sv[3])
    {
        // 创建一个矩阵
        Eigen::MatrixXd mat(Parameters::LATTICE_SIZE, Parameters::LATTICE_SIZE);

        // 填充矩阵
        for (int i = 0; i < Parameters::LATTICE_SIZE; ++i)
        {
            for (int j = 0; j < Parameters::LATTICE_SIZE; ++j)
            {
                mat(i, j) = 0;
                for (int k = 0; k < Parameters::LATTICE_SIZE; ++k)
                {
                    if (spins[i * Parameters::LATTICE_SIZE + k] == spins[j * Parameters::LATTICE_SIZE + k])
                    {
                        mat(i, j) += 1;
                    }
                    else
                    {
                        mat(i, j) -= (1.0 / (1.0 * (Parameters::Q - 1)));
                    }
                }
            }
        }

        // // 打印复数相位矩阵
        // std::cout << "复数相位矩阵 (matA):" << std::endl;
        // std::cout << matA << std::endl
        //           << std::endl;

        // // 创建一个复数矩阵
        // Eigen::MatrixXcd mat(Parameters::LATTICE_SIZE, Parameters::LATTICE_SIZE);

        // // 填充矩阵
        // for (int i = 0; i < Parameters::LATTICE_SIZE; ++i)
        // {
        //     for (int j = 0; j < Parameters::LATTICE_SIZE; ++j)
        //     {
        //         double theta = (2 * M_PI * spins[i * Parameters::LATTICE_SIZE + j]) / Parameters::Q;
        //         mat(i, j) = std::exp(std::complex<double>(0, theta)); // e^(i * theta)
        //     }
        // }

        // Eigen::MatrixXcd mat_self = mat * mat.adjoint(); // 矩阵乘以它的复共轭
        // Eigen::MatrixXd mat_self_real = mat_self.real(); // 只取实部

        // // std::cout << "自旋一致性矩阵 (real part of mat_self):" << std::endl;
        // // std::cout << mat_self_real << std::endl
        // //           << std::endl;

        // 计算特征值，开平方之后就是奇异值
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen_solver(mat); // 如果矩阵自伴随

        Eigen::VectorXd eigenvalues = eigen_solver.eigenvalues().real();         // 只取实部
        Eigen::VectorXd singularValues = eigenvalues.cwiseMax(0).array().sqrt(); // 开方，负值处理

        // 计算最大奇异值
        double norm_2 = singularValues.maxCoeff();

        // Step 2: 归一化奇异值
        double sum_singular_values = singularValues.sum();                    // 所有奇异值之和
        Eigen::VectorXd probabilities = singularValues / sum_singular_values; // 归一化

        // Step 3: 计算香农熵
        double shannon_entropy = 0.0;
        for (int i = 0; i < probabilities.size(); ++i)
        {
            if (probabilities(i) > 0)
            { // 忽略p_i为0的情况，避免log(0)
                shannon_entropy -= probabilities(i) * std::log(probabilities(i));
            }
        }

        // 存储结果
        sv[0] = norm_2 / Parameters::LATTICE_SIZE;
        sv[1] = sum_singular_values / Parameters::LATTICE_SIZE;
        sv[2] = shannon_entropy;
    }

}