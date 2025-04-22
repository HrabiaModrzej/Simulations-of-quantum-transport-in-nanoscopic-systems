#include <iostream>
#include <cmath>
#include <fstream>
#include "matrix.h"
#include <iomanip>
#include <eigen3/Eigen/Eigenvalues>

double gaussian(double x, double y, double x_c, double y_c, double alpha_x, double alpha_y)
{
    double x_part = pow(alpha_x * M_PI, -0.25) * exp(-0.5 * pow(x - x_c, 2) / alpha_x);
    double y_part = pow(alpha_y * M_PI, -0.25) * exp(-0.5 * pow(y - y_c, 2) / alpha_y);
    return x_part * y_part;
}

void set_S_matrix(matrix &S, double alpha_x, double alpha_y, double Delta, int N)
{
    for (int k = 0; k < S.sizeX; k++)
    {
        for (int l = k; l < S.sizeY; l++)
        {
            double x_k = int(k / N) * Delta;
            double y_k = k % N * Delta;

            double x_l = int(l / N) * Delta;
            double y_l = l % N * Delta;

            S[k][l] = S[l][k] = exp(-0.25 * pow(x_k - x_l, 2.) / alpha_x - 0.25 * pow(y_k - y_l, 2.) / alpha_y);
        }
    }
}

void set_H_matrix(matrix &H, matrix &S, double m_star, double alpha_x, double alpha_y, double Delta, int N, double omega_x, double omega_y)
{
    // printf("wx %10f   wy %10f", omega_x, omega_y);
    for (int k = 0; k < H.sizeX; k++)
    {
        for (int l = k; l < H.sizeY; l++)
        {
            double x_k = int(k / N) * Delta;
            double y_k = k % N * Delta;

            double x_l = int(l / N) * Delta;
            double y_l = l % N * Delta;

            double K_kl = -0.5 / m_star * ((pow(x_k - x_l, 2.) - 2 * alpha_x) / pow(2 * alpha_x, 2) + (pow(y_k - y_l, 2.) - 2 * alpha_y) / pow(2 * alpha_y, 2)) * S[k][l];
            double V_kl = 0.5 * m_star * (pow(omega_x, 2.) * (pow(x_k + x_l, 2.) + 2. * alpha_x) * 0.25 + pow(omega_y, 2.) * (pow(y_k + y_l, 2.) + 2. * alpha_y) * 0.25) * S[k][l];

            H[k][l] = H[l][k] = K_kl + V_kl;
        }
    }
}

int main()
{
    double hartree_to_meV = 27211.6;
    double bohr_radius = 0.0529;

    int N = 9;   // basis set grid
    int n = 101; // plotting map grid
    double a = 4 / bohr_radius;
    double Delta = 2 * a / (N - 1);
    double omega_x = 80. / hartree_to_meV;
    double omega_y = 200. / hartree_to_meV;
    double m_star = 0.24;

    double alpha_x = 1. / (m_star * omega_x);
    double alpha_y = 1. / (m_star * omega_y);

    matrix grid(N, N);
    matrix phi(n, n);

    double Dx = 2 * a / (n - 1);

    grid.set_all_elements(0.0);
    phi.set_all_elements(0.0);

    // plot gaussian base functions
    for (int k : {0, 8, 9})
    {
        for (int i = 0; i < phi.sizeX; i++)
        {
            for (int j = 0; j < phi.sizeY; j++)
            {
                double i_k = int(k / N);
                double j_k = k % N;
                phi[i][j] = gaussian(i * Dx, j * Dx, i_k * Delta, j_k * Delta, alpha_x, alpha_y);
            }
        }
        phi.save_to_file("data/gaussian_" + std::to_string(k) + ".dat");
    }

    // set matrices for generalized eigenproblem
    matrix S(N * N, N * N);
    set_S_matrix(S, alpha_x, alpha_y, Delta, N);

    matrix H(N * N, N * N);
    set_H_matrix(H, S, m_star, alpha_x, alpha_y, Delta, N, omega_x, omega_y);

    // transform into Eigen::Matrix
    Eigen::MatrixXd H_matrix(H.sizeX, H.sizeY);
    Eigen::MatrixXd S_matrix(S.sizeX, S.sizeY);
    for (int i = 0; i < H.sizeX; i++)
    {
        for (int j = 0; j < H.sizeY; j++)
        {
            H_matrix(i, j) = H[i][j];
            S_matrix(i, j) = S[i][j];
        }
    }

    // class eigensolver
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> ES(H_matrix, S_matrix);

    // Eigen State output
    for (int which_state : {0, 1, 2, 3, 4, 5, 6})
    {
        std::cout << "Eigen value number: " << which_state << "\n"
                  << ES.eigenvalues()[which_state] * hartree_to_meV << std::endl;

        // Wavefunction matrix for plotting
        matrix eigen_state(n, n);
        eigen_state.set_all_elements(0.0);
        for (int i = 0; i < eigen_state.sizeX; i++)
        {
            for (int j = 0; j < eigen_state.sizeY; j++)
            {
                double x = -a + i * Dx;
                double y = -a + j * Dx;
                for (int gauss_index = 0; gauss_index < H.sizeX; gauss_index++)
                {
                    double i_k = int(gauss_index / N);
                    double j_k = gauss_index % N;
                    eigen_state[i][j] += ES.eigenvectors().col(which_state)[gauss_index] * gaussian(x, y, i_k * Delta, j_k * Delta, alpha_x, alpha_y);
                }
                eigen_state[i][j] = pow(eigen_state[i][j],2.);
            }
        }

        eigen_state.save_to_file("data/state_" + std::to_string(which_state) + ".dat");
    }

    // task 4
    std::ofstream outfile("E_omega.dat");
    for (omega_x = 20.0 / hartree_to_meV; omega_x <= 500 / hartree_to_meV; omega_x += 10 / hartree_to_meV)
    {
        alpha_x = 1. / (m_star * omega_x);
        set_S_matrix(S, alpha_x, alpha_y, Delta, N);
        set_H_matrix(H, S, m_star, alpha_x, alpha_y, Delta, N, omega_x, omega_y);
        for (int i = 0; i < H.sizeX; i++)
        {
            for (int j = 0; j < H.sizeY; j++)
            {
                H_matrix(i, j) = H[i][j];
                S_matrix(i, j) = S[i][j];
            }
        }
        outfile << omega_x*hartree_to_meV << " ";

        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> NewES(H_matrix, S_matrix);
        for (int state = 0; state < 10; state++)
        {
            outfile << NewES.eigenvalues()[state]*hartree_to_meV << " ";
        }
        outfile << "\n";
    }
    outfile.close();

    // Task 5

    omega_x = 80/hartree_to_meV;
    omega_y = 360/hartree_to_meV;
    alpha_x = 1. / (m_star * omega_x);
    alpha_y = 1. / (m_star * omega_y);
    set_S_matrix(S, alpha_x, alpha_y, Delta, N);
    set_H_matrix(H, S, m_star, alpha_x, alpha_y, Delta, N, omega_x, omega_y);
    for (int i = 0; i < H.sizeX; i++)
    {
        for (int j = 0; j < H.sizeY; j++)
        {
            H_matrix(i, j) = H[i][j];
            S_matrix(i, j) = S[i][j];
        }
    }
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> NewES(H_matrix, S_matrix);

    for (int which_state : {0, 1, 2, 3, 4, 5})
    {
        std::cout << "Eigen value number: " << which_state << "\n"
                  << NewES.eigenvalues()[which_state] * hartree_to_meV << std::endl;

        // Wavefunction matrix for plotting
        matrix eigen_state(n, n);
        eigen_state.set_all_elements(0.0);
        for (int i = 0; i < eigen_state.sizeX; i++)
        {
            for (int j = 0; j < eigen_state.sizeY; j++)
            {
                double x = -a + i * Dx;
                double y = -a + j * Dx;
                for (int gauss_index = 0; gauss_index < H.sizeX; gauss_index++)
                {
                    double i_k = int(gauss_index / N);
                    double j_k = gauss_index % N;
                    eigen_state[i][j] += NewES.eigenvectors().col(which_state)[gauss_index] * gaussian(x, y, i_k * Delta, j_k * Delta, alpha_x, alpha_y);
                }
                eigen_state[i][j] = pow(eigen_state[i][j],2.);
            }
        }

        eigen_state.save_to_file("data/y_state_" + std::to_string(which_state) + ".dat");
    }

    return 0;
}
