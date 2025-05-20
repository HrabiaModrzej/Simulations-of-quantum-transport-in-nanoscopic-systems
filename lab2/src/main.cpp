#include <iostream>
#include <cmath>
#include <fstream>
#include "matrix.h"
#include <complex>
#include <vector>
#include <eigen3/Eigen/Eigenvalues>

void progress(double percent)
{
    int length = 50;
    int filled = int(percent * length) + 1;
    int empty = length - filled;
    std::string out = "[";
    for (int i = 0; i < filled; i++)
    {
        out += "█";
    }
    for (int i = 0; i < empty; i++)
    {
        out += " ";
    }
    out += "]  " + std::to_string(percent * 100) + "%";
    std::cout << "\r" + out << std::flush;
}

class params
{
public:
    params();
    ~params();
    int N = 1e3;
    double m_star_1 = 0.063;
    double m_star_2 = 0.063;
    double hartree_to_meV = 27211.6;
    double bohr_radius = 0.0529;
    double L = 15.0 / bohr_radius;
    double d_AlGaAs = 5.0 / bohr_radius;
    double d_GaAs = 3.0 / bohr_radius;

    double m_star_2_x(double x)
    {
        if (x < 0.0 || x > 1.0)
        {
            std::cerr << "The ration in Al_xGa_1-xAs isn't in range [0,1]\n";
        }
        m_star_2 = 0.063 + 0.083 * x;
        return m_star_2;
    }
};
params::params() {}
params::~params() {}

matrix<complex<double>> multiply_matrices(matrix<complex<double>> &A, matrix<complex<double>> &B)
{
    matrix<complex<double>> out(B.sizeX, A.sizeY);
    out.set_all_elements(-10000.0);

    if (A.sizeX != B.sizeY)
    {
        std::cerr << "Multiplyed matrices: wrong size\n";
        return out;
    }

    for (int i = 0; i < out.sizeX; i++)
    {
        for (int j = 0; j < out.sizeY; j++)
        {
            out[i][j] = 0.0;
            for (int k = 0; k < A.sizeX; k++)
            {
                out[i][j] += A[k][j] * B[i][k];
            }
        }
    }
    return out;
}

complex<double> k(double m_star_n, double E, double U_n)
{
    if (E - U_n < 0.0)
    {
        return 1i * sqrt(2.0 * m_star_n * (U_n - E));
    }
    else
    {
        return sqrt(2.0 * m_star_n * (E - U_n));
    }
}

matrix<complex<double>> M_n(double m_star_n, double m_star_np1, double U_n, double U_np1, double z_n, double E)
{
    matrix<complex<double>> out(2, 2);
    complex<double> kn = k(m_star_n, E, U_n);
    complex<double> knp1 = k(m_star_np1, E, U_np1);
    complex<double> frac = (knp1 * m_star_n) / (kn * m_star_np1);
    out[0][0] = 0.5 * (1.0 + frac) * exp(+1i * (knp1 - kn) * z_n);
    out[1][0] = 0.5 * (1.0 - frac) * exp(-1i * (knp1 + kn) * z_n);
    out[0][1] = 0.5 * (1.0 - frac) * exp(+1i * (knp1 + kn) * z_n);
    out[1][1] = 0.5 * (1.0 + frac) * exp(-1i * (knp1 - kn) * z_n);

    return out;
}

void set_U_m_single_barrier(params p, vector<double> &U, vector<double> &m)
{
    double len_of_nothing = (p.L - p.d_AlGaAs) / 2;
    int st_interface = len_of_nothing / (p.L / p.N);
    int nd_interface = st_interface + p.d_AlGaAs / (p.L / p.N);

    for (int i = 0; i < st_interface; i++)
    {
        U[i] = 0;
        m[i] = p.m_star_1;
    }
    for (int i = st_interface; i < nd_interface; i++)
    {
        U[i] = 0.27 * 1000 / p.hartree_to_meV;
        m[i] = p.m_star_2_x(0.3);
    }
    for (int i = nd_interface; i < p.N; i++)
    {
        U[i] = 0;
        m[i] = p.m_star_1;
    }
}

void set_U_m_double_barrier(params p, vector<double> &U, vector<double> &m)
{
    int st_interface = int((p.L - 2 * p.d_AlGaAs - p.d_GaAs) * p.N / (2 * p.L));
    int nd_interface = st_interface + int(p.d_AlGaAs * p.N / p.L);
    int rd_interface = nd_interface + int(p.d_GaAs * p.N / p.L);
    int th_interface = rd_interface + int(p.d_AlGaAs * p.N / p.L);

    for (int i = 0; i < st_interface; i++)
    {
        U[i] = 0;
        m[i] = p.m_star_1;
    }
    for (int i = st_interface; i < nd_interface; i++)
    {
        U[i] = 0.27 * 1000 / p.hartree_to_meV;
        m[i] = p.m_star_2_x(0.3);
    }
    for (int i = nd_interface; i < rd_interface; i++)
    {
        U[i] = 0;
        m[i] = p.m_star_1;
    }
    for (int i = rd_interface; i < th_interface; i++)
    {
        U[i] = 0.27 * 1000 / p.hartree_to_meV;
        m[i] = p.m_star_2_x(0.3);
    }
    for (int i = th_interface; i < p.N; i++)
    {
        U[i] = 0;
        m[i] = p.m_star_1;
    }
}

void set_U_double_barrier(params p, vector<double> &U)
{
    int st_interface = int((p.L - 2 * p.d_AlGaAs - p.d_GaAs) * p.N / (2 * p.L));
    int nd_interface = st_interface + int(p.d_AlGaAs * p.N / p.L);
    int rd_interface = nd_interface + int(p.d_GaAs * p.N / p.L);
    int th_interface = rd_interface + int(p.d_AlGaAs * p.N / p.L);

    for (int i = 0; i < st_interface; i++)
    {
        U[i] = 0;
    }
    for (int i = st_interface; i < nd_interface; i++)
    {
        U[i] = 0.27 * 1000 / p.hartree_to_meV;
    }
    for (int i = nd_interface; i < rd_interface; i++)
    {
        U[i] = 0;
    }
    for (int i = rd_interface; i < th_interface; i++)
    {
        U[i] = 0.27 * 1000 / p.hartree_to_meV;
    }
    for (int i = th_interface; i < p.N; i++)
    {
        U[i] = 0;
    }
}

void set_z_vec(vector<double> &z, double z_min, double z_max, int Steps)
{
    for (int i = 0; i < Steps; i++)
    {
        z[i] = z_min + (z_max - z_min) * i / (Steps - 1);
    }
}

matrix<complex<double>> get_M_from1_to_N(params p, vector<double> &U, vector<double> &m, vector<double> &z, double E)
{
    if (U.size() != m.size() || m.size() != z.size() || U.size() != z.size())
    {
        std::cerr << "Incoherent vector lengths\n";
        return -1;
    }
    matrix<complex<double>> out(2, 2);
    out.identity_matrix();

    matrix<complex<double>> temp(2, 2);
    for (int i = 1; i < int(U.size() - 1); i++)
    {
        temp = M_n(m[i], m[i + 1], U[i], U[i + 1], z[i], E);
        out = multiply_matrices(out, temp);
    }

    return out;
}

double Transmitance(params p, matrix<complex<double>> M_1_N, vector<double> &m, vector<double> &U, double E)
{
    int last_index = U.size() - 1;
    complex<double> I = k(m[last_index], E, U[last_index]) * m[0] / (k(m[0], E, U[0]) * m[last_index]);
    complex<double> II = 1.0 / pow(fabs(M_1_N[0][0]), 2.);
    return (I * II).real();
}

double Reflectance(matrix<complex<double>> M_1_N)
{
    return pow(fabs(M_1_N[0][1]), 2.0) / pow(fabs(M_1_N[0][0]), 2.0);
}

void set_U_V_bias(params p, vector<double> &U, double V_bias)
{
    set_U_double_barrier(p, U);
    for (int i = 0; i < p.N; i++)
    {
        U[i] -= V_bias * i / (p.N - 1);
    }
}

double Tsu_Esaki(params p, vector<double> &U, vector<double> &m, vector<double> &z, double V_bias)
{
    // V_bias*e in atomic units → e=1 → V_bias*e[meV] = V_bias*1/hartree_to_meV [ha]

    double integral = 0.;
    double kB = 8.61733e-2 / p.hartree_to_meV; // meV/K → ha/K
    double T = 77.;
    double mu_s = 87. / p.hartree_to_meV;
    double mu_d = 87. / p.hartree_to_meV;
    double frac = p.m_star_1 * kB * T / (2 * M_PI * M_PI);

    double dE = 1 / p.hartree_to_meV;
    set_U_V_bias(p, U, V_bias);

    for (double E = 1.0 / p.hartree_to_meV; E < 2 * mu_s; E += dE)
    {
        matrix<complex<double>> M_1_N = get_M_from1_to_N(p, U, m, z, E);
        double ln_in_integral = std::log((1.0 + std::exp((mu_s - E) / (kB * T))) / (1 + std::exp((mu_d - V_bias - E) / (kB * T))));
        /*
        Integral approximated as ∫f(x)dx = Σ Dx * f(x) in formalism of simple Newton-Cotes quadrature: 0 degree
        */
        integral += dE * Transmitance(p, M_1_N, m, U, E) * ln_in_integral;
    }
    return integral * frac;
}

void write_U_m_dist(params p, vector<double> &z, vector<double> &U, vector<double> &m, std::string file_name = "spacial_dist.dat")
{
    std::ofstream outfile("./../data/"+file_name);
    outfile << "# index | z | U | m\n";
    for (int i = 0; i < p.N; i++)
    {
        outfile << i << " " << z[i] * p.bohr_radius << " " << U[i] * p.hartree_to_meV / 1000 << " " << m[i] << "\n";
    }

    outfile.close();
}

void calculate_T_R_for_system(params p, vector<double> &z, vector<double> &U, vector<double> &m, std::string file_name = "T_R_E.dat")
{
    std::ofstream outfile("./../data/"+file_name);
    for (double E = 0.0; E < 1000 / p.hartree_to_meV; E += 1 / p.hartree_to_meV)
    {
        matrix<complex<double>> M_1_N(get_M_from1_to_N(p, U, m, z, E));
        outfile << E * p.hartree_to_meV / 1000 << " " << Transmitance(p, M_1_N, m, U, E) << " " << Reflectance(M_1_N) << "\n";
    }
    outfile.close();
}

void calculate_current_voltage_char_Tsu_Esaki(params p, vector<double> &z, vector<double> &U, vector<double> &m, std::string file_name = "I_V_char.dat")
{
    std::cout << "Start of current-voltage chatacteristic calculation\n";
    std::ofstream outfile("./../data/"+file_name);
    for (double V_bias = 2.0 / p.hartree_to_meV; V_bias <= 500 / p.hartree_to_meV; V_bias += 2 / p.hartree_to_meV)
    {
        outfile << V_bias * p.hartree_to_meV/1000 << " " << Tsu_Esaki(p, U, m, z, V_bias) << "\n";
        progress(V_bias / (500 / p.hartree_to_meV));
    }
    outfile.close();
}

class QPC_params
{
public:
    QPC_params() {};
    ~QPC_params() {};
    params p;
    double epsilon = 13.6;
    double gates_dist = 3 / p.bohr_radius;
    double V_gate = 4000 / p.hartree_to_meV;
    double width = 50 / p.bohr_radius;
    double length = 100 / p.bohr_radius;
    double gate_start = 0.3 * length;
    double gate_end = 0.7 * length;
    double gate_spacing = 0.6 * width;
    int xpoints = 100;
    int ypoints = int(xpoints * width / length);
    void set_V_gate(double gate_voltage) { V_gate = gate_voltage; }
    void set_xpoints(int points) { xpoints = points; }
    void set_ypoints(int points) { ypoints = points; }
};

double f(double u, double v, QPC_params QPC)
{
    return 1. / (2 * M_PI * QPC.epsilon) * std::atan(u * v / (QPC.gates_dist * std::sqrt(pow(QPC.gates_dist, 2.0) + pow(u, 2.0) + pow(v, 2.0))));
}

double QPC_potential(double x, double y, QPC_params QPC)
{
    double t = -1 * QPC.gate_spacing / 2;
    double b = -1 * QPC.width * 2;
    double l = QPC.gate_start;
    double r = QPC.gate_end;

    double f1 = f(x - l, y - b, QPC);
    double f2 = f(x - l, t - y, QPC);
    double f3 = f(r - x, y - b, QPC);
    double f4 = f(r - x, t - y, QPC);

    double gate_1 = QPC.V_gate * (f1 + f2 + f3 + f4);

    b = QPC.gate_spacing / 2;
    t = QPC.width * 2;

    f1 = f(x - l, y - b, QPC);
    f2 = f(x - l, t - y, QPC);
    f3 = f(r - x, y - b, QPC);
    f4 = f(r - x, t - y, QPC);

    double gate_2 = QPC.V_gate * (f1 + f2 + f3 + f4);

    return gate_1 + gate_2;
}

void plot_QPC_potential(params p, QPC_params QPC)
{
    matrix<double> qpcPot(QPC.xpoints, QPC.ypoints);
    double x, y;
    for (int i = 0; i < QPC.xpoints; i++)
    {
        x = double(i) / (QPC.xpoints - 1) * QPC.length;
        for (int j = 0; j < QPC.ypoints; j++)
        {
            y = double(j) / (QPC.ypoints - 1) * QPC.width - QPC.width / 2;
            qpcPot[i][j] = QPC_potential(x, y, QPC) * p.hartree_to_meV / 1000;
            // if(qpcPot[i][j] > 0.25){printf("%6f  %6f  %10f \n", x, y, qpcPot[i][j]);}
        }
    }
    qpcPot.plot_heatmap("QPCpot.png", true, "");
}

vector<Eigen::VectorXcd> get_eigen_energies(params p, QPC_params QPC, bool show_progress = true)
{
    vector<Eigen::VectorXcd> Eigenvals_matrix;
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(QPC.ypoints, QPC.ypoints);
    Eigen::MatrixXd one = Eigen::MatrixXd::Identity(QPC.ypoints, QPC.ypoints);
    double dy = QPC.width / (QPC.ypoints - 1);
    double alpha = 1 / (2 * p.m_star_1 * dy * dy);
    if (show_progress)
        std::cout << "\nStarting the effective potential calculations:\n";
    double x, y;

    for (int i = 0; i < QPC.xpoints; i++)
    {
        H = Eigen::MatrixXd::Zero(QPC.ypoints, QPC.ypoints);
        x = double(i) / (QPC.xpoints - 1) * QPC.length;
        for (int j = 0; j < QPC.ypoints; j++)
        {
            y = double(j) / (QPC.ypoints - 1) * QPC.width - QPC.width / 2;
            if (j == 0)
            {
                H(j, j) = 2 * alpha + QPC_potential(x, y, QPC);
                H(j, j + 1) = -alpha;
            }
            else if (j == QPC.ypoints - 1)
            {
                H(j, j) = 2 * alpha + QPC_potential(x, y, QPC);
                H(j, j - 1) = -alpha;
            }
            else
            {
                H(j, j) = 2 * alpha + QPC_potential(x, y, QPC);
                H(j, j - 1) = -alpha;
                H(j, j + 1) = -alpha;
            }
            if (show_progress)
                progress(double(i * QPC.ypoints + j) / double(QPC.xpoints * QPC.ypoints));
        }
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> ES(H, one);
        Eigenvals_matrix.push_back(ES.eigenvalues());
    }
    return Eigenvals_matrix;
}

double get_conductance_Landauer(params p, QPC_params QPC, vector<Eigen::VectorXcd> &Eigenvals_matrix, double E, int n_included_states)
{
    std::vector<double> x_vec(QPC.xpoints);
    std::vector<double> U(QPC.xpoints);
    std::vector<double> m(QPC.xpoints);
    double conduct_temp = 0;

    for (int state = 0; state < n_included_states; state++)
    {
        for (int i = 0; i < QPC.xpoints; i++)
        {
            U[i] = Eigenvals_matrix[i][state].real();
            x_vec[i] = QPC.length / (QPC.xpoints - 1) * i;
            m[i] = p.m_star_1;
        }
        matrix<complex<double>> M_1_N = get_M_from1_to_N(p, U, m, x_vec, E);
        conduct_temp += Transmitance(p, M_1_N, m, U, E);
    }
    return conduct_temp;
}

void plot_effective_potential_for_n_states(params p, QPC_params QPC, int n_states, vector<Eigen::VectorXcd> &Eigenvals_matrix)
{
    std::ofstream outfile("./../data/En_x.dat");
    for (int i = 0; i < QPC.xpoints; i++)
    {
        double x = double(i) / (QPC.xpoints - 1) * QPC.length;
        outfile << x * p.bohr_radius << " ";
        for (int state = 0; state < n_states; state++)
        {
            outfile << Eigenvals_matrix[i][state].real() * p.hartree_to_meV / 1000 << " ";
        }
        outfile << "\n";
    }
    outfile.close();
}

void plot_QPC_effective_potential(params p, QPC_params QPC, vector<Eigen::VectorXcd> &Eigenvals_matrix)
{
    matrix<double> qpcPot(QPC.xpoints, QPC.ypoints);
    for (int i = 0; i < QPC.xpoints; i++)
    {
        for (int j = 0; j < QPC.ypoints; j++)
        {
            qpcPot[i][j] = Eigenvals_matrix[i][j].real() * p.hartree_to_meV / 1000;
            // if(qpcPot[i][j] > 0.25){printf("%6f  %6f  %10f \n", x, y, qpcPot[i][j]);}
        }
    }
    qpcPot.plot_heatmap("QPC_effpot.png", false, "");
}

void calculate_conductance_vs_electron_energy(params p, std::string OutFileName, QPC_params QPC, vector<Eigen::VectorXcd> &Eigenvals_matrix)
{
    std::ofstream outfile;
    outfile.open("./../data/"+OutFileName);
    std::cout << "\nStarting conductance calculations:\n";
    for (double E = 0 / p.hartree_to_meV; E < 200 / p.hartree_to_meV; E += 1 / p.hartree_to_meV)
    {
        double cond = get_conductance_Landauer(p, QPC, Eigenvals_matrix, E, 10);
        outfile << E * p.hartree_to_meV / 1000 << " " << cond << "\n";
    }
    outfile.close();
}

void calculate_conductance_vs_gate_voltage(params p, QPC_params QPC, double ElectronEnergy, std::string OutFileName, int included_states)
{
    std::vector<Eigen::VectorXcd> Eigenvals_matrix;
    std::ofstream outfile;
    outfile.open("./../data/"+OutFileName);
    std::cout << "\nStarting conductance as a gate voltage function calculations\nIncident electron energy: " << std::round(ElectronEnergy * p.hartree_to_meV) << " meV\n";
    for (double V = 0; V <= 25000 / p.hartree_to_meV; V += 200 / p.hartree_to_meV)
    {
        QPC.set_V_gate(V);
        Eigenvals_matrix = get_eigen_energies(p, QPC, false);
        double cond = get_conductance_Landauer(p, QPC, Eigenvals_matrix, ElectronEnergy, included_states);
        outfile << V * p.hartree_to_meV / 1000 << " " << cond << "\n";
        progress(V / (25000 / p.hartree_to_meV));
    }
}

int main()
{
    params p;
    QPC_params QPC;

    vector<double> U(p.N);
    vector<double> m(p.N);
    vector<double> z(p.N);

    /* ==============================
    sets the spacial vector, it's not necessery, however it's convinient to not worry about indexes and spacial steps
    ==============================*/
    set_z_vec(z, 0, p.L, p.N);

    // set_U_m_single_barrier(p, U, m);
    set_U_m_double_barrier(p, U, m);

    set_U_V_bias(p, U, 50 / p.hartree_to_meV);

    /*==============================
    writes into a file spacial distribution of potential and effective masses for set system
    ==============================*/
    write_U_m_dist(p, z, U, m, "double_barier_check.dat");

    /* ==============================
    calculates transmitance and reflectance for set system using transfer matrix
    method and writes an output into a file
    ==============================*/
    // calculate_T_R_for_system(p, z, U, m, "T_R_E_double_barier.dat");

    /*==============================
    calculates current-voltage characteristic using Tsu-Esaki formula, and writes an output into a file
    ==============================*/
    // calculate_current_voltage_char_Tsu_Esaki(p, z, U, m, "I_V_char.dat");

    // plot_QPC_potential(p, QPC);

    // vector<Eigen::VectorXcd> Eigenvals_matrix = get_eigen_energies(p, QPC);
    // plot_effective_potential_for_n_states(p, QPC, 5, Eigenvals_matrix);
    // calculate_conductance_vs_electron_energy(p, "cond_vs_E.dat", QPC, Eigenvals_matrix);

    // double ElectronEnergy = 50 / p.hartree_to_meV;
    // int included_states = 8;
    // calculate_conductance_vs_gate_voltage(p, QPC, ElectronEnergy, "cond_vs_V_1.dat", included_states);
    // ElectronEnergy = 100 / p.hartree_to_meV;
    // calculate_conductance_vs_gate_voltage(p, QPC, ElectronEnergy, "cond_vs_V_2.dat", included_states);

    std::cout << std::endl;
    return 0;
}