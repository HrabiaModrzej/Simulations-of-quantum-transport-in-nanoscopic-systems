#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;

#ifndef MATRIX_H_
#define MATRIX_H_

template <class val_type>
class matrix; // wstępna deklaracja

template <class val_type>
std::ostream &operator<<(std::ostream &, const matrix<val_type> &); // <-- DODAJ TO

template <class val_type>
class matrix
{

public:
    val_type **tab;
    int sizeX;
    int sizeY;
    matrix(int, int);
    matrix(int);
    matrix() = delete;
    ~matrix();
    matrix<val_type> &operator=(const matrix<val_type> &);
    matrix(const matrix<val_type> &);
    friend ostream &operator<< <>(ostream &, const matrix<val_type> &);
    val_type *operator[](int);
    void operator*=(val_type);
    void operator/=(val_type);
    void save_to_file(std::string);
    void plot_heatmap(std::string, bool, std::string);
    void set_all_elements(val_type);
    matrix<val_type> operator+(matrix<val_type> &);
    matrix<val_type> operator*(val_type);
    void identity_matrix();
    // matrix<val_type> operator*(martix<val_type> &);
};

#endif

template <class val_type>
matrix<val_type>::matrix(int x, int y)
{
    sizeX = x;
    sizeY = y;
    tab = new val_type *[sizeX];
    for (int i = 0; i < sizeX; i++)
    {
        tab[i] = new val_type[sizeY];
    }
}

template <class val_type>
matrix<val_type>::matrix(int a)
{
    sizeX = a;
    sizeY = a;
    tab = new val_type *[sizeX];
    for (int i = 0; i < sizeX; i++)
    {
        tab[i] = new val_type[sizeY];
    }
}

template <class val_type>
matrix<val_type>::~matrix()
{
    for (int i = 0; i < sizeX; i++)
    {
        delete[] tab[i];
    }
    delete[] tab;
}
template <class val_type>
matrix<val_type> &matrix<val_type>::operator=(const matrix<val_type> &Matrix)
{
    if (sizeX == Matrix.sizeX && sizeY == Matrix.sizeY)
    {
        for (int i = 0; i < sizeX; i++)
        {
            for (int j = 0; j < sizeY; j++)
            {
                tab[i][j] = Matrix.tab[i][j];
            }
        }
    }
    return *this;
}

template <class val_type>
matrix<val_type>::matrix(const matrix<val_type> &Matrix)
{
    sizeX = Matrix.sizeX;
    sizeY = Matrix.sizeY;
    tab = new val_type *[sizeX];
    for (int i = 0; i < sizeX; i++)
    {
        tab[i] = new val_type[sizeY];
    }

    *this = Matrix;
}
template <class val_type>
std::ostream &operator<<(std::ostream &os, const matrix<val_type> &M)
{
    for (int j = 0; j < M.sizeY; j++)
    {
        os << "| ";
        for (int i = 0; i < M.sizeX; i++)
        {
            os << M.tab[i][j] << " ";
        }
        os << "|" << endl;
    }
    return os;
}

template <class val_type>
val_type *matrix<val_type>::operator[](int i)
{
    return tab[i];
}

template <class val_type>
void matrix<val_type>::operator*=(val_type val)
{
    for (int i = 0; i < sizeX; i++)
    {
        for (int j = 0; j < sizeY; j++)
        {
            tab[i][j] *= val;
        }
    }
}

template <class val_type>
void matrix<val_type>::operator/=(val_type val)
{
    *this *= (1 / val);
}

template <class val_type>
void matrix<val_type>::save_to_file(std::string path)
{
    ofstream file(path);
    for (int i = 0; i < sizeX; i++)
    {
        for (int j = 0; j < sizeY; j++)
        {
            file << i << " " << j << " " << tab[i][j] << endl;
        }
        file << endl;
    }
    file.close();
}

template <class val_type>
void matrix<val_type>::plot_heatmap(std::string path, bool save_temp, std::string additional_line)
{
    save_to_file("./temp_data.dat");
    ofstream file("./temp_plot_script.gnuplot");
    file << "set terminal png size 700,700 font \"Latin Modern Math,18\" background rgb 'white'\n";
    file << "set size ratio -1\n";
    file << "set view map\n";
    file << "unset xtics\n";
    file << "unset ytics\n";
    file << "set out '" + path + "'\n";
    file << additional_line;
    file << "splot './temp_data.dat' u 1:2:3 w pm3d notitle";

    file.close();
    system("gnuplot temp_plot_script.gnuplot");

    if (save_temp == false)
    {
        system("rm temp*");
    }
}

template <class val_type>
void matrix<val_type>::set_all_elements(val_type val)
{
    for (int i = 0; i < sizeX; i++)
    {
        for (int j = 0; j < sizeY; j++)
        {
            tab[i][j] = val;
        }
    }
}

template <class val_type>
matrix<val_type> matrix<val_type>::operator+(matrix<val_type> &other_matrix)
{
    if (other_matrix.sizeX != sizeX || other_matrix.sizeY != sizeY)
    {
        std::cerr << "Przy dodawaniu nierówny rozmiar macierzy!!!";
        return matrix(1, 1);
    }

    matrix out(sizeX, sizeY);

    for (int i = 0; i < sizeX; i++)
    {
        for (int j = 0; j < sizeY; j++)
        {
            out[i][j] = other_matrix[i][j] + tab[i][j];
        }
    }
    return out;
}

template <class val_type>
matrix<val_type> matrix<val_type>::operator*(val_type val)
{
    matrix out(*this);
    out *= val;
    return out;
}

template <class val_type>
void matrix<val_type>::identity_matrix()
{
    if (this->sizeX != this->sizeY)
    {
        std::runtime_error("This matrix isn't a square matrix.\n");
    }
    for (int i = 0; i < this->sizeX; i++)
    {
        this->tab[i][i] = 1;
    }
}

// template <class val_type>
// matrix<val_type> matrix<val_type>::operator*(matrix<val_type> &M)
// {
//     matrix<complex<double>> out(this->sizeX, M.sizeY);
//     out.set_all_elements(-10000.0);

//     if (this->sizeX != M.sizeY)
//     {
//         std::cerr << "Multiplyed matrices: wrong size\n";
//         return out;
//     }

//     for (int i = 0; i < out.sizeX; i++)
//     {
//         for (int j = 0; j < out.sizeY; j++)
//         {
//             out[i][j] = 0.0;
//             for (int k = 0; k < this->sizeX; k++)
//             {
//                 out[i][j] += this->tab[k][j] * M[i][k];
//             }
//         }
//     }
//     return out;
// }