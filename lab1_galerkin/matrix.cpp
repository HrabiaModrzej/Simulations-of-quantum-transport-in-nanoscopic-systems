#include <cmath>
#include <iostream>
#include <fstream>
#include "matrix.h"

using namespace std;

matrix::matrix(int x, int y)
{
    sizeX = x;
    sizeY = y;
    tab = new double *[sizeX];
    for (int i = 0; i < sizeX; i++)
    {
        tab[i] = new double[sizeY];
    }
}

matrix::matrix(int a)
{
    sizeX = a;
    sizeY = a;
    tab = new double *[sizeX];
    for (int i = 0; i < sizeX; i++)
    {
        tab[i] = new double[sizeY];
    }
}

matrix::~matrix()
{
    for (int i = 0; i < sizeX; i++)
    {
        delete[] tab[i];
    }
    delete[] tab;
}

matrix &matrix::operator=(const matrix &Matrix)
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

matrix::matrix(const matrix &Matrix)
{
    sizeX = Matrix.sizeX;
    sizeY = Matrix.sizeY;
    tab = new double *[sizeX];
    for (int i = 0; i < sizeX; i++)
    {
        tab[i] = new double[sizeY];
    }

    *this = Matrix;
}

ostream &operator<<(ostream &os, const matrix &M)
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

double *matrix::operator[](int i)
{
    return tab[i];
}

void matrix::operator*=(double val)
{
    for (int i = 0; i < sizeX; i++)
    {
        for (int j = 0; j < sizeY; j++)
        {
            tab[i][j] *= val;
        }
    }
}

void matrix::operator/=(double val)
{
    *this *= (1 / val);
}

void matrix::save_to_file(std::string path)
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

void matrix::plot_heatmap(std::string path = "./heatmap.png", bool save_temp = false, std::string additional_line = "")
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

void matrix::set_all_elements(double val = 0.0)
{
    for (int i = 0; i < sizeX; i++)
    {
        for (int j = 0; j < sizeY; j++)
        {
            tab[i][j] = val;
        }
    }
}

matrix matrix::operator+(matrix &other_matrix)
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

matrix matrix::operator*(double val){
    matrix out(*this);
    out *= val;
    return out;
}