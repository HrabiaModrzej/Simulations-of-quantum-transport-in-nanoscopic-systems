#include<cmath>
#include<iostream>
#include<fstream>

using namespace std;

#ifndef MATRIX_H_
#define MATRIX_H_

class matrix{

    public:
        double** tab;    
        int sizeX;
        int sizeY;
        matrix(int,int);
        matrix(int);
        matrix()=delete;
        ~matrix();
        matrix& operator=(const matrix&);
        matrix(const matrix&);
        friend ostream& operator<<(ostream&, const matrix&);
        double* operator[](int);
        void operator*=(double);
        void operator/=(double);
        void save_to_file(std::string);
        void plot_heatmap(std::string,bool,std::string);
        void set_all_elements(double);
        matrix operator+(matrix&);
        matrix operator*(double);

};





#endif
