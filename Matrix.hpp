//
//  Matrix.hpp
//  Computational_Methods_Final
//
//  Created by 黄佳铭 on 2019/3/14.
//  Copyright © 2019 Jiaming Huang. All rights reserved.
//

#ifndef Matrix_hpp
#define Matrix_hpp

#include <stdio.h>
#include <iostream>
#include <cmath>
using namespace std;

class Matrix {
private:
    double ** matrix;
    int nrow;
    int ncol;
public:
    Matrix(double value,int nrow, int ncol);
    Matrix(double ** matrix, int nrow, int ncol);
    Matrix(double * matrix, int nrow, int ncol);
    ~Matrix();
    double elementAt(int row, int col);
    void setValue(double value, int row, int col);
    shared_ptr<Matrix> subMatrix(int startrow, int endrow, int startcol, int endcol);
    shared_ptr<Matrix> Merge(int axis, shared_ptr<Matrix> m);
    shared_ptr<Matrix> operator+(double n);
    shared_ptr<Matrix> operator-(double n);
    shared_ptr<Matrix> max(double n);
    shared_ptr<Matrix> min(double n);
    shared_ptr<Matrix> operator+(shared_ptr<Matrix> m);
    shared_ptr<Matrix>operator-(shared_ptr<Matrix> m);
    shared_ptr<Matrix>multi(shared_ptr<Matrix> m);
    shared_ptr<Matrix>operator*(shared_ptr<Matrix> m);
    shared_ptr<Matrix>operator*(double n);
    shared_ptr<Matrix>cumsum(int axis);
    double getDerterminant(double ** old, int numrow, int numcol);
    shared_ptr<Matrix>getAdjMatrix();
    double getAdj(double ** old, int row, int col, int numrow, int numcol);
    shared_ptr<Matrix>Trans();
    shared_ptr<Matrix>Inverse();
    void Print();
};
#endif /* Matrix_hpp */
