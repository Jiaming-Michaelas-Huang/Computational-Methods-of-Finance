//
//  Matrix.cpp
//  Computational_Methods_Final
//
//  Created by 黄佳铭 on 2019/3/14.
//  Copyright © 2019 Jiaming Huang. All rights reserved.
//

#include "Matrix.hpp"
#include <iostream>
#include <cmath>
using namespace std;

Matrix::Matrix(double value, int nrow, int ncol){
    this->nrow = nrow;
    this->ncol = ncol;
    this->matrix = new double* [nrow];
    for (int i = 0; i < nrow; i++) {
        this->matrix[i] = new double[ncol];
        for (int j = 0; j < ncol; j++) {
            this->matrix[i][j] = value;
        }
    }
}

Matrix::Matrix(double ** m, int nrow, int ncol){
    this->nrow = nrow;
    this->ncol = ncol;
    this->matrix = new double* [nrow];
    for (int i = 0; i < nrow; i++) {
        this->matrix[i] = new double[ncol];
        for (int j = 0; j < ncol; j++) {
            this->matrix[i][j] = m[i][j];
        }
    }
}

Matrix::Matrix(double * m, int nrow,int ncol){
    this->nrow = nrow;
    this->ncol = ncol;
    this->matrix = new double* [nrow];
    for (int i = 0; i < nrow; i++) {
        this->matrix[i] = new double[ncol];
        for (int j = 0; j < ncol; j++) {
            this->matrix[i][j] = m[i*ncol+j];
        }
    }
}

Matrix::~Matrix(){
    for (int i = 0; i<this->nrow; i++) {
        delete [] matrix[i];
    }
    delete []this->matrix;
    //cout<<"Deleted"<<endl;
}

double Matrix::elementAt(int row, int col){
    return this->matrix[row][col];
}

void Matrix::setValue(double value, int row, int col){
    this->matrix[row][col] = value;
}

shared_ptr<Matrix> Matrix::subMatrix(int startrow, int endrow, int startcol, int endcol){
    double ** result = new double * [endrow-startrow];
    int rowindex = 0;
    for (int i = startrow; i < endrow; i++) {
        result[rowindex] = new double[endcol-startcol];
        int colindex = 0;
        for (int j = startcol; j < endcol; j++) {
            result[rowindex][colindex] = this->matrix[i][j];
            colindex++;
        }
        rowindex++;
    }
    shared_ptr<Matrix> re = make_shared<Matrix>(result,endrow-startrow, endcol-startcol);
    for (int i = 0; i < endrow-startrow; i++) {
        delete [] result[i];
    }
    delete [] result;
    return re;
}

shared_ptr<Matrix> Matrix::Merge(int axis, shared_ptr<Matrix> m){
    if (axis == 0) {
        double ** result = new double * [this->nrow];
        for (int i = 0; i < this->nrow; i++) {
            result[i] = new double[this->ncol+m->ncol];
            for (int j = 0; j < this->ncol; j++) {
                result[i][j] = this->matrix[i][j];
            }
            for (int j = 0; j < m->ncol; j++) {
                result[i][this->ncol+j] = m->matrix[i][j];
            }
        }
        shared_ptr<Matrix> re = make_shared<Matrix>(result, this->nrow, this->ncol+m->ncol);
        for (int i = 0; i < this->nrow; i++) {
            delete [] result[i];
        }
        delete [] result;
        return re;
    }
    else{
        double ** result = new double * [this->nrow+m->nrow];
        for (int i = 0; i < this->nrow; i++) {
            result[i] = new double[this->ncol];
            for (int j = 0; j < this->ncol; j++) {
                result[i][j] = this->matrix[i][j];
            }
        }
        for (int i = 0; i < m->nrow; i++) {
            result[i+this->nrow] = new double[this->ncol];
            for (int j = 0; j < this->ncol; j++) {
                result[i+this->nrow][j] = m->matrix[i][j];
            }
        }
        shared_ptr<Matrix> re = make_shared<Matrix>(result, this->nrow+m->nrow, this->ncol);
        for (int i = 0; i < this->nrow+m->nrow; i++) {
            delete [] result[i];
        }
        delete [] result;
        return re;
    }
}

shared_ptr<Matrix> Matrix::operator+(shared_ptr<Matrix> m){
    double ** result = new double * [this->nrow];
    for (int i = 0; i < this->nrow; i++) {
        result[i] = new double[this->ncol];
    }
    if (this->nrow!=m->nrow||this->ncol!=m->ncol) {
        cout<<"Two Matrix cannot be added"<<endl;
    }
    else{
        for (int i = 0; i < this->nrow; i++) {
            for (int j = 0; j<this->ncol; j++) {
                result[i][j] = this->matrix[i][j]+m->matrix[i][j];
            }
        }
    }
    shared_ptr<Matrix> re = make_shared<Matrix>(result, this->nrow, this->ncol);
    for (int i = 0; i < this->nrow; i++) {
        delete [] result[i];
    }
    delete [] result;
    return re;
}

shared_ptr<Matrix> Matrix::operator-(shared_ptr<Matrix> m){
    double ** result = new double * [this->nrow];
    for (int i = 0; i < this->nrow; i++) {
        result[i] = new double[this->ncol];
    }
    if (this->nrow!=m->nrow||this->ncol!=m->ncol) {
        cout<<"Two Matrix cannot be added"<<endl;
    }
    else{
        for (int i = 0; i < this->nrow; i++) {
            for (int j = 0; j<this->ncol; j++) {
                result[i][j] = this->matrix[i][j]-m->matrix[i][j];
            }
        }
    }
    shared_ptr<Matrix> re = make_shared<Matrix>(result, this->nrow, this->ncol);
    for (int i = 0; i < this->nrow; i++) {
        delete [] result[i];
    }
    delete [] result;
    return re;
}

shared_ptr<Matrix> Matrix::operator*(shared_ptr<Matrix> m){
    double ** result = new double * [this->nrow];
    for (int i = 0; i < nrow; i++) {
        result[i] = new double[m->ncol];
    }
    if (this->ncol != m->nrow) {
        cout<<"Two Matrix cannot be multiplied"<<endl;
    }
    else{
        for (int i = 0; i < this->nrow; i++) {
            for (int j = 0; j < m->ncol; j++) {
                double sum = 0;
                for (int k = 0; k < this->ncol; k++) {
                    sum += this->matrix[i][k]*m->matrix[k][j];
                }
                result[i][j] = sum;
            }
        }
    }
    shared_ptr<Matrix> re = make_shared<Matrix>(result,this->nrow, m->ncol);
    for (int i = 0; i < this->nrow; i++) {
        delete [] result[i];
    }
    delete [] result;
    return re;
}

shared_ptr<Matrix> Matrix::operator+(double n){
    double ** result = new double * [this->nrow];
    for (int i = 0; i < this->nrow; i++) {
        result[i] = new double[this->ncol];
    }
    for (int i = 0; i<this->nrow; i++) {
        for (int j = 0; j < this->ncol; j++) {
            result[i][j] = n + this->matrix[i][j];
        }
    }
    shared_ptr<Matrix> re = make_shared<Matrix>(result, this->nrow, this->ncol);
    for (int i = 0; i < this->nrow; i++) {
        delete [] result[i];
    }
    delete [] result;
    return re;
}

shared_ptr<Matrix> Matrix::operator-(double n){
    double ** result = new double * [this->nrow];
    for (int i = 0; i < this->nrow; i++) {
        result[i] = new double[this->ncol];
    }
    for (int i = 0; i<this->nrow; i++) {
        for (int j = 0; j < this->ncol; j++) {
            result[i][j] = this->matrix[i][j] - n;
        }
    }
    shared_ptr<Matrix> re = make_shared<Matrix>(result, this->nrow, this->ncol);
    for (int i = 0; i < this->nrow; i++) {
        delete [] result[i];
    }
    delete [] result;
    return re;
}

shared_ptr<Matrix> Matrix::operator*(double n){
    double ** result = new double * [this->nrow];
    for (int i = 0; i < this->nrow; i++) {
        result[i] = new double[this->ncol];
    }
    for (int i = 0; i<this->nrow; i++) {
        for (int j = 0; j < this->ncol; j++) {
            result[i][j] = n * this->matrix[i][j];
        }
    }
    shared_ptr<Matrix> re = make_shared<Matrix>(result, this->nrow, this->ncol);
    for (int i = 0; i < this->nrow; i++) {
        delete [] result[i];
    }
    delete [] result;
    return re;
}

shared_ptr<Matrix> Matrix::max(double n){
    double ** result = new double * [this->nrow];
    for (int i = 0; i < this->nrow; i++) {
        result[i] = new double[this->ncol];
    }
    for (int i = 0; i<this->nrow; i++) {
        for (int j = 0; j < this->ncol; j++) {
            result[i][j] = this->matrix[i][j]>n?this->matrix[i][j]:n;
        }
    }
    shared_ptr<Matrix> re = make_shared<Matrix>(result, this->nrow, this->ncol);
    for (int i = 0; i < this->nrow; i++) {
        delete [] result[i];
    }
    delete [] result;
    return re;
}

shared_ptr<Matrix> Matrix::min(double n){
    double ** result = new double * [this->nrow];
    for (int i = 0; i < this->nrow; i++) {
        result[i] = new double[this->ncol];
    }
    for (int i = 0; i<this->nrow; i++) {
        for (int j = 0; j < this->ncol; j++) {
            result[i][j] = this->matrix[i][j]<n?this->matrix[i][j]:n;
        }
    }
    shared_ptr<Matrix> re = make_shared<Matrix>(result, this->nrow, this->ncol);
    for (int i = 0; i < this->nrow; i++) {
        delete [] result[i];
    }
    delete [] result;
    return re;
}

shared_ptr<Matrix> Matrix::cumsum(int axis){
    if (axis == 0) {
        double * result = new double[this->nrow];
        for (int i = 0; i < this->nrow; i++) {
            double sum = 0;
            for (int j = 0; j < this->ncol; j++) {
                sum += this->matrix[i][j];
            }
            result[i] = sum;
        }
        shared_ptr<Matrix> re = make_shared<Matrix>(result, this->nrow,1);
        delete [] result;
        return re;
    }
    else{
        double * result = new double[this->ncol];
        for (int i = 0; i < this->ncol; i++) {
            double sum = 0;
            for (int j = 0; j < this->nrow; j++) {
                sum += this->matrix[j][i];
            }
            result[i] = sum;
        }
        shared_ptr<Matrix> re = make_shared<Matrix>(result, 1,this->ncol);
        delete [] result;
        return re;
    }
}

double Matrix::getAdj(double ** old,int row, int col, int numrow, int numcol){
    double ** result = new double * [numrow-1];
    for (int i = 0; i < numrow-1; i++) {
        result[i] = new double[numcol-1];
    }
    int i1 = 0;
    int j1 = 0;
    for (int i = 0; i < numrow; i++) {
        if(i == row)continue;
        j1 = 0;
        for (int j = 0; j < numcol; j++) {
            if (j == col) {
                continue;
            }
            else{
                result[i1][j1] = old[i][j];
                j1++;
            }
        }
        i1++;
    }
    double re = (row+col)%2 == 0?getDerterminant(result, numrow-1, numcol-1):-1*getDerterminant(result,numrow-1,numcol-1);
    for (int i = 0; i < numrow-1; i++) {
        delete [] result[i];
    }
    delete [] result;
    return re;
}

shared_ptr<Matrix> Matrix::getAdjMatrix(){
    double ** result = new double * [this->nrow];
    for (int i = 0; i < nrow; i++) {
        result[i] = new double[this->ncol];
    }
    for (int i = 0; i<this->nrow; i++) {
        for (int j = 0; j<this->ncol; j++) {
            result[i][j] = getAdj(this->matrix,i, j, this->nrow, this->ncol);
        }
    }
    shared_ptr<Matrix> re = make_shared<Matrix>(result, this->nrow, this->ncol);
    for (int i = 0; i < this->nrow; i++) {
        delete [] result[i];
    }
    delete [] result;
    return re;
}

double Matrix::getDerterminant(double ** old, int numrow, int numcol){
    double result = 0;
    if (numrow!=numcol) {
        cout<<"Matrix's determinant cannot be calculated"<<endl;
        return 0;
    }
    else{
        if (numcol>1) {
            for (int j = 0; j < numcol; j++) {
                if (old[0][j] == 0) {
                    continue;
                }
                double adj = getAdj(old, 0, j, numrow,numcol);
                result += old[0][j]*adj;
            }
        }
        else{
            result = old[0][0];
        }
        return result;
    }
}

shared_ptr<Matrix> Matrix::Trans(){
    double ** result = new double * [this->ncol];
    for (int i = 0; i < nrow; i++) {
        result[i] = new double[this->nrow];
    }
    for (int i = 0; i < this->nrow; i++) {
        for (int j = 0; j < this->ncol; j++) {
            result[j][i] = this->matrix[i][j];
        }
    }
    shared_ptr<Matrix> re = make_shared<Matrix>(result, this->ncol,this->nrow);
    for (int i = 0; i < this->ncol; i++) {
        delete [] result[i];
    }
    delete [] result;
    return re;
}

shared_ptr<Matrix> Matrix::Inverse(){
    
    if(nrow!=ncol)return NULL;
    
    int n = this->nrow*nrow;
    int s = this->nrow;
    
    double a[n];
    
    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) {
            a[j+i*nrow] = this->matrix[i][j];
        }
    }
    
    double *L = NULL;
    double *U = NULL;
    double tmp=0;
    
    
    L = new double[n];
    U = new double[n];
    
    for (int i = 0; i < s; i++)
    {
        for (int j = 0; j < s; j++)
        {
            if (i == j)
                L[i*s+j] = 1;
            if (i < j)
                L[i*s+j] = 0;
            if (i > j)
                U[i*s+j] = 0;
            
            U[0*s+j] = a[0*s+j];
            L[i*s+0] = a[i*s+0] / U[0*s+0];
        }
    }
    
    for (int k = 1; k < s; k++)
    {
        
        for (int j = k; j < s; j++)
        {
            tmp = 0;
            for (int m = 0; m < k; m++)
            {
                tmp += L[k*s+m] * U[m*s+j];
            }
            
            U[k*s+j] = a[k*s+j] - tmp;
        }
        
        for (int i = k+1; i < s; i++)
        {
            tmp = 0;
            for (int m = 0; m < k; m++)
            {
                tmp += L[i*s+m] * U[m*s+k];
            }
            
            L[i*s+k] = ( a[i*s+k] - tmp ) / U[k*s+k];
        }
    }
    double *u = new double[s*s];
    double *r = new double[s*s];
    double ss;
    for (int i = 0; i < s; i++) {
        for (int j = 0; j < s; j++) {
            u[i*s+j] = 0;
            r[i*s+j] = 0;
        }
    }
    for (int i=0;i<s;i++) /*求矩阵U的逆 */
    {
        u[i*s+i]=1/U[i*s+i];//对角元素的值，直接取倒数
        for (int k=i-1;k>=0;k--)
        {
            ss=0;
            for (int j=k+1;j<=i;j++)
                ss=ss+U[k*s+j]*u[j*s+i];
            u[k*s+i]=-ss/U[k*s+k];//迭代计算，按列倒序依次得到每一个值，
        }
    }
    for (int i=0;i<s;i++) //求矩阵L的逆
    {
        r[i*s+i]=1; //对角元素的值，直接取倒数，这里为1
        for (int k=i+1;k<s;k++)
        {
            for (int j=i;j<=k-1;j++)
                r[k*s+i]=r[k*s+i]-L[k*s+j]*r[j*s+i];   //迭代计算，按列顺序依次得到每一个值
        }
    }
    
    shared_ptr<Matrix> matrixIL = make_shared<Matrix>(r,s,s);
    shared_ptr<Matrix> matrixIU = make_shared<Matrix>(u,s,s);
    matrixIL->Print();
    matrixIU->Print();
    return (*matrixIU)*matrixIL;
}

void Matrix::Print(){
    for (int i = 0; i < this->nrow; i++) {
        for (int j = 0; j < this->ncol; j++) {
            cout<<this->matrix[i][j];
            cout<<" ";
        }
        cout<<endl;
    }
    cout<<endl;
}

shared_ptr<Matrix> Matrix::multi(shared_ptr<Matrix> m){
    double ** result = new double * [this->nrow];
    for (int i = 0; i < this->nrow; i++) {
        result[i] = new double[this->ncol];
    }
    if (this->nrow!=m->nrow||this->ncol!=m->ncol) {
        cout<<"Two Matrix cannot be multiplied"<<endl;
    }
    else{
        for (int i = 0; i < this->nrow; i++) {
            for (int j = 0; j<this->ncol; j++) {
                result[i][j] = this->matrix[i][j]*m->matrix[i][j];
            }
        }
    }
    shared_ptr<Matrix> re = make_shared<Matrix>(result, this->nrow, this->ncol);
    for (int i = 0; i < this->nrow; i++) {
        delete [] result[i];
    }
    delete []result;
    return re;
}
