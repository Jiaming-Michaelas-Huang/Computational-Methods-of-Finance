//
//  Utils.cpp
//  Computational_Methods_Final
//
//  Created by 黄佳铭 on 2019/3/13.
//  Copyright © 2019 Jiaming Huang. All rights reserved.
//

#include "Utils.hpp"
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>

using namespace std;

double Utils::prob(double N, double * distribution, int size){
    double p = 0;
    for (int i = 0; i < size; i++) {
        if(distribution[i] < N)p++;
    }
    return p/size;
}

int Utils::factorial(int n){
    int sum = 1;
    for (int i = n; i > 0; i--) {
        sum = sum * i;
    }
    return sum;
}

double Utils::mean(double *arr, int size){
    double sum = 0;
    for (int i = 0; i < size; i++) {
        sum = sum+arr[i];
    }
    return sum/size;
}

double Utils::Std(double *arr, int size){
    double m = mean(arr, size);
    double sum = 0;
    for (int i = 0; i < size; i++) {
        sum = sum + (arr[i] - m) * (arr[i] - m);
    }
    return sum/size;
}

double Utils::corr(double *arr1, double *arr2, int size){
    double mean1 = mean(arr1,size);
    double mean2 = mean(arr2, size);
    double sum1 = 0, sum2 = 0, sum3 = 0;
    for (int i = 0; i < size; i++) {
        sum1 = sum1 + (arr1[i] - mean1)*(arr2[i] - mean2);
        sum2 = sum2 + (arr1[i] - mean1)*(arr1[i] - mean1);
        sum3 = sum3 + (arr2[i] - mean2)*(arr2[i] - mean2);
    }
    return sum1 / (sqrt(sum2)*sqrt(sum3));
}

double Utils::cov(double *arr1, double *arr2, int size){
    double mean1 = mean(arr1,size);
    double mean2 = mean(arr2, size);
    double sum = 0;
    for (int i = 0; i < size; i++) {
        sum = sum + (arr1[i] - mean1)*(arr2[i] - mean2);
    }
    return sum/size;
}

void Utils::tocsv(double *array, string filename, int size){
    ofstream ofile;
    ofile.open("/Users/jiaminghuang/Desktop/Computational-Methods-Project-1/data/hjm/"+filename+".csv");
    for (int i = 0; i < size; i++)
    {
        ofile << array[i] << ",\n";
    }
    ofile << "\n";
    ofile.close();
}

double Utils::pnorm(double quantile){
    const double d1 = 0.0498673470;
    const double d2 = 0.0211410061;
    const double d3 = 0.0032776263;
    const double d4 = 0.0000380036;
    const double d5 = 0.0000488906;
    const double d6 = 0.0000053830;
    
    long double probability = 0.0;
    long double absq = (quantile >=0)?quantile:-quantile;
    long double temp =  1 + (d1 * absq) + (d2 * absq * absq) + (d3 * absq * absq * absq) + (d4 * pow(absq,4))
    + (d5 * pow(absq,5)) + (d6 * pow(absq, 6));
    probability = 1 - 0.5 * pow(temp, -16);
    
    return (quantile >= 0)? probability : (1- probability);
}

double fun(double x)
{
    return 2*x;
}

double Utils::Integral(double a, double b, int n){
    double w;
    double k=(b-a)/n;
    double s=0.0;
    for(int i=1;i<=n;i++)
    {
        w=fun(a+(i-1)*k)*k;
        s=s+w;
    }
    return s;
}

double Utils::Derivative(double (*f)(double x, int n), double x, int n){
    if(n==1){
        double result = (f(x+0.001,n)-f(x,n))/0.001;
        return result;
    }
    else{
        double result = (Derivative(f, x+0.001, n-1)-Derivative(f, x, n-1))/0.001;
        return result;
    }
}
