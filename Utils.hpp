//
//  Utils.hpp
//  Computational_Methods_Final
//
//  Created by 黄佳铭 on 2019/3/13.
//  Copyright © 2019 Jiaming Huang. All rights reserved.
//

#ifndef Utils_hpp
#define Utils_hpp

#include <stdio.h>
#include <iostream>
#include <cmath>

using namespace std;

class Utils{
public:
    static double prob(double N, double * distribution, int size);
    static int factorial(int n);
    static double mean(double * arr, int size);
    static double Std(double * arr, int size);
    static double corr(double * arr1, double * arr2, int size);
    static double cov(double * arr1, double * arr2, int size);
    static void tocsv(double * array, string filename, int size);
    static double pnorm(double quantile);
    static double Integral(double a, double b, int n);
    static double Derivative(double(*f)(double x,int n), double x, int n);
};
#endif /* Utils_hpp */
