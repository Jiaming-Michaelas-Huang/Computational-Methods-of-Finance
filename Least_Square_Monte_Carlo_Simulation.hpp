//
//  Least_Square_Monte_Carlo_Simulation.hpp
//  Computational_Methods_Final
//
//  Created by 黄佳铭 on 2019/3/14.
//  Copyright © 2019 Jiaming Huang. All rights reserved.
//

#ifndef Least_Square_Monte_Carlo_Simulation_hpp
#define Least_Square_Monte_Carlo_Simulation_hpp

#include <stdio.h>
#include <iostream>
#include "Matrix.hpp"

using namespace std;

class Least_Square_Monte_Carlo_Simulation{
public:
    static double Hermite(double x,int n);
    static double Lagueree(double x,int n);
    static double Monomial(double x,int n);
    static shared_ptr<Matrix> MultiStockPathsSimulation(double S0, double T, int n, double mu, double sigma, int seed, int size);
    static shared_ptr<Matrix> LeastSquaredRegression(shared_ptr<Matrix> St,shared_ptr<Matrix> Yt, string type, int k, int size, int n);
    static double EuroForwardStartPricing(double t, double S0, double T, int n, double mu, double sigma, int size);
};
#endif /* Least_Square_Monte_Carlo_Simulation_hpp */
