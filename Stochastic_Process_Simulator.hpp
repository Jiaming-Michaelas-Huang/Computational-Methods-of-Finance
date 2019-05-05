//
//  Stochastic_Process_Simulator.hpp
//  Computational_Methods_Final
//
//  Created by 黄佳铭 on 2019/3/14.
//  Copyright © 2019 Jiaming Huang. All rights reserved.
//

#ifndef Stochastic_Process_Simulator_hpp
#define Stochastic_Process_Simulator_hpp

#include <stdio.h>
class Stochastic_Process_Simulation{
public:
    static double * WinierPathSimulation(double T, int n, int seed);
    static double * WinierProcessSimulation(double T, int n, int size, int seed);
    static double * StockPathSimulation(double S0, double T, int n, double mu, double sigma, double * winierpath);
    static double * StockPathSimulationWithJump(double S0, double T, int n, double mu, double sigma,double gamma, double * winierpath, double * poisson);
    static double * StochasticPathSimulationEulerSchema(double x0, double T, int n, double * winier,double (*a)(double x),double (*b)(double x));
    static double * StochasticPathSimulationMilshteinSchema(double x0, double T, int n, double * winier,double (*a)(double x),double (*b)(double x));
    static double * HestonModelReflectionSchema(double *w1, double *w2, double T, int n,double r, double s0, double v0, double sigma, double alpha, double beta);
    static double * HestonModelPartialTruncationSchema(double *w1, double *w2, double T, int n,double r, double s0, double v0, double sigma, double alpha, double beta);
    static double * HestonModelFullTruncationSchema(double *w1, double *w2, double T, int n,double r, double s0, double v0, double sigma, double alpha, double beta);
};
#endif /* Stochastic_Process_Simulator_hpp */
