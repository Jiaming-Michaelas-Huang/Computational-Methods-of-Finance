//
//  GreeksSimulation.hpp
//  Computational_Methods_Final
//
//  Created by 黄佳铭 on 2019/3/14.
//  Copyright © 2019 Jiaming Huang. All rights reserved.
//

#ifndef GreeksSimulation_hpp
#define GreeksSimulation_hpp

#include <stdio.h>
class GreeksSimulation{
public:
    static double CallDelta(double T, double S0, double r, double mu, double sigma, double K, int size);
    static double PutDelta(double T, double S0, double r, double mu, double sigma, double K, int size);
    static double CallGamma(double T, double S0, double r, double mu, double sigma, double K, int size);
    static double PutGamma(double T, double S0, double r, double mu, double sigma, double K, int size);
    static double CallTheta(double T, double S0, double r, double mu, double sigma, double K, int size);
    static double PutTheta(double T, double S0, double r, double mu, double sigma, double K, int size);
    static double CallVega(double T, double S0, double r, double mu, double sigma, double K, int size);
    static double PutVega(double T, double S0, double r, double mu, double sigma, double K, int size);
    static double CallRho(double T, double S0, double r, double mu, double sigma, double K, int size);
    static double PutRho(double T, double S0, double r, double mu, double sigma, double K, int size);
};
#endif /* GreeksSimulation_hpp */
