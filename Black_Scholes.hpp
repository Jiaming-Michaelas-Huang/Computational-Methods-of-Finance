//
//  Black_Scholes.hpp
//  Computational_Methods_Final
//
//  Created by 黄佳铭 on 2019/3/14.
//  Copyright © 2019 Jiaming Huang. All rights reserved.
//

#ifndef Black_Scholes_hpp
#define Black_Scholes_hpp

#include <stdio.h>
class Black_Scholes{
public:
    static double EuroCallPriceByBS(double S0, double r, double sigma, double T, double K);
    static double EuroPutPriceByBS(double S0, double r, double sigma, double T, double K);
    static double CallDeltaByBS(double S0, double r, double sigma, double T, double K);
    static double CallGammaByBS(double S0, double r, double sigma, double T, double K);
    static double CallThetaByBS(double S0, double r, double sigma, double T, double K);
    static double CallVegaByBS(double S0, double r, double sigma, double T, double K);
    static double CallRhoByBS(double S0, double r, double sigma, double T, double K);
    static double PutDeltaByBS(double S0, double r, double sigma, double T, double K);
    static double PutGammaByBS(double S0, double r, double sigma, double T, double K);
    static double PutThetaByBS(double S0, double r, double sigma, double T, double K);
    static double PutVegaByBS(double S0, double r, double sigma, double T, double K);
    static double PutRhoByBS(double S0, double r, double sigma, double T, double K);
};
#endif /* Black_Scholes_hpp */
