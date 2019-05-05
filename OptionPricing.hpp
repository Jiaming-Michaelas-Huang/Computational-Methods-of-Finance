//
//  OptionPricing.hpp
//  Computational_Methods_Final
//
//  Created by 黄佳铭 on 2019/3/13.
//  Copyright © 2019 Jiaming Huang. All rights reserved.
//

#ifndef OptionPricing_hpp
#define OptionPricing_hpp

#include <stdio.h>
#include <cmath>
class OptionPricing{
public:
    static double EuroCallOptionPricingByMonteCarlo(double T, double S0, double r,double mu, double sigma, double K, int size);
    static double EuroCallOptionPricingByStochasticProcessSimulation(double T, double S0, double r,double mu, double sigma, double K, int size,int n);
};
#endif /* OptionPricing_hpp */
