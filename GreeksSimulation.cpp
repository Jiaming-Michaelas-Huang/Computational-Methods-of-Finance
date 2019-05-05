//
//  GreeksSimulation.cpp
//  Computational_Methods_Final
//
//  Created by 黄佳铭 on 2019/3/14.
//  Copyright © 2019 Jiaming Huang. All rights reserved.
//

#include "GreeksSimulation.hpp"
#include "OptionPricing.hpp"
#include <iostream>
#include <cmath>

using namespace std;

double GreeksSimulation::CallDelta(double T, double S0, double r, double mu, double sigma, double K, int size){
    return (OptionPricing::EuroCallOptionPricingByMonteCarlo(T, S0+0.1, r, mu, sigma, K, size) - OptionPricing::EuroCallOptionPricingByMonteCarlo(T, S0, r, mu, sigma, K, size))/0.1;
}

double GreeksSimulation::CallGamma(double T, double S0, double r, double mu, double sigma, double K, int size){
    return (OptionPricing::EuroCallOptionPricingByMonteCarlo(T, S0+0.2, r, mu, sigma, K, size) - 2*OptionPricing::EuroCallOptionPricingByMonteCarlo(T, S0+0.1, r, mu, sigma, K, size))/0.01;
}

double GreeksSimulation::CallTheta(double T, double S0, double r, double mu, double sigma, double K, int size){
    return (OptionPricing::EuroCallOptionPricingByMonteCarlo(T+0.1, S0, r, mu, sigma, K, size) - OptionPricing::EuroCallOptionPricingByMonteCarlo(T, S0, r, mu, sigma, K, size))/0.1;
}

double GreeksSimulation::CallVega(double T, double S0, double r, double mu, double sigma, double K, int size){
    return (OptionPricing::EuroCallOptionPricingByMonteCarlo(T, S0, r, mu, sigma+0.1, K, size) - OptionPricing::EuroCallOptionPricingByMonteCarlo(T, S0, r, mu, sigma, K, size))/0.1;
}

double GreeksSimulation::CallRho(double T, double S0, double r, double mu, double sigma, double K, int size){
    return (OptionPricing::EuroCallOptionPricingByMonteCarlo(T, S0, r+0.1, mu, sigma, K, size) - OptionPricing::EuroCallOptionPricingByMonteCarlo(T, S0, r, mu, sigma, K, size))/0.1;
}
