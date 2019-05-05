//
//  OptionPricing.cpp
//  Computational_Methods_Final
//
//  Created by 黄佳铭 on 2019/3/13.
//  Copyright © 2019 Jiaming Huang. All rights reserved.
//

#include "OptionPricing.hpp"
#include "Monte_Carlo_Simulator.hpp"
#include "Stochastic_Process_Simulator.hpp"
#include "Utils.hpp"
#include <iostream>
#include <random>
#include <cmath>


using namespace std;

double OptionPricing::EuroCallOptionPricingByMonteCarlo(double T, double S0, double r, double mu, double sigma, double K, int size){
    int seed = rand();
    double * brownie = Monte_Carlo_Simulator::BrownieMotionSimulation(T, seed, size);
    double * stock_simulation = Monte_Carlo_Simulator::StockPriceSimulation(brownie ,S0, T, mu, sigma,size);
    double * price = new double[size];
    for (int i = 0; i < size; i++) {
        double payoff = stock_simulation[i]-K>0?stock_simulation[i]-K:0;
        price[i] = payoff*exp(-1*r*T);
    }
    double p = Utils::mean(price, size);
    delete [] stock_simulation;
    delete [] price;
    return p;
}

double OptionPricing::EuroCallOptionPricingByStochasticProcessSimulation(double T, double S0, double r, double mu, double sigma, double K, int size, int n){
    int seed = rand();
    double * winier = Stochastic_Process_Simulation::WinierProcessSimulation(T, n, size, seed);
    double * stock_simulation = Monte_Carlo_Simulator::StockPriceSimulation(winier ,S0, T, mu, sigma,size);
    double * price = new double[size];
    for (int i = 0; i < size; i++) {
        double payoff = stock_simulation[i]-K>0?stock_simulation[i]-K:0;
        price[i] = payoff*exp(-1*r*T);
    }
    double p = Utils::mean(price, size);
    delete [] stock_simulation;
    delete [] price;
    return p;
}
