//
//  Monte_Carlo_Simulator.cpp
//  Computational_Methods_Final
//
//  Created by 黄佳铭 on 2019/3/13.
//  Copyright © 2019 Jiaming Huang. All rights reserved.
//

#include "Monte_Carlo_Simulator.hpp"
#include "Random_Number_Generator.hpp"
#include "Utils.hpp"
#include <iostream>
#include <cmath>

using namespace std;

double * Monte_Carlo_Simulator::BrownieMotionSimulation(double T, int seed, int size){
    double * normal_samples = Random_Number_Generator::Normal_Distribution_Generator_Box_Muller(size, seed, "LGM");
    double * brownie_simulation = new double[size];
    for (int i = 0; i < size; i++) {
        brownie_simulation[i] = sqrt(T) * normal_samples[i];
    }
    delete [] normal_samples;
    return brownie_simulation;
}

double * Monte_Carlo_Simulator::StockPriceSimulation(double * brownie ,double S0 ,double T, double mu, double sigma, int size){
    double * S = new double[size];
    for (int i = 0; i < size; i++) {
        S[i] = S0*exp((mu-0.5*sigma*sigma)*T+sigma*brownie[i]);
    }
    delete []brownie;
    return S;
}


