//
//  Monte_Carlo_Simulator.hpp
//  Computational_Methods_Final
//
//  Created by 黄佳铭 on 2019/3/13.
//  Copyright © 2019 Jiaming Huang. All rights reserved.
//

#ifndef Monte_Carlo_Simulator_hpp
#define Monte_Carlo_Simulator_hpp

#include <stdio.h>
#include <cmath>

using namespace std;

class Monte_Carlo_Simulator{
public:
    static double * BrownieMotionSimulation(double T,int seed, int size);
    static double * StockPriceSimulation(double * brownie ,double S0, double T, double mu, double sigma, int size);
};
#endif /* Monte_Carlo_Simulator_hpp */
