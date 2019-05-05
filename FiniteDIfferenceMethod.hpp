//
//  FiniteDIfferenceMethod.hpp
//  Computational_Methods_Final
//
//  Created by 黄佳铭 on 2019/3/21.
//  Copyright © 2019 Jiaming Huang. All rights reserved.
//

#ifndef FiniteDIfferenceMethod_hpp
#define FiniteDIfferenceMethod_hpp

#include <stdio.h>
#include <iostream>
#include "Matrix.hpp"

using namespace std;

class FiniteDifferenceMethod{
private:
    double S0;
    double r;
    double sigma;
    double delta_t;
    double delta_x;
    double delta_s;
    double maturity;
    double K;
    double paths;
    double N;
    string type;
public:
    FiniteDifferenceMethod(double S0, double r, double sigma, double delta_t,double delta_x,double delta_s, double maturity, double K,double paths, double N,string type);
    shared_ptr<Matrix> xsimulation();
    shared_ptr<Matrix> ssimulation();
    double ExplicitFiniteDifferenceMethodForX();
    double ImpliciteFiniteDifferenceMethodForX();
    double CrankNicolsonFiniteDifferenceMethodForX();
    double ExplicitFiniteDifferenceMethodForS();
    double ImpliciteFiniteDifferenceMethodForS();
    double CrankNicolsonFiniteDifferenceMethodForS();
};
#endif /* FiniteDIfferenceMethod_hpp */
