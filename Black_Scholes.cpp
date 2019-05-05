//
//  Black_Scholes.cpp
//  Computational_Methods_Final
//
//  Created by 黄佳铭 on 2019/3/14.
//  Copyright © 2019 Jiaming Huang. All rights reserved.
//

#include "Black_Scholes.hpp"
#include "Utils.hpp"
#include <iostream>
#include <cmath>

using namespace std;

double Black_Scholes::EuroCallPriceByBS(double S0, double r, double sigma, double T, double K){
    double d1 = (log(S0/K)+(r+0.5*sigma*sigma)*T)/(sigma*sqrt(T));
    double d2 = d1 - sigma*sqrt(T);
    double c = S0 * Utils::pnorm(d1)-K*exp(-1*r*T)*Utils::pnorm(d2);
    return c;
}

double Black_Scholes::EuroPutPriceByBS(double S0, double r, double sigma, double T, double K){
    double d1 = (log(S0/K)+(r+0.5*sigma*sigma)*T)/(sigma*sqrt(T));
    double d2 = d1 - sigma*sqrt(T);
    double p = K*exp(-1*r*T)*Utils::pnorm(-1*d2) - S0 * Utils::pnorm(-1*d1);
    return p;
}

double Black_Scholes::CallDeltaByBS(double S0, double r, double sigma, double T, double K){
    double d1 = (log(S0/K)+(r+0.5*sigma*sigma)*T)/(sigma*sqrt(T));
    return Utils::pnorm(d1);
}

double Black_Scholes::PutDeltaByBS(double S0, double r, double sigma, double T, double K){
    double d1 = (log(S0/K)+(r+0.5*sigma*sigma)*T)/(sigma*sqrt(T));
    return Utils::pnorm(d1)-1;
}

double Black_Scholes::CallGammaByBS(double S0, double r, double sigma, double T, double K){
    double d1 = (log(S0/K)+(r+0.5*sigma*sigma)*T)/(sigma*sqrt(T));
    return (Utils::pnorm(d1))/(S0*sigma*sqrt(T));
}

double Black_Scholes::PutGammaByBS(double S0, double r, double sigma, double T, double K){
    double d1 = (log(S0/K)+(r+0.5*sigma*sigma)*T)/(sigma*sqrt(T));
    return (Utils::pnorm(d1))/(S0*sigma*sqrt(T));
}

double Black_Scholes::CallThetaByBS(double S0, double r, double sigma, double T, double K){
    double d1 = (log(S0/K)+(r+0.5*sigma*sigma)*T)/(sigma*sqrt(T));
    double d2 = d1 - sigma*sqrt(T);
    return (-1*S0*sigma*Utils::pnorm(d1))/(2*sqrt(T)) - r*K*exp(-1*r*T)*Utils::pnorm(d2);
}

double Black_Scholes::PutThetaByBS(double S0, double r, double sigma, double T, double K){
    double d1 = (log(S0/K)+(r+0.5*sigma*sigma)*T)/(sigma*sqrt(T));
    double d2 = d1 - sigma*sqrt(T);
    return (-1*S0*sigma*Utils::pnorm(d1))/(2*sqrt(T)) + r*K*exp(-1*r*T)*Utils::pnorm(d2);
}

double Black_Scholes::CallVegaByBS(double S0, double r, double sigma, double T, double K){
    double d1 = (log(S0/K)+(r+0.5*sigma*sigma)*T)/(sigma*sqrt(T));
    return S0*sqrt(T)*Utils::pnorm(d1);
}

double Black_Scholes::PutVegaByBS(double S0, double r, double sigma, double T, double K){
    double d1 = (log(S0/K)+(r+0.5*sigma*sigma)*T)/(sigma*sqrt(T));
    return S0*sqrt(T)*Utils::pnorm(d1);
}

double Black_Scholes::CallRhoByBS(double S0, double r, double sigma, double T, double K){
    double d1 = (log(S0/K)+(r+0.5*sigma*sigma)*T)/(sigma*sqrt(T));
    double d2 = d1 - sigma*sqrt(T);
    return K*T*exp(-1*r*T)*Utils::pnorm(d2);
}

double Black_Scholes::PutRhoByBS(double S0, double r, double sigma, double T, double K){
    double d1 = (log(S0/K)+(r+0.5*sigma*sigma)*T)/(sigma*sqrt(T));
    double d2 = d1 - sigma*sqrt(T);
    return -1*K*T*exp(-1*r*T)*Utils::pnorm(-1*d2);
}


