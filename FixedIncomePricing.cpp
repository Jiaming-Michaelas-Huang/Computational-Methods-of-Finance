//
//  FixedIncomePricing.cpp
//  Computational_Methods_Final
//
//  Created by 黄佳铭 on 2019/3/21.
//  Copyright © 2019 Jiaming Huang. All rights reserved.
//

#include "FixedIncomePricing.hpp"
#include <iostream>
#include <cmath>
#include <random>
#include "Matrix.hpp"
#include "Utils.hpp"
#include "Random_Number_Generator.hpp"
#include "Monte_Carlo_Simulator.hpp"

using namespace std;

FixedIncomePricing::FixedIncomePricing(double r0, double r_bar, double sigma, double kai){
    this->r0 = r0;
    this->r_bar = r_bar;
    this->sigma = sigma;
    this->kai = kai;
}

shared_ptr<Matrix> FixedIncomePricing::Vasicek_Model(int seed, double T, double delta_t, int paths){
    double ** result = new double*[paths];
    for (int i = 0; i < paths; i++) {
        double * z1 = Random_Number_Generator::Normal_Distribution_Generator_Box_Muller(T/delta_t,seed+i, "LGM");
        result[i] = new double[(int)(T/delta_t)];
        double rt = r0;
        result[i][0] = rt;
        for (int j = 1; j < T/delta_t; j++) {
            rt += kai*(r_bar - rt)*delta_t + sigma*sqrt(delta_t)*z1[j];
            //rt = rt>0?rt:0;
            result[i][j] = rt;
        }
        delete []z1;
    }
    shared_ptr<Matrix> re = make_shared<Matrix>(result,paths,T/delta_t);
    for (int i = 0; i < paths; i++) {
        delete [] result[i];
    }
    delete []result;
    return re;
}

shared_ptr<Matrix> FixedIncomePricing::Final_Model(int seed, double T, double delta_t, int paths, double alpha, double beta, double gamma){
    double ** result = new double*[paths];
    for (int i = 0; i < paths; i++) {
        double * z1 = Random_Number_Generator::Normal_Distribution_Generator_Box_Muller(T/delta_t,seed+i, "LGM");
        result[i] = new double[(int)(T/delta_t)];
        double rt = r0;
        result[i][0] = rt;
        for (int j = 1; j < T/delta_t; j++) {
            rt += (alpha+beta*rt)*delta_t + sigma*pow(rt, gamma)*sqrt(delta_t)*z1[j];
            rt = rt>0?rt:0;
            result[i][j] = rt;
        }
        delete []z1;
    }
    shared_ptr<Matrix> re = make_shared<Matrix>(result,paths,T/delta_t);
    for (int i = 0; i < paths; i++) {
        delete [] result[i];
    }
    delete []result;
    return re;
}

shared_ptr<Matrix> FixedIncomePricing::CIR_Model(int seed, double T, double delta_t, int paths){
    double ** result = new double*[paths];
    for (int i = 0; i < paths; i++) {
        double * z1 = Random_Number_Generator::Normal_Distribution_Generator_Box_Muller(T/delta_t,seed+i, "LGM");
        result[i] = new double[(int)(T/delta_t)];
        double rt = r0;
        result[i][0] = r0;
        for (int j = 1; j < T/delta_t; j++) {
            rt += kai*(r_bar - rt)*delta_t + sigma*sqrt(abs(rt))*sqrt(delta_t)*z1[j];
            result[i][j] = rt;
        }
        delete []z1;
    }
    shared_ptr<Matrix> re = make_shared<Matrix>(result,paths,T/delta_t);
    for (int i = 0; i < paths; i++) {
        delete [] result[i];
    }
    delete []result;
    return re;
}

shared_ptr<Matrix>FixedIncomePricing::Gpp_Model(int seed, double T, double delta_t, int paths,double roh,double a, double b,double sig, double nah, double phi){
    double ** result = new double*[paths];
    for (int i = 0; i < paths; i++) {
        double * w1 = Monte_Carlo_Simulator::BrownieMotionSimulation(delta_t, 19951203+i, T/delta_t);
        double * w2 = Monte_Carlo_Simulator::BrownieMotionSimulation(delta_t, 19961219+i, T/delta_t);
        double ** w = Random_Number_Generator::Bivariate_Normal_Distribution_Generator(T/delta_t, w1, w2, roh);
        result[i] = new double[(int)(T/delta_t)];
        double rt = 0;
        double xt = 0;
        double yt = 0;
        for (int j = 0; j < T/delta_t; j++) {
            double dxt = -a*xt*delta_t+sig*w[0][j];
            double dyt = -b*yt*delta_t+nah*w[1][j];
            xt = xt+dxt;
            yt = yt+dyt;
            rt = phi+xt+yt;
            result[i][j] = rt;
        }
        delete [] w[0];
        delete [] w[1];
        delete [] w;
    }
    shared_ptr<Matrix> re = make_shared<Matrix>(result,paths,T/delta_t);
    for (int i = 0; i < paths; i++) {
        delete [] result[i];
    }
    delete []result;
    return re;
}

double FixedIncomePricing::PureDiscountBondPricing(double FaceValue, double T, double delta_t, int paths,string type){
    int seed = rand();
    shared_ptr<Matrix> simulated_r = NULL;
    if(type == "Vasicek") simulated_r = Vasicek_Model(seed,T,delta_t,paths);
    if(type == "CIR") simulated_r = CIR_Model(seed,T,delta_t,paths);
    if(type == "Final") simulated_r = Final_Model(seed,T,delta_t,paths,0.36,-5.86,2);
    if(type == "GPP")simulated_r = Gpp_Model(seed, T, delta_t, paths, 0.7, 0.1, 0.3, 0.03, 0.08, 0.03);
    double *price = new double[paths];
    for (int i = 0; i < paths; i++) {
        double R = 0;
        for (int j = 0; j < T/delta_t; j++) {
            R += delta_t* simulated_r->elementAt(i, j);
        }
        price[i] = FaceValue * exp(-1*R);
    }
    return Utils::mean(price, paths);
}

double FixedIncomePricing::CoponPayingBondPricing(double FaceValue, double frequency, double Coupon, double T, double delta_t, int paths,string type){
    int seed = rand();
    shared_ptr<Matrix> simulated_r = NULL;
    if(type == "Vasicek") simulated_r = Vasicek_Model(seed,T,delta_t,paths);
    if(type == "CIR") simulated_r = CIR_Model(seed,T,delta_t,paths);
    if(type == "Final") simulated_r = Final_Model(seed,T,delta_t,paths,0.36,-5.86,2);
    if(type == "GPP")simulated_r = Gpp_Model(seed, T, delta_t, paths, 0.7, 0.1, 0.3, 0.03, 0.08, 0.03);
    double *prices = new double[paths];
    for (int i = 0; i < paths; i++) {
        double price = 0;
        double R = 0;
        for (int j = 0; j < T/delta_t; j++) {
            R += delta_t* simulated_r->elementAt(i, j);
            if(((int)(1000*(j+1)*delta_t))%((int)(1000*frequency)) == 0){
                if(j!=0)price += Coupon*exp(-1*R);
            }
        }
        price += (FaceValue)*exp(-1*R);
        prices[i] = price;
    }
    return Utils::mean(prices, paths);
}

double FixedIncomePricing::EuroCallOptionOnPureDiscountedBondPricing(double FaceValue, double T, double delta_t, int paths,double Ti,double K, string type){
    int seed = rand();
    shared_ptr<Matrix> simulated_r = NULL;
    if(type == "Vasicek") simulated_r = Vasicek_Model(seed,Ti,delta_t,paths);
    if(type == "CIR") simulated_r = CIR_Model(seed,Ti,delta_t,paths);
    if(type == "Final") simulated_r = Final_Model(seed,T,delta_t,paths,0.36,-5.86,2);
    if(type == "GPP")simulated_r = Gpp_Model(seed, Ti, delta_t, paths, 0.7, 0.1, 0.3, 0.03, 0.08, 0.03);
    double *callprice = new double[paths];
    for (int i = 0; i < paths; i++) {
        double R = 0;
        double bondprice = 0;
        for (int j = T/delta_t; j < (Ti)/delta_t; j++) {
            R += delta_t* simulated_r->elementAt(i, j);
        }
        bondprice = FaceValue * exp(-1*R);
        R = 0;
        for (int j = 0; j < T/delta_t; j++) {
            R += delta_t* simulated_r->elementAt(i, j);
        }
        callprice[i] = (bondprice-K>0?bondprice-K:0)*exp(-1*R);
        cout<<callprice[i]<<endl;
    }
    return Utils::mean(callprice, paths);
}

double FixedIncomePricing::EuroPutOptionOnPureDiscountedBondPricing(double FaceValue, double T, double delta_t, int paths, double Ti, double K, string type){
    int seed = rand();
    shared_ptr<Matrix> simulated_r = NULL;
    if(type == "Vasicek") simulated_r = Vasicek_Model(seed,Ti,delta_t,paths);
    if(type == "CIR") simulated_r = CIR_Model(seed,Ti,delta_t,paths);
    if(type == "Final") simulated_r = Final_Model(seed,T,delta_t,paths,0.36,-5.86,2);
    if(type == "GPP")simulated_r = Gpp_Model(seed, Ti, delta_t, paths, 0.7, 0.1, 0.3, 0.03, 0.08, 0.03);
    double *putprice = new double[paths];
    for (int i = 0; i < paths; i++) {
        double R = 0;
        double bondprice = 0;
        for (int j = T/delta_t; j < (Ti)/delta_t; j++) {
            R += delta_t* simulated_r->elementAt(i, j);
        }
        bondprice = FaceValue * exp(-1*R);
        R = 0;
        for (int j = 0; j < T/delta_t; j++) {
            R += delta_t* simulated_r->elementAt(i, j);
        }
        putprice[i] = (K-bondprice>0?K-bondprice:0)*exp(-1*R);
        //cout<<putprice[i]<<endl;
    }
    return Utils::mean(putprice, paths);
}

double FixedIncomePricing::EuroCallOptionOnCouponPayingBondPricing(double FaceValue, double frequency, double Coupon, double T, double delta_t, int paths, double Ti, double K,string type){
    int seed = rand();
    shared_ptr<Matrix> simulated_r = NULL;
    if(type == "Vasicek") simulated_r = Vasicek_Model(seed,Ti,delta_t,paths);
    if(type == "CIR") simulated_r = CIR_Model(seed,Ti,delta_t,paths);
    if(type == "Final") simulated_r = Final_Model(seed,T,delta_t,paths,0.36,-5.86,2);
    if(type == "GPP")simulated_r = Gpp_Model(seed, Ti, delta_t, paths, 0.7, 0.1, 0.3, 0.03, 0.08, 0.03);
    double *callprice = new double[paths];
    for (int i = 0; i < paths; i++) {
        double bondprice = 0;
        double R = 0;
        for (int j = T/delta_t; j < Ti/delta_t; j++) {
            R += delta_t* simulated_r->elementAt(i, j);
            if(((int)(1000*(j+1)*delta_t))%((int)(1000*frequency)) == 0){
                bondprice += Coupon*exp(-1*R);
            }
        }
        bondprice += (FaceValue)*exp(-1*R);
        callprice[i] = (bondprice-K>0?bondprice-K:0)*exp(-1*T*simulated_r->elementAt(i, T/delta_t-1));
    }
    return Utils::mean(callprice, paths);
}

double FixedIncomePricing::EuroCallOptionByVasicekFormula(double FaceValue, double T,double Ti,double K, double delta_t, int paths){
    int seed = rand();
    double Bt = 1/kai*(1-exp(-kai*(Ti-T)));
    double At = exp((r_bar-sigma*sigma/(2*kai*kai))*(Bt-Ti+T)-sigma*sigma*Bt*Bt/(4*kai));
    shared_ptr<Matrix> rt = Vasicek_Model(seed, T+delta_t, delta_t, paths);
    double payoffs = 0;
    double price = 0;
    for (int i = 0; i < paths; i++) {
        double Pt = FaceValue*At*exp(-Bt*rt->elementAt(i, (int)(T/delta_t)));
        payoffs = Pt-K>0?Pt-K:0;
        double R = 0;
        for (int j = 0; j < T/delta_t; j++) {
            R += delta_t* rt->elementAt(i, j);
        }
        price += payoffs*exp(-1*R);
    }
    return price/paths;
}

double FixedIncomePricing::EuroCallOptionByCIRFormula(double FaceValue, double T, double Ti, double K){
    double h1 = sqrt(kai*kai+2*sigma*sigma);
    double h2 = (h1+kai)/2;
    double h3 = 2*kai*r_bar/(sigma*sigma);
    double A0s = pow((h1*exp(h2*(Ti-0))/(h1+h2*(exp(h1*(Ti-0))-1))), h3);
    double B0s = (exp(h1*(Ti-0))-1)/(h1+h2*(exp(h1*(Ti-0))-1));
    double P0s = A0s*exp(-1*B0s*r0);
    double phi = sqrt(kai*kai+2*sigma*sigma);
    double tau = 2*phi/(sigma*sigma*(exp(phi*(T-0))-1));
    double dao = (kai+phi)/(sigma*sigma);
    double Ats = pow((h1*exp(h2*(Ti-T))/(h1+h2*(exp(h1*(Ti-T))-1))), h3);
    double Bts = (exp(h1*(Ti-T))-1)/(h1+h2*(exp(h1*(Ti-T))-1));
    double r_star = log(Ats/K)/Bts;
    double chisquared1 = 0.2893;
    double chisquared2 = 0.2864;
    double A0t = pow((h1*exp(h2*(T-0))/(h1+h2*(exp(h1*(T-0))-1))), h3);
    double B0t = (exp(h1*(T-0))-1)/(h1+h2*(exp(h1*(T-0))-1));
    double P0t = A0t*exp(-1*B0t*r0);
    return FaceValue*P0s*chisquared1-K*P0t*chisquared2;
}


