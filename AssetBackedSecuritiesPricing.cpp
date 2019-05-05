//
//  AssetBackedSecuritiesPricing.cpp
//  Computational_Methods_Final
//
//  Created by 黄佳铭 on 2019/3/21.
//  Copyright © 2019 Jiaming Huang. All rights reserved.
//

#include "AssetBackedSecuritiesPricing.hpp"
#include "Random_Number_Generator.hpp"
#include "FixedIncomePricing.hpp"
#include "Matrix.hpp"
#include "Utils.hpp"
#include <iostream>
#include <cmath>

using namespace std;

AssetBackedSecuritiesPricing::AssetBackedSecuritiesPricing(double PrincipleValue, double maturity, double WAC, double r0, double r_bar, double kai, double sigma,double delta_t, int paths){
    this->PrincipleValue = PrincipleValue;
    this->maturity = maturity;
    this->WAC = WAC;
    this->r0 = r0;
    this->r_bar = r_bar;
    this->kai = kai;
    this->sigma = sigma;
    this->delta_t = delta_t;
    this->paths = paths;
}

double AssetBackedSecuritiesPricing::NumerixPrepaymentModel(double R, double r_tm1_10, double t, double PV_tm1, double PV0){
    double SY[12] = {0.94,0.76,0.74,0.95,0.98,0.92,0.98,1.10,1.18,1.22,1.23,0.98};
    double RIt = 0.28+0.14*atan(-8.57+430*(R-r_tm1_10));
    double BUt = 0.3+0.7*PV_tm1/PV0;
    double SGt = 1>t/30.0?t/30.0:1;
    int SYt_index = ((int)(t)%12)-1;
    double SYt = SY[SYt_index];
    return RIt*BUt*SGt*SYt;
}

double AssetBackedSecuritiesPricing::PrepaymentPricing(){
    int seed = rand();
    double PV0 = PrincipleValue;
    double PVt = PrincipleValue;
    double P0 = 0;
    FixedIncomePricing * fip = new FixedIncomePricing(r0,r_bar,sigma,kai);
    shared_ptr<Matrix> discount_rate = fip->CIR_Model(seed, maturity+10, delta_t, paths);
    for (int i = 0; i < paths; i++) {
        double R = 0;
        PVt = PrincipleValue;
        for (int j = 0; j < maturity/delta_t; j++) {
            R += delta_t*discount_rate->elementAt(i, j);
            int t = j + 1;
            if (t%30 == 0) {
                int month = t/30;
                int T = 10*360;
                double R_temp = 0;
                double wac = WAC/12;
                for(int k = 0; k < T; k++)
                {
                    double dis = 0;
                    for (int k1 = 0; k1<paths; k1++) {
                        dis += discount_rate->elementAt(i, j+k);
                    }
                    R_temp += delta_t*(dis/paths);
                }
                double r_tm1_10 = (1/((double)120))*R_temp;
                double CRPt = NumerixPrepaymentModel(WAC, r_tm1_10, month, PVt, PV0);
                //double ct = PVt*wac/(1-pow(1+wac, month-1-maturity*12))+(PVt-PVt*wac*(1/(1-pow(1+wac, month-1-maturity*12))-1))*(1-pow(1-CRPt, 1.0/12.0));
                double TPPt = PVt*wac*(1/(1-pow(1+wac, month-1-maturity*12))-1)+(PVt-PVt*wac*(1/(1-pow(1+wac, month-1-maturity*12))-1))*(1-pow(1-CRPt, 1.0/12.0));
                PVt = PVt - TPPt;
                double ct = TPPt+PVt*wac;
                if (t == 10800) {
                    int baba = 0;
                }
                P0 += ct*exp(-R);
            }
        }
    }
    return P0/paths;
}

double AssetBackedSecuritiesPricing::OASPrepaymentPricing(double deviation){
    int seed = rand();
    double PV0 = PrincipleValue;
    double PVt = PrincipleValue;
    double P0 = 0;
    FixedIncomePricing * fip = new FixedIncomePricing(r0,r_bar,sigma,kai);
    shared_ptr<Matrix> discount_rate = fip->CIR_Model(seed, maturity+10+delta_t, delta_t, paths);
    double *R10 = new double [(int)(maturity/delta_t)+1];
    for (int i = 0; i < maturity/delta_t+1; i++) {
        double R_temp=0;
        for(int k = 0; k < 120; k++)
        {
            double dis = 0;
            for (int k1 = 0; k1<paths; k1++) {
                dis += discount_rate->elementAt(k1, i+k);
            }
            R_temp += delta_t*(dis/paths);
        }
        R10[i] = R_temp/120;
    }
    for (int i = 0; i < paths; i++) {
        double R = 0;
        PVt = PrincipleValue;
        for (int j = 1; j < maturity/delta_t+1; j++) {
            R += delta_t*(discount_rate->elementAt(i, j)+deviation);
            int t = j;
            int month = t;
            int T = 120;
            double wac = WAC/12;
            double r_tm1_10 = (-1/((double)120))*R10[j];
            double CRPt = NumerixPrepaymentModel(WAC, r_tm1_10, month, PVt, PV0);
                //double ct = PVt*wac/(1-pow(1+wac, month-1-maturity*12))+(PVt-PVt*wac*(1/(1-pow(1+wac, month-1-maturity*12))-1))*(1-pow(1-CRPt, 1.0/12.0));
            double TPPt = PVt*wac*(1.0/(1.0-pow(1.0+wac, month-1.0-maturity*12))-1.0)+(PVt-PVt*wac*(1.0/(1.0-pow(1.0+wac, month-1.0-maturity*12))-1.0))*(1-pow(1-CRPt, 1.0/12.0));
            PVt = PVt - TPPt;
            double ct = TPPt+PVt*wac;
            P0 += ct*exp(-R);
        }
    }
    return P0/paths;
}
