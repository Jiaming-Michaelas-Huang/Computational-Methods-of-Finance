//
//  main.cpp
//  Computational_Methods_Final
//
//  Created by 黄佳铭 on 2019/3/13.
//  Copyright © 2019 Jiaming Huang. All rights reserved.
//

#include <iostream>
#include <cmath>
#include "Random_Number_Generator.hpp"
#include "Utils.hpp"
#include "AssetBackedSecuritiesPricing.hpp"
#include "FixedIncomePricing.hpp"
#include "Least_Square_Monte_Carlo_Simulation.hpp"
#include "Stochastic_Process_Simulator.hpp"
#include "Monte_Carlo_Simulator.hpp"
#include "FixedIncomePricing.hpp"

using namespace std;

int main(int argc, const char * argv[]) {
   //Q1
    int k = 0;
    double * payoff = new double[10000];
    int seed = rand();
    for (int i = 0; i < 10000; i++) {
        double s = 0;
        double * u = Random_Number_Generator::Uniform_Distribution_Generator(1000, seed+i, "LGM");
        for (int j = 0; j < 1000; j++) {
            s += u[j];
            if(s>1.1){k = j+1;break;}
        }
        payoff[i] = 4.54-k>0?4.54-k:0;
    }
    cout<<"Ansewer to Q1 is:"<<Utils::mean(payoff, 10000)<<endl;
    
    //Q2
    double T = 2;
    double n = 100;
    double v0 = 0.06;
    double S0 = 20;
    double sigma = 0.25;
    double r = 0.05;
    double alpha = 0.45;
    double beta = -5.105;
    double rho[3] = {-0.75,0,0.75};
    double * result = new double[3];
    for (int rh = 0; rh < 3; rh++) {
        double rho1 = rho[rh];
        int seed = rand();
        double payoff = 0;
        for (int i = 0; i < 10000; i++) {
            double * w1 = Monte_Carlo_Simulator::BrownieMotionSimulation(T/n, seed+i, n);
            double * w2 = Monte_Carlo_Simulator::BrownieMotionSimulation(T/n, seed+30000+i, n);
            double ** w = Random_Number_Generator::Bivariate_Normal_Distribution_Generator(n, w1, w2, rho1);
            double * S = Stochastic_Process_Simulation::HestonModelFullTruncationSchema(w[0], w[1], T, n, r, S0, v0, sigma, alpha, beta);
            payoff+=(S[(int)n]-Utils::mean(S, n+1)>0)?(S[(int)n]-Utils::mean(S, n+1)):0;
        }
        result[rh] = (payoff/10000)*exp(-r*T);
    }
    
    cout<<"When rho is -0.75, price is: "<<result[0]<<endl;
    cout<<"When rho is 0, price is: "<<result[1]<<endl;
    cout<<"When rho is 0.75, price is: "<<result[2]<<endl;
    
    //Q3
    S0 = 100;
    r = 0.05;
    sigma = 0.35;
    T = 5;
    double K = 100;
    n = 100;
    double * Lt = new double[n+1];
    Lt[0] = 50;
    double * Ut = new double[n+1];
    Ut[0] = 150;
    for (int i = 1; i < n+1; i++) {
        double time = i*T/n;
        Lt[i] = 50*exp(0.138629*time);
        Ut[i] = 200-Lt[i];
    }
    double *payoffs = new double[10000];
    seed = rand();
    int num_exercise = 0;
    int Lexercise = 0;
    for (int i = 0; i < 10000; i++) {
        double * w = Monte_Carlo_Simulator::BrownieMotionSimulation(T/n,seed+i ,n);
        double * S3 = Stochastic_Process_Simulation::StockPathSimulation(S0, T, n, r, sigma, w);
        int L1U2 = 0;
        double exetime = 0;
        double payoff = 0;
        for (int j = 0; j < n; j++) {
            if (S3[j]<=Lt[j]) {
                L1U2 = 1;
                exetime = j;
                break;
            }
            if (S3[j]>=Ut[j]) {
                L1U2 = 2;
                exetime = j;
                break;
            }
        }
        if(L1U2 == 1){payoff = K-S3[(int)exetime];num_exercise++;Lexercise++;}
        if(L1U2 == 2){payoff = S3[(int)exetime] - K;num_exercise++;}
        payoffs[i] = payoff*exp(-r*T);
    }
    cout<<"The price is: "<<Utils::mean(payoffs, 10000)<<endl;
    cout<<"The probability is: "<<((double)Lexercise)/((double)num_exercise)<<endl;
    //Q4
    double gamma = 2;
    alpha = 0.36;
    beta = -5.86;
    sigma = 0.36;
    double r0 = 0.05;
    K = 9800;
    double FaceValue = 10000;
    double Ti = 1;
    T = 0.5;
    double delta_t = Ti/n;
    FixedIncomePricing * fip = new FixedIncomePricing(r0,0,sigma,0);
    double price4 = fip->EuroPutOptionOnPureDiscountedBondPricing(FaceValue, T, delta_t, 10000, Ti, K, "Final");
    cout<<"The price is: "<<price4<<endl;
    
    //Q5
    S0 = 6000;
    double E0 = 0.0096;
    double rho5 = -0.25;
    r = 0.05;
    double rf = 0.04;
    double sigma1 = 0.1;
    double sigma2 = 0.15;
    gamma = -0.04;
    double lambda = 1.5;
    K  = 60;
    T = 1;
    seed = rand();
    double payoff5 = 0;
    for (int i = 0; i < 10000; i++) {
        double * w1 = Monte_Carlo_Simulator::BrownieMotionSimulation(T/n, seed+i, n);
        double * w2 = Monte_Carlo_Simulator::BrownieMotionSimulation(T/n, seed + 30000+i, n);
        double ** w = Random_Number_Generator::Bivariate_Normal_Distribution_Generator(n, w1, w2, rho5);
        double * Jt = Random_Number_Generator::Poisson_Distribution_Generator(n+1, seed+i, "LGM", lambda);
        double *S5 = Stochastic_Process_Simulation::StockPathSimulationWithJump(S0, T, n, r, sigma1, gamma, w[0], Jt);
        double *E = Stochastic_Process_Simulation::StockPathSimulation(E0, T, n, r-rf, sigma2, w[1]);
        payoff5 += ((S5[(int)n]*E[(int)n]-K>0)?(S5[(int)n]*E[(int)n]-K):0)*exp(-r*T);
    }
    cout<<"The price is: "<<payoff5/10000<<endl;
    
    return 0;
    
}
