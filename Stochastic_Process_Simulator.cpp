//
//  Stochastic_Process_Simulator.cpp
//  Computational_Methods_Final
//
//  Created by 黄佳铭 on 2019/3/14.
//  Copyright © 2019 Jiaming Huang. All rights reserved.
//

#include "Stochastic_Process_Simulator.hpp"
#include "Monte_Carlo_Simulator.hpp"
#include <cmath>
#include <iostream>

using namespace std;

double * Stochastic_Process_Simulation::WinierPathSimulation(double T, int n, int seed){
    double * dw = Monte_Carlo_Simulator::BrownieMotionSimulation(T/n, seed, n);
    return dw;
}

double * Stochastic_Process_Simulation::WinierProcessSimulation(double T, int n, int size, int seed){
    double * winierprocess = new double[size];
    for (int i = 0; i < size; i++) {
        double * winierpath = WinierPathSimulation(T, n, seed+i);
        winierprocess[i] = winierpath[n-1];
        delete [] winierpath;
    }
    return winierprocess;
}

double * Stochastic_Process_Simulation::StockPathSimulation(double S0, double T, int n, double mu, double sigma, double *winierpath){
    double * stockpath = new double[n+1];
    double st = S0;
    stockpath[0] = S0;
    for (int i = 1; i < n+1; i++) {
        stockpath[i] = st + mu*st*T/n + sigma*st*winierpath[i-1];
        //stockpath[i] = st*exp((mu-0.5*sigma*sigma)*T/n+sigma*winierpath[i-1]);
        st = stockpath[i];
    }
    delete [] winierpath;
    return stockpath;
}

double * Stochastic_Process_Simulation::StockPathSimulationWithJump(double S0, double T, int n, double mu, double sigma, double gamma, double *winierpath, double *poisson){
    double * stockpath = new double[n+1];
    double st = S0;
    stockpath[0] = S0;
    for (int i = 1; i < n+1; i++) {
        double delta_j = poisson[i]-poisson[i-1];
        stockpath[i] = st + mu*st*T/n + sigma*st*winierpath[i-1] + gamma*st*(delta_j);
        //stockpath[i] = st*exp((mu-0.5*sigma*sigma)*T/n+sigma*winierpath[i-1]);
        st = stockpath[i];
    }
    delete [] winierpath;
    return stockpath;
}

double * Stochastic_Process_Simulation::StochasticPathSimulationEulerSchema(double x0, double T, int n,double * winier,double (*a)(double), double (*b)(double)){
    double * result = new double[n];
    double lastresult = x0;
    for (int i = 0; i < n; i++) {
        result[i] = lastresult + a(lastresult)*T/n + b(lastresult)*winier[i];
        lastresult = result[i];
    }
    delete [] winier;
    return result;
}

double * Stochastic_Process_Simulation::StochasticPathSimulationMilshteinSchema(double x0, double T, int n,double * winier,double (*a)(double), double (*b)(double)){
    double * result = new double[n];
    double lastresult = x0;
    for (int i = 0; i < n; i++) {
        result[i] = lastresult + a(lastresult)*T/n + b(lastresult)*winier[i] + 0.5*b(lastresult)*((b(lastresult+0.01)-b(lastresult))/0.01)*(winier[i]*winier[i]-T/n);
        lastresult = result[i];
    }
    delete [] winier;
    return result;
}

double * Stochastic_Process_Simulation::HestonModelReflectionSchema(double *w1, double *w2, double T, int n, double r, double s0, double v0, double sigma, double alpha, double beta){
    double delta = T/n;
    double *v_reflect = new double[n+1];
    double *s_reflect = new double[n+1];
    v_reflect[0] = v0;
    s_reflect[0] = s0;
    double last_s = s0;
    double last_v = v0;
    for(int i = 1; i < n+1; i++){
        s_reflect[i] = last_s + r * last_s * delta + sqrt(abs(last_v)) * last_s * w1[i-1];
        v_reflect[i] = abs(last_v) + alpha * (beta * abs(last_v)) * delta + sigma * sqrt(abs(last_v))* w2[i-1];
        last_s = s_reflect[i];
        last_v = v_reflect[i];
    }
    delete [] v_reflect;
    return s_reflect;
}

double * Stochastic_Process_Simulation::HestonModelPartialTruncationSchema(double *w1, double *w2, double T, int n, double r, double s0, double v0, double sigma, double alpha, double beta){
    double delta = T/n;
    double *v_reflect = new double[n];
    double *s_reflect = new double[n];
    double last_s = s0;
    double last_v = v0;
    for(int i = 0; i < n; i++){
        s_reflect[i] = last_s + r * last_s * delta + sqrt(max(0.0, last_v)) * last_s * w1[i];
        v_reflect[i] = last_v + alpha * (beta - last_v) * delta + sigma * sqrt(max(0.0, last_v))* w2[i];
        last_s = s_reflect[i];
        last_v = v_reflect[i];
    }
    delete [] v_reflect;
    return s_reflect;
}

double * Stochastic_Process_Simulation::HestonModelFullTruncationSchema(double *w1, double *w2, double T, int n, double r, double s0, double v0, double sigma, double alpha, double beta){
    double delta = T/n;
    double *v_reflect = new double[n+1];
    double *s_reflect = new double[n+1];
    v_reflect[0] = v0;
    s_reflect[0] = s0;
    double last_s = s0;
    double last_v = v0;
    for(int i = 1; i < n+1; i++){
        s_reflect[i] = last_s + r * last_s * delta + sqrt(max(0.0, last_v)) * last_s * w1[i-1];
        v_reflect[i] = last_v + (alpha + (beta * max(0.0, last_v))) * delta + sigma * sqrt(max(0.0, last_v))* w2[i-1];
        last_s = s_reflect[i];
        last_v = v_reflect[i];
    }
    delete [] v_reflect;
    return s_reflect;
}


