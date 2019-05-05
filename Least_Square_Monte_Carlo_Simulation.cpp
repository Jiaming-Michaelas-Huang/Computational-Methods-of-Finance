//
//  Least_Square_Monte_Carlo_Simulation.cpp
//  Computational_Methods_Final
//
//  Created by 黄佳铭 on 2019/3/14.
//  Copyright © 2019 Jiaming Huang. All rights reserved.
//

#include "Least_Square_Monte_Carlo_Simulation.hpp"
#include "Stochastic_Process_Simulator.hpp"
#include "Matrix.hpp"
#include "Utils.hpp"
#include <iostream>
#include <cmath>

using namespace std;

double Least_Square_Monte_Carlo_Simulation::Hermite(double x,int n){
    if (n==4) return 8*x*x*x-12*x;
    if(n==3)return 4*x*x-2;
    if(n==2)return 2*x;
    else return 1;
}

double Least_Square_Monte_Carlo_Simulation::Lagueree(double x, int n){
    if (n==4) return exp(-0.5*x)*(1-3*x+3*x*x/2-x*x*x/6);
    if(n==3)return exp(-0.5*x)*(1-2*x+x*x/2);
    if(n==2)return exp(-0.5*x)*(1-x);
    else return exp(-0.5*x);
}

double Least_Square_Monte_Carlo_Simulation::Monomial(double x, int n){
    if (n==4) return x*x*x;
    if(n==3)return x*x;
    if(n==2)return x;
    else return 1;
}

shared_ptr<Matrix> Least_Square_Monte_Carlo_Simulation::MultiStockPathsSimulation(double S0, double T, int n, double mu, double sigma, int seed, int size){
    double ** result = new double*[size];
    for (int i = 0; i < size; i++) {
        double * winierpath = Stochastic_Process_Simulation::WinierPathSimulation(T, n, seed+i);
        cout<< Utils::mean(winierpath, n)<<endl;
        cout<< Utils::Std(winierpath, n)<<endl;
        result[i] = Stochastic_Process_Simulation::StockPathSimulation(S0, T, n, mu, sigma, winierpath);
    }
    shared_ptr<Matrix> re = make_shared<Matrix>(result,size,n+1);
    for (int i = 0; i < size; i++) {
        delete [] result[i];
    }
    delete []result;
    return re;
}

shared_ptr<Matrix> Least_Square_Monte_Carlo_Simulation::LeastSquaredRegression(shared_ptr<Matrix> St, shared_ptr<Matrix> Yt, string type, int k, int size, int n){
    shared_ptr<Matrix> A = make_shared<Matrix>(0.0, k, k);
    shared_ptr<Matrix> b = make_shared<Matrix>(0.0, k, 1);
    for (int i = 0; i < k; i++) {
        double sum2 = 0;
        for (int j = 0; j < k; j++) {
            double sum = 0;
            for (int p = 0; p < size; p++) {
                if(type == "Hermite")sum+=Hermite(St->elementAt(p, 0), i+1)*Hermite(St->elementAt(p, 0), j+1);
                if(type == "Monomial")sum+=Monomial(St->elementAt(p, 0), i+1)*Monomial(St->elementAt(p, 0), j+1);
                if(type == "Lagueree")sum+=Lagueree(St->elementAt(p, 0), i+1)*Lagueree(St->elementAt(p, 0), j+1);
            }
            A->setValue(sum, i, j);
        }
        for (int t = 0; t < size; t++) {
            if(type == "Hermite")sum2+=Yt->elementAt(t, 0)*Hermite(St->elementAt(t, 0), i+1);
            if(type == "Monomial")sum2+=Yt->elementAt(t, 0)*Monomial(St->elementAt(t, 0), i+1);
            if(type == "Lagueree")sum2+=Yt->elementAt(t, 0)*Lagueree(St->elementAt(t, 0), i+1);
        }
        b->setValue(sum2, i, 0);
    }
    shared_ptr<Matrix> a = (*(A->Inverse()))*b;
    shared_ptr<Matrix> result = make_shared<Matrix>(0.0, size, 1);
    for (int i = 0; i < size; i++) {
        double ecv = 0;
        for (int j = 0; j < k; j++) {
            if(type == "Hermite")ecv+=a->elementAt(j, 0)*Hermite(St->elementAt(i, 0), j+1);
            if(type == "Monomial")ecv+=a->elementAt(j, 0)*Monomial(St->elementAt(i, 0), j+1);
            if(type == "Lagueree")ecv+=a->elementAt(j, 0)*Lagueree(St->elementAt(i, 0), j+1);
        }
        result->setValue(ecv, i, 0);
    }
    return result;
}

double Least_Square_Monte_Carlo_Simulation::EuroForwardStartPricing(double t, double S0, double T, int n, double mu, double sigma,int size){
    int seed = rand();
    shared_ptr<Matrix> sim = Least_Square_Monte_Carlo_Simulation::MultiStockPathsSimulation(S0, T, n, mu, sigma,seed, size);
    sim->Print();
    double * ST = new double[size];
    for (int i = 0; i < size; i++) {
        ST[i] = sim->elementAt(i, n);
    }
    cout<<Utils::mean(ST, size);
    double sum = 0;
    for (int i= 0; i < size; i++) {
        sum += (sim->elementAt(i, t/T*n)-sim->elementAt(i, n))>0?(sim->elementAt(i, t/T*n)-sim->elementAt(i, n)):0;
    }
    return (sum/size)*exp(-1*mu*T);
}
