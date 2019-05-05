//
//  FiniteDIfferenceMethod.cpp
//  Computational_Methods_Final
//
//  Created by 黄佳铭 on 2019/3/21.
//  Copyright © 2019 Jiaming Huang. All rights reserved.
//

#include "FiniteDIfferenceMethod.hpp"
#include <iostream>
#include <cmath>
#include "Matrix.hpp"
#include "Utils.hpp"


using namespace std;

FiniteDifferenceMethod::FiniteDifferenceMethod(double S0, double r, double sigma, double delta_t,double delta_x,double delta_s, double maturity, double K,double paths, double N,string type){
    this->S0 = S0;
    this->r = r;
    this->sigma = sigma;
    this->delta_t = delta_t;
    this->delta_x = delta_x;
    this->delta_s = delta_s;
    this->maturity = maturity;
    this->K = K;
    this->paths = paths;
    this->N = N;
    this->type = type;
}

shared_ptr<Matrix> FiniteDifferenceMethod::xsimulation(){
    double *multi = new double[2*(int)paths+1];
    double x0 = log(S0);
    for (int i = 0; i < paths; i++) {
        multi[i] = (paths-i)*delta_x+x0;
    }
    multi[(int)paths] = x0;
    for (int i = paths+1; i < 2*(int)paths+1; i++) {
        multi[i] = (paths-i)*delta_x+x0;
    }
    shared_ptr<Matrix> Xt = make_shared<Matrix>(multi,2*paths+1,1);
    delete [] multi;
    return Xt;
}

shared_ptr<Matrix> FiniteDifferenceMethod::ssimulation(){
    double *multi = new double[(int)paths+1];
    for (int i = 0; i <= paths; i++) {
        multi[i] = (paths-i)*delta_s+0;
    }
    shared_ptr<Matrix> St = make_shared<Matrix>(multi,paths+1,1);
    St->Print();
    delete [] multi;
    return St;
}

double FiniteDifferenceMethod::ExplicitFiniteDifferenceMethodForX(){
    double Nforpaths = 2*paths+1;
    shared_ptr<Matrix> x = xsimulation();
    x->Print();
    shared_ptr<Matrix> s = make_shared<Matrix>(0,Nforpaths,1);
    for (int i = 0; i < Nforpaths; i++) {
        s->setValue(exp(x->elementAt(i, 0)), i, 0);
    }
    s->Print();
    shared_ptr<Matrix> A = make_shared<Matrix>(0,Nforpaths,Nforpaths);
    //construct matrix A
    double Pu = delta_t*(sigma*sigma/(2*delta_x*delta_x)+(r-sigma*sigma/2)/(2*delta_x));
    double Pm = 1-delta_t*sigma*sigma/(delta_x*delta_x)-r*delta_t;
    double Pd = delta_t*(sigma*sigma/(2*delta_x*delta_x)-(r-sigma*sigma/2)/(2*delta_x));
    for (int i = 0; i < Nforpaths; i++) {
        if (i == 0) {
            A->setValue(Pu, 0, 0);
            A->setValue(Pm, 0, 1);
            A->setValue(Pd, 0, 2);
        }
        else{
            if (i == Nforpaths-1) {
                A->setValue(Pu, i, i-2);
                A->setValue(Pm, i, i-1);
                A->setValue(Pd, i, i);
            }
            else{
                A->setValue(Pu, i, i-1);
                A->setValue(Pm, i, i);
                A->setValue(Pd, i, i+1);
            }
        }
    }
    A->Print();
    //end constructing matrix A
    //Calculation of option price
    shared_ptr<Matrix> Ct = NULL;
    for (int i = N-1; i >= 0; i--) {
        cout<<i<<endl;
        if (i == N-1) {
            if(type == "EuroCall"||type == "AmeriCall")Ct = (*s - K)->max(0);
            if(type == "EuroPut"||type == "AmeriPut")Ct = ((*((*s)*(-1)))+K)->max(0);
        }
        else{
            shared_ptr<Matrix> F = Ct;
            shared_ptr<Matrix> B = make_shared<Matrix>(0,Nforpaths,1);
            B->setValue(s->elementAt(0, 0)-s->elementAt(1, 0), 0, 0);
            Ct = (*((*A)*F))+B;
            Ct->Print();
            if (type == "AmeriCall"||type == "AmeriPut") {
                for (int j = 0; j < Nforpaths; j++) {
                    double execercise_value = 0;
                    if(type == "AmeriCall")execercise_value = s->elementAt(j, 0)-K>0?s->elementAt(j, 0)-K:0;
                    if(type == "AmeriPut")execercise_value = K-s->elementAt(j, 0)>0?K-s->elementAt(j, 0):0;
                    double continuous_value = Ct->elementAt(j, 0);
                    double value = execercise_value>=continuous_value?execercise_value:continuous_value;
                    Ct->setValue(value, j, 0);
                }
            }
        }
    }
    double price = 0;
    for (int i = 0; i < Nforpaths; i++) {
        if(i == paths)price = Ct->elementAt(i, 0);
    }
    return price;
}

double FiniteDifferenceMethod::ImpliciteFiniteDifferenceMethodForX(){
    double Nforpaths = 2*paths+1;
    shared_ptr<Matrix> x = xsimulation();
    shared_ptr<Matrix> s = make_shared<Matrix>(0,Nforpaths,1);
    for (int i = 0; i < Nforpaths; i++) {
        s->setValue(exp(x->elementAt(i, 0)), i, 0);
    }
    shared_ptr<Matrix> A = make_shared<Matrix>(0,Nforpaths,Nforpaths);
    //construct matrix A
    double v = r-sigma*sigma/2;
    double Pu = -0.5*delta_t*(sigma*sigma/(delta_x*delta_x)+(v)/(delta_x));
    double Pm = 1+delta_t*sigma*sigma/(delta_x*delta_x)+r*delta_t;
    double Pd = -0.5*delta_t*(sigma*sigma/(delta_x*delta_x)-(v)/(delta_x));
    for (int i = 0; i < Nforpaths; i++) {
        if (i == 0) {
            A->setValue(1, 0, 0);
            A->setValue(-1, 0, 1);
        }
        else{
            if (i == Nforpaths-1) {
                A->setValue(1, i, i-1);
                A->setValue(-1, i, i);
            }
            else{
                A->setValue(Pu, i, i-1);
                A->setValue(Pm, i, i);
                A->setValue(Pd, i, i+1);
            }
        }
    }
    //end constructing matrix A
    //Calculation of option price
    shared_ptr<Matrix> Ct = NULL;
    shared_ptr<Matrix> IA = A->Inverse();
    for (int i = N-1; i >= 0; i--) {
        cout<<i<<endl;
        if (i == N-1) {
            if(type == "EuroCall"||type == "AmeriCall")Ct = (*s - K)->max(0);
            if(type == "EuroPut"||type == "AmeriPut")Ct = ((*((*s)*(-1)))+K)->max(0);
        }
        else{
            // constructing matrix B
            shared_ptr<Matrix> B = Ct;
            B->setValue(s->elementAt(0, 0)-s->elementAt(1, 0), 0, 0);
            B->setValue(0, Nforpaths-1, 0);
            Ct = (*(IA))*B;
            if (type == "AmeriCall"||type == "AmeriPut") {
                for (int j = 0; j < Nforpaths; j++) {
                    double execercise_value = 0;
                    if(type == "AmeriCall")execercise_value = s->elementAt(j, 0)-K>0?s->elementAt(j, 0)-K:0;
                    if(type == "AmeriPut")execercise_value = K-s->elementAt(j, 0)>0?K-s->elementAt(j, 0):0;
                    double continuous_value = Ct->elementAt(j, 0);
                    double value = execercise_value>=continuous_value?execercise_value:continuous_value;
                    Ct->setValue(value, j, 0);
                }
            }
        }
    }
    double price = 0;
    for (int i = 0; i < Nforpaths; i++) {
        if(i == paths)price = Ct->elementAt(i, 0);
    }
    return price;
}

double FiniteDifferenceMethod::CrankNicolsonFiniteDifferenceMethodForX(){
    double Nforpaths = 2*paths+1;
    shared_ptr<Matrix> x = xsimulation();
    shared_ptr<Matrix> s = make_shared<Matrix>(0,Nforpaths,1);
    for (int i = 0; i < Nforpaths; i++) {
        s->setValue(exp(x->elementAt(i, 0)), i, 0);
    }
    shared_ptr<Matrix> A = make_shared<Matrix>(0,Nforpaths,Nforpaths);
    //construct matrix A
    double v = r-sigma*sigma/2;
    double Pu = -0.5*delta_t*(sigma*sigma/(delta_x*delta_x)+(v)/(delta_x));
    double Pm = 1+delta_t*sigma*sigma/(delta_x*delta_x)+r*delta_t;
    double Pd = -0.5*delta_t*(sigma*sigma/(delta_x*delta_x)-(v)/(delta_x));
    for (int i = 0; i < Nforpaths; i++) {
        if (i == 0) {
            A->setValue(1, 0, 0);
            A->setValue(-1, 0, 1);
        }
        else{
            if (i == Nforpaths-1) {
                A->setValue(1, i, i-1);
                A->setValue(-1, i, i);
            }
            else{
                A->setValue(Pu, i, i-1);
                A->setValue(Pm, i, i);
                A->setValue(Pd, i, i+1);
            }
        }
    }
    //end constructing matrix A
    //Calculation of option price
    shared_ptr<Matrix> Ct = NULL;
    shared_ptr<Matrix> IA = A->Inverse();
    for (int i = N-1; i >= 0; i--) {
        cout<<i<<endl;
        if (i == N-1) {
            if(type == "EuroCall"||type == "AmeriCall")Ct = (*s - K)->max(0);
            if(type == "EuroPut"||type == "AmeriPut")Ct = ((*((*s)*(-1)))+K)->max(0);
        }
        else{
            
            // constructing matrix B
            shared_ptr<Matrix> B = make_shared<Matrix>(0,Nforpaths,1);
            for (int i = 0; i < Nforpaths; i++) {
                if (i == 0) {
                    B->setValue(s->elementAt(0, 0)-s->elementAt(1, 0), 0, 0);
                }
                else{
                    if (i == Nforpaths-1) {
                        B->setValue(0, i, 0);
                    }
                    else{
                        double bvalue = -1*Pu*Ct->elementAt(i-1, 0)+(-1)*(Pm-2)*Ct->elementAt(i, 0)-Pd*Ct->elementAt(i+1, 0);
                        B->setValue(bvalue, i, 0);
                    }
                }
            }
            
            Ct = (*(IA))*B;
            if (type == "AmeriCall"||type == "AmeriPut") {
                for (int j = 0; j < Nforpaths; j++) {
                    double execercise_value = 0;
                    if(type == "AmeriCall")execercise_value = s->elementAt(j, 0)-K>0?s->elementAt(j, 0)-K:0;
                    if(type == "AmeriPut")execercise_value = K-s->elementAt(j, 0)>0?K-s->elementAt(j, 0):0;
                    double continuous_value = Ct->elementAt(j, 0);
                    double value = execercise_value>=continuous_value?execercise_value:continuous_value;
                    Ct->setValue(value, j, 0);
                }
            }
        }
    }
    double price = 0;
    for (int i = 0; i < Nforpaths; i++) {
        if(i == paths)price = Ct->elementAt(i, 0);
    }
    return price;
}

double FiniteDifferenceMethod::ExplicitFiniteDifferenceMethodForS(){
    shared_ptr<Matrix> s = ssimulation();
    s->Print();
    shared_ptr<Matrix> A = make_shared<Matrix>(0,paths+1,paths+1);
    //construct matrix A
    for (int i = 0; i <= paths; i++) {
        if (i == 0) {
            double Pu = delta_t*(r*(paths-1)/2+sigma*sigma*(paths-1)*(paths-1)/2);
            double Pm = 1-delta_t*(sigma*sigma*(paths-1)*(paths-1)+r);
            double Pd = delta_t*(-r*(paths-1)/2+sigma*sigma*(paths-1)*(paths-1)/2);
            A->setValue(Pu, 0, 0);
            A->setValue(Pm, 0, 1);
            A->setValue(Pd, 0, 2);
        }
        else{
            if (i == paths) {
                double Pu = delta_t*(r*(1)/2+sigma*sigma*(1)*(1)/2);
                double Pm = 1-delta_t*(sigma*sigma*(1)*(1)+r);
                double Pd = delta_t*(-r*(1)/2+sigma*sigma*(1)*(1)/2);
                A->setValue(Pu, i, i-2);
                A->setValue(Pm, i, i-1);
                A->setValue(Pd, i, i);
            }
            else{
                int j = paths-i;
                double Pu = delta_t*(r*(j)/2+sigma*sigma*(j)*(j)/2);
                double Pm = 1-delta_t*(sigma*sigma*(j)*(j)+r);
                double Pd = delta_t*(-r*(j)/2+sigma*sigma*(j)*(j)/2);
                A->setValue(Pu, i, i-1);
                A->setValue(Pm, i, i);
                A->setValue(Pd, i, i+1);
            }
        }
    }
    A->Print();
    //end constructing matrix A
    //Calculation of option price
    shared_ptr<Matrix> Ct = NULL;
    for (int i = N-1; i >= 0; i--) {
        if (i == N-1) {
            if(type == "EuroCall"||type == "AmeriCall")Ct = (*s - K)->max(0);
            if(type == "EuroPut"||type == "AmeriPut")Ct = ((*((*s)*(-1)))+K)->max(0);
        }
        else{
            shared_ptr<Matrix> F = Ct;
            shared_ptr<Matrix> B = make_shared<Matrix>(0,paths+1,1);
            B->setValue(s->elementAt(0, 0)-s->elementAt(1, 0), 0, 0);
            Ct = (*((*A)*F))+B;
            if (type == "AmeriCall"||type == "AmeriPut") {
                for (int j = 0; j <= paths; j++) {
                    double execercise_value = 0;
                    if(type == "AmeriCall")execercise_value = s->elementAt(j, 0)-K>0?s->elementAt(j, 0)-K:0;
                    if(type == "AmeriPut")execercise_value = K-s->elementAt(j, 0)>0?K-s->elementAt(j, 0):0;
                    double continuous_value = Ct->elementAt(j, 0);
                    double value = execercise_value>=continuous_value?execercise_value:continuous_value;
                    Ct->setValue(value, j, 0);
                }
            }
        }
        Ct->Print();
    }
    shared_ptr<Matrix> F = Ct;
    shared_ptr<Matrix> B = make_shared<Matrix>(0,paths+1,1);
    B->setValue(s->elementAt(0, 0)-s->elementAt(1, 0), 0, 0);
    Ct = (*((*A)*F))+B;
    double price = 0;
    for (int i = 0; i <= paths; i++) {
        if(i == int(paths/2))price = Ct->elementAt(i, 0);
    }
    return price;
}

double FiniteDifferenceMethod::ImpliciteFiniteDifferenceMethodForS(){
    shared_ptr<Matrix> s = ssimulation();
    shared_ptr<Matrix> A = make_shared<Matrix>(0,paths+1,paths+1);
    //construct matrix A
    for (int i = 0; i <= paths; i++) {
        if (i == 0) {
            A->setValue(1, 0, 0);
            A->setValue(-1, 0, 1);
        }
        else{
            if (i == paths) {
                A->setValue(1, i, i-1);
                A->setValue(-1, i, i);
            }
            else{
                double j = paths-i;
                double Pu = -0.5*delta_t*(sigma*sigma*j*j+r*j);
                double Pm = 1+delta_t*sigma*sigma*j*j+delta_t*r;
                double Pd = 0.5*delta_t*(r*j-sigma*sigma*j*j);
                A->setValue(Pu, i, i-1);
                A->setValue(Pm, i, i);
                A->setValue(Pd, i, i+1);
            }
        }
    }
    A->Print();
    //end constructing matrix A
    //Calculation of option price
    shared_ptr<Matrix> IA = A->Inverse();
    shared_ptr<Matrix> Ct = NULL;
    for (int i = N-1; i >= 0; i--) {
        if (i == N-1) {
            if(type == "EuroCall"||type == "AmeriCall")Ct = (*s - K)->max(0);
            if(type == "EuroPut"||type == "AmeriPut")Ct = ((*((*s)*(-1)))+K)->max(0);
        }
        else{
            // constructing matrix B
            shared_ptr<Matrix> B = Ct;
            B->setValue(s->elementAt(0, 0)-s->elementAt(1, 0), 0, 0);
            B->setValue(0, paths, 0);
            Ct = (*(IA))*B;
            if (type == "AmeriCall"||type == "AmeriPut") {
                for (int j = 0; j <= paths; j++) {
                    double execercise_value = 0;
                    if(type == "AmeriCall")execercise_value = s->elementAt(j, 0)-K>0?s->elementAt(j, 0)-K:0;
                    if(type == "AmeriPut")execercise_value = K-s->elementAt(j, 0)>0?K-s->elementAt(j, 0):0;
                    double continuous_value = Ct->elementAt(j, 0);
                    double value = execercise_value>=continuous_value?execercise_value:continuous_value;
                    Ct->setValue(value, j, 0);
                }
            }
        }
    }
    double price = 0;
    for (int i = 0; i <= paths; i++) {
        if(i == int(paths/2))price = Ct->elementAt(i, 0);
    }
    return price;
}

double FiniteDifferenceMethod::CrankNicolsonFiniteDifferenceMethodForS(){
    shared_ptr<Matrix> s = ssimulation();
    shared_ptr<Matrix> A = make_shared<Matrix>(0,paths+1,paths+1);
    //construct matrix A
    for (int i = 0; i <= paths; i++) {
        if (i == 0) {
            A->setValue(1, 0, 0);
            A->setValue(-1, 0, 1);
        }
        else{
            if (i == paths) {
                A->setValue(1, i, i-1);
                A->setValue(-1, i, i);
            }
            else{
                double j = paths-i;
                double Pu = sigma*sigma*j*j*delta_t+r*delta_t*j;
                double Pm = -1*(4+2*sigma*sigma*j*j*delta_t+2*r*delta_t);
                double Pd = sigma*sigma*j*j*delta_t-r*j*delta_t;
                A->setValue(Pu, i, i-1);
                A->setValue(Pm, i, i);
                A->setValue(Pd, i, i+1);
            }
        }
    }
    //end constructing matrix A
    //Calculation of option price
    shared_ptr<Matrix> IA = A->Inverse();
    shared_ptr<Matrix> Ct = NULL;
    for (int i = N-1; i >= 0; i--) {
        if (i == N-1) {
            if(type == "EuroCall"||type == "AmeriCall")Ct = (*s - K)->max(0);
            if(type == "EuroPut"||type == "AmeriPut")Ct = ((*((*s)*(-1)))+K)->max(0);
        }
        else{
            
            // constructing matrix B
            shared_ptr<Matrix> B = make_shared<Matrix>(0,paths+1,1);
            for (int i = 0; i <= paths; i++) {
                if (i == 0) {
                    B->setValue(s->elementAt(0, 0)-s->elementAt(1, 0), 0, 0);
                }
                else{
                    if (i == paths) {
                        B->setValue(0, i, 0);
                    }
                    else{
                        double j = paths-i;
                        double Pu = sigma*sigma*j*j*delta_t+r*delta_t*j;
                        double Pm = -1*(4+2*sigma*sigma*j*j*delta_t+2*r*delta_t);
                        double Pd = sigma*sigma*j*j*delta_t-r*j*delta_t;
                        double bvalue = -1*Pu*Ct->elementAt(i-1, 0)+(-1)*(Pm+8)*Ct->elementAt(i, 0)-Pd*Ct->elementAt(i+1, 0);
                        B->setValue(bvalue, i, 0);
                    }
                }
            }
            
            Ct = (*(IA))*B;
            if (type == "AmeriCall"||type == "AmeriPut") {
                for (int j = 0; j <= paths; j++) {
                    double execercise_value = 0;
                    if(type == "AmeriCall")execercise_value = s->elementAt(j, 0)-K>0?s->elementAt(j, 0)-K:0;
                    if(type == "AmeriPut")execercise_value = K-s->elementAt(j, 0)>0?K-s->elementAt(j, 0):0;
                    double continuous_value = Ct->elementAt(j, 0);
                    double value = execercise_value>=continuous_value?execercise_value:continuous_value;
                    Ct->setValue(value, j, 0);
                }
            }
        }
    }
    double price = 0;
    for (int i = 0; i <= paths; i++) {
        if(i == int(paths/2))price = Ct->elementAt(i, 0);
    }
    return price;
}
