//
//  Random_Number_Generator.cpp
//  Computational_Methods_Final
//
//  Created by 黄佳铭 on 2019/3/13.
//  Copyright © 2019 Jiaming Huang. All rights reserved.
//

#include "Random_Number_Generator.hpp"
#include <iostream>
#include <cmath>

using namespace std;

unsigned long int * Random_Number_Generator::LGM_Generator(int size, int seed){
    unsigned long int a = pow(7, 5);
    unsigned long int m = pow(2, 31)-1;
    unsigned long int b = 0;
    unsigned long int * pseudorandom_samples = new unsigned long int[size];
    unsigned long int x0= seed;
    pseudorandom_samples[0]  = x0;
    //pseudorandom number generator
    for (int i = 1; i<size; i++) {
        pseudorandom_samples[i] = (pseudorandom_samples[i-1] * a + b) % m;
    }
    return pseudorandom_samples;
}

unsigned long int * Random_Number_Generator::RANDU_Generator(int size, int seed){
    unsigned long int a = pow(2, 16)+3;
    unsigned long int m = pow(2, 31);
    unsigned long int b = 0;
    unsigned long int * pseudorandom_samples = new unsigned long int[size];
    unsigned long int x0= seed;
    pseudorandom_samples[0]  = x0;
    //pseudorandom number generator
    for (int i = 1; i<size; i++) {
        pseudorandom_samples[i] = (pseudorandom_samples[i-1] * a + b) % m;
    }
    return pseudorandom_samples;
}

unsigned long int * Random_Number_Generator::APPLE_Generator(int size, int seed){
    unsigned long int a = pow(5, 13);
    unsigned long int m = pow(2, 35);
    unsigned long int b = 0;
    unsigned long int * pseudorandom_samples = new unsigned long int[size];
    unsigned long int x0= seed;
    pseudorandom_samples[0]  = x0;
    //pseudorandom number generator
    for (int i = 1; i<size; i++) {
        pseudorandom_samples[i] = (pseudorandom_samples[i-1] * a + b) % m;
    }
    return pseudorandom_samples;
}

unsigned long int * Random_Number_Generator::Turbo_Pascal_Generator(int size, int seed){
    unsigned long int a = 134775813;
    unsigned long int m = pow(2, 32);
    unsigned long int b = 1;
    unsigned long int * pseudorandom_samples = new unsigned long int[size];
    unsigned long int x0= seed;
    pseudorandom_samples[0]  = x0;
    //pseudorandom number generator
    for (int i = 1; i<size; i++) {
        pseudorandom_samples[i] = (pseudorandom_samples[i-1] * a + b) % m;
    }
    return pseudorandom_samples;
}

unsigned long int * Random_Number_Generator::Wu_1997_Generator(int size, int seed){
    unsigned long int a = pow(2, 19)-1;
    unsigned long int m = pow(2, 61)-1;
    unsigned long int b = 0;
    unsigned long int * pseudorandom_samples = new unsigned long int[size];
    unsigned long int x0= seed;
    pseudorandom_samples[0]  = x0;
    //pseudorandom number generator
    for (int i = 1; i<size; i++) {
        pseudorandom_samples[i] = (pseudorandom_samples[i-1] * a + b) % m;
    }
    return pseudorandom_samples;
}

double * Random_Number_Generator::Uniform_Distribution_Generator(int size, int seed, string generator_type){
    unsigned long int * pseudorandom_samples = NULL;
    unsigned long int m = 0;
    if(generator_type == "LGM"){pseudorandom_samples = LGM_Generator(size, seed);m = pow(2, 31)-1;}
    if(generator_type == "RANDU"){pseudorandom_samples = RANDU_Generator(size, seed);m = pow(2, 31);}
    if(generator_type == "APPLE"){pseudorandom_samples = APPLE_Generator(size, seed);m = pow(2, 35);}
    if(generator_type == "Turbo-Pascal"){pseudorandom_samples = Turbo_Pascal_Generator(size, seed);m = pow(2, 32);}
    if(generator_type == "Wu-1997"){pseudorandom_samples = Wu_1997_Generator(size, seed);m = pow(2, 61)-1;}
    double * uniform_samples = new double[size];
    for (int i = 0; i < size; i++) {
        uniform_samples[i] = (1.0*pseudorandom_samples[i])/(1.0*m);
    }
    delete []pseudorandom_samples;
    return uniform_samples;
}

double * Random_Number_Generator::Bernoulli_Distribution_Genrator(int size, int seed, string generator_type, double p){
    double * uniform_samples = Uniform_Distribution_Generator(size, seed, generator_type);
    double * bernoulli_samples = new double[size];
    for (int i = 0; i < size; i++) {
        if(uniform_samples[i]<p)bernoulli_samples[i] = 1;
        else bernoulli_samples[i] = 0;
    }
    delete []uniform_samples;
    return bernoulli_samples;
}

double * Random_Number_Generator::Binomial_Distribution_Genrator(int size, int seed, string generator_type, double p, int n){
    double * binomial_samples = new double[size];
    for (int i = 0; i < size; i++) binomial_samples[i] = 0;
    double * bernoulli_samples = Bernoulli_Distribution_Genrator(size*n, seed, generator_type, p);
    for (int j = 0; j < n; j++) {
        for (int k = 0; k < size; k++) {
            binomial_samples[k] = binomial_samples[k] + bernoulli_samples[k+j*size];
        }
    }
    delete []bernoulli_samples;
    return binomial_samples;
}

double * Random_Number_Generator::Poisson_Distribution_Generator(int size, int seed, string generator_type, double lambda){
    double * poisson_samples = new double[size];
    double * uniform_samples = Uniform_Distribution_Generator(size, seed, generator_type);
    for (int i = 0; i < size; i++) {
        double k = 0;
        double x = exp(-1*lambda);
        double z = x;
        while (uniform_samples[i]>=x) {
            z = (lambda/(k+1))*z;
            x = x + z;
            k = k+1;
        }
        poisson_samples[i] = k;
    }
    delete []uniform_samples;
    return poisson_samples;
}

double * Random_Number_Generator::Exponential_Distribution_Generator(int size, int seed, string generator_type, double lambda){
    double * exponential_samples = new double[size];
    double * uniform_samples = Uniform_Distribution_Generator(size, seed, generator_type);
    for (int i = 0; i < size; i++) {
        exponential_samples[i] = -1*lambda*log(1-uniform_samples[i]);
    }
    delete []uniform_samples;
    return exponential_samples;
}

double * Random_Number_Generator::Gamma_Distribution_Generator(int size, int seed, string generator_type, double lambda, int n){
    double * gamma_samples = new double[size];
    for (int i = 0; i < size; i++) gamma_samples[i] = 0;
    double * exponential_samples = Exponential_Distribution_Generator(size*n, seed, generator_type, lambda);
    for (int j = 0; j < n; j++) {
        for (int k = 0; k < size; k++) {
            gamma_samples[k] = gamma_samples[k] + exponential_samples[k+j*size];
        }
    }
    delete []exponential_samples;
    return gamma_samples;
}

double * Random_Number_Generator::Logistic_Distribution_Generator(int size, int seed, string generator_type, double a, double b){
    double * logistic_samples = new double[size];
    double * uniform_samples = Uniform_Distribution_Generator(size, seed, generator_type);
    for (int i = 0; i < size; i++) {
        logistic_samples[i] = a + b * log(uniform_samples[i]/(1-uniform_samples[i]));
    }
    delete []uniform_samples;
    return logistic_samples;
}

double * Random_Number_Generator::Normal_Distribution_Generator_Box_Muller(int size, int seed, string generator_type){
    double * normal_samples = new double[size];
    double * u1 = Uniform_Distribution_Generator(size, seed, generator_type);
    double * u2 = Uniform_Distribution_Generator(size, seed+1, generator_type);
    int t = 0;
    for(int i = 0 ; i < size; i+=2){
        normal_samples[i] = sqrt((-2*log(u1[t])))*cos(2*M_PI*u2[t]);
        normal_samples[i+1] = sqrt((-2*log(u1[t])))*sin(2*M_PI*u2[t]);
        t += 1;
    }
    delete []u1;
    delete []u2;
    return normal_samples;
}

double * Random_Number_Generator::Normal_Distribution_Generator_Pollar_Masaglia(int size, int seed, string generator_type){
    double * normal_samples = new double[size];
    double * u1 = Uniform_Distribution_Generator(size, seed, generator_type);
    double * u2 = Uniform_Distribution_Generator(size, seed, generator_type);
    double v1,v2,w;
    int t = 0;
    for(int i = 0 ; i < size; i++){
        v1 = 2*u1[i]-1;
        v2 = 2*u2[i]-1;
        w = v1*v1+v2*v2;
        if (w <= 1) {
            normal_samples[t] = v1 * sqrt((-2*log(w))/w);
            normal_samples[t+1] = v2 * sqrt((-2*log(w))/w);
            t += 2;
        }
        if (t >= size) {
            break;
        }
    }
    return normal_samples;
}

double * Random_Number_Generator::Discrete_Distribution_Genrator(int size, int seed,string generator_type, double *probs, double *values, int n){
    double * discrete_samples = new double[size];
    double * uniform_samples = Uniform_Distribution_Generator(size, seed, generator_type);
    for (int i = 0; i < size; i++) {
        double low = 0;
        double high = probs[0];
        for (int j = 0; j < n; j++) {
            if (uniform_samples[i]>=low && uniform_samples[i]<high) {
                discrete_samples[i] = values[j];
                break;
            }
            low = low + probs[j];
            high = high + probs[j+1];
        }
    }
    return discrete_samples;
}

double ** Random_Number_Generator::Bivariate_Normal_Distribution_Generator(int size, double *z1, double *z2, double rho){
    double ** binormal_samples = new double *[2];
    binormal_samples[0] = new double[size];
    binormal_samples[1] = new double[size];
    for (int i = 0; i < size; i++) {
        binormal_samples[0][i] = z1[i];
        binormal_samples[1][i] = rho*z1[i] + sqrt(1-rho*rho)*z2[i];
    }
    return binormal_samples;
}

double * Random_Number_Generator::Halton_Sequence_Generator(int base, int size){
    double *seq = new double[size];
    for(int i=0; i<size; i++){
        seq[i] = 0;
    }
    
    int NumBits = 1 + ceil(log(size)/log(base));
    
    double vetBase[NumBits];
    double workVet[NumBits];
    for(int i=0; i< NumBits; i++){
        vetBase[i] = pow(base, -(i+1));
        workVet[i] = 0;
    }
    
    for(int i=1; i<=size; i++){
        int num = i;
        int counter = 0;
        while (num > 0){
            workVet[counter] = num % base;
            num /= base;
            counter++;
        }
        seq[i-1] = 0;
        for(int i = 0; i<NumBits; i++){
            seq[i-1] += workVet[i] * vetBase[i];
        }
    }
    
    return seq;
}
