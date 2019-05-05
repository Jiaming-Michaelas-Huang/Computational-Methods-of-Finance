//
//  Random_Number_Generator.hpp
//  Computational_Methods_Final
//
//  Created by 黄佳铭 on 2019/3/13.
//  Copyright © 2019 Jiaming Huang. All rights reserved.
//

#ifndef Random_Number_Generator_hpp
#define Random_Number_Generator_hpp

#include <stdio.h>
#include <iostream>
#include <cmath>

using namespace std;

class Random_Number_Generator{
public:
    static unsigned long int * LGM_Generator(int size, int seed);
    static unsigned long int * RANDU_Generator(int size, int seed);
    static unsigned long int * APPLE_Generator(int size, int seed);
    static unsigned long int * Turbo_Pascal_Generator(int size, int seed);
    static unsigned long int * Wu_1997_Generator(int size, int seed);
    static double * Uniform_Distribution_Generator(int size, int seed, string generator_type);
    static double * Bernoulli_Distribution_Genrator(int size, int seed, string generator_type, double p);
    static double * Binomial_Distribution_Genrator(int size, int seed, string generator_type, double p, int n);
    static double * Poisson_Distribution_Generator(int size, int seed, string generator_type, double lambda);
    static double * Exponential_Distribution_Generator(int size, int seed, string generator_type, double lambda);
    static double * Gamma_Distribution_Generator(int size, int seed, string generator_type, double lambda, int n);
    static double * Logistic_Distribution_Generator(int size, int seed, string generator_type, double a, double b);
    static double * Normal_Distribution_Generator_Box_Muller(int size, int seed, string generator_type);
    static double * Normal_Distribution_Generator_Pollar_Masaglia(int size, int seed, string generator_type);
    static double * Discrete_Distribution_Genrator(int size,int seed,string generator_type, double * probs, double * values, int n);
    static double ** Bivariate_Normal_Distribution_Generator(int size, double * z1, double * z2, double rho);
    static double * Halton_Sequence_Generator(int base, int size);
};
#endif /* Random_Number_Generator_hpp */
