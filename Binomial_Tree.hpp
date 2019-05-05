//
//  Binomial_Tree.hpp
//  Computational_Methods_Final
//
//  Created by 黄佳铭 on 2019/3/14.
//  Copyright © 2019 Jiaming Huang. All rights reserved.
//

#ifndef Binomial_Tree_hpp
#define Binomial_Tree_hpp

#include <stdio.h>
#include "Binomial_Node.hpp"
class Binomial_Tree{
private:
    Binomial_Node * root;
    double S0;
    double p;
    double up;
    double down;
    double r;
    double K;
    double T;
    int n;
public:
    Binomial_Tree(double S0, double p, double up, double down, double r, double K, double T, int n);
    double Binomial_Pricing(double (*payoff)(double x, double k));
};
#endif /* Binomial_Tree_hpp */
