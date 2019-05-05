//
//  Trinomial_Tree.hpp
//  Computational_Methods_Final
//
//  Created by 黄佳铭 on 2019/3/14.
//  Copyright © 2019 Jiaming Huang. All rights reserved.
//

#ifndef Trinomial_Tree_hpp
#define Trinomial_Tree_hpp

#include <stdio.h>
#include "Trinomial_Node.hpp"
class Trinomial_Tree{
private:
    Trinomial_Node * root;
    double S0;
    double p;
    double m;
    double up;
    double down;
    double r;
    double K;
    double T;
    int n;
public:
    Trinomial_Tree(double S0, double p, double m, double up, double down, double r, double K, double T, int n);
    double Trinomial_Pricing(double (*payoff)(double x, double k));
};
#endif /* Trinomial_Tree_hpp */
