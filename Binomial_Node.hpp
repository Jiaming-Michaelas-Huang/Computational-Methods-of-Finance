//
//  Binomial_Node.hpp
//  Computational_Methods_Final
//
//  Created by 黄佳铭 on 2019/3/14.
//  Copyright © 2019 Jiaming Huang. All rights reserved.
//

#ifndef Binomial_Node_hpp
#define Binomial_Node_hpp

#include <stdio.h>
class Binomial_Node{
private:
    double S;
    double C;
    double p;
    double u;
    double d;
    double time;
    Binomial_Node * up;
    Binomial_Node * down;
public:
    Binomial_Node(double time, double S, double p, double u, double d);
    void Set_Up(Binomial_Node * node);
    void Set_Down(Binomial_Node * node);
    void Set_OptionPrice(double (*payoff)(double x, double k), double r, double delta_t,double K);
    Binomial_Node * Get_Up();
    Binomial_Node * Get_Down();
    double Get_OptionPrice();
    double Get_Stock_Value();
    double Get_Time();
    double Get_Prob();
};
#endif /* Binomial_Node_hpp */
