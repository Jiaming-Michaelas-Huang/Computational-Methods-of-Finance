//
//  Trinomial_Node.hpp
//  Computational_Methods_Final
//
//  Created by 黄佳铭 on 2019/3/14.
//  Copyright © 2019 Jiaming Huang. All rights reserved.
//

#ifndef Trinomial_Node_hpp
#define Trinomial_Node_hpp

#include <stdio.h>
class Trinomial_Node{
private:
    double S;
    double C;
    double p;
    double m;
    double u;
    double d;
    double time;
    Trinomial_Node * up;
    Trinomial_Node * mid;
    Trinomial_Node * down;
public:
    Trinomial_Node(double time, double S, double p,double m, double u, double d);
    void Set_Up(Trinomial_Node * node);
    void Set_Mid(Trinomial_Node * node);
    void Set_Down(Trinomial_Node * node);
    void Set_OptionPrice(double (*payoff)(double x, double k), double r, double delta_t,double K);
    Trinomial_Node * Get_Up();
    Trinomial_Node * Get_Mid();
    Trinomial_Node * Get_Down();
    double Get_OptionPrice();
    double Get_Stock_Value();
    double Get_Time();
    double Get_Prob();
    double Get_Mrob();
};
#endif /* Trinomial_Node_hpp */
