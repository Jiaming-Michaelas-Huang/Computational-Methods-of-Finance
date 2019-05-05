//
//  Binomial_Node.cpp
//  Computational_Methods_Final
//
//  Created by 黄佳铭 on 2019/3/14.
//  Copyright © 2019 Jiaming Huang. All rights reserved.
//

#include "Binomial_Node.hpp"
#include <iostream>
#include <cmath>

using namespace std;

Binomial_Node::Binomial_Node(double time, double S, double p, double u, double d){
    this->time = time;
    this->C = -100;
    this->S = S;
    this->p = p;
    this->u = u;
    this->d = d;
}

void Binomial_Node::Set_Up(Binomial_Node *node){
    this->up = node;
}

void Binomial_Node::Set_Down(Binomial_Node *node){
    this->down = node;
}

Binomial_Node *Binomial_Node::Get_Up(){
    return this->up;
}

Binomial_Node *Binomial_Node::Get_Down(){
    return down;
}

double Binomial_Node::Get_Prob(){
    return this->p;
}

double Binomial_Node::Get_Time(){
    return this->time;
}

double Binomial_Node::Get_Stock_Value(){
    return this->S;
}

double Binomial_Node::Get_OptionPrice(){
    return C;
}

void Binomial_Node::Set_OptionPrice(double (*payoff)(double x, double k), double r, double delta_t,double K){
    if (this->up == NULL && this->down == NULL) {
        this->C = payoff(this->S,K);
    }
    else{
        if(this->up->Get_OptionPrice()<0)this->up->Set_OptionPrice(payoff, r, delta_t, K);
        if(this->down->Get_OptionPrice()<0)this->down->Set_OptionPrice(payoff, r, delta_t, K);
        this->C = exp(-1*r*delta_t)*(this->Get_Prob()*this->up->Get_OptionPrice() + (1-this->Get_Prob())*this->down->Get_OptionPrice());
    }
}
