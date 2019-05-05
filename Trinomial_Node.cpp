//
//  Trinomial_Node.cpp
//  Computational_Methods_Final
//
//  Created by 黄佳铭 on 2019/3/14.
//  Copyright © 2019 Jiaming Huang. All rights reserved.
//

#include "Trinomial_Node.hpp"
#include <iostream>
#include <cmath>

using namespace std;

Trinomial_Node::Trinomial_Node(double time, double S, double p,double m, double u, double d){
    this->time = time;
    this->C = -100;
    this->S = S;
    this->p = p;
    this->m = m;
    this->u = u;
    this->d = d;
}

void Trinomial_Node::Set_Up(Trinomial_Node *node){
    this->up = node;
}

void Trinomial_Node::Set_Mid(Trinomial_Node *node){
    this->mid = node;
}

void Trinomial_Node::Set_Down(Trinomial_Node *node){
    this->down = node;
}

Trinomial_Node *Trinomial_Node::Get_Up(){
    return this->up;
}

Trinomial_Node * Trinomial_Node::Get_Mid(){
    return this->mid;
}

Trinomial_Node *Trinomial_Node::Get_Down(){
    return down;
}

double Trinomial_Node::Get_Prob(){
    return this->p;
}

double Trinomial_Node::Get_Mrob(){
    return this->m;
}

double Trinomial_Node::Get_Time(){
    return this->time;
}

double Trinomial_Node::Get_Stock_Value(){
    return this->S;
}

double Trinomial_Node::Get_OptionPrice(){
    return C;
}

void Trinomial_Node::Set_OptionPrice(double (*payoff)(double x, double k), double r, double delta_t,double K){
    if (this->up == NULL && this->down == NULL&&this->mid == NULL) {
        this->C = payoff(this->S,K);
    }
    else{
        if(this->up->Get_OptionPrice()<0)this->up->Set_OptionPrice(payoff, r, delta_t, K);
        if(this->mid->Get_OptionPrice()<0)this->mid->Set_OptionPrice(payoff, r, delta_t, K);
        if(this->down->Get_OptionPrice()<0)this->down->Set_OptionPrice(payoff, r, delta_t, K);
        this->C = exp(-1*r*delta_t)*(this->Get_Prob()*this->up->Get_OptionPrice() +(this->Get_Mrob()*this->mid->Get_OptionPrice()) + (1-this->Get_Prob()-this->Get_Mrob())*this->down->Get_OptionPrice());
    }
}

