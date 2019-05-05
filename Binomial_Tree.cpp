//
//  Binomial_Tree.cpp
//  Computational_Methods_Final
//
//  Created by 黄佳铭 on 2019/3/14.
//  Copyright © 2019 Jiaming Huang. All rights reserved.
//

#include "Binomial_Tree.hpp"
#include <iostream>
#include <cmath>
#include "Binomial_Node.hpp"

using namespace std;

Binomial_Tree::Binomial_Tree(double S0, double p, double up, double down, double r, double K, double T, int n){
    this->S0 = S0;
    this->p = p;
    this->up = up;
    this->down = down;
    this->r = r;
    this->K = K;
    this->T = T;
    this->n = n;
    this->root = new Binomial_Node(0, S0, p, up, down);
    int tt = (1+n)*n/2;
    Binomial_Node * nodes[tt];
    for (int i = 0; i < n; i++) {
        if (i == 0) {
            nodes[i] = root;
        }
        else{
            int count = 0;
            for (int j = (1+i)*i/2; j < (2+i)*(1+i)/2; j++) {
                if(((i*(i-1)/2)+count)>=(1+i)*i/2)break;
                if(j == (1+i)*i/2){
                    Binomial_Node * upnode = new Binomial_Node(i+1, up * nodes[(i*(i-1)/2)+count]->Get_Stock_Value(),p, up, down);
                    nodes[(i*(i-1)/2)+count]->Set_Up(upnode);
                }
                else nodes[(i*(i-1)/2)+count]->Set_Up(nodes[(i*(i-1)/2)+count-1]->Get_Down());
                Binomial_Node * downnode = new Binomial_Node(i+1, down*nodes[(i*(i-1)/2)+count]->Get_Stock_Value(),p, up, down);
                nodes[(i*(i-1)/2)+count]->Set_Down(downnode);
                nodes[j] = nodes[(i*(i-1)/2)+count]->Get_Up();
                nodes[j+1] = nodes[(i*(i-1)/2)+count]->Get_Down();
                count++;
            }
        }
    }
}

double Binomial_Tree::Binomial_Pricing(double (*payoff)(double, double)){
    this->root->Set_OptionPrice(payoff, this->r, T/n, K);
    return this->root->Get_OptionPrice();
}
