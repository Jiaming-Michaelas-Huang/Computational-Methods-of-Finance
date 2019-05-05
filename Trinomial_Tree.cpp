//
//  Trinomial_Tree.cpp
//  Computational_Methods_Final
//
//  Created by 黄佳铭 on 2019/3/14.
//  Copyright © 2019 Jiaming Huang. All rights reserved.
//

#include "Trinomial_Tree.hpp"
#include "Trinomial_Node.hpp"
#include <iostream>
#include <cmath>

using namespace std;

Trinomial_Tree::Trinomial_Tree(double S0, double p, double m, double up, double down, double r, double K, double T, int n){
    this->S0 = S0;
    this->p = p;
    this->m = m;
    this->up = up;
    this->down = down;
    this->r = r;
    this->K = K;
    this->T = T;
    this->n = n;
    this->root = new Trinomial_Node(0, S0, p, m, up, down);
    int tt = (1+2*n-1)*n/2;
    Trinomial_Node * nodes[tt];
    for (int i = 0; i < n; i++) {
        if (i == 0) {
            nodes[i] = root;
        }
        else{
            int count = 0;
            for (int j = (1+2*i-1)*i/2; j < (1+2*(1+i)-1)*(1+i)/2; j++) {
                if ((((1+2*(i-1)-1)*(i-1)/2)+count)>=((1+2*i-1)*i/2)) {
                    break;
                }
                if(j == (1+2*i-1)*i/2){
                    Trinomial_Node * upnode = new Trinomial_Node(i+1,up*nodes[((1+2*(i-1)-1)*(i-1)/2)+count]->Get_Stock_Value(), p, m, up, down);
                    Trinomial_Node * midnode = new Trinomial_Node(i+1, nodes[((1+2*(i-1)-1)*(i-1)/2)+count]->Get_Stock_Value(), p, m, up, down);
                    nodes[((1+2*(i-1)-1)*(i-1)/2)+count]->Set_Up(upnode);
                    nodes[((1+2*(i-1)-1)*(i-1)/2)+count]->Set_Mid(midnode);
                }
                else {
                    nodes[((1+2*(i-1)-1)*(i-1)/2)+count]->Set_Up(nodes[((1+2*(i-1)-1)*(i-1)/2)+count-1]->Get_Mid());
                    nodes[((1+2*(i-1)-1)*(i-1)/2)+count]->Set_Mid(nodes[((1+2*(i-1)-1)*(i-1)/2)+count-1]->Get_Down());
                }
                Trinomial_Node * downnode = new Trinomial_Node(i+1, down*nodes[((1+2*(i-1)-1)*(i-1)/2)+count]->Get_Stock_Value(), p, m, up, down);
                nodes[((1+2*(i-1)-1)*(i-1)/2)+count]->Set_Down(downnode);
                nodes[j] = nodes[((1+2*(i-1)-1)*(i-1)/2)+count]->Get_Up();
                nodes[j+1] = nodes[((1+2*(i-1)-1)*(i-1)/2)+count]->Get_Mid();
                nodes[j+2] = nodes[((1+2*(i-1)-1)*(i-1)/2)+count]->Get_Down();
                count++;
            }
        }
    }
}
