//
//  FixedIncomePricing.hpp
//  Computational_Methods_Final
//
//  Created by 黄佳铭 on 2019/3/21.
//  Copyright © 2019 Jiaming Huang. All rights reserved.
//

#ifndef FixedIncomePricing_hpp
#define FixedIncomePricing_hpp

#include <stdio.h>
#include "Matrix.hpp"
class FixedIncomePricing{
private:
    double r0;
    double r_bar;
    double sigma;
    double kai;
public:
    FixedIncomePricing(double r0, double r_bar, double sigma, double kai);
    ~FixedIncomePricing();
    shared_ptr<Matrix> Vasicek_Model(int seed, double T, double delta_t, int paths);
    shared_ptr<Matrix> Final_Model(int seed, double T,double delta_t, int paths, double alpha, double beta, double gamma);
    shared_ptr<Matrix> CIR_Model(int seed, double T, double delta_t, int paths);
    shared_ptr<Matrix>Gpp_Model(int seed, double T, double delta_t, int paths, double roh,double a, double b, double sig, double nah, double phi);
    double PureDiscountBondPricing(double FaceValue, double T, double delta_t, int paths,string type);
    double CoponPayingBondPricing(double FaceValue, double frequency, double Coupon, double T, double delta_t, int paths,string type);
    double EuroCallOptionOnPureDiscountedBondPricing(double FaceValue, double T, double delta_t, int paths,double Ti,double K,string type);
    double EuroPutOptionOnPureDiscountedBondPricing(double FaceValue, double T, double delta_t, int paths,double Ti,double K,string type);
    double EuroCallOptionOnCouponPayingBondPricing(double FaceValue, double frequency, double Coupon, double T, double delta_t, int paths, double Ti, double K,string type);
    double EuroCallOptionByVasicekFormula(double FaceValue, double T,double Ti,double K, double delta_t, int paths);
    double EuroCallOptionByCIRFormula(double FaceValue, double T, double Ti, double K);
};
#endif /* FixedIncomePricing_hpp */
