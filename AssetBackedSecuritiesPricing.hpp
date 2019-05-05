//
//  AssetBackedSecuritiesPricing.hpp
//  Computational_Methods_Final
//
//  Created by 黄佳铭 on 2019/3/21.
//  Copyright © 2019 Jiaming Huang. All rights reserved.
//

#ifndef AssetBackedSecuritiesPricing_hpp
#define AssetBackedSecuritiesPricing_hpp

#include <stdio.h>
class AssetBackedSecuritiesPricing{
private:
    double PrincipleValue;
    double maturity;
    double WAC;
    double r0;
    double r_bar;
    double kai;
    double sigma;
    double delta_t;
    int paths;
public:
    AssetBackedSecuritiesPricing(double PrincipleValue, double maturity, double WAC, double r0, double r_bar, double kai, double sigma, double delta_t,int paths);
    double NumerixPrepaymentModel(double R, double r_tm1_10, double t, double PV_tm1, double PV0);
    double PrepaymentPricing();
    double OASPrepaymentPricing(double deviation);
};
#endif /* AssetBackedSecuritiesPricing_hpp */
