/*
 * @Author:Lyn
 * @Date: 2022-05-31 18:56:32
 * @LastEditors: error: git config user.name && git config user.email & please set dead value or install git
 * @LastEditTime: 2022-11-22 20:46:24
 * @FilePath: \test1\Fx.h
 * @Citrine
 */

#ifndef __FX_H__
#define __FX_H__
#include <iostream>
#include <cmath>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

class Fx
{
    private:
        VectorXd Yold;
        VectorXd Ynew;
        VectorXd temp;

    public:
        Fx(VectorXd& arr, VectorXd& brr);
        ~Fx();

        VectorXd fx(int k);
};

#endif