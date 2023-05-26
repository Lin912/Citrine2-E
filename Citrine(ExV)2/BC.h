/*
 * @Author: Lyn
 * @Date: 2022-06-16 18:13:22
 * @LastEditors: error: git config user.name && git config user.email & please set dead value or install git
 * @LastEditTime: 2022-11-22 20:46:37
 * @FilePath: \test1\BC.h
 * @Citrine
 */
#ifndef __BC_H__
#define __BC_H__
#include <iostream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

class BC
{
    private:
            VectorXd Yold;
            VectorXd Ynew;

    public:
            BC(VectorXd &arr, VectorXd &brr);
            ~BC();

            VectorXd yold();
            VectorXd ynew();
};

#endif