/*
 * @Author: Lyn
 * @Date: 2022-06-27 19:58:59
 * @LastEditors: error: git config user.name && git config user.email & please set dead value or install git
 * @LastEditTime: 2022-11-22 20:46:15
 * @FilePath: \test1tt\Jacobian.h
 * @Citrine
 */
#ifndef __Jacobian_H__
#define __Jacobian_H__

#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "MNQ.h"

class Jacobian
{
    private:
        VectorXd Yold;
        VectorXd Ynew;
        MatrixXd temp;

    public:
        Jacobian(VectorXd& arr, VectorXd& brr);
        ~Jacobian();

        int judg(double a); //判断符号正负
        
        MatrixXd jacobian(int k);

};



#endif