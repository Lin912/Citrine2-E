/*
 * @Author: error: git config user.name && git config user.email & please set dead value or install git
 * @Date: 2022-10-25 21:09:02
 * @LastEditors: error: git config user.name && git config user.email & please set dead value or install git
 * @LastEditTime: 2022-11-30 22:50:56
 * @FilePath: \Siano2.1\Add.h
 * @Citrine
 */
#ifndef __ADD_H__
#define __ADD_H__
#include <iostream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

class Add
{
    private:
            VectorXd Yold;
            VectorXd Ynew;

    public:
            Add(VectorXd &arr, VectorXd &brr);
            ~Add();

            VectorXd Addyold(int k);
            VectorXd Addynew(int k);
};

#endif