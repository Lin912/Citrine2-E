/*
 * @Author: Lyn
 * @Date: 2022-05-31 16:39:53
 * @LastEditors: error: git config user.name && git config user.email & please set dead value or install git
 * @LastEditTime: 2022-11-22 20:46:12
 * @FilePath: \test1\Load.cpp
 * @Citrine
 */
#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "Load.h"

Load::Load(VectorXd& arr, VectorXd& brr)
{
    Yold = arr;
    Ynew = brr;

}

Load::~Load()
{

}

VectorXd Load::LF(int k)
{
    Fx A(Yold, Ynew);
    return A.fx(k);


}

MatrixXd Load::LJ(int k)
{
    Jacobian A(Yold, Ynew);
    return A.jacobian(k);
    
}