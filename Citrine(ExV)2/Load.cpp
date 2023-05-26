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

    // cout << "Read in Yold And Ynew !!" << endl;
}

Load::~Load()
{
    // cout << "Load done! Time to get Iteration !!" << endl;

}

VectorXd Load::LF(int k)
{
    Fx A(Yold, Ynew);
    return A.fx(k);//读入两个向量并设置f(x)


}

MatrixXd Load::LJ(int k)
{
    Jacobian A(Yold, Ynew);
    return A.jacobian(k);//读入两个向量并设置Jacobi
    
}