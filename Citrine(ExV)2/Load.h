/*
 * @Author: Lyn
 * @Date: 2022-05-26 21:59:38
 * @LastEditors: error: git config user.name && git config user.email & please set dead value or install git
 * @LastEditTime: 2022-11-22 20:46:07
 * @FilePath: \test1\Load.h
 * @Citrine
 */
#ifndef __LOAD_H__
#define __LOAD_H__
#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "Fx.h"
#include "Jacobian.h"

using namespace std;
using namespace Eigen;

class Load
{
private:
    VectorXd Yold;//用于存储临时Yold
    VectorXd Ynew;//用于存储临时Ynew


public:
    Load(VectorXd& arr, VectorXd& brr);
    ~Load();

    VectorXd LF(int k);
    MatrixXd LJ(int k);
  
    
};


#endif