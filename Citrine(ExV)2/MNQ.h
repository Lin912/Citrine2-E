/*
 * @Author: Lyn
 * @Date: 2022-06-01 21:49:17
 * @LastEditors: error: git config user.name && git config user.email & please set dead value or install git
 * @LastEditTime: 2022-11-22 21:51:25
 * @FilePath: \test1\MNQ.h
 * @Citrine
 */
#ifndef __MNQ_H__
#define __MNQ_H__
#include <iostream>
#include <fstream>
#include <cmath>
#include <Eigen/Dense>
#include "EXf.h"

using namespace std;
using namespace Eigen;

class MNQ
{
    private:
        VectorXd Yold;//用于存储临时Yold
        VectorXd Ynew;

    public:
        MNQ(VectorXd& arr, VectorXd& brr);
        ~MNQ();

        MatrixXd Mnew();
        MatrixXd Mold();
        MatrixXd Nold();
        MatrixXd Nnew();
        VectorXd Qold();
        VectorXd Qnew();

        void savetxt(Eigen::MatrixXd mat, string filename);

};



#endif
