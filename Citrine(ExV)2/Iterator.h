/*
 * @Author: error: Lyn
 * @Date: 2022-05-26 21:05:14
 * @LastEditors: error: git config user.name && git config user.email & please set dead value or install git
 * @LastEditTime: 2022-11-22 20:46:21
 * @FilePath: \test1\Iterator.h
 * @Citrine
 */
#ifndef __ITERATOR_H__
#define __ITERATOR_H__
#include <iostream>
#include <fstream>
#include <cmath>
#include <Eigen/Dense>
#include "Load.h"
#include "BC.h"
#include "Add.h"


using namespace std;
using namespace Eigen;

class Iterator
{
        private:
            int times;//Iteartion Times
		    double Error;//Convergence Error

            VectorXd fx;//类存储fx
            MatrixXd jac;//类存储jac

            VectorXd Yold;//临时存储Yold
            VectorXd Ynew;//临时存储Ynew

        public:
            Iterator (VectorXd& arr, VectorXd& brr, int a, double b);
            ~Iterator ();

            void begin (int k);
            VectorXd out();

            void savetxt(Eigen::MatrixXd mat, string filename);
};





#endif