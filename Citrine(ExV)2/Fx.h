
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