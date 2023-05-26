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
        VectorXd Yold;
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
