/*
 * @Author: Lyn
 * @Date: 2022-06-01 21:49:23
 * @LastEditors: error: git config user.name && git config user.email & please set dead value or install git
 * @LastEditTime: 2022-11-30 21:06:59
 * @FilePath: \test1\MNQ.cpp
 * @Citrine
 */
#include <vector>
#include "MNQ.h"
#include "EXf.h"


MNQ::MNQ(VectorXd& arr, VectorXd& brr)
{
    Yold = arr;
    Ynew = brr;

}

MNQ::~MNQ()
{
    
}

MatrixXd MNQ::Mold()
{
    EXf Ev;
    vector<double> arr; 
    arr = Ev.Val(1, "./Data./PhyChar.csv");
	double A = arr[0];
	double rho = arr[1];
	double d0 = arr[2];
    double E = arr[3];
	double I = arr[4];


    //数学特性
    double pi = arr[11];
    double g = arr.back();

    //质量特性
    double M = arr[5];
    double ma = arr[6];
    double Cdt = arr[7];
    double Cdn = arr[8];
    double Cdb = arr[9];

    MatrixXd AA(500, 500);
    MatrixXd temp(500, 500);
    temp = AA.Zero(500, 500);

    for(int i = 0; i < 50; i++)
    {
        temp(10*i + 0, 10*i + 0) = M;
        temp(10*i + 0, 10*i + 6) = -M*Yold(i*10 + 1)*cos(Yold(i*10 + 7));
        temp(10*i + 0, 10*i + 7) = M*Yold(i*10 + 2);

        temp(10*i + 1, 10*i + 1) = M + ma;
        temp(10*i + 1, 10*i + 6) = M*(Yold(i*10 + 0)*cos(Yold(i*10 + 7)) + Yold(i*10 + 2)*sin(Yold(i*10 + 7)));

        temp(10*i + 2, 10*i + 2) = M + ma;
        temp(10*i + 2, 10*i + 6) = -M*Yold(i*10 + 1)*sin(Yold(i*10 + 7));
        temp(10*i + 2, 10*i + 7) = -M*Yold(i*10 + 0);

        temp(10*i + 5, 10*i + 3) = 1/(A*E);

        temp(10*i + 6, 10*i + 6) = cos(Yold(i*10 + 7))*(Yold(i*10 + 3)/(A*E) + 1);

        temp(10*i + 7, 10*i + 7) = Yold(i*10 + 3)/(A*E) + 1;
    }

    return temp;
}

MatrixXd MNQ::Mnew()
{

    EXf Ev;
    vector<double> arr; 
    arr = Ev.Val(1, "./Data./PhyChar.csv");
	double A = arr[0];
	double rho = arr[1];
	double d0 = arr[2];
    double E = arr[3];
	double I = arr[4];
   

    //数学特性
    double pi = arr[11];
    double g = arr.back();

    //质量特性
    double M = arr[5];
    double ma = arr[6];
    double Cdt = arr[7];
    double Cdn = arr[8];
    double Cdb = arr[9];


    MatrixXd AA(500, 500);
    MatrixXd temp(500, 500);
    temp = AA.Zero(500, 500);

    for(int i = 0; i < 50; i++)
    {

        temp(10*i + 0, 10*i + 0) = M;
        temp(10*i + 0, 10*i + 6) = -M*Ynew(i*10 + 1)*cos(Ynew(i*10 + 7));
        temp(10*i + 0, 10*i + 7) = M*Ynew(i*10 + 2);

        temp(10*i + 1, 10*i + 1) = M + ma;
        temp(10*i + 1, 10*i + 6) = M*(Ynew(i*10 + 0)*cos(Ynew(i*10 + 7)) + Ynew(i*10 + 2)*sin(Ynew(i*10 + 7)));

        temp(10*i + 2, 10*i + 2) = M + ma;
        temp(10*i + 2, 10*i + 6) = -M*Ynew(i*10 + 1)*sin(Ynew(i*10 + 7));
        temp(10*i + 2, 10*i + 7) = -M*Ynew(i*10 + 0);

        temp(10*i + 5, 10*i + 3) = 1/(A*E);

        temp(10*i + 6, 10*i + 6) = cos(Ynew(i*10 + 7))*(Ynew(i*10 + 3)/(A*E) + 1);

        temp(10*i + 7, 10*i + 7) = Ynew(i*10 + 3)/(A*E) + 1;
    }


    return temp;

}

MatrixXd MNQ::Nold()
{

    EXf Ev;
    vector<double> arr; 
    arr = Ev.Val(1, "./Data./PhyChar.csv");
	double A = arr[0];
	double rho = arr[1];
	double d0 = arr[2];
    double E = arr[3];
	double I = arr[4];

    
    double pi = arr[11];
    double g = arr.back();


    double M = arr[5];
    double ma = arr[6];
    double Cdt = arr[7];
    double Cdn = arr[8];
    double Cdb = arr[9];
    double w0 = arr[12];

    MatrixXd AA(500, 500);
    MatrixXd temp(500, 500);
    temp = AA.Zero(500, 500);

    for(int i = 0; i < 50; i++)
    {
        temp(i*10 + 0, i*10 + 3) = -1.0;
        temp(i*10 + 1, i*10 + 4) = -1.0;
        temp(i*10 + 2, i*10 + 5) = -1.0;
        temp(i*10 + 3, i*10 + 8) = E*I;
        temp(i*10 + 4, i*10 + 9) = E*I;
        temp(i*10 + 5, i*10 + 0) = -1.0;
        temp(i*10 + 6, i*10 + 1) = -1.0;
        temp(i*10 + 7, i*10 + 2) = 1.0;
        temp(i*10 + 8, i*10 + 7) = -1.0;
        temp(i*10 + 9, i*10 + 6) = -cos(Yold(i*10 + 7));
    }

    return temp;

}

MatrixXd MNQ::Nnew()
{

    EXf Ev;
    vector<double> arr; 
    arr = Ev.Val(1, "./Data./PhyChar.csv");
	double A = arr[0];
	double rho = arr[1];
	double d0 = arr[2];
    double E = arr[3];
	double I = arr[4];

    double pi = arr[11];
    double g = arr.back();

    double M = arr[5];
    double ma = arr[6];
    double Cdt = arr[7];
    double Cdn = arr[8];
    double Cdb = arr[9];
    double w0 = arr[12];

    MatrixXd AA(500, 500);
    MatrixXd temp(500, 500);
    temp = AA.Zero(500, 500);

    for(int i = 0; i < 50; i++)
    {   
        temp(i*10 + 0, i*10 + 3) = -1.0;
        temp(i*10 + 1, i*10 + 4) = -1.0;
        temp(i*10 + 2, i*10 + 5) = -1.0;
        temp(i*10 + 3, i*10 + 8) = E*I;
        temp(i*10 + 4, i*10 + 9) = E*I;
        temp(i*10 + 5, i*10 + 0) = -1.0;
        temp(i*10 + 6, i*10 + 1) = -1.0;
        temp(i*10 + 7, i*10 + 2) = 1.0;
        temp(i*10 + 8, i*10 + 7) = -1.0;
        temp(i*10 + 9, i*10 + 6) = -cos(Ynew(i*10 + 7));
    }

    return temp;
}

VectorXd MNQ::Qold()
{
   
    EXf Ev;
    vector<double> arr; 
    arr = Ev.Val(1, "./Data./PhyChar.csv");
	double A = arr[0];
	double rho = arr[1];
	double d0 = arr[2];
    double E = arr[3];
	double I = arr[4];
    double pi = arr[11];
    double g = arr.back();
    double M = arr[5];
    double ma = arr[6];
    double Cdt = arr[7];
    double Cdn = arr[8];
    double Cdb = arr[9];
    double w0 = arr[12];
    VectorXd temp(500);
    for(int i = 0; i < 50; i++)
    {
        temp(i*10 + 0) = Yold(i*10 + 9)*Yold(i*10 + 4) - Yold(i*10 + 8)*Yold(i*10 + 5) + w0*cos(Yold(i*10 + 7))*cos(Yold(i*10 + 6)) + 0.5*Cdt*d0*rho*Yold(i*10 + 0)*pi*fabs(Yold(i*10 + 0))*sqrt(Yold(i*10 + 3)/(A*E) + 1);
        temp(i*10 + 1) = -w0*sin(Yold(i*10 + 6)) - Yold(i*10 + 9)*Yold(i*10 + 3) - Yold(i*10 + 9)*Yold(i*10 + 5)*tan(Yold(i*10 + 7)) + 0.5*Cdn*d0*rho*Yold(i*10 + 1)*sqrt(pow(Yold(i*10 + 1),2) + pow(Yold(i*10 + 2),2))*sqrt(Yold(i*10 + 3)/(A*E) + 1);
        temp(i*10 + 2) =  Yold(i*10 + 8)*Yold(i*10 + 3) + Yold(i*10 + 9)*Yold(i*10 + 4)*tan(Yold(i*10 + 7)) - w0*cos(Yold(i*10 + 6))*sin(Yold(i*10 + 7)) + 0.5*Cdb*d0*rho*Yold(i*10 + 2)*sqrt(pow(Yold(i*10 + 1),2) + pow(Yold(i*10 + 2),2))*sqrt(Yold(i*10 + 3)/(A*E) + 1);
        temp(i*10 + 3) = E*I*pow(Yold(i*10 + 9),2)*tan(Yold(i*10 + 7)) - Yold(i*10 + 5)*pow((Yold(i*10 + 3)/(A*E) + 1),3);
        temp(i*10 + 4) = Yold(i*10 + 4)*pow((Yold(i*10 + 3)/(A*E) + 1),3) - E*I*Yold(i*10 + 8)*Yold(i*10 + 9)*tan(Yold(i*10 + 7));
        temp(i*10 + 5) = Yold(i*10 + 9)*Yold(i*10 + 1) - Yold(i*10 + 8)*Yold(i*10 + 2);
        temp(i*10 + 6) = -Yold(i*10 + 9)*(Yold(i*10 + 0) + Yold(i*10 + 2)*tan(Yold(i*10 + 7)));
        temp(i*10 + 7) = -Yold(i*10 + 8)*Yold(i*10 + 0) - Yold(i*10 + 9)*Yold(i*10 + 1)*tan(Yold(i*10 + 7));
        temp(i*10 + 8) = Yold(i*10 + 8);
        temp(i*10 + 9) = Yold(i*10 + 9);
    }
    return temp;   
}

VectorXd MNQ::Qnew()
{

    EXf Ev;
    vector<double> arr; 
    arr = Ev.Val(1, "./Data./PhyChar.csv");
	double A = arr[0];
	double rho = arr[1];
	double d0 = arr[2];
    double E = arr[3];
	double I = arr[4];
    double pi = arr[11];
    double g = arr.back();
    double M = arr[5];
    double ma = arr[6];
    double Cdt = arr[7];
    double Cdn = arr[8];
    double Cdb = arr[9];
    double w0 = arr[12];
    VectorXd temp(500);
   for(int i = 0; i < 50; i++)
    {
        temp(i*10 + 0) = Ynew(i*10 + 9)*Ynew(i*10 + 4) - Ynew(i*10 + 8)*Ynew(i*10 + 5) + w0*cos(Ynew(i*10 + 7))*cos(Ynew(i*10 + 6)) + 0.5*Cdt*d0*rho*Ynew(i*10 + 0)*pi*fabs(Ynew(i*10 + 0))*sqrt(Ynew(i*10 + 3)/(A*E) + 1);
        temp(i*10 + 1) = -w0*sin(Ynew(i*10 + 6)) - Ynew(i*10 + 9)*Ynew(i*10 + 3) - Ynew(i*10 + 9)*Ynew(i*10 + 5)*tan(Ynew(i*10 + 7)) + 0.5*Cdn*d0*rho*Ynew(i*10 + 1)*sqrt(pow(Ynew(i*10 + 1),2) + pow(Ynew(i*10 + 2),2))*sqrt(Ynew(i*10 + 3)/(A*E) + 1);
        temp(i*10 + 2) =  Ynew(i*10 + 8)*Ynew(i*10 + 3) + Ynew(i*10 + 9)*Ynew(i*10 + 4)*tan(Ynew(i*10 + 7)) - w0*cos(Ynew(i*10 + 6))*sin(Ynew(i*10 + 7)) + 0.5*Cdb*d0*rho*Ynew(i*10 + 2)*sqrt(pow(Ynew(i*10 + 1),2) + pow(Ynew(i*10 + 2),2))*sqrt(Ynew(i*10 + 3)/(A*E) + 1);
        temp(i*10 + 3) = E*I*pow(Ynew(i*10 + 9),2)*tan(Ynew(i*10 + 7)) - Ynew(i*10 + 5)*pow((Ynew(i*10 + 3)/(A*E) + 1),3);
        temp(i*10 + 4) = Ynew(i*10 + 4)*pow((Ynew(i*10 + 3)/(A*E) + 1),3) - E*I*Ynew(i*10 + 8)*Ynew(i*10 + 9)*tan(Ynew(i*10 + 7));
        temp(i*10 + 5) = Ynew(i*10 + 9)*Ynew(i*10 + 1) - Ynew(i*10 + 8)*Ynew(i*10 + 2);
        temp(i*10 + 6) = -Ynew(i*10 + 9)*(Ynew(i*10 + 0) + Ynew(i*10 + 2)*tan(Ynew(i*10 + 7)));
        temp(i*10 + 7) = -Ynew(i*10 + 8)*Ynew(i*10 + 0) - Ynew(i*10 + 9)*Ynew(i*10 + 1)*tan(Ynew(i*10 + 7));
        temp(i*10 + 8) = Ynew(i*10 + 8);
        temp(i*10 + 9) = Ynew(i*10 + 9);
    }
    return temp;

    
}

void MNQ::savetxt(Eigen::MatrixXd mat, string filename)
{
    ofstream outfile(filename, ios::trunc);
    outfile << mat;
    outfile.close();
}