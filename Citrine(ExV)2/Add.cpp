#include <iostream>
#include <Eigen/Dense>
#include "Add.h"
#include "EXf.h"

Add::Add(VectorXd &arr, VectorXd &brr)
{
    Yold = arr;
    Ynew = brr;

}

Add::~Add()
{

}


VectorXd Add::Addyold(int k)
{
    EXf a;
    vector<double> arr; 
    arr = a.Val(k, "./Data./input.csv");
	double V1 = arr[0];
	double V2 = arr[1];
	double V3 = arr[2];
    double Sd = arr[8];
	double G = arr.back();

    EXf c;
    vector<double> crr;
    crr = c.Val(1, "./Data./PhyChar.csv");
    double ma = crr[6];
    double Cdt = crr[7];
    double Cdn = crr[8];
    double Cdb = crr[9];
    double rho = crr[1];
    double d0 = crr[2];
    double pi = crr[11];



    VectorXd temp(500);

    VectorXd point00(3);
    point00(0) = V1*cos(Yold(7))*cos(Yold(6)) + V2*sin(Yold(6))*cos(Yold(7)) + V3*sin(Yold(7));            
    point00(1) = V2*cos(Yold(6)) - V1*sin(Yold(6));            
    point00(2) = - V1*sin(Yold(7))*cos(Yold(6)) - V2*sin(Yold(7))*sin(Yold(6)) + V3*cos(Yold(7));             

    VectorXd point01(2);
    point01(0) = 0;              
    point01(1) = 0;              

    VectorXd point02(3);
    point02(0) = - G*cos(Yold(497))*cos(Yold(496));           
    point02(1) = + G*sin(Yold(496)) - 0.5*rho*Cdn*Sd*Yold(491)*sqrt(pow(Yold(491),2)+pow(Yold(492),2));           
    point02(2) = + G*cos(Yold(496))*sin(Yold(497)) - 0.5*rho*Cdn*Sd*Yold(492)*sqrt(pow(Yold(491),2)+pow(Yold(492),2));          

    VectorXd point03(2);
    point03(0) = 0;             
    point03(1) = 0;            

    temp.segment(0,3) = point00;
    temp.segment(3,5) = Yold.segment(3, 5);
    temp.segment(8,2) = point01;


    temp.segment(10, 480) = Yold.segment(10, 480);

    temp.segment(490, 3) = Yold.segment(490, 3);
    temp.segment(493, 3) = point02;
    temp.segment(496, 2) = Yold.segment(496, 2);
    temp.tail(2) = point03;

    return temp;
}

VectorXd Add::Addynew(int k)
{
    EXf a;
    vector<double> arr; 
    arr = a.Val(k, "./Data./input.csv");
	double V1 = arr[0];
	double V2 = arr[1];
	double V3 = arr[2];
    double Sd = arr[8];
	double G = arr.back();

    EXf c;
    vector<double> crr;
    crr = c.Val(1, "./Data./PhyChar.csv");
    double ma = crr[6];
    double Cdt = crr[7];
    double Cdn = crr[8];
    double Cdb = crr[9];
    double rho = crr[1];
    double d0 = crr[2];
    double pi = crr[11];


    VectorXd temp(500);

    VectorXd point00(3);
    point00(0) = + V1*cos(Ynew(7))*cos(Ynew(6)) + V2*sin(Ynew(6))*cos(Ynew(7)) + V3*sin(Ynew(7));           
    point00(1) = + V2*cos(Ynew(6)) - V1*sin(Ynew(6));            //速度V边界条件
    point00(2) = - V1*sin(Ynew(7))*cos(Ynew(6)) - V2*sin(Ynew(7))*sin(Ynew(6)) + V3*cos(Ynew(7));              

    VectorXd point01(2);
    point01(0) = 0;            
    point01(1) = 0;           

    VectorXd point02(3);
    point02(0) = - G*cos(Ynew(497))*cos(Ynew(496));          
    point02(1) = + G*sin(Ynew(496)) - 0.5*rho*Cdn*Sd*Ynew(491)*sqrt(pow(Ynew(491),2)+pow(Ynew(492),2));        
    point02(2) = + G*cos(Ynew(496))*sin(Ynew(497)) - 0.5*rho*Cdn*Sd*Ynew(492)*sqrt(pow(Ynew(491),2)+pow(Ynew(492),2));            


    VectorXd point03(2);
    point03(0) = 0;              
    point03(1) = 0;             

    temp.segment(0,3) = point00;
    temp.segment(3,5) = Ynew.segment(3, 5);
    temp.segment(8,2) = point01;


    temp.segment(10, 480) = Ynew.segment(10, 480);

    temp.segment(490, 3) = Ynew.segment(490, 3);
    temp.segment(493, 3) = point02;
    temp.segment(496, 2) = Ynew.segment(496, 2);
    temp.tail(2) = point03;


    return temp;

}