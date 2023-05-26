/*
 * @Author: Lyn
 * @Date: 2022-06-27 19:59:06
 * @LastEditors: error: git config user.name && git config user.email & please set dead value or install git
 * @LastEditTime: 2022-12-09 18:02:19
 * @FilePath: \test1tt\Jacobian.cpp
 * @Citrine
 */
#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "Jacobian.h"
#include "MNQ.h"
#include <vector>
#include "EXf.h"


Jacobian::Jacobian(VectorXd& arr, VectorXd& brr)
{
    Ynew = arr;
    Yold = brr;

}

Jacobian::~Jacobian()
{
    // cout << "Jacobian is read in !!" << endl;

}

MatrixXd Jacobian::jacobian(int k)
{

    //控制步长
    EXf b;
    vector<double> brr;
    brr = b.Val(0, "./Data./delta.csv");
    double deltaS = brr[0];
    double deltaT = brr.back();


    EXf Ev;
    vector<double> crr; 
    crr = Ev.Val(1, "./Data./PhyChar.csv");
	double A = crr[0];//截面面积 A
	double rho = crr[1];//介质密度
	double d0 = crr[2];//带缆直径 d0
    double E = crr[3];//弹性模量 E
	double I = crr[4];
    // double MNQ::Ip = 2.80e-7;
    // double MNQ::G = 8.11e9;

    //数学特性
    double pi = crr[11];
    double g = crr.back();

    //质量特性
    double M = crr[5];//质量 M = Rhocable*A
    double ma = crr[6];//附加质量 ma = rho*A*(Cm)  此时Cm = 1.00
    double Cdt = crr[7];
    double Cdn = crr[8];
    double Cdb = crr[9];//待修改
    double w0 = crr[12];


    EXf a;
    vector<double> arr; 
    arr = a.Val(k, "./Data./input.csv");
	double V1 = arr[0];             //Velocity of water
	double V2 = arr[1];
	double V3 = arr[2];
    double Sd = arr[8];
	double G = arr.back();


    MatrixXd AA(500, 500);
    MatrixXd temp(500, 500);
    temp = AA.Zero(500, 500);
    
    //BC01
    temp(0, 0) = 1;
    temp(0, 1) = 0;
    temp(0, 2) = 0;
    temp(0, 3) = 0;
    temp(0, 4) = 0;
    temp(0, 5) = 0;
    temp(0, 6) = -V2*cos(Ynew(6))*cos(Ynew(7)) + V1*sin(Ynew(6))*cos(Ynew(7));
    temp(0, 7) = V1*cos(Ynew(6))*sin(Ynew(7)) - V3*cos(Ynew(7)) + V2*sin(Ynew(6))*sin(Ynew(7));
    temp(0, 8) = 0;
    temp(0, 9) = 0;

    temp(1, 0) = 0;
    temp(1, 1) = 1;
    temp(1, 2) = 0;
    temp(1, 3) = 0;
    temp(1, 4) = 0;
    temp(1, 5) = 0;
    temp(1, 6) = V1*cos(Ynew(6)) + V2*sin(Ynew(6));
    temp(1, 7) = 0;
    temp(1, 8) = 0;
    temp(1, 9) = 0;

    temp(2, 0) = 0;
    temp(2, 1) = 0;
    temp(2, 2) = 1;
    temp(2, 3) = 0;
    temp(2, 4) = 0;
    temp(2, 5) = 0;
    temp(2, 6) = V2*cos(Ynew(6))*sin(Ynew(7)) - V1*sin(Ynew(6))*sin(Ynew(7));
    temp(2, 7) = V3*sin(Ynew(7)) + V1*cos(Ynew(6))*cos(Ynew(7)) + V2*cos(Ynew(7))*sin(Ynew(6));
    temp(2, 8) = 0;
    temp(2, 9) = 0;

    temp(3, 0) = 0;
    temp(3, 1) = 0;
    temp(3, 2) = 0;
    temp(3, 3) = 0;
    temp(3, 4) = 0;
    temp(3, 5) = 0;
    temp(3, 6) = 0;
    temp(3, 7) = 0;
    temp(3, 8) = 1;
    temp(3, 9) = 0;

    temp(4, 0) = 0;
    temp(4, 1) = 0;
    temp(4, 2) = 0;
    temp(4, 3) = 0;
    temp(4, 4) = 0;
    temp(4, 5) = 0;
    temp(4, 6) = 0;
    temp(4, 7) = 0;
    temp(4, 8) = 0;
    temp(4, 9) = 1;
    

    //BC49
    temp(495, 490) = 0;
    temp(495, 491) = 0;
    temp(495, 492) = 0;
    temp(495, 493) = 1;
    temp(495, 494) = 0;
    temp(495, 495) = 0;
    temp(495, 496) = -G*sin(Ynew(496))*cos(Ynew(497));
    temp(495, 497) = -G*cos(Ynew(496))*sin(Ynew(497));
    temp(495, 498) = 0;
    temp(495, 499) = 0;

    temp(496, 490) = 0;
    temp(496, 491) = (Sd*Cdn*rho*sqrt(pow(Ynew(491),2) + pow(Ynew(492),2)))/2 + (Sd*Cdn*rho*pow(Ynew(491),2))/(2*sqrt(pow(Ynew(491),2) + pow(Ynew(492),2)));
    temp(496, 492) = (Sd*Cdn*rho*Ynew(491)*Ynew(492))/(2*sqrt(pow(Ynew(491),2) + pow(Ynew(492),2)));
    temp(496, 493) = 0;
    temp(496, 494) = 1;
    temp(496, 495) = 0;
    temp(496, 496) = -G*cos(Ynew(496));
    temp(496, 497) = 0;
    temp(496, 498) = 0;
    temp(496, 499) = 0;

    temp(497, 490) = 0;
    temp(497, 491) = (Sd*Cdn*rho*Ynew(491)*Ynew(492))/(2*sqrt(pow(Ynew(491),2) + pow(Ynew(492),2)));
    temp(497, 492) = (Sd*Cdn*rho*(2*sqrt(pow(Ynew(491),2) + pow(Ynew(492),2))))/2 + (Sd*Cdn*rho*pow(Ynew(492),2))/(2*(2*sqrt(pow(Ynew(491),2) + pow(Ynew(492),2))));
    temp(497, 493) = 0;
    temp(497, 494) = 0;
    temp(497, 495) = 1;
    temp(497, 496) = G*sin(Ynew(496))*sin(Ynew(497));
    temp(497, 497) = -G*cos(Ynew(497))*cos(Ynew(496));
    temp(497, 498) = 0;
    temp(497, 499) = 0;

    temp(498, 490) = 0;
    temp(498, 491) = 0;
    temp(498, 492) = 0;
    temp(498, 493) = 0;
    temp(498, 494) = 0;
    temp(498, 495) = 0;
    temp(498, 496) = 0;
    temp(498, 497) = 0;
    temp(498, 498) = 1;
    temp(498, 499) = 0;

    temp(499, 490) = 0;
    temp(499, 491) = 0;
    temp(499, 492) = 0;
    temp(499, 493) = 0;
    temp(499, 494) = 0;
    temp(499, 495) = 0;
    temp(499, 496) = 0;
    temp(499, 497) = 0;
    temp(499, 498) = 0;
    temp(499, 499) = 1;   


    

    //Block01    
    for(int i = 0; i < 49; i++)
    {
        temp(i*10 + 5, i*10 + 0) = (2*M)*deltaS + (0.5*Cdt*d0*rho*pi*fabs(Ynew(i*10 + 0))*sqrt(Ynew(i*10 + 3)/(A*E) + 1) + 0.5*Cdt*d0*rho*Ynew(i*10 + 0)*pi*judg(Ynew(i*10 + 0))*sqrt(Ynew(i*10 + 3)/(A*E) + 1))*deltaS*deltaT;
        temp(i*10 + 6, i*10 + 0) = (M*cos(Ynew(i*10 + 7))*(Yold(i*10 + 6) - Ynew(i*10 + 6)))*deltaS;
        temp(i*10 + 7, i*10 + 0) = (M*(Yold(i*10 + 7) - Ynew(i*10 + 7)))*deltaS;
        temp(i*10 + 10, i*10 + 0) = 2*deltaT;
        temp(i*10 + 11, i*10 + 0) = -Ynew(i*10 + 9)*deltaT*deltaS;
        temp(i*10 + 12, i*10 + 0) = -Ynew(i*10 + 8)*deltaT*deltaS;

        temp(i*10 + 5, i*10 + 1) = (M*cos(Ynew(i*10 + 7))*(Yold(i*10 + 6) - Ynew(i*10 + 6)))*deltaS;
        temp(i*10 + 6, i*10 + 1) = (2*M + 2*ma)*deltaS + deltaS*deltaT*(0.5*Cdn*d0*rho*sqrt(pow(Ynew(i*10 + 1),2) + pow(Ynew(i*10 + 2),2))*sqrt(Ynew(i*10 + 3)/(A*E) + 1) + (0.5*Cdn*d0*rho*pow(Ynew(i*10 + 1),2)*sqrt(Ynew(i*10 + 3)/(A*E) + 1))/sqrt(pow(Ynew(i*10 + 1),2) + pow(Ynew(i*10 + 2),2)));
        temp(i*10 + 7, i*10 + 1) = (M*sin(Ynew(i*10 + 7))*(Yold(i*10 + 6) - Ynew(i*10 + 6)))*deltaS + (0.5*Cdb*d0*rho*deltaS*deltaT*Ynew(i*10 + 1)*Ynew(i*10 + 2)*sqrt(Ynew(i*10 + 3)/(A*E) + 1))/sqrt(pow(Ynew(i*10 + 1),2) + pow(Ynew(i*10 + 2),2));
        temp(i*10 + 10, i*10 + 1) = Ynew(i*10 + 9)*deltaS*deltaT;
        temp(i*10 + 11, i*10 + 1) = 2*deltaT;
        temp(i*10 + 12, i*10 + 1) = -Ynew(i*10 + 9)*tan(Ynew(i*10 + 7))*deltaS*deltaT;

        temp(i*10 + 5, i*10 + 2) = -(M*(Yold(i*10 + 7) - Ynew(i*10 + 7)))*deltaS;
        temp(i*10 + 6, i*10 + 2) = (M*sin(Ynew(i*10 + 7))*(Yold(i*10 + 6) - Ynew(i*10 + 6)))*deltaS + (0.5*Cdn*d0*deltaS*deltaT*rho*Ynew(i*10 + 1)*Ynew(i*10 + 2)*sqrt(Ynew(i*10 + 3)/(A*E) + 1))/sqrt(pow(Ynew(i*10 + 1),2) + pow(Ynew(i*10 + 2),2));
        temp(i*10 + 7, i*10 + 2) = (2*M + 2*ma)*deltaS + deltaS*deltaT*(0.5*Cdb*d0*rho*sqrt(pow(Ynew(i*10 + 1),2) + pow(Ynew(i*10 + 2),2))*sqrt(Ynew(i*10 + 3)/(A*E) + 1) + (0.5*Cdb*d0*rho*pow(Ynew(i*10 + 2),2)*sqrt(Ynew(i*10 + 3)/(A*E) + 1))/sqrt(pow(Ynew(i*10 + 1),2) + pow(Ynew(i*10 + 2),2)));
        temp(i*10 + 10, i*10 + 2) = deltaS*deltaT*(-Ynew(i*10 + 8));
        temp(i*10 + 11, i*10 + 2) = deltaS*deltaT*(-Ynew(i*10 + 9)*tan(Ynew(i*10 + 7)));
        temp(i*10 + 12, i*10 + 2) = -2*deltaT;

        temp(i*10 + 5, i*10 + 3) = 2*deltaT + (0.25*Cdt*d0*deltaS*deltaT*rho*Ynew(i*10 + 0)*pi*fabs(Ynew(i*10 + 0)))/(A*E*sqrt(Ynew(i*10 + 3)/(A*E) + 1));
        temp(i*10 + 6, i*10 + 3) = -deltaS*deltaT*(-(0.25*Cdn*d0*rho*Ynew(i*10 + 1)*sqrt(pow(Ynew(i*10 + 1),2) + pow(Ynew(i*10 + 2),2)))/(A*E*sqrt(Ynew(i*10 + 3)/(A*E) + 1)) + Ynew(i*10 + 9));
        temp(i*10 + 7, i*10 + 3) = deltaS*deltaT*(Ynew(i*10 + 8) + (0.25*Cdb*d0*rho*Ynew(i*10 + 2)*sqrt(pow(Ynew(i*10 + 1),2) + pow(Ynew(i*10 + 2),2)))/(A*E*sqrt(Ynew(i*10 + 3)/(A*E) + 1)));
        temp(i*10 + 8, i*10 + 3) = -(3*deltaS*deltaT*Ynew(i*10 + 5)*pow((Ynew(i*10 + 3)/(A*E) + 1),2))/(A*E);
        temp(i*10 + 9, i*10 + 3) = (3*deltaS*deltaT*Ynew(i*10 + 4)*pow((Ynew(i*10 + 3)/(A*E) + 1),2))/(A*E);
        temp(i*10 + 10, i*10 + 3) = 2*deltaS/(A*E);
        temp(i*10 + 11, i*10 + 3) = -((cos(Ynew(i*10 + 7))*(Yold(i*10 + 6) - Ynew(i*10 + 6)))/(A*E*deltaT))*deltaS*deltaT;
        temp(i*10 + 12, i*10 + 3) = -((Yold(i*10 + 7) - Ynew(i*10 + 7))/(A*E*deltaT))*deltaS*deltaT;

        temp(i*10 + 5, i*10 + 4) = Ynew(i*10 + 9)*deltaS*deltaT;
        temp(i*10 + 6, i*10 + 4) = 2*deltaT;
        temp(i*10 + 7, i*10 + 4) = Ynew(i*10 + 9)*tan(Ynew(i*10 + 7))*deltaS*deltaT;
        temp(i*10 + 9, i*10 + 4) = pow((Ynew(i*10 + 3)/(A*E) + 1),3)*deltaS*deltaT;

        temp(i*10 + 5, i*10 + 5) = -Ynew(i*10 + 8)*deltaS*deltaT;
        temp(i*10 + 6, i*10 + 5) = -Ynew(i*10 + 9)*tan(Ynew(i*10 + 7))*deltaS*deltaT;
        temp(i*10 + 7, i*10 + 5) = 2*deltaT;
        temp(i*10 + 8, i*10 + 5) = -pow((Ynew(i*10 + 3)/(A*E) + 1),3)*deltaS*deltaT;

        temp(i*10 + 5, i*10 + 6) = -w0*sin(Ynew(i*10 + 6))*cos(Ynew(i*10 + 7))*deltaS*deltaT - (M*Ynew(i*10 + 1)*cos(Ynew(i*10 + 7)) + M*Yold(i*10 + 1)*cos(Yold(i*10 + 7)))*deltaS;;
        temp(i*10 + 6, i*10 + 6) = - (M*(Ynew(i*10 + 0)*cos(Ynew(i*10 + 7)) + Ynew(i*10 + 2)*sin(Ynew(i*10 + 7))) + M*(Yold(i*10 + 0)*cos(Yold(i*10 + 7)) + Yold(i*10 + 2)*sin(Yold(i*10 + 7))))*deltaS-w0*cos(Ynew(i*10 + 6))*deltaS*deltaT;
        temp(i*10 + 7, i*10 + 6) = w0*sin(Ynew(i*10 + 6))*sin(Ynew(i*10 + 7))*deltaS*deltaT - (M*Ynew(i*10 + 1)*sin(Ynew(i*10 + 7)) + M*Yold(i*10 + 1)*sin(Yold(i*10 + 7)))*deltaS;
        temp(i*10 + 11, i*10 + 6) = (cos(Ynew(i*10 + 7))*(Ynew(i*10 + 3)/(A*E) + 1) + cos(Yold(i*10 + 7))*(Yold(i*10 + 3)/(A*E) + 1))*deltaS;
        temp(i*10 + 14, i*10 + 6) = (cos(Ynew(i*10 + 7)) + cos(Ynew(i*10 + 17)))*deltaT;

        temp(i*10 + 5, i*10 + 7) = (M*Ynew(i*10 + 2) + M*Yold(i*10 + 2) - M*Ynew(i*10 + 1)*sin(Ynew(i*10 + 7))*(Yold(i*10 + 6) - Ynew(i*10 + 6)))*deltaS-w0*deltaS*deltaT*cos(Ynew(i*10 + 6))*sin(Ynew(i*10 + 7));
        temp(i*10 + 6, i*10 + 7) = (M*(Yold(i*10 + 6) - Ynew(i*10 + 6))*(Ynew(i*10 + 2)*cos(Ynew(i*10 + 7)) - Ynew(i*10 + 0)*sin(Ynew(i*10 + 7))))*deltaS - Ynew(i*10 + 9)*Ynew(i*10 + 5)*deltaS*deltaT*(pow(tan(Ynew(i*10 + 7)),2) + 1);
        temp(i*10 + 7, i*10 + 7) = (Ynew(i*10 + 9)*Ynew(i*10 + 4)*(pow(tan(Ynew(i*10 + 7)),2) + 1))*deltaS*deltaT - (M*Ynew(i*10 + 0) + M*Yold(i*10 + 0) - M*Ynew(i*10 + 1)*cos(Ynew(i*10 + 7))*(Yold(i*10 + 6) - Ynew(i*10 + 6)))*deltaS - (w0*cos(Ynew(i*10 + 7))*cos(Ynew(i*10 + 6)))*deltaS*deltaT;
        temp(i*10 + 8, i*10 + 7) = E*I*pow(Ynew(i*10 + 9),2)*(pow(tan(Ynew(i*10 + 7)),2) + 1)*deltaS*deltaT;
        temp(i*10 + 9, i*10 + 7) = -E*I*Ynew(i*10 + 8)*Ynew(i*10 + 9)*(pow(tan(Ynew(i*10 + 7)),2) + 1)*deltaS*deltaT;
        temp(i*10 + 11, i*10 + 7) = (sin(Ynew(i*10 + 7))*(Ynew(i*10 + 3)/(A*E) + 1)*(Yold(i*10 + 6) - Ynew(i*10 + 6)))*deltaS-Ynew(i*10 + 9)*Ynew(i*10 + 2)*(pow(tan(Ynew(i*10 + 7)),2) + 1)*deltaS*deltaT;
        temp(i*10 + 12, i*10 + 7) = (Ynew(i*10 + 3)/(A*E) + Yold(i*10 + 3)/(A*E) + 2)*deltaS -Ynew(i*10 + 9)*Ynew(i*10 + 1)*(pow(tan(Ynew(i*10 + 7)),2) + 1)*deltaS*deltaT;
        temp(i*10 + 13, i*10 + 7) = 2*deltaT;
        temp(i*10 + 14, i*10 + 7) = -(sin(Ynew(i*10 + 7))*(Ynew(i*10 + 6) - Ynew(i*10 + 16)))*deltaT;

        temp(i*10 + 5, i*10 + 8) = -Ynew(i*10 + 5)*deltaS*deltaT;
        temp(i*10 + 7, i*10 + 8) = Ynew(i*10 + 3)*deltaS*deltaT;
        temp(i*10 + 8, i*10 + 8) = -(2*E*I)*deltaT;
        temp(i*10 + 9, i*10 + 8) = -E*I*Ynew(i*10 + 9)*tan(Ynew(i*10 + 7))*deltaS*deltaT;
        temp(i*10 + 10, i*10 + 8) = -Ynew(i*10 + 2)*deltaS*deltaT;
        temp(i*10 + 12, i*10 + 8) = -Ynew(i*10 + 0)*deltaS*deltaT;
        temp(i*10 + 13, i*10 + 8) = deltaS*deltaT;

        temp(i*10 + 5, i*10 + 9) = Ynew(i*10 + 4)*deltaS*deltaT;
        temp(i*10 + 6, i*10 + 9) = (- Ynew(i*10 + 3) - Ynew(i*10 + 5)*tan(Ynew(i*10 + 7)))*deltaS*deltaT;
        temp(i*10 + 7, i*10 + 9) = Ynew(i*10 + 4)*tan(Ynew(i*10 + 7))*deltaS*deltaT;
        temp(i*10 + 8, i*10 + 9) = 2*E*I*Ynew(i*10 + 9)*tan(Ynew(i*10 + 7))*deltaS*deltaT;
        temp(i*10 + 9, i*10 + 9) = - (2*E*I)*deltaT - E*I*Ynew(i*10 + 8)*tan(Ynew(i*10 + 7))*deltaS*deltaT;
        temp(i*10 + 10, i*10 + 9) = Ynew(i*10 + 1)*deltaS*deltaT;
        temp(i*10 + 11, i*10 + 9) = (- Ynew(i*10 + 0) - Ynew(i*10 + 2)*tan(Ynew(i*10 + 7)))*deltaS*deltaT;
        temp(i*10 + 12, i*10 + 9) = -Ynew(i*10 + 1)*tan(Ynew(i*10 + 7))*deltaS*deltaT;
        temp(i*10 + 14, i*10 + 9) = deltaS*deltaT;

    }

    //Block02
    for(int i = 0; i < 49; i++)
    {
        temp(i*10 + 5, i*10 + 10) = (2*M)*deltaS + (0.5*Cdt*d0*rho*pi*fabs(Ynew(i*10 + 10))*sqrt(Ynew(i*10 + 13)/(A*E) + 1) + 0.5*Cdt*d0*rho*Ynew(i*10 + 10)*pi*judg(Ynew(i*10 + 10))*sqrt(Ynew(i*10 + 13)/(A*E) + 1))*deltaS*deltaT;
        temp(i*10 + 6, i*10 + 10) = (M*cos(Ynew(i*10 + 17))*(Yold(i*10 + 16) - Ynew(i*10 + 16)))*deltaS;
        temp(i*10 + 7, i*10 + 10) = (M*(Yold(i*10 + 17) - Ynew(i*10 + 17)))*deltaS;
        temp(i*10 + 10, i*10 + 10) = -2*deltaT;
        temp(i*10 + 11, i*10 + 10) = -Ynew(i*10 + 19)*deltaS*deltaT;
        temp(i*10 + 12, i*10 + 10) = -Ynew(i*10 + 18)*deltaS*deltaT;

        temp(i*10 + 5, i*10 + 11) = (M*cos(Ynew(i*10 + 17))*(Yold(i*10 + 16) - Ynew(i*10 + 16)))*deltaS; 
        temp(i*10 + 6, i*10 + 11) = (2*M + 2*ma)*deltaS + (0.5*Cdn*d0*rho*sqrt(pow(Ynew(i*10 + 11),2) + pow(Ynew(i*10 + 12),2))*sqrt(Ynew(i*10 + 13)/(A*E) + 1) + (0.5*Cdn*d0*rho*pow(Ynew(i*10 + 11),2)*sqrt(Ynew(i*10 + 13)/(A*E) + 1))/sqrt(pow(Ynew(i*10 + 11),2) + pow(Ynew(i*10 + 12),2)))*deltaS*deltaT;
        temp(i*10 + 7, i*10 + 11) = (M*sin(Ynew(i*10 + 17))*(Yold(i*10 + 16) - Ynew(i*10 + 16)))*deltaS + (0.5*Cdb*deltaS*deltaT*d0*rho*Ynew(i*10 + 11)*Ynew(i*10 + 12)*sqrt(Ynew(i*10 + 13)/(A*E) + 1))/sqrt(pow(Ynew(i*10 + 11),2) + pow(Ynew(i*10 + 12),2));
        temp(i*10 + 10, i*10 + 11) = Ynew(i*10 + 19)*deltaS*deltaT;
        temp(i*10 + 11, i*10 + 11) = -2*deltaT;
        temp(i*10 + 12, i*10 + 11) = -Ynew(i*10 + 19)*tan(Ynew(i*10 + 17))*deltaS*deltaT;

        temp(i*10 + 5, i*10 + 12) = -(M*(Yold(i*10 + 17) - Ynew(i*10 + 17)))*deltaS;
        temp(i*10 + 6, i*10 + 12) = (M*sin(Ynew(i*10 + 17))*(Yold(i*10 + 16) - Ynew(i*10 + 16)))*deltaS + (0.5*Cdn*d0*deltaS*deltaT*rho*Ynew(i*10 + 11)*Ynew(i*10 + 12)*sqrt(Ynew(i*10 + 13)/(A*E) + 1))/sqrt(pow(Ynew(i*10 + 11),2) + pow(Ynew(i*10 + 12),2));
        temp(i*10 + 7, i*10 + 12) = (2*M + 2*ma)*deltaS + (0.5*Cdb*d0*rho*sqrt(pow(Ynew(i*10 + 11),2) + pow(Ynew(i*10 + 12),2))*sqrt(Ynew(i*10 + 13)/(A*E) + 1) + (0.5*Cdb*d0*rho*pow(Ynew(i*10 + 12),2)*sqrt(Ynew(i*10 + 13)/(A*E) + 1))/sqrt(pow(Ynew(i*10 + 11),2) + pow(Ynew(i*10 + 12),2)))*deltaS*deltaT;
        temp(i*10 + 10, i*10 + 12) = -Ynew(i*10 + 18)*deltaS*deltaT;
        temp(i*10 + 11, i*10 + 12) = -Ynew(i*10 + 19)*tan(Ynew(i*10 + 17))*deltaS*deltaT;
        temp(i*10 + 12, i*10 + 12) = 2*deltaT;

        temp(i*10 + 5, i*10 + 13) =  (0.25*Cdt*d0*deltaS*deltaT*rho*Ynew(i*10 + 10)*pi*fabs(Ynew(i*10 + 10)))/(A*E*sqrt(Ynew(i*10 + 13)/(A*E) + 1)) - 2*deltaT;
        temp(i*10 + 6, i*10 + 13) = ((0.25*Cdn*d0*rho*Ynew(i*10 + 11)*sqrt(pow(Ynew(i*10 + 11),2) + pow(Ynew(i*10 + 12),2)))/(A*E*sqrt(Ynew(i*10 + 13)/(A*E) + 1)) - Ynew(i*10 + 19))*deltaS*deltaT;
        temp(i*10 + 7, i*10 + 13) = (Ynew(i*10 + 18) + (0.25*Cdb*d0*rho*Ynew(i*10 + 12)*sqrt(pow(Ynew(i*10 + 11),2) + pow(Ynew(i*10 + 12),2)))/(A*E*sqrt(Ynew(i*10 + 13)/(A*E) + 1)))*deltaS*deltaT;
        temp(i*10 + 8, i*10 + 13) = -(3*Ynew(i*10 + 15)*deltaS*deltaT*pow((Ynew(i*10 + 13)/(A*E) + 1),2))/(A*E);
        temp(i*10 + 9, i*10 + 13) = (3*Ynew(i*10 + 14)*deltaS*deltaT*pow((Ynew(i*10 + 13)/(A*E) + 1),2))/(A*E);
        temp(i*10 + 10, i*10 + 13) = 2*deltaS/(A*E);
        temp(i*10 + 11, i*10 + 13) =-((cos(Ynew(i*10 + 17))*(Yold(i*10 + 16) - Ynew(i*10 + 16)))*deltaS)/(A*E);
        temp(i*10 + 12, i*10 + 13) =-((Yold(i*10 + 17) - Ynew(i*10 + 17))*deltaS)/(A*E);

        temp(i*10 + 5, i*10 + 14) = Ynew(i*10 + 19)*deltaT;
        temp(i*10 + 6, i*10 + 14) = -2*deltaT;
        temp(i*10 + 7, i*10 + 14) = Ynew(i*10 + 19)*tan(Ynew(i*10 + 17))*deltaS*deltaT;
        temp(i*10 + 9, i*10 + 14) = pow((Ynew(i*10 + 13)/(A*E) + 1),3)*deltaS*deltaT;

        temp(i*10 + 5, i*10 + 15) = -Ynew(i*10 + 18)*deltaS*deltaT;
        temp(i*10 + 6, i*10 + 15) = -Ynew(i*10 + 19)*tan(Ynew(i*10 + 17))*deltaS*deltaT;
        temp(i*10 + 7, i*10 + 15) = -2*deltaT;
        temp(i*10 + 8, i*10 + 15) = -pow((Ynew(i*10 + 13)/(A*E) + 1),3)*deltaS*deltaT;
        

        temp(i*10 + 5, i*10 + 16) = -w0*sin(Ynew(i*10 + 16))*cos(Ynew(i*10 + 17))*deltaS*deltaT- (M*Ynew(i*10 + 11)*cos(Ynew(i*10 + 17)) + M*Yold(i*10 + 11)*cos(Yold(i*10 + 17)))*deltaS;
        temp(i*10 + 6, i*10 + 16) = - (M*(Ynew(i*10 + 10)*cos(Ynew(i*10 + 17)) + Ynew(i*10 + 12)*sin(Ynew(i*10 + 17))) + M*(Yold(i*10 + 10)*cos(Yold(i*10 + 17)) + Yold(i*10 + 12)*sin(Yold(i*10 + 17))))*deltaS -w0*cos(Ynew(i*10 + 16))*deltaS*deltaT;
        temp(i*10 + 7, i*10 + 16) = w0*sin(Ynew(i*10 + 16))*sin(Ynew(i*10 + 17))*deltaS*deltaT- (M*Ynew(i*10 + 11)*sin(Ynew(i*10 + 17)) + M*Yold(i*10 + 11)*sin(Yold(i*10 + 17)))*deltaS;
        temp(i*10 + 11, i*10 + 16) = (cos(Ynew(i*10 + 17))*(Ynew(i*10 + 13)/(A*E) + 1) + cos(Yold(i*10 + 17))*(Yold(i*10 + 13)/(A*E) + 1))*deltaS;
        temp(i*10 + 14, i*10 + 16) = -(cos(Ynew(i*10 + 7)) + cos(Ynew(i*10 + 17)))*deltaT;
        

        temp(i*10 + 5, i*10 + 17) = (M*Ynew(i*10 + 12) + M*Yold(i*10 + 12) - M*Ynew(i*10 + 11)*sin(Ynew(i*10 + 17))*(Yold(i*10 + 16) - Ynew(i*10 + 16)))*deltaS-w0*cos(Ynew(i*10 + 16))*sin(Ynew(i*10 + 17))*deltaS*deltaT;
        temp(i*10 + 6, i*10 + 17) = (M*(Yold(i*10 + 16) - Ynew(i*10 + 16))*(Ynew(i*10 + 12)*cos(Ynew(i*10 + 17)) - Ynew(i*10 + 10)*sin(Ynew(i*10 + 17))))*deltaS-(Ynew(i*10 + 19)*Ynew(i*10 + 15)*(pow(tan(Ynew(i*10 + 17)),2) + 1))*deltaS*deltaT;
        temp(i*10 + 7, i*10 + 17) = (Ynew(i*10 + 19)*Ynew(i*10 + 14)*(pow(tan(Ynew(i*10 + 17)),2) + 1))*deltaS*deltaT- (M*Ynew(i*10 + 10) + M*Yold(i*10 + 10) - M*Ynew(i*10 + 11)*cos(Ynew(i*10 + 17))*(Yold(i*10 + 16) - Ynew(i*10 + 16)))*deltaS - (w0*cos(Ynew(i*10 + 17))*cos(Ynew(i*10 + 16)))*deltaS*deltaT;
        temp(i*10 + 8, i*10 + 17) = E*I*pow(Ynew(i*10 + 19),2)*(pow(tan(Ynew(i*10 + 17)),2) + 1)*deltaS*deltaT;
        temp(i*10 + 9, i*10 + 17) = -E*I*Ynew(i*10 + 18)*Ynew(i*10 + 19)*(pow(tan(Ynew(i*10 + 17)),2) + 1)*deltaS*deltaT;
        temp(i*10 + 11, i*10 + 17) = (sin(Ynew(i*10 + 17))*(Ynew(i*10 + 13)/(A*E) + 1)*(Yold(i*10 + 16) - Ynew(i*10 + 16)))*deltaS-(Ynew(i*10 + 19)*Ynew(i*10 + 12)*(pow(tan(Ynew(i*10 + 17)),2) + 1))*deltaS*deltaT;
        temp(i*10 + 12, i*10 + 17) = (Ynew(i*10 + 13)/(A*E) + Yold(i*10 + 13)/(A*E) + 2)*deltaS -(Ynew(i*10 + 19)*Ynew(i*10 + 11)*(pow(tan(Ynew(i*10 + 17)),2) + 1))*deltaS*deltaT;
        temp(i*10 + 13, i*10 + 17) = -2*deltaT;
        temp(i*10 + 14, i*10 + 17) = -(sin(Ynew(i*10 + 17))*(Ynew(i*10 + 6) - Ynew(i*10 + 16)))*deltaT;

        temp(i*10 + 5, i*10 + 18) = -Ynew(i*10 + 15)*deltaS*deltaT;
        temp(i*10 + 7, i*10 + 18) = Ynew(i*10 + 13)*deltaS*deltaT;
        temp(i*10 + 8, i*10 + 18) = (2*E*I)*deltaT;
        temp(i*10 + 9, i*10 + 18) = -E*I*Ynew(i*10 + 19)*tan(Ynew(i*10 + 17))*deltaS*deltaT;
        temp(i*10 + 10, i*10 + 18) = -Ynew(i*10 + 12)*deltaS*deltaT;
        temp(i*10 + 12, i*10 + 18) = -Ynew(i*10 + 10)*deltaS*deltaT;
        temp(i*10 + 13, i*10 + 18) = deltaS*deltaT;

        temp(i*10 + 5, i*10 + 19) = Ynew(i*10 + 14)*deltaS*deltaT;
        temp(i*10 + 6, i*10 + 19) = (- Ynew(i*10 + 13) - Ynew(i*10 + 15)*tan(Ynew(i*10 + 17)))*deltaS*deltaT;
        temp(i*10 + 7, i*10 + 19) = (Ynew(i*10 + 14)*tan(Ynew(i*10 + 17)))*deltaS*deltaT;
        temp(i*10 + 8, i*10 + 19) = 2*E*I*Ynew(i*10 + 19)*tan(Ynew(i*10 + 17))*deltaS*deltaT;
        temp(i*10 + 9, i*10 + 19) = (2*E*I)*deltaT - (E*I*Ynew(i*10 + 18)*tan(Ynew(i*10 + 17)))*deltaS*deltaT;
        temp(i*10 + 10, i*10 + 19) = Ynew(i*10 + 11)*deltaS*deltaT;
        temp(i*10 + 11, i*10 + 19) = (- Ynew(i*10 + 10) - Ynew(i*10 + 12)*tan(Ynew(i*10 + 17)))*deltaS*deltaT;
        temp(i*10 + 12, i*10 + 19) = -Ynew(i*10 + 11)*tan(Ynew(i*10 + 17))*deltaS*deltaT;
        temp(i*10 + 14, i*10 + 19) = deltaS*deltaT;

    }


    //Orgin JacobianMatrix 

    return temp;

}

int Jacobian::judg(double a)
{
    if(a > 0)
    {
        return 1;
    }
    else
    {
        if(a < 0)
        {
            return -1;
        }
        else
        {
            return 0;
        }
    }
}