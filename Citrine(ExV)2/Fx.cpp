/*
 * @Author: Lyn
 * @Date: 2022-05-31 18:56:38
 * @LastEditors: error: git config user.name && git config user.email & please set dead value or install git
 * @LastEditTime: 2022-12-09 18:02:49
 * @FilePath: \test1\Fx.cpp
 * @Citrine
 */
#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "Fx.h"
#include "MNQ.h"
#include "EXf.h"
#include <vector>



Fx::Fx(VectorXd& arr, VectorXd& brr)
{
    Yold = arr;
    Ynew = brr;

}

Fx::~Fx()
{
    //  cout << "Function is read in!!" << endl;

}


VectorXd Fx::fx(int k)
{
    EXf a;
    vector<double> arr; 
    arr = a.Val(k, "./Data./input.csv");
	double V1 = arr[0];
	double V2 = arr[1];
	double V3 = arr[2];
    double Sd = arr[8];
	double G = arr.back();

    EXf b;
    vector<double> brr;
    brr = b.Val(0, "./Data./delta.csv");
    double deltaS = brr[0];
    double deltaT = brr.back();
    
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


    VectorXd Y0old(10);
    VectorXd Y1old(10);
    VectorXd Y2old(10);
    VectorXd Y3old(10);
    VectorXd Y4old(10);
    VectorXd Y5old(10);
    VectorXd Y6old(10);
    VectorXd Y7old(10);
    VectorXd Y8old(10);
    VectorXd Y9old(10);
    VectorXd Y10old(10);
    VectorXd Y11old(10);
    VectorXd Y12old(10);
    VectorXd Y13old(10);
    VectorXd Y14old(10);
    VectorXd Y15old(10);
    VectorXd Y16old(10);
    VectorXd Y17old(10);
    VectorXd Y18old(10);
    VectorXd Y19old(10);
    VectorXd Y20old(10);
    VectorXd Y21old(10);
    VectorXd Y22old(10);
    VectorXd Y23old(10);
    VectorXd Y24old(10);
    VectorXd Y25old(10);
    VectorXd Y26old(10);
    VectorXd Y27old(10);
    VectorXd Y28old(10);
    VectorXd Y29old(10);
    VectorXd Y30old(10);
    VectorXd Y31old(10);
    VectorXd Y32old(10);
    VectorXd Y33old(10);
    VectorXd Y34old(10);
    VectorXd Y35old(10);
    VectorXd Y36old(10);
    VectorXd Y37old(10);
    VectorXd Y38old(10);
    VectorXd Y39old(10);
    VectorXd Y40old(10);
    VectorXd Y41old(10);
    VectorXd Y42old(10);
    VectorXd Y43old(10);
    VectorXd Y44old(10);
    VectorXd Y45old(10);
    VectorXd Y46old(10);
    VectorXd Y47old(10);
    VectorXd Y48old(10);
    VectorXd Y49old(10);


    VectorXd Y0new(10);
    VectorXd Y1new(10);
    VectorXd Y2new(10);
    VectorXd Y3new(10);
    VectorXd Y4new(10);
    VectorXd Y5new(10);
    VectorXd Y6new(10);
    VectorXd Y7new(10);
    VectorXd Y8new(10);
    VectorXd Y9new(10);
    VectorXd Y10new(10);
    VectorXd Y11new(10);
    VectorXd Y12new(10);
    VectorXd Y13new(10);
    VectorXd Y14new(10);
    VectorXd Y15new(10);
    VectorXd Y16new(10);
    VectorXd Y17new(10);
    VectorXd Y18new(10);
    VectorXd Y19new(10);
    VectorXd Y20new(10);
    VectorXd Y21new(10);
    VectorXd Y22new(10);
    VectorXd Y23new(10);
    VectorXd Y24new(10);
    VectorXd Y25new(10);
    VectorXd Y26new(10);
    VectorXd Y27new(10);
    VectorXd Y28new(10);
    VectorXd Y29new(10);
    VectorXd Y30new(10);
    VectorXd Y31new(10);
    VectorXd Y32new(10);
    VectorXd Y33new(10);
    VectorXd Y34new(10);
    VectorXd Y35new(10);
    VectorXd Y36new(10);
    VectorXd Y37new(10);
    VectorXd Y38new(10);
    VectorXd Y39new(10);
    VectorXd Y40new(10);
    VectorXd Y41new(10);
    VectorXd Y42new(10);
    VectorXd Y43new(10);
    VectorXd Y44new(10);
    VectorXd Y45new(10);
    VectorXd Y46new(10);
    VectorXd Y47new(10);
    VectorXd Y48new(10);
    VectorXd Y49new(10);


    Y0old = Yold.head(10);
    Y1old = Yold.segment(10, 10);
    Y2old = Yold.segment(20, 10);
    Y3old = Yold.segment(30, 10);
    Y4old = Yold.segment(40, 10);
    Y5old = Yold.segment(50, 10);
    Y6old = Yold.segment(60, 10);
    Y7old = Yold.segment(70, 10);
    Y8old = Yold.segment(80, 10);
    Y9old = Yold.segment(90, 10);
    Y10old = Yold.segment(100, 10);
    Y11old = Yold.segment(110, 10);
    Y12old = Yold.segment(120, 10);
    Y13old = Yold.segment(130, 10);
    Y14old = Yold.segment(140, 10);
    Y15old = Yold.segment(150, 10);
    Y16old = Yold.segment(160, 10);
    Y17old = Yold.segment(170, 10);
    Y18old = Yold.segment(180, 10);
    Y19old = Yold.segment(190, 10);
    Y20old = Yold.segment(200, 10);
    Y21old = Yold.segment(210, 10);
    Y22old = Yold.segment(220, 10);
    Y23old = Yold.segment(230, 10);
    Y24old = Yold.segment(240, 10);
    Y25old = Yold.segment(250, 10);
    Y26old = Yold.segment(260, 10);
    Y27old = Yold.segment(270, 10);
    Y28old = Yold.segment(280, 10);
    Y29old = Yold.segment(290, 10);
    Y30old = Yold.segment(300, 10);
    Y31old = Yold.segment(310, 10);
    Y32old = Yold.segment(320, 10);
    Y33old = Yold.segment(330, 10);
    Y34old = Yold.segment(340, 10);
    Y35old = Yold.segment(350, 10);
    Y36old = Yold.segment(360, 10);
    Y37old = Yold.segment(370, 10);
    Y38old = Yold.segment(380, 10);
    Y39old = Yold.segment(390, 10);
    Y40old = Yold.segment(400, 10);
    Y41old = Yold.segment(410, 10);
    Y42old = Yold.segment(420, 10);
    Y43old = Yold.segment(430, 10);
    Y44old = Yold.segment(440, 10);
    Y45old = Yold.segment(450, 10);
    Y46old = Yold.segment(460, 10);
    Y47old = Yold.segment(470, 10);
    Y48old = Yold.segment(480, 10);
    Y49old = Yold.tail(10);

    Y0new = Ynew.head(10);
    Y1new = Ynew.segment(10, 10);
    Y2new = Ynew.segment(20, 10);
    Y3new = Ynew.segment(30, 10);
    Y4new = Ynew.segment(40, 10);
    Y5new = Ynew.segment(50, 10);
    Y6new = Ynew.segment(60, 10);
    Y7new = Ynew.segment(70, 10);
    Y8new = Ynew.segment(80, 10);
    Y9new = Ynew.segment(90, 10);
    Y10new = Ynew.segment(100, 10);
    Y11new = Ynew.segment(110, 10);
    Y12new = Ynew.segment(120, 10);
    Y13new = Ynew.segment(130, 10);
    Y14new = Ynew.segment(140, 10);
    Y15new = Ynew.segment(150, 10);
    Y16new = Ynew.segment(160, 10);
    Y17new = Ynew.segment(170, 10);
    Y18new = Ynew.segment(180, 10);
    Y19new = Ynew.segment(190, 10);
    Y20new = Ynew.segment(200, 10);
    Y21new = Ynew.segment(210, 10);
    Y22new = Ynew.segment(220, 10);
    Y23new = Ynew.segment(230, 10);
    Y24new = Ynew.segment(240, 10);
    Y25new = Ynew.segment(250, 10);
    Y26new = Ynew.segment(260, 10);
    Y27new = Ynew.segment(270, 10);
    Y28new = Ynew.segment(280, 10);
    Y29new = Ynew.segment(290, 10);
    Y30new = Ynew.segment(300, 10);
    Y31new = Ynew.segment(310, 10);
    Y32new = Ynew.segment(320, 10);
    Y33new = Ynew.segment(330, 10);
    Y34new = Ynew.segment(340, 10);
    Y35new = Ynew.segment(350, 10);
    Y36new = Ynew.segment(360, 10);
    Y37new = Ynew.segment(370, 10);
    Y38new = Ynew.segment(380, 10);
    Y39new = Ynew.segment(390, 10);
    Y40new = Ynew.segment(400, 10);
    Y41new = Ynew.segment(410, 10);
    Y42new = Ynew.segment(420, 10);
    Y43new = Ynew.segment(430, 10);
    Y44new = Ynew.segment(440, 10);
    Y45new = Ynew.segment(450, 10);
    Y46new = Ynew.segment(460, 10);
    Y47new = Ynew.segment(470, 10);
    Y48new = Ynew.segment(480, 10);
    Y49new = Ynew.tail(10);


    MNQ AA(Yold, Ynew);
    MatrixXd M0old(10, 10);
    MatrixXd M1old(10, 10);
    MatrixXd M2old(10, 10);
    MatrixXd M3old(10, 10);
    MatrixXd M4old(10, 10);
    MatrixXd M5old(10, 10);
    MatrixXd M6old(10, 10);
    MatrixXd M7old(10, 10);
    MatrixXd M8old(10, 10);
    MatrixXd M9old(10, 10);
    MatrixXd M10old(10, 10);
    MatrixXd M11old(10, 10);
    MatrixXd M12old(10, 10);
    MatrixXd M13old(10, 10);
    MatrixXd M14old(10, 10);
    MatrixXd M15old(10, 10);
    MatrixXd M16old(10, 10);
    MatrixXd M17old(10, 10);
    MatrixXd M18old(10, 10);
    MatrixXd M19old(10, 10);
    MatrixXd M20old(10, 10);
    MatrixXd M21old(10, 10);
    MatrixXd M22old(10, 10);
    MatrixXd M23old(10, 10);
    MatrixXd M24old(10, 10);
    MatrixXd M25old(10, 10);
    MatrixXd M26old(10, 10);
    MatrixXd M27old(10, 10);
    MatrixXd M28old(10, 10);
    MatrixXd M29old(10, 10);
    MatrixXd M30old(10, 10);
    MatrixXd M31old(10, 10);
    MatrixXd M32old(10, 10);
    MatrixXd M33old(10, 10);
    MatrixXd M34old(10, 10);
    MatrixXd M35old(10, 10);
    MatrixXd M36old(10, 10);
    MatrixXd M37old(10, 10);
    MatrixXd M38old(10, 10);
    MatrixXd M39old(10, 10);
    MatrixXd M40old(10, 10);
    MatrixXd M41old(10, 10);
    MatrixXd M42old(10, 10);
    MatrixXd M43old(10, 10);
    MatrixXd M44old(10, 10);
    MatrixXd M45old(10, 10);
    MatrixXd M46old(10, 10);
    MatrixXd M47old(10, 10);
    MatrixXd M48old(10, 10);
    MatrixXd M49old(10, 10);



    MatrixXd M0new(10, 10);
    MatrixXd M1new(10, 10);
    MatrixXd M2new(10, 10);
    MatrixXd M3new(10, 10);
    MatrixXd M4new(10, 10);
    MatrixXd M5new(10, 10);
    MatrixXd M6new(10, 10);
    MatrixXd M7new(10, 10);
    MatrixXd M8new(10, 10);
    MatrixXd M9new(10, 10);
    MatrixXd M10new(10, 10);
    MatrixXd M11new(10, 10);
    MatrixXd M12new(10, 10);
    MatrixXd M13new(10, 10);
    MatrixXd M14new(10, 10);
    MatrixXd M15new(10, 10);
    MatrixXd M16new(10, 10);
    MatrixXd M17new(10, 10);
    MatrixXd M18new(10, 10);
    MatrixXd M19new(10, 10);
    MatrixXd M20new(10, 10);
    MatrixXd M21new(10, 10);
    MatrixXd M22new(10, 10);
    MatrixXd M23new(10, 10);
    MatrixXd M24new(10, 10);
    MatrixXd M25new(10, 10);
    MatrixXd M26new(10, 10);
    MatrixXd M27new(10, 10);
    MatrixXd M28new(10, 10);
    MatrixXd M29new(10, 10);
    MatrixXd M30new(10, 10);
    MatrixXd M31new(10, 10);
    MatrixXd M32new(10, 10);
    MatrixXd M33new(10, 10);
    MatrixXd M34new(10, 10);
    MatrixXd M35new(10, 10);
    MatrixXd M36new(10, 10);
    MatrixXd M37new(10, 10);
    MatrixXd M38new(10, 10);
    MatrixXd M39new(10, 10);
    MatrixXd M40new(10, 10);
    MatrixXd M41new(10, 10);
    MatrixXd M42new(10, 10);
    MatrixXd M43new(10, 10);
    MatrixXd M44new(10, 10);
    MatrixXd M45new(10, 10);
    MatrixXd M46new(10, 10);
    MatrixXd M47new(10, 10);
    MatrixXd M48new(10, 10);
    MatrixXd M49new(10, 10);//矩阵分块


    M0old = AA.Mold().block(0, 0, 10, 10);//分块赋值
    M1old = AA.Mold().block(10, 10, 10, 10);
    M2old = AA.Mold().block(20, 20, 10, 10);
    M3old = AA.Mold().block(30, 30, 10, 10);
    M4old = AA.Mold().block(40, 40, 10, 10);
    M5old = AA.Mold().block(50, 50, 10, 10);
    M6old = AA.Mold().block(60, 60, 10, 10);
    M7old = AA.Mold().block(70, 70, 10, 10);
    M8old = AA.Mold().block(80, 80, 10, 10);
    M9old = AA.Mold().block(90, 90, 10, 10);
    M10old = AA.Mold().block(100, 100, 10, 10);//分块赋值
    M11old = AA.Mold().block(110, 110, 10, 10);
    M12old = AA.Mold().block(120, 120, 10, 10);
    M13old = AA.Mold().block(130, 130, 10, 10);
    M14old = AA.Mold().block(140, 140, 10, 10);
    M15old = AA.Mold().block(150, 150, 10, 10);
    M16old = AA.Mold().block(160, 160, 10, 10);
    M17old = AA.Mold().block(170, 170, 10, 10);
    M18old = AA.Mold().block(180, 180, 10, 10);
    M19old = AA.Mold().block(190, 190, 10, 10);
    M20old = AA.Mold().block(200, 200, 10, 10);//分块赋值
    M21old = AA.Mold().block(210, 210, 10, 10);
    M22old = AA.Mold().block(220, 220, 10, 10);
    M23old = AA.Mold().block(230, 230, 10, 10);
    M24old = AA.Mold().block(240, 240, 10, 10);
    M25old = AA.Mold().block(250, 250, 10, 10);
    M26old = AA.Mold().block(260, 260, 10, 10);
    M27old = AA.Mold().block(270, 270, 10, 10);
    M28old = AA.Mold().block(280, 280, 10, 10);
    M29old = AA.Mold().block(290, 290, 10, 10);
    M30old = AA.Mold().block(300, 300, 10, 10);//分块赋值
    M31old = AA.Mold().block(310, 310, 10, 10);
    M32old = AA.Mold().block(320, 320, 10, 10);
    M33old = AA.Mold().block(330, 330, 10, 10);
    M34old = AA.Mold().block(340, 340, 10, 10);
    M35old = AA.Mold().block(350, 350, 10, 10);
    M36old = AA.Mold().block(360, 360, 10, 10);
    M37old = AA.Mold().block(370, 370, 10, 10);
    M38old = AA.Mold().block(380, 380, 10, 10);
    M39old = AA.Mold().block(390, 390, 10, 10);
    M40old = AA.Mold().block(400, 400, 10, 10);//分块赋值
    M41old = AA.Mold().block(410, 410, 10, 10);
    M42old = AA.Mold().block(420, 420, 10, 10);
    M43old = AA.Mold().block(430, 430, 10, 10);
    M44old = AA.Mold().block(440, 440, 10, 10);
    M45old = AA.Mold().block(450, 450, 10, 10);
    M46old = AA.Mold().block(460, 460, 10, 10);
    M47old = AA.Mold().block(470, 470, 10, 10);
    M48old = AA.Mold().block(480, 480, 10, 10);
    M49old = AA.Mold().block(490, 490, 10, 10);


    M0new = AA.Mnew().block(0, 0, 10, 10);//分块赋值
    M1new = AA.Mnew().block(10, 10, 10, 10);
    M2new = AA.Mnew().block(20, 20, 10, 10);
    M3new = AA.Mnew().block(30, 30, 10, 10);
    M4new = AA.Mnew().block(40, 40, 10, 10);
    M5new = AA.Mnew().block(50, 50, 10, 10);
    M6new = AA.Mnew().block(60, 60, 10, 10);
    M7new = AA.Mnew().block(70, 70, 10, 10);
    M8new = AA.Mnew().block(80, 80, 10, 10);
    M9new = AA.Mnew().block(90, 90, 10, 10);
    M10new = AA.Mnew().block(100, 100, 10, 10);//分块赋值
    M11new = AA.Mnew().block(110, 110, 10, 10);
    M12new = AA.Mnew().block(120, 120, 10, 10);
    M13new = AA.Mnew().block(130, 130, 10, 10);
    M14new = AA.Mnew().block(140, 140, 10, 10);
    M15new = AA.Mnew().block(150, 150, 10, 10);
    M16new = AA.Mnew().block(160, 160, 10, 10);
    M17new = AA.Mnew().block(170, 170, 10, 10);
    M18new = AA.Mnew().block(180, 180, 10, 10);
    M19new = AA.Mnew().block(190, 190, 10, 10);
    M20new = AA.Mnew().block(200, 200, 10, 10);//分块赋值
    M21new = AA.Mnew().block(210, 210, 10, 10);
    M22new = AA.Mnew().block(220, 220, 10, 10);
    M23new = AA.Mnew().block(230, 230, 10, 10);
    M24new = AA.Mnew().block(240, 240, 10, 10);
    M25new = AA.Mnew().block(250, 250, 10, 10);
    M26new = AA.Mnew().block(260, 260, 10, 10);
    M27new = AA.Mnew().block(270, 270, 10, 10);
    M28new = AA.Mnew().block(280, 280, 10, 10);
    M29new = AA.Mnew().block(290, 290, 10, 10);
    M30new = AA.Mnew().block(300, 300, 10, 10);//分块赋值
    M31new = AA.Mnew().block(310, 310, 10, 10);
    M32new = AA.Mnew().block(320, 320, 10, 10);
    M33new = AA.Mnew().block(330, 330, 10, 10);
    M34new = AA.Mnew().block(340, 340, 10, 10);
    M35new = AA.Mnew().block(350, 350, 10, 10);
    M36new = AA.Mnew().block(360, 360, 10, 10);
    M37new = AA.Mnew().block(370, 370, 10, 10);
    M38new = AA.Mnew().block(380, 380, 10, 10);
    M39new = AA.Mnew().block(390, 390, 10, 10);
    M40new = AA.Mnew().block(400, 400, 10, 10);//分块赋值
    M41new = AA.Mnew().block(410, 410, 10, 10);
    M42new = AA.Mnew().block(420, 420, 10, 10);
    M43new = AA.Mnew().block(430, 430, 10, 10);
    M44new = AA.Mnew().block(440, 440, 10, 10);
    M45new = AA.Mnew().block(450, 450, 10, 10);
    M46new = AA.Mnew().block(460, 460, 10, 10);
    M47new = AA.Mnew().block(470, 470, 10, 10);
    M48new = AA.Mnew().block(480, 480, 10, 10);
    M49new = AA.Mnew().block(490, 490, 10, 10);



    MatrixXd N0old(10, 10);
    MatrixXd N1old(10, 10);
    MatrixXd N2old(10, 10);
    MatrixXd N3old(10, 10);
    MatrixXd N4old(10, 10);
    MatrixXd N5old(10, 10);
    MatrixXd N6old(10, 10);
    MatrixXd N7old(10, 10);
    MatrixXd N8old(10, 10);
    MatrixXd N9old(10, 10);
    MatrixXd N10old(10, 10);
    MatrixXd N11old(10, 10);
    MatrixXd N12old(10, 10);
    MatrixXd N13old(10, 10);
    MatrixXd N14old(10, 10);
    MatrixXd N15old(10, 10);
    MatrixXd N16old(10, 10);
    MatrixXd N17old(10, 10);
    MatrixXd N18old(10, 10);
    MatrixXd N19old(10, 10);
    MatrixXd N20old(10, 10);
    MatrixXd N21old(10, 10);
    MatrixXd N22old(10, 10);
    MatrixXd N23old(10, 10);
    MatrixXd N24old(10, 10);
    MatrixXd N25old(10, 10);
    MatrixXd N26old(10, 10);
    MatrixXd N27old(10, 10);
    MatrixXd N28old(10, 10);
    MatrixXd N29old(10, 10);
    MatrixXd N30old(10, 10);
    MatrixXd N31old(10, 10);
    MatrixXd N32old(10, 10);
    MatrixXd N33old(10, 10);
    MatrixXd N34old(10, 10);
    MatrixXd N35old(10, 10);
    MatrixXd N36old(10, 10);
    MatrixXd N37old(10, 10);
    MatrixXd N38old(10, 10);
    MatrixXd N39old(10, 10);
    MatrixXd N40old(10, 10);
    MatrixXd N41old(10, 10);
    MatrixXd N42old(10, 10);
    MatrixXd N43old(10, 10);
    MatrixXd N44old(10, 10);
    MatrixXd N45old(10, 10);
    MatrixXd N46old(10, 10);
    MatrixXd N47old(10, 10);
    MatrixXd N48old(10, 10);
    MatrixXd N49old(10, 10);


    MatrixXd N0new(10, 10);
    MatrixXd N1new(10, 10);
    MatrixXd N2new(10, 10);
    MatrixXd N3new(10, 10);
    MatrixXd N4new(10, 10);
    MatrixXd N5new(10, 10);
    MatrixXd N6new(10, 10);
    MatrixXd N7new(10, 10);
    MatrixXd N8new(10, 10);
    MatrixXd N9new(10, 10);//矩阵分块
    MatrixXd N10new(10, 10);
    MatrixXd N11new(10, 10);
    MatrixXd N12new(10, 10);
    MatrixXd N13new(10, 10);
    MatrixXd N14new(10, 10);
    MatrixXd N15new(10, 10);
    MatrixXd N16new(10, 10);
    MatrixXd N17new(10, 10);
    MatrixXd N18new(10, 10);
    MatrixXd N19new(10, 10);//矩阵分块
    MatrixXd N20new(10, 10);
    MatrixXd N21new(10, 10);
    MatrixXd N22new(10, 10);
    MatrixXd N23new(10, 10);
    MatrixXd N24new(10, 10);
    MatrixXd N25new(10, 10);
    MatrixXd N26new(10, 10);
    MatrixXd N27new(10, 10);
    MatrixXd N28new(10, 10);
    MatrixXd N29new(10, 10);//矩阵分块
    MatrixXd N30new(10, 10);
    MatrixXd N31new(10, 10);
    MatrixXd N32new(10, 10);
    MatrixXd N33new(10, 10);
    MatrixXd N34new(10, 10);
    MatrixXd N35new(10, 10);
    MatrixXd N36new(10, 10);
    MatrixXd N37new(10, 10);
    MatrixXd N38new(10, 10);
    MatrixXd N39new(10, 10);//矩阵分块
    MatrixXd N40new(10, 10);
    MatrixXd N41new(10, 10);
    MatrixXd N42new(10, 10);
    MatrixXd N43new(10, 10);
    MatrixXd N44new(10, 10);
    MatrixXd N45new(10, 10);
    MatrixXd N46new(10, 10);
    MatrixXd N47new(10, 10);
    MatrixXd N48new(10, 10);
    MatrixXd N49new(10, 10);//矩阵分块


    N0old = AA.Nold().block(0, 0, 10, 10);//分块赋值
    N1old = AA.Nold().block(10, 10, 10, 10);
    N2old = AA.Nold().block(20, 20, 10, 10);
    N3old = AA.Nold().block(30, 30, 10, 10);
    N4old = AA.Nold().block(40, 40, 10, 10);
    N5old = AA.Nold().block(50, 50, 10, 10);
    N6old = AA.Nold().block(60, 60, 10, 10);
    N7old = AA.Nold().block(70, 70, 10, 10);
    N8old = AA.Nold().block(80, 80, 10, 10);
    N9old = AA.Nold().block(90, 90, 10, 10);
    N10old = AA.Nold().block(100, 100, 10, 10);//分块赋值
    N11old = AA.Nold().block(110, 110, 10, 10);
    N12old = AA.Nold().block(120, 120, 10, 10);
    N13old = AA.Nold().block(130, 130, 10, 10);
    N14old = AA.Nold().block(140, 140, 10, 10);
    N15old = AA.Nold().block(150, 150, 10, 10);
    N16old = AA.Nold().block(160, 160, 10, 10);
    N17old = AA.Nold().block(170, 170, 10, 10);
    N18old = AA.Nold().block(180, 180, 10, 10);
    N19old = AA.Nold().block(190, 190, 10, 10);
    N20old = AA.Nold().block(200, 200, 10, 10);//分块赋值
    N21old = AA.Nold().block(210, 210, 10, 10);
    N22old = AA.Nold().block(220, 220, 10, 10);
    N23old = AA.Nold().block(230, 230, 10, 10);
    N24old = AA.Nold().block(240, 240, 10, 10);
    N25old = AA.Nold().block(250, 250, 10, 10);
    N26old = AA.Nold().block(260, 260, 10, 10);
    N27old = AA.Nold().block(270, 270, 10, 10);
    N28old = AA.Nold().block(280, 280, 10, 10);
    N29old = AA.Nold().block(290, 290, 10, 10);
    N30old = AA.Nold().block(300, 300, 10, 10);//分块赋值
    N31old = AA.Nold().block(310, 310, 10, 10);
    N32old = AA.Nold().block(320, 320, 10, 10);
    N33old = AA.Nold().block(330, 330, 10, 10);
    N34old = AA.Nold().block(340, 340, 10, 10);
    N35old = AA.Nold().block(350, 350, 10, 10);
    N36old = AA.Nold().block(360, 360, 10, 10);
    N37old = AA.Nold().block(370, 370, 10, 10);
    N38old = AA.Nold().block(380, 380, 10, 10);
    N39old = AA.Nold().block(390, 390, 10, 10);
    N40old = AA.Nold().block(400, 400, 10, 10);//分块赋值
    N41old = AA.Nold().block(410, 410, 10, 10);
    N42old = AA.Nold().block(420, 420, 10, 10);
    N43old = AA.Nold().block(430, 430, 10, 10);
    N44old = AA.Nold().block(440, 440, 10, 10);
    N45old = AA.Nold().block(450, 450, 10, 10);
    N46old = AA.Nold().block(460, 460, 10, 10);
    N47old = AA.Nold().block(470, 470, 10, 10);
    N48old = AA.Nold().block(480, 480, 10, 10);
    N49old = AA.Nold().block(490, 490, 10, 10);


    N0new = AA.Nnew().block(0, 0, 10, 10);//分块赋值
    N1new = AA.Nnew().block(10, 10, 10, 10);
    N2new = AA.Nnew().block(20, 20, 10, 10);
    N3new = AA.Nnew().block(30, 30, 10, 10);
    N4new = AA.Nnew().block(40, 40, 10, 10);
    N5new = AA.Nnew().block(50, 50, 10, 10);
    N6new = AA.Nnew().block(60, 60, 10, 10);
    N7new = AA.Nnew().block(70, 70, 10, 10);
    N8new = AA.Nnew().block(80, 80, 10, 10);
    N9new = AA.Nnew().block(90, 90, 10, 10);
    N10new = AA.Nnew().block(100, 100, 10, 10);//分块赋值
    N11new = AA.Nnew().block(110, 110, 10, 10);
    N12new = AA.Nnew().block(120, 120, 10, 10);
    N13new = AA.Nnew().block(130, 130, 10, 10);
    N14new = AA.Nnew().block(140, 140, 10, 10);
    N15new = AA.Nnew().block(150, 150, 10, 10);
    N16new = AA.Nnew().block(160, 160, 10, 10);
    N17new = AA.Nnew().block(170, 170, 10, 10);
    N18new = AA.Nnew().block(180, 180, 10, 10);
    N19new = AA.Nnew().block(190, 190, 10, 10);
    N20new = AA.Nnew().block(200, 200, 10, 10);//分块赋值
    N21new = AA.Nnew().block(210, 210, 10, 10);
    N22new = AA.Nnew().block(220, 220, 10, 10);
    N23new = AA.Nnew().block(230, 230, 10, 10);
    N24new = AA.Nnew().block(240, 240, 10, 10);
    N25new = AA.Nnew().block(250, 250, 10, 10);
    N26new = AA.Nnew().block(260, 260, 10, 10);
    N27new = AA.Nnew().block(270, 270, 10, 10);
    N28new = AA.Nnew().block(280, 280, 10, 10);
    N29new = AA.Nnew().block(290, 290, 10, 10);
    N30new = AA.Nnew().block(300, 300, 10, 10);//分块赋值
    N31new = AA.Nnew().block(310, 310, 10, 10);
    N32new = AA.Nnew().block(320, 320, 10, 10);
    N33new = AA.Nnew().block(330, 330, 10, 10);
    N34new = AA.Nnew().block(340, 340, 10, 10);
    N35new = AA.Nnew().block(350, 350, 10, 10);
    N36new = AA.Nnew().block(360, 360, 10, 10);
    N37new = AA.Nnew().block(370, 370, 10, 10);
    N38new = AA.Nnew().block(380, 380, 10, 10);
    N39new = AA.Nnew().block(390, 390, 10, 10);
    N40new = AA.Nnew().block(400, 400, 10, 10);//分块赋值
    N41new = AA.Nnew().block(410, 410, 10, 10);
    N42new = AA.Nnew().block(420, 420, 10, 10);
    N43new = AA.Nnew().block(430, 430, 10, 10);
    N44new = AA.Nnew().block(440, 440, 10, 10);
    N45new = AA.Nnew().block(450, 450, 10, 10);
    N46new = AA.Nnew().block(460, 460, 10, 10);
    N47new = AA.Nnew().block(470, 470, 10, 10);
    N48new = AA.Nnew().block(480, 480, 10, 10);
    N49new = AA.Nnew().block(490, 490, 10, 10);


    VectorXd Q0old(10);
    VectorXd Q1old(10);
    VectorXd Q2old(10);
    VectorXd Q3old(10);
    VectorXd Q4old(10);
    VectorXd Q5old(10);
    VectorXd Q6old(10);
    VectorXd Q7old(10);
    VectorXd Q8old(10);
    VectorXd Q9old(10);
    VectorXd Q10old(10);
    VectorXd Q11old(10);
    VectorXd Q12old(10);
    VectorXd Q13old(10);
    VectorXd Q14old(10);
    VectorXd Q15old(10);
    VectorXd Q16old(10);
    VectorXd Q17old(10);
    VectorXd Q18old(10);
    VectorXd Q19old(10);
    VectorXd Q20old(10);
    VectorXd Q21old(10);
    VectorXd Q22old(10);
    VectorXd Q23old(10);
    VectorXd Q24old(10);
    VectorXd Q25old(10);
    VectorXd Q26old(10);
    VectorXd Q27old(10);
    VectorXd Q28old(10);
    VectorXd Q29old(10);
    VectorXd Q30old(10);
    VectorXd Q31old(10);
    VectorXd Q32old(10);
    VectorXd Q33old(10);
    VectorXd Q34old(10);
    VectorXd Q35old(10);
    VectorXd Q36old(10);
    VectorXd Q37old(10);
    VectorXd Q38old(10);
    VectorXd Q39old(10);
    VectorXd Q40old(10);
    VectorXd Q41old(10);
    VectorXd Q42old(10);
    VectorXd Q43old(10);
    VectorXd Q44old(10);
    VectorXd Q45old(10);
    VectorXd Q46old(10);
    VectorXd Q47old(10);
    VectorXd Q48old(10);
    VectorXd Q49old(10);

    
    VectorXd Q0new(10);
    VectorXd Q1new(10);
    VectorXd Q2new(10);
    VectorXd Q3new(10);
    VectorXd Q4new(10);
    VectorXd Q5new(10);
    VectorXd Q6new(10);
    VectorXd Q7new(10);
    VectorXd Q8new(10);
    VectorXd Q9new(10);
    VectorXd Q10new(10);
    VectorXd Q11new(10);
    VectorXd Q12new(10);
    VectorXd Q13new(10);
    VectorXd Q14new(10);
    VectorXd Q15new(10);
    VectorXd Q16new(10);
    VectorXd Q17new(10);
    VectorXd Q18new(10);
    VectorXd Q19new(10);
    VectorXd Q20new(10);
    VectorXd Q21new(10);
    VectorXd Q22new(10);
    VectorXd Q23new(10);
    VectorXd Q24new(10);
    VectorXd Q25new(10);
    VectorXd Q26new(10);
    VectorXd Q27new(10);
    VectorXd Q28new(10);
    VectorXd Q29new(10);
    VectorXd Q30new(10);
    VectorXd Q31new(10);
    VectorXd Q32new(10);
    VectorXd Q33new(10);
    VectorXd Q34new(10);
    VectorXd Q35new(10);
    VectorXd Q36new(10);
    VectorXd Q37new(10);
    VectorXd Q38new(10);
    VectorXd Q39new(10);
    VectorXd Q40new(10);
    VectorXd Q41new(10);
    VectorXd Q42new(10);
    VectorXd Q43new(10);
    VectorXd Q44new(10);
    VectorXd Q45new(10);
    VectorXd Q46new(10);
    VectorXd Q47new(10);
    VectorXd Q48new(10);
    VectorXd Q49new(10);


    Q0old = AA.Qold().head(10);//分块赋值
    Q1old = AA.Qold().segment(10, 10);
    Q2old = AA.Qold().segment(20, 10);
    Q3old = AA.Qold().segment(30, 10);
    Q4old = AA.Qold().segment(40, 10);
    Q5old = AA.Qold().segment(50, 10);
    Q6old = AA.Qold().segment(60, 10);
    Q7old = AA.Qold().segment(70, 10);
    Q8old = AA.Qold().segment(80, 10);
    Q9old = AA.Qold().segment(90, 10);
    Q10old = AA.Qold().segment(100, 10);
    Q11old = AA.Qold().segment(110, 10);
    Q12old = AA.Qold().segment(120, 10);
    Q13old = AA.Qold().segment(130, 10);
    Q14old = AA.Qold().segment(140, 10);
    Q15old = AA.Qold().segment(150, 10);
    Q16old = AA.Qold().segment(160, 10);
    Q17old = AA.Qold().segment(170, 10);
    Q18old = AA.Qold().segment(180, 10);
    Q19old = AA.Qold().segment(190, 10);
    Q20old = AA.Qold().segment(200, 10);
    Q21old = AA.Qold().segment(210, 10);
    Q22old = AA.Qold().segment(220, 10);
    Q23old = AA.Qold().segment(230, 10);
    Q24old = AA.Qold().segment(240, 10);
    Q25old = AA.Qold().segment(250, 10);
    Q26old = AA.Qold().segment(260, 10);
    Q27old = AA.Qold().segment(270, 10);
    Q28old = AA.Qold().segment(280, 10);
    Q29old = AA.Qold().segment(290, 10);
    Q30old = AA.Qold().segment(300, 10);
    Q31old = AA.Qold().segment(310, 10);
    Q32old = AA.Qold().segment(320, 10);
    Q33old = AA.Qold().segment(330, 10);
    Q34old = AA.Qold().segment(340, 10);
    Q35old = AA.Qold().segment(350, 10);
    Q36old = AA.Qold().segment(360, 10);
    Q37old = AA.Qold().segment(370, 10);
    Q38old = AA.Qold().segment(380, 10);
    Q39old = AA.Qold().segment(390, 10);
    Q40old = AA.Qold().segment(400, 10);
    Q41old = AA.Qold().segment(410, 10);
    Q42old = AA.Qold().segment(420, 10);
    Q43old = AA.Qold().segment(430, 10);
    Q44old = AA.Qold().segment(440, 10);
    Q45old = AA.Qold().segment(450, 10);
    Q46old = AA.Qold().segment(460, 10);
    Q47old = AA.Qold().segment(470, 10);
    Q48old = AA.Qold().segment(480, 10);
    Q49old = AA.Qold().tail(10);

    Q0new = AA.Qnew().head(10);//分块赋值
    Q1new = AA.Qnew().segment(10, 10);
    Q2new = AA.Qnew().segment(20, 10);
    Q3new = AA.Qnew().segment(30, 10);
    Q4new = AA.Qnew().segment(40, 10);
    Q5new = AA.Qnew().segment(50, 10);
    Q6new = AA.Qnew().segment(60, 10);
    Q7new = AA.Qnew().segment(70, 10);
    Q8new = AA.Qnew().segment(80, 10);
    Q9new = AA.Qnew().segment(90, 10);
    Q10new = AA.Qnew().segment(100, 10);
    Q11new = AA.Qnew().segment(110, 10);
    Q12new = AA.Qnew().segment(120, 10);
    Q13new = AA.Qnew().segment(130, 10);
    Q14new = AA.Qnew().segment(140, 10);
    Q15new = AA.Qnew().segment(150, 10);
    Q16new = AA.Qnew().segment(160, 10);
    Q17new = AA.Qnew().segment(170, 10);
    Q18new = AA.Qnew().segment(180, 10);
    Q19new = AA.Qnew().segment(190, 10);
    Q20new = AA.Qnew().segment(200, 10);
    Q21new = AA.Qnew().segment(210, 10);
    Q22new = AA.Qnew().segment(220, 10);
    Q23new = AA.Qnew().segment(230, 10);
    Q24new = AA.Qnew().segment(240, 10);
    Q25new = AA.Qnew().segment(250, 10);
    Q26new = AA.Qnew().segment(260, 10);
    Q27new = AA.Qnew().segment(270, 10);
    Q28new = AA.Qnew().segment(280, 10);
    Q29new = AA.Qnew().segment(290, 10);
    Q30new = AA.Qnew().segment(300, 10);
    Q31new = AA.Qnew().segment(310, 10);
    Q32new = AA.Qnew().segment(320, 10);
    Q33new = AA.Qnew().segment(330, 10);
    Q34new = AA.Qnew().segment(340, 10);
    Q35new = AA.Qnew().segment(350, 10);
    Q36new = AA.Qnew().segment(360, 10);
    Q37new = AA.Qnew().segment(370, 10);
    Q38new = AA.Qnew().segment(380, 10);
    Q39new = AA.Qnew().segment(390, 10);
    Q40new = AA.Qnew().segment(400, 10);
    Q41new = AA.Qnew().segment(410, 10);
    Q42new = AA.Qnew().segment(420, 10);
    Q43new = AA.Qnew().segment(430, 10);
    Q44new = AA.Qnew().segment(440, 10);
    Q45new = AA.Qnew().segment(450, 10);
    Q46new = AA.Qnew().segment(460, 10);
    Q47new = AA.Qnew().segment(470, 10);
    Q48new = AA.Qnew().segment(480, 10);
    Q49new = AA.Qnew().tail(10);


    VectorXd temp0(10);//0节点赋值
    VectorXd temp1(10);//1节点赋值
    VectorXd temp2(10);//2节点赋值
    VectorXd temp3(10);
    VectorXd temp4(10);
    VectorXd temp5(10);
    VectorXd temp6(10);
    VectorXd temp7(10);
    VectorXd temp8(10);
    VectorXd temp9(10);
    VectorXd temp10(10);//0节点赋值
    VectorXd temp11(10);//1节点赋值
    VectorXd temp12(10);//2节点赋值
    VectorXd temp13(10);
    VectorXd temp14(10);
    VectorXd temp15(10);
    VectorXd temp16(10);
    VectorXd temp17(10);
    VectorXd temp18(10);
    VectorXd temp19(10);
    VectorXd temp20(10);//0节点赋值
    VectorXd temp21(10);//1节点赋值
    VectorXd temp22(10);//2节点赋值
    VectorXd temp23(10);
    VectorXd temp24(10);
    VectorXd temp25(10);
    VectorXd temp26(10);
    VectorXd temp27(10);
    VectorXd temp28(10);
    VectorXd temp29(10);
    VectorXd temp30(10);//0节点赋值
    VectorXd temp31(10);//1节点赋值
    VectorXd temp32(10);//2节点赋值
    VectorXd temp33(10);
    VectorXd temp34(10);
    VectorXd temp35(10);
    VectorXd temp36(10);
    VectorXd temp37(10);
    VectorXd temp38(10);
    VectorXd temp39(10);
    VectorXd temp40(10);//0节点赋值
    VectorXd temp41(10);//1节点赋值
    VectorXd temp42(10);//2节点赋值
    VectorXd temp43(10);
    VectorXd temp44(10);
    VectorXd temp45(10);
    VectorXd temp46(10);
    VectorXd temp47(10);
    VectorXd temp48(10);
    VectorXd temp49(10);


    //temp(n) n为中间节点号
    temp0 = ((N1new + N0new)*(Y1new - Y0new) * deltaT) + ((N1old + N0old)*(Y1old - Y0old) * deltaT) + ((M1new + M1old)*(Y1new - Y1old) * deltaS) + ((M0new + M0old)*(Y0new - Y0old) * deltaS) + (Q0old + Q1old + Q0new + Q1new) * (deltaT*deltaS);
    temp1 = ((N2new + N1new)*(Y2new - Y1new) * deltaT) + ((N2old + N1old)*(Y2old - Y1old) * deltaT) + ((M2new + M2old)*(Y2new - Y2old) * deltaS) + ((M1new + M1old)*(Y1new - Y1old) * deltaS) + (Q1old + Q2old + Q1new + Q2new) * (deltaT*deltaS);
    temp2 = ((N3new + N2new)*(Y3new - Y2new) * deltaT) + ((N3old + N2old)*(Y3old - Y2old) * deltaT) + ((M3new + M3old)*(Y3new - Y3old) * deltaS) + ((M2new + M2old)*(Y2new - Y2old) * deltaS) + (Q2old + Q3old + Q2new + Q3new) * (deltaT*deltaS);
    temp3 = ((N4new + N3new)*(Y4new - Y3new) * deltaT) + ((N4old + N3old)*(Y4old - Y3old) * deltaT) + ((M4new + M4old)*(Y4new - Y4old) * deltaS) + ((M3new + M3old)*(Y3new - Y3old) * deltaS) + (Q3old + Q4old + Q3new + Q4new) * (deltaT*deltaS);
    temp4 = ((N5new + N4new)*(Y5new - Y4new) * deltaT) + ((N5old + N4old)*(Y5old - Y4old) * deltaT) + ((M5new + M5old)*(Y5new - Y5old) * deltaS) + ((M4new + M4old)*(Y4new - Y4old) * deltaS) + (Q4old + Q5old + Q4new + Q5new) * (deltaT*deltaS);
    temp5 = ((N6new + N5new)*(Y6new - Y5new) * deltaT) + ((N6old + N5old)*(Y6old - Y5old) * deltaT) + ((M6new + M6old)*(Y6new - Y6old) * deltaS) + ((M5new + M5old)*(Y5new - Y5old) * deltaS) + (Q5old + Q6old + Q5new + Q6new) * (deltaT*deltaS);
    temp6 = ((N7new + N6new)*(Y7new - Y6new) * deltaT) + ((N7old + N6old)*(Y7old - Y6old) * deltaT) + ((M7new + M7old)*(Y7new - Y7old) * deltaS) + ((M6new + M6old)*(Y6new - Y6old) * deltaS) + (Q6old + Q7old + Q6new + Q7new) * (deltaT*deltaS);
    temp7 = ((N8new + N7new)*(Y8new - Y7new) * deltaT) + ((N8old + N7old)*(Y8old - Y7old) * deltaT) + ((M8new + M8old)*(Y8new - Y8old) * deltaS) + ((M7new + M7old)*(Y7new - Y7old) * deltaS) + (Q7old + Q8old + Q7new + Q8new) * (deltaT*deltaS);
    temp8 = ((N9new + N8new)*(Y9new - Y8new) * deltaT) + ((N9old + N8old)*(Y9old - Y8old) * deltaT) + ((M9new + M9old)*(Y9new - Y9old) * deltaS) + ((M8new + M8old)*(Y8new - Y8old) * deltaS) + (Q8old + Q9old + Q8new + Q9new) * (deltaT*deltaS);
    temp9 = ((N10new + N9new)*(Y10new - Y9new) * deltaT) + ((N10old + N9old)*(Y10old - Y9old) * deltaT) + ((M10new + M10old)*(Y10new - Y10old) * deltaS) + ((M9new + M9old)*(Y9new - Y9old) * deltaS) + (Q9old + Q10old + Q9new + Q10new) * (deltaT*deltaS);

    temp10 = ((N11new + N10new)*(Y11new - Y10new) * deltaT) + ((N11old + N10old)*(Y11old - Y10old) * deltaT) + ((M11new + M11old)*(Y11new - Y11old) * deltaS) + ((M10new + M10old)*(Y10new - Y10old) * deltaS) + (Q10old + Q11old + Q10new + Q11new) * (deltaT*deltaS);
    temp11 = ((N12new + N11new)*(Y12new - Y11new) * deltaT) + ((N12old + N11old)*(Y12old - Y11old) * deltaT) + ((M12new + M12old)*(Y12new - Y12old) * deltaS) + ((M11new + M11old)*(Y11new - Y11old) * deltaS) + (Q11old + Q12old + Q11new + Q12new) * (deltaT*deltaS);
    temp12 = ((N13new + N12new)*(Y13new - Y12new) * deltaT) + ((N13old + N12old)*(Y13old - Y12old) * deltaT) + ((M13new + M13old)*(Y13new - Y13old) * deltaS) + ((M12new + M12old)*(Y12new - Y12old) * deltaS) + (Q12old + Q13old + Q12new + Q13new) * (deltaT*deltaS);
    temp13 = ((N14new + N13new)*(Y14new - Y13new) * deltaT) + ((N14old + N13old)*(Y14old - Y13old) * deltaT) + ((M14new + M14old)*(Y14new - Y14old) * deltaS) + ((M13new + M13old)*(Y13new - Y13old) * deltaS) + (Q13old + Q14old + Q13new + Q14new) * (deltaT*deltaS);
    temp14 = ((N15new + N14new)*(Y15new - Y14new) * deltaT) + ((N15old + N14old)*(Y15old - Y14old) * deltaT) + ((M15new + M15old)*(Y15new - Y15old) * deltaS) + ((M14new + M14old)*(Y14new - Y14old) * deltaS) + (Q14old + Q15old + Q14new + Q15new) * (deltaT*deltaS);
    temp15 = ((N16new + N15new)*(Y16new - Y15new) * deltaT) + ((N16old + N15old)*(Y16old - Y15old) * deltaT) + ((M16new + M16old)*(Y16new - Y16old) * deltaS) + ((M15new + M15old)*(Y15new - Y15old) * deltaS) + (Q15old + Q16old + Q15new + Q16new) * (deltaT*deltaS);
    temp16 = ((N17new + N16new)*(Y17new - Y16new) * deltaT) + ((N17old + N16old)*(Y17old - Y16old) * deltaT) + ((M17new + M17old)*(Y17new - Y17old) * deltaS) + ((M16new + M16old)*(Y16new - Y16old) * deltaS) + (Q16old + Q17old + Q16new + Q17new) * (deltaT*deltaS);
    temp17 = ((N18new + N17new)*(Y18new - Y17new) * deltaT) + ((N18old + N17old)*(Y18old - Y17old) * deltaT) + ((M18new + M18old)*(Y18new - Y18old) * deltaS) + ((M17new + M17old)*(Y17new - Y17old) * deltaS) + (Q17old + Q18old + Q17new + Q18new) * (deltaT*deltaS);
    temp18 = ((N19new + N18new)*(Y19new - Y18new) * deltaT) + ((N19old + N18old)*(Y19old - Y18old) * deltaT) + ((M19new + M19old)*(Y19new - Y19old) * deltaS) + ((M18new + M18old)*(Y18new - Y18old) * deltaS) + (Q18old + Q19old + Q18new + Q19new) * (deltaT*deltaS);
    temp19 = ((N20new + N19new)*(Y20new - Y19new) * deltaT) + ((N20old + N19old)*(Y20old - Y19old) * deltaT) + ((M20new + M20old)*(Y20new - Y20old) * deltaS) + ((M19new + M19old)*(Y19new - Y19old) * deltaS) + (Q19old + Q20old + Q19new + Q20new) * (deltaT*deltaS);

    temp20 = ((N21new + N20new)*(Y21new - Y20new) * deltaT) + ((N21old + N20old)*(Y21old - Y20old) * deltaT) + ((M21new + M21old)*(Y21new - Y21old) * deltaS) + ((M20new + M20old)*(Y20new - Y20old) * deltaS) + (Q20old + Q21old + Q20new + Q21new) * (deltaT*deltaS);
    temp21 = ((N22new + N21new)*(Y22new - Y21new) * deltaT) + ((N22old + N21old)*(Y22old - Y21old) * deltaT) + ((M22new + M22old)*(Y22new - Y22old) * deltaS) + ((M21new + M21old)*(Y21new - Y21old) * deltaS) + (Q21old + Q22old + Q21new + Q22new) * (deltaT*deltaS);
    temp22 = ((N23new + N22new)*(Y23new - Y22new) * deltaT) + ((N23old + N22old)*(Y23old - Y22old) * deltaT) + ((M23new + M23old)*(Y23new - Y23old) * deltaS) + ((M22new + M22old)*(Y22new - Y22old) * deltaS) + (Q22old + Q23old + Q22new + Q23new) * (deltaT*deltaS);
    temp23 = ((N24new + N23new)*(Y24new - Y23new) * deltaT) + ((N24old + N23old)*(Y24old - Y23old) * deltaT) + ((M24new + M24old)*(Y24new - Y24old) * deltaS) + ((M23new + M23old)*(Y23new - Y23old) * deltaS) + (Q23old + Q24old + Q23new + Q24new) * (deltaT*deltaS);
    temp24 = ((N25new + N24new)*(Y25new - Y24new) * deltaT) + ((N25old + N24old)*(Y25old - Y24old) * deltaT) + ((M25new + M25old)*(Y25new - Y25old) * deltaS) + ((M24new + M24old)*(Y24new - Y24old) * deltaS) + (Q24old + Q25old + Q24new + Q25new) * (deltaT*deltaS);
    temp25 = ((N26new + N25new)*(Y26new - Y25new) * deltaT) + ((N26old + N25old)*(Y26old - Y25old) * deltaT) + ((M26new + M26old)*(Y26new - Y26old) * deltaS) + ((M25new + M25old)*(Y25new - Y25old) * deltaS) + (Q25old + Q26old + Q25new + Q26new) * (deltaT*deltaS);
    temp26 = ((N27new + N26new)*(Y27new - Y26new) * deltaT) + ((N27old + N26old)*(Y27old - Y26old) * deltaT) + ((M27new + M27old)*(Y27new - Y27old) * deltaS) + ((M26new + M26old)*(Y26new - Y26old) * deltaS) + (Q26old + Q27old + Q26new + Q27new) * (deltaT*deltaS);
    temp27 = ((N28new + N27new)*(Y28new - Y27new) * deltaT) + ((N28old + N27old)*(Y28old - Y27old) * deltaT) + ((M28new + M28old)*(Y28new - Y28old) * deltaS) + ((M27new + M27old)*(Y27new - Y27old) * deltaS) + (Q27old + Q28old + Q27new + Q28new) * (deltaT*deltaS);
    temp28 = ((N29new + N28new)*(Y29new - Y28new) * deltaT) + ((N29old + N28old)*(Y29old - Y28old) * deltaT) + ((M29new + M29old)*(Y29new - Y29old) * deltaS) + ((M28new + M28old)*(Y28new - Y28old) * deltaS) + (Q28old + Q29old + Q28new + Q29new) * (deltaT*deltaS);
    temp29 = ((N30new + N29new)*(Y30new - Y29new) * deltaT) + ((N30old + N29old)*(Y30old - Y29old) * deltaT) + ((M30new + M30old)*(Y30new - Y30old) * deltaS) + ((M29new + M29old)*(Y29new - Y29old) * deltaS) + (Q29old + Q30old + Q29new + Q30new) * (deltaT*deltaS);

    temp30 = ((N31new + N30new)*(Y31new - Y30new) * deltaT) + ((N31old + N30old)*(Y31old - Y30old) * deltaT) + ((M31new + M31old)*(Y31new - Y31old) * deltaS) + ((M30new + M30old)*(Y30new - Y30old) * deltaS) + (Q30old + Q31old + Q30new + Q31new) * (deltaT*deltaS);
    temp31 = ((N32new + N31new)*(Y32new - Y31new) * deltaT) + ((N32old + N31old)*(Y32old - Y31old) * deltaT) + ((M32new + M32old)*(Y32new - Y32old) * deltaS) + ((M31new + M31old)*(Y31new - Y31old) * deltaS) + (Q31old + Q32old + Q31new + Q32new) * (deltaT*deltaS);
    temp32 = ((N33new + N32new)*(Y33new - Y32new) * deltaT) + ((N33old + N32old)*(Y33old - Y32old) * deltaT) + ((M33new + M33old)*(Y33new - Y33old) * deltaS) + ((M32new + M32old)*(Y32new - Y32old) * deltaS) + (Q32old + Q33old + Q32new + Q33new) * (deltaT*deltaS);
    temp33 = ((N34new + N33new)*(Y34new - Y33new) * deltaT) + ((N34old + N33old)*(Y34old - Y33old) * deltaT) + ((M34new + M34old)*(Y34new - Y34old) * deltaS) + ((M33new + M33old)*(Y33new - Y33old) * deltaS) + (Q33old + Q34old + Q33new + Q34new) * (deltaT*deltaS);
    temp34 = ((N35new + N34new)*(Y35new - Y34new) * deltaT) + ((N35old + N34old)*(Y35old - Y34old) * deltaT) + ((M35new + M35old)*(Y35new - Y35old) * deltaS) + ((M34new + M34old)*(Y34new - Y34old) * deltaS) + (Q34old + Q35old + Q34new + Q35new) * (deltaT*deltaS);
    temp35 = ((N36new + N35new)*(Y36new - Y35new) * deltaT) + ((N36old + N35old)*(Y36old - Y35old) * deltaT) + ((M36new + M36old)*(Y36new - Y36old) * deltaS) + ((M35new + M35old)*(Y35new - Y35old) * deltaS) + (Q35old + Q36old + Q35new + Q36new) * (deltaT*deltaS);
    temp36 = ((N37new + N36new)*(Y37new - Y36new) * deltaT) + ((N37old + N36old)*(Y37old - Y36old) * deltaT) + ((M37new + M37old)*(Y37new - Y37old) * deltaS) + ((M36new + M36old)*(Y36new - Y36old) * deltaS) + (Q36old + Q37old + Q36new + Q37new) * (deltaT*deltaS);
    temp37 = ((N38new + N37new)*(Y38new - Y37new) * deltaT) + ((N38old + N37old)*(Y38old - Y37old) * deltaT) + ((M38new + M38old)*(Y38new - Y38old) * deltaS) + ((M37new + M37old)*(Y37new - Y37old) * deltaS) + (Q37old + Q38old + Q37new + Q38new) * (deltaT*deltaS);
    temp38 = ((N39new + N38new)*(Y39new - Y38new) * deltaT) + ((N39old + N38old)*(Y39old - Y38old) * deltaT) + ((M39new + M39old)*(Y39new - Y39old) * deltaS) + ((M38new + M38old)*(Y38new - Y38old) * deltaS) + (Q38old + Q39old + Q38new + Q39new) * (deltaT*deltaS);
    temp39 = ((N40new + N39new)*(Y40new - Y39new) * deltaT) + ((N40old + N39old)*(Y40old - Y39old) * deltaT) + ((M40new + M40old)*(Y40new - Y40old) * deltaS) + ((M39new + M39old)*(Y39new - Y39old) * deltaS) + (Q39old + Q40old + Q39new + Q40new) * (deltaT*deltaS);

    temp40 = ((N41new + N40new)*(Y41new - Y40new) * deltaT) + ((N41old + N40old)*(Y41old - Y40old) * deltaT) + ((M41new + M41old)*(Y41new - Y41old) * deltaS) + ((M40new + M40old)*(Y40new - Y40old) * deltaS) + (Q40old + Q41old + Q40new + Q41new) * (deltaT*deltaS);
    temp41 = ((N42new + N41new)*(Y42new - Y41new) * deltaT) + ((N42old + N41old)*(Y42old - Y41old) * deltaT) + ((M42new + M42old)*(Y42new - Y42old) * deltaS) + ((M41new + M41old)*(Y41new - Y41old) * deltaS) + (Q41old + Q42old + Q41new + Q42new) * (deltaT*deltaS);
    temp42 = ((N43new + N42new)*(Y43new - Y42new) * deltaT) + ((N43old + N42old)*(Y43old - Y42old) * deltaT) + ((M43new + M43old)*(Y43new - Y43old) * deltaS) + ((M42new + M42old)*(Y42new - Y42old) * deltaS) + (Q42old + Q43old + Q42new + Q43new) * (deltaT*deltaS);
    temp43 = ((N44new + N43new)*(Y44new - Y43new) * deltaT) + ((N44old + N43old)*(Y44old - Y43old) * deltaT) + ((M44new + M44old)*(Y44new - Y44old) * deltaS) + ((M43new + M43old)*(Y43new - Y43old) * deltaS) + (Q43old + Q44old + Q43new + Q44new) * (deltaT*deltaS);
    temp44 = ((N45new + N44new)*(Y45new - Y44new) * deltaT) + ((N45old + N44old)*(Y45old - Y44old) * deltaT) + ((M45new + M45old)*(Y45new - Y45old) * deltaS) + ((M44new + M44old)*(Y44new - Y44old) * deltaS) + (Q44old + Q45old + Q44new + Q45new) * (deltaT*deltaS);
    temp45 = ((N46new + N45new)*(Y46new - Y45new) * deltaT) + ((N46old + N45old)*(Y46old - Y45old) * deltaT) + ((M46new + M46old)*(Y46new - Y46old) * deltaS) + ((M45new + M45old)*(Y45new - Y45old) * deltaS) + (Q45old + Q46old + Q45new + Q46new) * (deltaT*deltaS);
    temp46 = ((N47new + N46new)*(Y47new - Y46new) * deltaT) + ((N47old + N46old)*(Y47old - Y46old) * deltaT) + ((M47new + M47old)*(Y47new - Y47old) * deltaS) + ((M46new + M46old)*(Y46new - Y46old) * deltaS) + (Q46old + Q47old + Q46new + Q47new) * (deltaT*deltaS);
    temp47 = ((N48new + N47new)*(Y48new - Y47new) * deltaT) + ((N48old + N47old)*(Y48old - Y47old) * deltaT) + ((M48new + M48old)*(Y48new - Y48old) * deltaS) + ((M47new + M47old)*(Y47new - Y47old) * deltaS) + (Q47old + Q48old + Q47new + Q48new) * (deltaT*deltaS);
    temp48 = ((N49new + N48new)*(Y49new - Y48new) * deltaT) + ((N49old + N48old)*(Y49old - Y48old) * deltaT) + ((M49new + M49old)*(Y49new - Y49old) * deltaS) + ((M48new + M48old)*(Y48new - Y48old) * deltaS) + (Q48old + Q49old + Q48new + Q49new) * (deltaT*deltaS);


    //BC 边界条件
    VectorXd BCtemp0(5);
    VectorXd BCtemp1(5);

    //第一点(u, v, w, O2mega, O3mega)(V1 轴向; V2，V3均为水平方向)
    BCtemp0(0) = Y0new(0) - V1*cos(Y0new(7))*cos(Y0new(6)) - V2*sin(Y0new(6))*cos(Y0new(7)) -V3*sin(Y0new(7));
    BCtemp0(1) = Y0new(1) - V2*cos(Y0new(6)) + V1*sin(Y0new(6));
    BCtemp0(2) = Y0new(2) + V1*sin(Y0new(7))*cos(Y0new(6)) + V2*sin(Y0new(7))*sin(Y0new(6)) - V3*cos(Y0new(7));
    BCtemp0(3) = Y0new(8);
    BCtemp0(4) = Y0new(9);


    //最后点(T, Sn, Sb, O2mega, O3mega)
    // BCtemp1(0) = Y49new(3) + G*cos(Y49new(7))*cos(Y49new(6)) + 0.5*pi*Rho*d0*Cdt*Y49new(0)*abs(Y49new(0));
    // BCtemp1(1) = Y49new(4) - G*sin(Y49new(6)) + 0.5*Rho*d0*Cdn*Y49new(1)*sqrt(pow(Y49new(1),2)+ pow(Y49new(2),2))+ ma*(Y49new(1)-Y49old(1))/deltaT;
    // BCtemp1(2) = Y49new(5) - G*cos(Y49new(6))*sin(Y49new(7)) + 0.5*Rho*d0*Cdb*Y49new(2)*sqrt(pow(Y49new(1),2)+ pow(Y49new(2),2))+ ma*(Y49new(2)-Y49old(2))/deltaT;
    // BCtemp1(3) = Y49new(8);
    // BCtemp1(4) = Y49new(9);


    BCtemp1(0) = Y49new(3) + G*cos(Y49new(7))*cos(Y49new(6));
    BCtemp1(1) = Y49new(4) - G*sin(Y49new(6)) + 0.5*rho*Cdn*Sd*Y49new(1)*sqrt(pow(Y49new(1),2)+pow(Y49new(2),2));
    BCtemp1(2) = Y49new(5) - G*cos(Y49new(6))*sin(Y49new(7))+ 0.5*rho*Cdn*Sd*Y49new(2)*sqrt(pow(Y49new(1),2)+pow(Y49new(2),2));
    BCtemp1(3) = Y49new(8);
    BCtemp1(4) = Y49new(9);



    VectorXd ret(500);
    ret.head(5) = BCtemp0;
    ret.segment(5, 10) = temp0;
    ret.segment(15, 10) = temp1;
    ret.segment(25, 10) = temp2;
    ret.segment(35, 10) = temp3;
    ret.segment(45, 10) = temp4;
    ret.segment(55, 10) = temp5;
    ret.segment(65, 10) = temp6;
    ret.segment(75, 10) = temp7;
    ret.segment(85, 10) = temp8;
    ret.segment(95, 10) = temp9;

    ret.segment(105, 10) = temp10;
    ret.segment(115, 10) = temp11;
    ret.segment(125, 10) = temp12;
    ret.segment(135, 10) = temp13;
    ret.segment(145, 10) = temp14;
    ret.segment(155, 10) = temp15;
    ret.segment(165, 10) = temp16;
    ret.segment(175, 10) = temp17;
    ret.segment(185, 10) = temp18;
    ret.segment(195, 10) = temp19;

    ret.segment(205, 10) = temp20;
    ret.segment(215, 10) = temp21;
    ret.segment(225, 10) = temp22;
    ret.segment(235, 10) = temp23;
    ret.segment(245, 10) = temp24;
    ret.segment(255, 10) = temp25;
    ret.segment(265, 10) = temp26;
    ret.segment(275, 10) = temp27;
    ret.segment(285, 10) = temp28;
    ret.segment(295, 10) = temp29;

    ret.segment(305, 10) = temp30;
    ret.segment(315, 10) = temp31;
    ret.segment(325, 10) = temp32;
    ret.segment(335, 10) = temp33;
    ret.segment(345, 10) = temp34;
    ret.segment(355, 10) = temp35;
    ret.segment(365, 10) = temp36;
    ret.segment(375, 10) = temp37;
    ret.segment(385, 10) = temp38;
    ret.segment(395, 10) = temp39;

    ret.segment(405, 10) = temp40;
    ret.segment(415, 10) = temp41;
    ret.segment(425, 10) = temp42;
    ret.segment(435, 10) = temp43;
    ret.segment(445, 10) = temp44;
    ret.segment(455, 10) = temp45;
    ret.segment(465, 10) = temp46;
    ret.segment(475, 10) = temp47;
    ret.segment(485, 10) = temp48;

    ret.tail(5) = BCtemp1;

    return ret;
}