/*
 * @Author: error: git config user.name && git config user.email & please set dead value or install git
 * @Date: 2022-11-14 15:28:20
 * @LastEditors: error: error: git config user.name & please set dead value or install git && error: git config user.email & please set dead value or install git & please set dead value or install git
 * @LastEditTime: 2023-03-02 19:56:41
 * @FilePath: \Siano3.4.1\csvmaker\csvmaker.cpp
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <cmath>

using namespace std;
using namespace Eigen;

int main()
{
    int Nodes = 50;
    int variable = 10;
    int TV = 500;                                                //总变量数TV= Nodes * variable
    int TimeStep = 1000;                                         //总时间步数
    double DelTime = 0.01;                                      //时间步长(真实时间步长)
    double pi = 3.1415926;

    MatrixXd a(TimeStep, TV);
    a.Zero(TimeStep, TV);


    //***********************************************************************//
    
    // for(int j = 0; j < TimeStep; j++)
    // {
    //     for(int i = 0; i < Nodes; i++)
    //     {
    //         a(j, i*10 + 0) = 0;    //u
    //         a(j, i*10 + 1) = 0.00000001;    //v
    //         a(j, i*10 + 2) = 0;    //w

    //         a(j, i*10 + 3) = 1.201;  //T
    //         a(j, i*10 + 4) = 0;    //Sn
    //         a(j, i*10 + 5) = 0;    //Sb

    //         a(j, i*10 + 6) = 0;    //Phi
    //         a(j, i*10 + 7) = 0;    //Theta

    //         a(j, i*10 + 8) = 0.1;    //O2mega
    //         a(j, i*10 + 9) = 0.1;    //O3mega
    //    }
    // }


    // ofstream dataFile;
    // dataFile.open(".././Data./out.csv", ios::out | ios::trunc);
    // for(int i = 0; i < TimeStep; i++)
    // {
    //     for(int j = 0; j < TV; j++)
    //     {
    //         dataFile << a(i, j) << ",";
    //     }
    //     dataFile << endl;
    // }
    // dataFile.close();

    //*********************************************************************//
    MatrixXd b(TimeStep, 11);
    b.Zero(TimeStep, 11);
    for(int i = 0; i < TimeStep; i++)
    {
        //01
        // if(i< 50)
        // {
        //     b(i, 1) = 0.05 * (i + 1);
        // }
        // else
        // {
        //     b(i, 1) = 2.5;
        // }

        //02
        b(i, 1) = -0.5 * 2*pi/1.0 *sin(2*pi/1.0 * (i+1) * DelTime);
        b(i, 2) = 0.5 * 2*pi/1.0 *cos(2*pi/1.0 * (i+1) * DelTime);


        b(i, 8) = 0.15;//Sd1
        b(i, 9) = 0.15;//Sd2
        b(i, 10) = -10;//边界条件负向，牵引需填写负值
    }

    ofstream dataFile;
    dataFile.open(".././Data./input.csv", ios::out | ios::trunc);
    for(int i = 0; i < TimeStep; i++)
    {
        for(int j = 0; j < 11; j++)
        {
            dataFile << b(i, j) << ",";
        }
        dataFile << endl;
    }
    dataFile.close();



    return 0;
}
