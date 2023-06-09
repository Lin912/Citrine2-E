/*
 * @Author: Lyn
 * @Date: 2022-05-26 17:54:40
 * @LastEditors: Lyn 18340802816@163.com
 * @LastEditTime: 2023-05-26 22:24:23
 * @FilePath: \test1\main.cpp
 */


#include <iostream>
#include <Eigen/Dense>
#include "Iterator.h"
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>


using namespace std;
using namespace Eigen;

void  getArrCol(Matrix<double,500,1> &arr_col,int index, string file)
{
  ifstream infile(file, ios::in);
    if (!infile)
    {exit(1);}
    string line;
    string field;
    for (int i = 0; i <=index; i++)
    {
      getline(infile, line);
    }
    stringstream sin(line);
    int colindex = 0;
    while (getline(sin, field, ','))
    {
      double dvalue = atof(field.c_str());
      arr_col(colindex, 0) = dvalue;
      colindex++;
    }
    infile.close();
}
int main()
{
    int Times = 200;                                          
    double Error = 5e-08;                                  
    int Nodes = 50;
    int variable = 10;
    int TV = 500;                                            
    int TimeStep = 1000;                                       
    double DelTime = 0.01;                               
    VectorXd a(TV);
    a.Zero(TV);
    MatrixXd Mcsv(TimeStep, TV);                                 
    Mcsv.Zero(TimeStep, TV);  
    
    VectorXd TransVal(TV);                                         

    Matrix<double, 500, 1> ans;                                                           
    for(int i = 0; i < TimeStep; i++)
    {
        cout << endl;
        cout << "Time Step : " << i << "     Now the real time is " << i*DelTime << "s";
        
      if(i > 0)
        {
          getArrCol(ans, i-1, "./Data./out.csv");
          for(int i = 0; i < TV; i++)
          {
              TransVal(i) = ans[i];
          }
        }
        else
        {
          TransVal = a;                                           
        }
        Iterator b(TransVal, TransVal, Times, Error);              
        b.begin(i);                                                 
        TransVal = b.out();                                            
        for(int j = 0; j < TV; j++)
        {
            Mcsv(i, j) = TransVal(j);
        }
        ofstream dataFile;
        dataFile.open("./Data./out.csv", ios::out | ios::trunc);
        for(int i = 0; i < TimeStep; i++)
        {
            for(int j = 0; j < TV; j++)
            {
                dataFile << Mcsv(i, j) << ",";
            }
            dataFile << endl;
        }
        dataFile.close();
        ofstream detaFile;
        detaFile.open("./Matrix./Proced./out.csv", ios::out | ios::trunc);
        for(int i = 0; i < TimeStep; i++)
        {
            for(int j = 0; j < TV; j++)
            {
                detaFile << Mcsv(i, j) << ",";
            }
            detaFile << endl;
        }
        detaFile.close();
    }   
    return 0;
}

