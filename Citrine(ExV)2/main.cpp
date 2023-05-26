/*
 * @Author: Lyn
 * @Date: 2022-05-26 17:54:40
 * @LastEditors: error: error: git config user.name & please set dead value or install git && error: git config user.email & please set dead value or install git & please set dead value or install git
 * @LastEditTime: 2023-02-28 16:09:18
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

void  getArrCol(Matrix<double,500,1> &arr_col,int index, string file) // index 读取第多少行的列向量 , file 待读取的文件名 
{
  ifstream infile(file, ios::in);
    if (!infile)
    {
      cout << "打开文件夹" + file + "失败" << endl;
      exit(1);
    }

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
   //迭代参数区
    /****************************************************************************************/
    int Times = 200;                                              //迭代次数
    double Error = 5e-08;                                      //迭代误差(在迭代器内调整)
    int Nodes = 50;
    int variable = 10;
    int TV = 500;                                                 //总变量数TV= Nodes * variable
    int TimeStep = 1000;                                          //总时间步数
    double DelTime = 0.01;                                        //时间步长(用于查看真实时间)
    /****************************************************************************************/


    //初值定义区
    /***************************************************************************************/
    VectorXd a(TV);
    a.Zero(TV);
    

    for(int i = 0; i < Nodes; i++)
    {
      a(i*10 + 0) = 0;    //u
      a(i*10 + 1) = 0.002;    //v
      a(i*10 + 2) = 0;    //w

      a(i*10 + 3) = 0;   //T
      a(i*10 + 4) = 0;    //Sn
      a(i*10 + 5) = 0;    //Sb

      a(i*10 + 6) = 0.00001;    //Phi
      a(i*10 + 7) = 0.00001;    //Theta

      a(i*10 + 8) = 0.00001;    //O2mega
      a(i*10 + 9) = 0.00001;    //O3mega
    }
      a(493) = 3480.0;


    /****************************************************************************************/

    //CSV初始化区
    /****************************************************************************************/
    MatrixXd Mcsv(TimeStep, TV);                                    //初始CSV空表格
    Mcsv.Zero(TimeStep, TV);  
    
    VectorXd TransVal(TV);                                          //计算容器

    Matrix<double, 500, 1> ans;                                     //传值容器                                
    /*****************************************************************************************/


    //主体区
    /*****************************************************************************************/
    for(int i = 0; i < TimeStep; i++)
    {
        cout << endl;
        cout << "Time Step : " << i << "     Now the real time is " << i*DelTime << "s";
        
        //读入部(读入上一个时间步结束的值)
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
          TransVal = a;                                             //初始值代入
        }
        

        //计算部
        Iterator b(TransVal, TransVal, Times, Error);                //迭代设置
        b.begin(i);                                                  //迭代开始(i控制读入流场速度值)
        TransVal = b.out();


        //输运部                                                   
        for(int j = 0; j < TV; j++)
        {
            Mcsv(i, j) = TransVal(j);
        }

        //存储部 
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

        ofstream detaFile;//数据处理
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
    /****************************************************************************************/

    }
    
    return 0;
}

