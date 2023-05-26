/*
 * @Author: error: git config user.name && git config user.email & please set dead value or install git
 * @Date: 2022-11-13 18:13:01
 * @LastEditors: error: git config user.name && git config user.email & please set dead value or install git
 * @LastEditTime: 2022-11-22 20:46:34
 * @FilePath: \Siano3.4.1\Exf.cpp
 * @Citrine
 */
#include<iostream>
#include<string>
#include<fstream>
#include<iomanip>
#include<sstream>
#include<Eigen/Dense>
#include<vector>
#include "EXf.h"

vector<double> EXf::Val(int index, string file)
{
    ifstream infile(file, ios::in);
    if (!infile)
    {
        cout << "Could not Open" + file + ".csv" << endl;
        exit(1);
    }
    vector<double> arrdata;
    string line;
    string field;
    for (int i = 0; i <= index; i++)
    {
        getline(infile, line);
    }
    stringstream sin(line);
    int colindex = 0;
    while (getline(sin, field, ','))
    {
        double dvalue = atof(field.c_str());
        arrdata.push_back(dvalue);
    }
    infile.close();

    return arrdata;
}