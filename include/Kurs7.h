#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <fstream>
#include <iomanip>  

using namespace std;

class Kurs7 {

private:

    double N = 10;

    double A = 0.03;
    double B = 0.04;
    double P_a = 0;
    double P_b = 10e6;
    double nu = 0.3;
    double E = 200e9;

    // 
    vector<double> progonka(vector<double> a, vector<double> b,
								vector<double> c, vector<double> r);
    double relativeError();
    double absoluteError();

    double lambda() {
        return nu * E / (1 + nu) / (1 - 2 * nu);
    }
    
    double mu() {
        return E / 2 / (1 + nu);
    }

    double ri(int i) {
        return A + (B - A) / (N-1) * i;
    }
    void printSystem(vector<double> &a, vector<double> &b, vector<double> &c,
                vector<double> &r);
    void outPutRes(vector<double> &res, ofstream &out);

public:
    // 
    void axisymmetric();
    void setSystemAxis(vector<double> &a, vector<double> &b, vector<double> &c,
                vector<double> &r);


};