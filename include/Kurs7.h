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

    double A = 3;
    double B = 4;
    double P_a = 10e5;
    double P_b = 10e7;
    double nu = 0.25;
    double E = 2e9;

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
        return A + (A - B) / N * i;
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