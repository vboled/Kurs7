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

    double N = 200;
    double M = 100;

    double A = 0.03;
    double B = 0.04;
    double P_a = 10e4;
    double P_b = 0;
    double nu = 0.3;
    double E = 2e8;
    double T = 1;
    double OMEGA = 6 * 10e-21;
    // 
    vector<double> progonka(vector<double> a, vector<double> b,
								vector<double> c, vector<double> r);
    vector<double> Gauss(vector<double> a, vector<double> b,
                                vector<double> c, vector<double> f);
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
    void outSigmaRR(vector<double> &res);
    void outSigmaFF(vector<double> &res);
    void outSigmaPols(vector<double> &res, ofstream &out);
    double exactSigmaRR(double r);
    double exactSigmaFF(double r);

public:
    // 
    void axisymmetric();
    void setSystemAxis(vector<double> &a, vector<double> &b, vector<double> &c,
                vector<double> &r);
        
    void pols();
    void setSystemPols(vector<double> &a, vector<double> &b, vector<double> &c,
                vector<double> &r, vector<double> &polsRR, vector<double> &polsFF);

};