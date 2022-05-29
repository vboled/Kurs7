#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <fstream>
#include <iomanip>  
#include <iterator>

using namespace std;

class Kurs7 {

private:

    double N = 10;
    double T = 50;
    double M = T * 10;

    double A = 0.03;
    double B = 0.04;
    double P_a = 10e5;
    double P_b = 0;
    double nu = 0.3;
    double E = 2e11;
    double OMEGA = 6 * 10e-27;
    // double OMEGA = 2.3 * 10e-27;
    // 
    vector<double> progonka(vector<double> a, vector<double> b,
								vector<double> c, vector<double> r);
    vector<double> Gauss(vector<double> a, vector<double> b,
                                vector<double> c, vector<double> f);
    vector<double> Gauss(vector<double> a, vector<double> b, vector<double> c,
                    vector<double> d, vector<double> e,vector<double> f);
    vector<double> progonka(vector<double> a,
                                 vector<double> b,
                                 vector<double> c,
                                 vector<double> d,
                                 vector<double> e,
                                 vector<double> f);

    double relativeError();
    double absoluteError();

    double lambda() {
        return nu * E / (1.0 + nu) / (1.0 - 2.0 * nu);
    }
    
    double mu() {
        return E / 2.0 / (1.0 + nu);
    }

    double ri(int i) {
        return A + 1.0 * (B - A) / (N-1) * i;
    }

    double riu(int i) {
        return A + 1.0 * (B - A) / (N) * i;
    }

    
    double Power(double num, int st) {
        double res = 1;
        for (int i = 0; i < st; i++) {
            res *= num;
        }
        return res;
    }
    
    vector<double> coefC(double right, double center, double left);
    vector<double> coefR(double right, double center, double left);
    vector<double> coefL(double right, double center, double left);
    
    void printSystem(vector<double> &a, vector<double> &b, vector<double> &c,
                vector<double> &r);
    void outPutRes(vector<double> &res, ofstream &out);
    void outPutResUpr(vector<double> &res, ofstream &out);
    void outSigmaRR(vector<double> &res);
    void outSigmaFF(vector<double> &res);
    void outSigmaRR2(vector<double> &res);
    void outSigmaFF2(vector<double> &res);
    void outSigmaPols(vector<double> &res, ofstream &out);
    double exactSigmaRR(double r);
    double exactSigmaFF(double r);
    double exactSigmaRRPols(double r);
    double exactSigmaFFPols(double r);
    double exacSolUpr(double r);
public:
    // 
    void axisymmetric();
    void setSystemAxis(vector<double> &a, vector<double> &b, vector<double> &c,
                vector<double> &r);
        
    void pols();
    void setSystemPols(vector<double> &a, vector<double> &b, vector<double> &c,
                vector<double> &r, vector<double> &polsRR, vector<double> &polsFF);

    void upr2();
    void pols2();
    void setSystemUpr2(vector<double> &a, vector<double> &b, vector<double> &c,
                vector<double> &d, vector<double> &e, vector<double> &r);
    void printSystem(vector<double> &a, vector<double> &b, vector<double> &c,
                vector<double> &d, vector<double> &e, vector<double> &r);
    void setSystemPols2(vector<double> &r, vector<double> &polsRR, vector<double> &polsFF);
    void outSigmaPols2(vector<double> &res, ofstream &out);
};