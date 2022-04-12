#include "../include/Kurs7.h"

void Kurs7::setSystemPols(vector<double> &a, vector<double> &b, vector<double> &c,
                vector<double> &r, vector<double> &polsRR, vector<double> &polsFF) {
    for (int i = 0; i < r.size(); i++) {
        r[i] = 0;
        // Из Математики: mu (-polsFF + polsRR) r1i + (lambda + mu) (polsFF + polsRR) ri
        if (i != 0) {
            r[i] += mu() * (-polsFF[i] + polsRR[i]) * ri(i-1) + (lambda() + mu()) * (polsFF[i] + polsRR[i]) * ri(i);
        }
        
        // Из Математики:   -((lambda + mu) (PolsFF + PolsRR) ri) + mu (PolsFF - PolsRR) ri1
        if (i != a.size() - 1) {
            r[i] += -((lambda() + mu()) * (polsFF[i+1] + polsRR[i+1]) * ri(i)) + mu() * (polsFF[i+1] - polsRR[i+1]) * ri(i+1);
        }
        if (i == 0)
            r[i] += (P_a * A);
        if (i == a.size() - 1) {
            r[i] += -(P_b * B);
            // cout << P_b * B << endl;
        }
    }
}

void Kurs7::outSigmaPols(vector<double> &res, ofstream &out) {
    for (int i = 1; i < N; i++) {
        out << std::fixed << res[i] << " ";
    }
    out << endl;
}

void Kurs7::pols() {
    double h = (B - A) / N;
    double tau = T / (M - 1);
    vector<double> a(N), b(N), c(N), r(N), res, polsFF(N), polsRR(N), sigmaRR(N), sigmaFF(N);

    ofstream out("pols.txt");
    ofstream outSigmaFF("sigmaFFPols.txt");
    ofstream outSigmaRR("sigmaRRPols.txt");
    ofstream outSet("polsSettings.txt");
    outSet << h << endl;
    outSet << A << endl;
    outSet << B << endl;
    outSet << P_a << endl;
    outSet << P_b << endl;
    outSet << E << endl;
    outSet << nu << endl;
    outSet << T << endl;
    outSet << M << endl;

    setSystemAxis(a, b, c, r); // задаем систему только с упругостью
    for (int k = 0; k < M; k++) {
        setSystemPols(a, b, c, r, polsRR, polsFF); // добавляем ползучесть
        // printSystem(a, b, c, r);
        res = progonka(a, b, c, r);

        for (int i = 1; i < N; i++) {

            double r = A + (i - 0.5) * h; // точка в между i-1 и i узлами
            // cout << r << endl;

            // Считаем компоненты упругой деформации
            double epsRR = (res[i] - res[i - 1]) / h;
            
            double epsFF = 0.5 * (res[i] + res[i-1]) / r;

            // Считаем компоненты напряжения
            sigmaRR[i] = (lambda() + 2.0 * mu()) * (epsRR - polsRR[i]) + lambda() * (epsFF - polsFF[i]);
            sigmaFF[i] = lambda() * (epsRR - polsRR[i]) + (lambda() + 2.0 * mu()) * (epsFF - polsFF[i]);
            
            // if (i == 3)
            //     cout << "sigmaRR = " << sigmaRR[i] << " sigmaFF = " << (epsRR - polsRR[i]) << endl;

            double sigmaRR_ = (sigmaRR[i] - sigmaFF[i]) / 2;
            double sigmaFF_ = (sigmaFF[i] - sigmaRR[i]) / 2;

            // Считаем компоненты деформации ползучести

            double tmp = 9.0 / 4.0 * OMEGA * tau * (sigmaRR_ * sigmaRR_ +  sigmaFF_ *  sigmaFF_);

            polsRR[i] += tmp * sigmaRR_;
            polsFF[i] += tmp * sigmaFF_;

            // polsRR[i] = 10e-4 * tau * (k + 1);
            // // cout << polsRR[i] << endl;
            // polsFF[i] = polsRR[i];
            // polsRR[i] = 10e-4 * 
            // if (i == 3)
            //     cout << "polsRR = " << polsRR[i] << " polsFF = " << polsFF[i] << endl;
            // cout << endl;
        }
        outSigmaPols(sigmaFF, outSigmaFF);
        outSigmaPols(sigmaRR, outSigmaRR);
        // cout << endl;
        // break;
    }
}