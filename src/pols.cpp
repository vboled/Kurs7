#include "../include/Kurs7.h"

void Kurs7::setSystemPols(vector<double> &a, vector<double> &b, vector<double> &c,
                vector<double> &r, vector<double> &polsRR, vector<double> &polsFF) {
    for (int i = 0; i < r.size(); i++) {
        // if (i != 0) {
        //     r[i] += ((A - B) * ((A + B) * (lambda() * (1 + lambda() + 2 * mu()) *polsFF[i] + 
        //             2 * (lambda() + mu()) * polsRR[i+1]) - 2 * lambda() * (lambda() * polsFF[i] + 2 * mu() * polsFF[i] + polsRR[i]) * ri(i-1)))
        //             /(2 * (ri(i-1) - ri(i)));
        // }
            
        // if (i != a.size() - 1) {
        //     r[i] += ((A - B) * ((A + B) * (-lambda() * (1 + lambda()+ 2 * mu()) * polsFF[i+1] - 
        //             2 * (lambda()+ mu()) * polsRR[i+1]) + 2 * lambda() * (lambda() * polsFF[i+1] + 2 * mu() * polsFF[i+1] + polsRR[i+1]) * ri(i+1)))
        //             /(2 * (ri(i) - ri(i+1)));
        // }

        // if (i != 0) {
        //     r[i] += -(1/2) * lambda() * (-1 + lambda() + 2 * mu()) * polsFF[i] * ri(i-1) + mu() * polsRR[i] * ri(i-1) + 
        //                 1/2 * lambda() * (1 + lambda() + 2 * mu()) * polsFF[i] * ri(i) + (lambda() + mu()) * polsRR[i] * ri(i);
        // }
            
        // if (i != a.size() - 1) {
        //     r[i] += 1/2 * (-lambda() * (1 + lambda() + 2 * mu()) * polsFF[i+1] * ri(i) - 2 *(lambda() + mu()) * polsRR[i+1] * ri(i) + 
        //             lambda() * (-1 + lambda() + 2 * mu()) * polsFF[i+1] * ri(i + 1) - 2 * mu() * polsRR[i + 1] * ri(i+1));
        //     if (i == 0)
        //         cout << "ri = " << r[i] << endl;
        // }

        if (i != 0) {
            r[i] += mu() * (-polsFF[i] + polsRR[i]) * ri(i-1) + (lambda() + mu()) * (polsFF[i] + polsRR[i]) * ri(i);
        }
            
        if (i != a.size() - 1) {
            r[i] += -((lambda() + mu()) * (polsFF[i+1] + polsRR[i+1]) * ri(i)) + mu() * (polsFF[i] - polsRR[i]) * ri(i+1);
        }
            
    }
}

void Kurs7::pols() {
    double h = (B - A) / N;
    double tau = T / (M - 1);
    vector<double> a(N), b(N), c(N), r(N), res, polsFF(N), polsRR(N);

    ofstream out("pols.txt");
    ofstream outSet("polsSettings.txt");
    outSet << h << endl;
    outSet << A << endl;
    outSet << B << endl;
    outSet << P_a << endl;
    outSet << P_b << endl;
    outSet << E << endl;
    outSet << nu << endl;
    outSet << T << endl;
    setSystemAxis(a, b, c, r);
    for (int i = 0; i < M; i++) {
        setSystemPols(a, b, c, r, polsRR, polsFF);
        printSystem(a, b, c, r);
        res = progonka(a, b, c, r);

        outPutRes(res, out);

        for (int i = 0; i < N; i++) {
            double r = A + (i - 0.5) * h;

            double epsFF = (lambda() + 2 * mu()) * (res[i] / h - res[i - 1] / h) + lambda() * 0.5 * 
                            (res[i] + res[i-1]) / r;
            double epsRR = lambda() * (res[i] / h - res[i - 1] / h) + (lambda() + 2 * mu()) * 0.5 * 
                            (res[i] + res[i-1]) / r;

            double sigmaRR = (lambda() + 2 * mu()) * (epsRR - polsRR[i]) + lambda() * (epsFF - polsFF[i]);
            double sigmaFF = lambda() * (lambda() + 2 * mu()) * (epsRR - polsRR[i]) + (lambda() + 2 * mu()) * (epsFF - polsFF[i]);

            polsRR[i] += 9.0 / 4.0 * OMEGA * tau * (sigmaRR * sigmaRR + sigmaFF * sigmaFF) * (sigmaRR - sigmaFF) / 2;
            polsFF[i] += 9.0 / 4.0 * OMEGA * tau * (sigmaRR * sigmaRR + sigmaFF * sigmaFF) * (sigmaFF - sigmaRR) / 2;
        }
        if (i == 2)
            break;
    }
    
    // res = Gauss(a, b, c, r);
}