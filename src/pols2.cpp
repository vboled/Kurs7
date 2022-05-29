#include "../include/Kurs7.h"


void Kurs7::pols2() {
    double h = (B - A) / (N);
    double tau = T / (M - 1);
    vector<double> a(N+1), b(N+1), c(N+1), d(N+1), e(N+1), r(N+1), res, polsFF(N+1), polsRR(N+1), sigmaRR(N+1), sigmaFF(N+1);

    ofstream out("pols2.txt");
    ofstream outSigmaFF("sigmaFFPols2.txt");
    // ofstream outSigmaFFAbsError("sigmaFFPolsAbsError2.txt");
    // ofstream outSigmaFFLError("sigmaFFPolsLError2.txt");
    ofstream outSigmaRR("sigmaRRPols2.txt");
    ofstream outSet("polsSettings2.txt");
    outSet << h << endl;
    outSet << A << endl;
    outSet << B << endl;
    outSet << P_a << endl;
    outSet << P_b << endl;
    outSet << E << endl;
    outSet << nu << endl;
    outSet << T << endl;
    outSet << M << endl;

    // double minAbs = 0, minL = 0;
    // double absTime = 0, lTime = 0;

    setSystemUpr2(a, b, c, d, e, r); // задаем систему только с упругостью
    for (int k = 0; k < M; k++) {
        setSystemPols2(r, polsRR, polsFF); // добавляем ползучесть
        // printSystem(a, b, c, d, e, r);
        res = Gauss(a, b, c, d, e, r);
        // outSigmaFFAbsError << tau * k * 1000000<< " ";
        // outSigmaFFLError << tau * k  * 1000000 << " ";
        // double absErrorFF = 0;
        // double sum = 0;
        for (int i = 1; i <= N; i++) {
            if (i % 2) {
                double r = riu(i);
                out << r << " ";

                vector<double> center = coefC(riu(i - 1), riu(i), riu(i + 1));
                double a = center[0], b = center[1];

                vector<double> right = coefR(riu(i - 1), riu(i), riu(i + 1));
                double ar = right[0], br = right[1];

                vector<double> left = coefL(riu(i - 1), riu(i), riu(i + 1));
                double al = left[0], bl = left[1];

                double epsFF = res[i] / r;
                double epsRR = res[i - 1] * (2 * ar * r + br) + res[i] * (2*a * r + b) + res[i + 1] * (2*al * r + bl);

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

                // double tmpMod = fabs(sigmaFF[i] - exactSigmaFFPols(r));
                // absErrorFF = max(absErrorFF, tmpMod);
                // sum += tmpMod * tmpMod;
            }
        }
        
        // double sqrt_tmp = sqrt(h * sum);

        // if (k == 0) {
        //     minAbs = absErrorFF;
        //     minL = sqrt_tmp;
        // } else {
        //     if (minAbs > absErrorFF) {
        //         minAbs = absErrorFF;
        //         absTime = tau * k;
        //     }
        //     if (minL > sqrt_tmp) {
        //         minL = sqrt_tmp;
        //         lTime = k;
        //     }
        // }
        // if (M - k < 4) {
        outSigmaPols2(sigmaFF, outSigmaFF);
        outSigmaPols2(sigmaRR, outSigmaRR);
        // }
        // outSigmaFFAbsError << absErrorFF << endl;
        // outSigmaFFLError << sqrt_tmp << endl;
        // if (k == 2)
        //     break;
    }
    // cout << "AbsError min = " << minAbs << ", time = " << std::fixed<< absTime << endl;
    // cout << "LError min = " << minL << ", time = " << lTime << endl; 
}

void Kurs7::setSystemPols2(vector<double> &r, vector<double> &polsRR, vector<double> &polsFF) {
    for (int i = 0; i < r.size(); i++) {
        r[i] = 0;
        if (i % 2) {
            vector<double> center = coefC(riu(i-1), riu(i), riu(i+1));
            double ar = center[0], br = center[1], cr = center[2];

            double r1i = riu(i - 1), rii = riu(i + 1), ri = riu(i);
            r[i] += -0.3333333333333333*((r1i - rii)*(3*cr*(2*mu()*polsFF[i] + lambda()*(polsFF[i] + polsRR[i])) + 3*br*(lambda() + mu())*(polsFF[i] + polsRR[i])*(r1i + rii) + 
            ar*(3*lambda()*(polsFF[i] + polsRR[i]) + 2*mu()*(polsFF[i] + 2*polsRR[i]))*(Power(r1i,2) + r1i*ri + Power(rii,2))));

        } else {
            if (i > 0) {
                vector<double> center = coefL(riu(i-2), riu(i-1), riu(i));
                double ar = center[0], br = center[1], cr = center[2];

                double r1i = riu(i - 2), rii = riu(i), ri = riu(i-1);
                r[i] += -0.3333333333333333*((r1i - rii)*(3*cr*(2*mu()*polsFF[i-1] + lambda()*(polsFF[i-1] + polsRR[i-1])) + 3*br*(lambda() + mu())*(polsFF[i-1] + polsRR[i-1])*(r1i + rii) + 
                ar*(3*lambda()*(polsFF[i-1] + polsRR[i-1]) + 2*mu()*(polsFF[i-1] + 2*polsRR[i-1]))*(Power(r1i,2) + r1i*ri + Power(rii,2))));
            }
            if (i < r.size()) {
                vector<double> center = coefR(riu(i), riu(i+1), riu(i+2));
                double ar = center[0], br = center[1], cr = center[2];

                double r1i = riu(i), rii = riu(i+2), ri = riu(i+1);
                r[i] += -0.3333333333333333*((r1i - rii)*(3*cr*(2*mu()*polsFF[i+1] + lambda()*(polsFF[i+1] + polsRR[i+1])) + 3*br*(lambda() + mu())*(polsFF[i+1] + polsRR[i+1])*(r1i + rii) + 
                ar*(3*lambda()*(polsFF[i+1] + polsRR[i+1]) + 2*mu()*(polsFF[i+1] + 2*polsRR[i+1]))*(Power(r1i,2) + r1i*ri + Power(rii,2))));
            }
        }
        if (i == 0)
            r[i] += (P_a * A);
        if (i == r.size() - 1) {
            r[i] += -(P_b * B);
        }
    }
}

void Kurs7::outSigmaPols2(vector<double> &res, ofstream &out) {
    for (int i = 1; i <= N; i++) {
        if (i%2)
            out << std::fixed << res[i] << " ";
    }
    out << endl;
}