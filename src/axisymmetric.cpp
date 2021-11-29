#include "../include/Kurs7.h"

void Kurs7::setSystemAxis(vector<double> &a, vector<double> &b, vector<double> &c,
                vector<double> &r) {
    double h = (A - B) / N;
    double lamPlusMu = lambda() + mu();

    for (int i = 0; i < a.size(); i++) {
        if (i == 0) {
            r[0] = (P_a * A);
            b[0] = lamPlusMu * ((ri(0) + ri(1)) / 2 / h + 
            (ri(1) - 3 * ri(0)) / 2 / h + log(ri(1) / ri(0))) -
            lambda() * (3 * ri(1) + ri(0)) / 2 / h;
            c[i] = lamPlusMu / pow(h, 2) * (ri(i) * ri(i + 1) * log(ri(i) / ri(i+1)));
        } else if (i == a.size() - 1){
            r[a.size() - 1] = (P_b * B);
            b[a.size() - 1] = lamPlusMu * ((ri(N) + ri(N-1)) / 2 / h +
            pow(ri(N), 2) * log(ri(N) / ri(N-1)) - 0.5) + lambda();
            a[i] =  lamPlusMu / pow(h, 2) * (ri(i) * ri(i - 1) * log(ri(i - 1) / ri(i)));
        } else {
            a[i] =  lamPlusMu / pow(h, 2) * (ri(i) * ri(i - 1) * log(ri(i - 1) / ri(i)));
            
            b[i] = lamPlusMu / pow(h, 2) * (h* (2 * ri(i) - ri(i-1) - ri(i+1)) / 2 +
            pow(ri(i+1), 2) * log(ri(i+1) / ri(i)) +
            pow(ri(i-1), 2) * log(ri(i-1)/ ri(i)));
            
            c[i] = lamPlusMu / pow(h, 2) * (ri(i) * ri(i + 1) * log(ri(i) / ri(i+1)));
        }
    }
}

void Kurs7::axisymmetric() {
    double h = (A - B) / N;
    vector<double> a(N), b(N), c(N), r(N), res;

    ofstream out("axis.txt");
    ofstream outSet("axisSettings.txt");
    outSet << h << endl;
    outSet << A << endl;
    outSet << B << endl;
    outSet << P_a << endl;
    outSet << P_b << endl;
    outSet << E << endl;
    outSet << nu << endl;
    setSystemAxis(a, b, c, r);
    printSystem(a, b, c, r);
    res = progonka(a, b, c, r);
    outPutRes(res, out);
}