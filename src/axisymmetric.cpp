#include "../include/Kurs7.h"

void Kurs7::setSystemAxis(vector<double> &a, vector<double> &b, vector<double> &c,
                vector<double> &r) {
    double h = (B - A) / N;
    double lamPlusMu = lambda() + 2 * mu();
    ofstream outa("a.txt");
    ofstream outb("b.txt");
    ofstream outc("c.txt");
    for (int i = 0; i < a.size(); i++) {
        if (i == 0) {
            r[0] = -(P_a * A);
            r[a.size() - 1] = -(P_b * B);
            b[0] = -lambda() - lamPlusMu + (lamPlusMu * pow(ri(i+1), 2) * log(ri(i+1) / ri(i)) / pow(ri(i) - ri(i+1), 2));
            c[0] = lamPlusMu * ri(i) * ri(i+1) * log(ri(i) / ri(i+1)) / pow(ri(i) - ri(i+1), 2);
        } else {
            a[i] = lamPlusMu * ri(i-1) * ri(i) *log(ri(i-1) / ri(i)) / pow(ri(i-1) - ri(i), 2);
            c[i] = lamPlusMu * ri(i) * ri(i+1) * log(ri(i) / ri(i+1)) / pow(ri(i) - ri(i+1), 2);
            b[i] = lambda() + lamPlusMu + (lamPlusMu * pow(ri(i-1), 2) * log(ri(i) / ri(i-1)) / pow(ri(i-1) - ri(i),2));
            if (i != a.size() - 1) 
                b[i] += -lambda() - lamPlusMu + (lamPlusMu * pow(ri(i+1), 2) * log(ri(i+1) / ri(i)) / pow(ri(i) - ri(i+1), 2)); 
        }
        if (i == a.size() - 1)
            c[i] = 0.0;
    }
}

double Kurs7::exactSigmaFF(double r) {
    return (P_a * A * A - P_b * B * B) / 
           (B * B - A * A) + (A * A * B * B) / 
           pow(r, 2) * (P_a - P_b) / (B * B - A * A);
}

double Kurs7::exactSigmaRR(double r) {
    return (P_a  * A * A - P_b * B * B) / 
    (B * B - A * A) - (A * A * B * B) / 
    pow(r, 2) * (P_a - P_b) / (B * B - A * A);
}

void Kurs7::outSigmaRR(vector<double> &res) {
    ofstream out("sigmaRR.txt");
    double error = 0;
    double h = (B - A) / N;
    for (int i = 1; i < N; i++) {
        double r = A + (i - 0.5) * h;
        out << r << " ";
        double tmp = (lambda() + 2 * mu()) * (res[i] / h - res[i - 1] / h) + lambda() * 0.5 * 
        (res[i] + res[i-1]) / r;
        out << tmp << "\n";
        if (abs(tmp - exactSigmaRR(r)) > error) 
            error = abs(tmp - exactSigmaRR(r));
    }
    cout << "sigmaRR Error = " << error << endl;
}

void Kurs7::outSigmaFF(vector<double> &res) {
    ofstream out("sigmaFF.txt");
    double h = (B - A) / N;
    double error = 0;
    for (int i = 1; i < N; i++) {
        double r = A + (i - 0.5) * h;
        out << r << " ";
        double tmp = lambda() * (res[i] / h - res[i - 1] / h) + (lambda() + 2 * mu()) * 0.5 * 
        (res[i] + res[i-1]) / r;
        out << tmp << "\n";
        if (abs(tmp - exactSigmaFF(r)) > error) 
            error = abs(tmp - exactSigmaFF(r));
    }
    cout << "sigmaFF Error = " << error << endl;
}

void Kurs7::axisymmetric() {
    double h = (B - A) / N;
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
    // printSystem(a, b, c, r);
    res = progonka(a, b, c, r);
    // res = Gauss(a, b, c, r);
    outPutRes(res, out);
    outSigmaRR(res);
    outSigmaFF(res);
}