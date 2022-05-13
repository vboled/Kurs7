 #include "../include/Kurs7.h"

void Kurs7::printSystem(vector<double> &a, vector<double> &b, vector<double> &c,
                vector<double> &d, vector<double> &e, vector<double> &r) {
    double h = (B - A) / N;
    for (size_t i = 0; i < c.size(); i++) {
        int j = 0, k = c.size() - 1;
        while (++j < i)
            cout << 0 << " ";
        if (i)
            cout << a[i] << " ";
        cout << b[i] << " ";
        if (i != k)
            cout << c[i] << " ";
        while (k-- > j)
            cout << 0 << " ";
        cout << "| " << r[i] << endl;
    }
    cout << endl;
}

void Kurs7::upr2() {
    double h = (B - A) / (N+1);
    vector<double> a(N+1), b(N+1), c(N+1), d(N+1), e(N+1), r(N+1), res;

    ofstream out("upr2.txt");
    ofstream outSet("upr2Settings.txt");
    outSet << h << endl;
    outSet << A << endl;
    outSet << B << endl;
    outSet << P_a << endl;
    outSet << P_b << endl;
    outSet << E << endl;
    outSet << nu << endl;
    setSystemUpr2(a, b, c, d, e, r);
    printSystem(a, b, c, d, e, r);
    // res = progonka(a, b, c, r);
    // res = Gauss(a, b, c, r);
    // outPutRes(res, out);
    // outSigmaRR(res);
    // outSigmaFF(res);
}

void Kurs7::setSystemUpr2(vector<double> &a, vector<double> &b, vector<double> &c, vector<double> &d, vector<double> &e, vector<double> &r) {
    double h = (B - A) / (N + 1);
    double lamPlusMu = lambda() + 2 * mu();
    ofstream outa("a.txt");
    ofstream outb("b.txt");
    ofstream outc("c.txt");
    for (int i = 0; i < a.size(); i++) {
        if (!(i % 2)) {
            if (i > 0) {
                r[i] = -2;
            } else {
                r[i] = -1;
            }
        }
        else 
            r[i] = 1;
    }
}