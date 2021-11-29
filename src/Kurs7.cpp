#include "../include/Kurs7.h"

void Kurs7::printSystem(vector<double> &a, vector<double> &b, vector<double> &c,
                vector<double> &r) {
    double step = (A - B) / N;
    for (size_t i = 0; i < c.size(); i++) {
        int j = 0, k = c.size() - 1;
        while (j++ < i)
            cout << "0 ";
        cout << a[i] << " " << b[i] << " " << c[i] << " ";
        while (k-- >= j)
            cout << "0 ";
        cout << "| " << r[i] << endl;
    }
}

void Kurs7::outPutRes(vector<double> &res, ofstream &out) {
    out << 0 << " ";
    for (int i = 0; i < res.size(); i++) {
        out << res[i] << " ";
    }
    out << endl;
}
