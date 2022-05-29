#include "../include/Kurs7.h"

void Kurs7::printSystem(vector<double> &a, vector<double> &b, vector<double> &c,
                vector<double> &r) {
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

void Kurs7::outPutResUpr(vector<double> &res, ofstream &out) {
    double step = (B - A) / (N);
    double error = 0;
    for (int i = 0; i < res.size(); i++) {
        out << A + step * (i) << " " << res[i] << endl;
        double tmp = fabs(exacSolUpr(A + step * (i)) - res[i]);
        if (tmp > error)
            error = tmp;
    }
    cout << "U error = " << error << endl;
}

void Kurs7::outPutRes(vector<double> &res, ofstream &out) {
    double step = (B - A) / (N-1);

    for (int i = 0; i < res.size(); i++) {
        out << A + step * (i) << " " << res[i] << endl;
        
    }
}
