#include "../include/Kurs7.h"

vector<double> Kurs7::progonka(vector<double> a,
                                 vector<double> b,
                                 vector<double> c,
                                 vector<double> f) {
    vector<double> res(a.size());
    double m;
	for (int i = 1; i < a.size(); i++)
	{
		m = a[i]/b[i-1];
		b[i] = b[i] - m*c[i-1];
		f[i] = f[i] - m*f[i-1];
	}
    
	res[a.size()-1] = f[a.size()-1]/b[a.size()-1];

	for (int i = a.size() - 2; i >= 0; i--)
    {
		res[i]=(f[i]-c[i]*res[i+1])/b[i];
    }        
    
    return res;
}

vector<double> Kurs7::progonka(vector<double> a,
                                 vector<double> b,
                                 vector<double> c,
                                 vector<double> d,
                                 vector<double> e,
                                 vector<double> f) {
    vector<double> res(a.size()), alpha(a.size()), beta(a.size()), gamma(a.size());

    for (int i = 0; i < a.size(); i++) {
        b[i] *= -1;
        d[i] *= -1;
    }
    
    double C = c[0], D = d[0], B = b[1], Q = c[1], T = b[2], R = 0, A = a[3], F = f[0], J = f[1], G = f[2], H = f[3], S = a[2];
    vector<double> kappa(a.size()), eta(a.size()), sigma(a.size());
    kappa[0] = 0; eta[0] = 1;

    double N = a.size() - 1;
    for (int i = 0; i <= N - 2; i++) {
        if (fabs(C) >= fabs(D) && fabs(C) >= fabs(e[i])) {
            alpha[i+1] = D / C;
            beta[i+1] = e[i] / C;
            gamma[i+1] = F / C;

            C = Q - B * alpha[i+1];
            D = d[i+1] - B * beta[i+1];
            F = J + B * gamma[i+1];

            B = T - S * alpha[i+1];
            Q = c[i+2] - S * beta[i+1];
            J = G - S * gamma[i+1];

            if (i != N - 2) {
                S = A - R * alpha[i+1];
                T = b[i+3] - R * beta[i+1];
                G = H + R * gamma[i+1];
            }
            if (i < N - 3) {
                R = 0; 
                A = a[i+4];
                H = f[i+4];
            }
            sigma[i+1] = kappa[i];
            kappa[i+1] = eta[i];
            eta[i+1] = i + 2;
        } else if (fabs(D) > fabs(C) && fabs(D) >= fabs(e[i])) {
            alpha[i+1] = C / D;
            beta[i+1] = -e[i] / D;
            gamma[i+1] = -F / D;

            C = -B + Q * alpha[i+1];
            D = d[i+1] + Q * beta[i+1];
            F = J - Q * gamma[i+1];

            B = -S + T * alpha[i+1];
            Q = c[i+2] + T * beta[i+1];
            J = G
             + T * gamma[i+1];
            if (i != N - 2) {
                S = -R + A * alpha[i+1];
                T = b[i+3] + A * beta[i+1];
                G = H - A * gamma[i+1];
            }
            if (i < N - 3) {
                R = 0; 
                A = a[i+4];
                H = f[i+4];
            }
            sigma[i+1] = eta[i];
            kappa[i+1] = kappa[i];
            eta[i+1] = i + 2;
        } else if (fabs(e[i]) > C && fabs(e[i]) > fabs(D)) {
            alpha[i+1] = D / e[i];
            beta[i+1] = C / e[i];
            gamma[i+1] = F / e[i];

            C = Q - d[i+1] * alpha[i+1];
            D = B - d[i+1] * beta[i+1];
            F = J + d[i+1] * gamma[i+1];

            B = T - c[i+2] * alpha[i+1];
            Q = S - c[i+2] * beta[i+1];
            J = G - c[i+2] * gamma[i+1];
            if (i != N - 2) {
                S = A - b[i+3] * alpha[i+1];
                T = R - b[i+3] * beta[i+1];
                G = H + b[i+3] * gamma[i+1];
            }
            if (i < N - 3) {
                R = -a[i+4] * alpha[i+1]; 
                A = -a[i+4] * beta[i+1];
                H = f[i+4] - a[i+4]*gamma[i+1];
            }
            sigma[i+1] = i+2;
            kappa[i+1] = eta[i];
            eta[i+1] = kappa[i];
        } else
            cout << "Error" << endl;
    }
    double gamma_last = 0;
    if (fabs(C) >= fabs(D)) {
        alpha[a.size() - 1] = D / C;
        gamma[a.size() - 1] = F / C;
        gamma_last = (J + B * gamma[a.size() - 1])/(Q - B * alpha[a.size() - 1]);
        sigma[a.size() - 1] = kappa[a.size() - 2];
        kappa[a.size() - 1] = eta[a.size() - 2];
    } else {
        alpha[a.size() - 1] = C / D;
        gamma[a.size() - 1] = -F / D;
        gamma_last = (J - Q * gamma[a.size() - 1])/(-B + Q * alpha[a.size() - 1]);
        sigma[a.size() - 1] = eta[a.size() - 2];
        kappa[a.size() - 1] = kappa[a.size() - 2];
    }
    res[kappa[N]] = gamma_last;
    res[sigma[N]] = alpha[N] * res[kappa[N]] + gamma[N];

    for (int i = N - 2; i >= 0; i--) {
        res[sigma[i+1]] = alpha[i + 1] * res[kappa[i+1]] - 
            beta[i+1] * res[eta[i+1]] + gamma[i+1];
    }
    return res;
}

vector<double> Kurs7::Gauss(vector<double> a,
                            vector<double> b,
                            vector<double> c,
                            vector<double> d,
                            vector<double> e,
                            vector<double> f) {

	vector<vector<double>> aa(a.size());
	vector<double> res(a.size());
	for (int i = 0; i < a.size(); i++) {
		aa[i].resize(a.size() + 1);
        if (i > 1)
            aa[i][i-2] = a[i];
        if (i > 0)
            aa[i][i-1] = b[i];
		aa[i][i] = c[i];
        if (i < a.size())
            aa[i][i+1] = d[i];
        if (i < a.size() - 1)
            aa[i][i+2] = e[i];    
		aa[i][a.size()] = f[i];
		// for (int j = 0; j < a.size(); j++) {
		// 	cout << aa[i][j] << " ";
		// }
		// cout << " | " <<aa[i][a.size()] << endl;
	}
    // cout << endl;
	double aaa;
	for (int k = 0; k < a.size(); k++) //Поиск максимального элемента в первом столбце
    {
        aaa = abs(aa[k][k]);
		double bb;
        int i = k;
        for(int m = k+1; m < a.size(); m++)
            if(abs(aa[m][k])>aaa)
            {
                i = m;
                aaa = abs(aa[m][k]);
            }
 
            if (aaa == 0)   //проверка на нулевой элемент
            {
                cout<<"Система не имеет решений"<<endl;
            }
 
            if (i != k)  //  перестановка i-ой строки, содержащей главный элемент k-ой строки
            {
                for (int j=k; j < a.size()+1; j++)
                {
                    bb = aa[k][j];
                    aa[k][j] = aa[i][j];
                    aa[i][j] = bb;
                }
            }
            aaa = aa[k][k];//преобразование k-ой строки (Вычисление масштабирующих множителей)
            aa[k][k] = 1;   
            for (int j=k+1;j<a.size()+1;j++) 
                aa[k][j] = aa[k][j]/aaa;
            for (i = k+1; i < a.size(); i++)//преобразование строк с помощью k-ой строки
            {
                bb = aa[i][k];
                aa[i][k] = 0;
                if (bb!=0)
                    for (int j=k+1; j< a.size()+1; j++)
                        aa[i][j]=aa[i][j]-bb*aa[k][j];
            }
    }
 
    for (int i=a.size()-1; i>=0; i--)   //Нахождение решений СЛАУ
    {
        res[i] = 0;
        aaa = aa[i][a.size()];
        for (int j = a.size(); j>i; j--) 
            aaa = aaa-aa[i][j]*res[j];
        res[i] = aaa;
    }
	return res;
}

vector<double> Kurs7::Gauss(vector<double> a,
                                 vector<double> b,
                                 vector<double> c,
                                 vector<double> f) {
	vector<vector<double>> aa(a.size());
	vector<double> res(a.size());
	for (int i = 0; i < a.size(); i++) {
		aa[i].resize(a.size() + 1);
		if (i != 0)
			aa[i][i-1] = a[i];
		if (i < a.size() - 1) 
			aa[i][i + 1] = c[i];
		aa[i][i] = b[i];
		aa[i][a.size()] = f[i];
		// for (int j = 0; j < a.size(); j++) {
		// 	cout << aa[i][j] << " ";
		// }
		// cout << " | " <<aa[i][a.size()] << endl;
	}
	double aaa;
	for (int k = 0; k < a.size(); k++) //Поиск максимального элемента в первом столбце
    {
        aaa = abs(aa[k][k]);
		double bb;
        int i = k;
        for(int m = k+1; m < a.size(); m++)
            if(abs(aa[m][k])>aaa)
            {
                i = m;
                aaa = abs(aa[m][k]);
            }
 
            if (aaa == 0)   //проверка на нулевой элемент
            {
                cout<<"Система не имеет решений"<<endl;
            }
 
            if (i != k)  //  перестановка i-ой строки, содержащей главный элемент k-ой строки
            {
                for (int j=k; j < a.size()+1; j++)
                {
                    bb = aa[k][j];
                    aa[k][j] = aa[i][j];
                    aa[i][j] = bb;
                }
            }
            aaa = aa[k][k];//преобразование k-ой строки (Вычисление масштабирующих множителей)
            aa[k][k] = 1;   
            for (int j=k+1;j<a.size()+1;j++) 
                aa[k][j] = aa[k][j]/aaa;
            for (i = k+1; i < a.size(); i++)//преобразование строк с помощью k-ой строки
            {
                bb = aa[i][k];
                aa[i][k] = 0;
                if (bb!=0)
                    for (int j=k+1; j< a.size()+1; j++)
                        aa[i][j]=aa[i][j]-bb*aa[k][j];
            }
    }
 
    for (int i=a.size()-1; i>=0; i--)   //Нахождение решений СЛАУ
    {
        res[i] = 0;
        aaa = aa[i][a.size()];
        for (int j = a.size(); j>i; j--) 
            aaa = aaa-aa[i][j]*res[j];
        res[i] = aaa;
    }
	return res;
}