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