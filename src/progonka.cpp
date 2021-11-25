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