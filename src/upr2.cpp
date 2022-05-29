 #include "../include/Kurs7.h"

void Kurs7::printSystem(vector<double> &a, vector<double> &b, vector<double> &c,
                vector<double> &d, vector<double> &e, vector<double> &r) {
    for (size_t i = 0; i < c.size(); i++) {
        int j = 1, k = c.size() - 1;
        while (++j < i)
            cout << 0 << " ";
        if (i > 1)
             cout << a[i] << " ";
        if (i > 0)
            cout << b[i] << " ";
        cout << c[i] << " ";
        if (i < k)
            cout << d[i] << " ";
        if (i < (k-1))
            cout << e[i] << " ";
        if (i == 1)
            j += 1;
        else if (i != 0)
            j += 2;
        while (k-- > j)
            cout << 0 << " ";
        cout << "| " << r[i] << endl;
    }
    cout << endl;
}

void Kurs7::upr2() {
    double h = (B - A) / (N);
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
    res = progonka(a, b, c, d, e, r);
    outPutResUpr(res, out);
    outSigmaRR2(res);
    outSigmaFF2(res);
}

vector<double> Kurs7::coefC(double right, double center, double left) {
    vector<double> res(3);

    res[0] = 1.0/((center - left) * (center - right));
    res[1] = -((left + right)/((center - left) * (center - right)));
    res[2] = (left * right)/((center - left) * (center - right));
    return res;
}

vector<double> Kurs7::coefR(double right, double center, double left) {
    vector<double> res(3);

    res[0] = -(1.0/((center - right) * (-left + right)));
    res[1] = -((-center - left)/((center - right) * (-left + right)));
    res[2] = -((center * left)/((center - right) * (-left + right)));
    return res;
}

vector<double> Kurs7::coefL(double right, double center, double left) {
    vector<double> res(3);

    res[0] =-(1.0/((center - left) * (left - right)));
    res[1] = -((-center - right)/((center - left) * (left - right)));
    res[2] = -((center * right)/((center - left) * (left - right)));
    return res;
}

void Kurs7::setSystemUpr2(vector<double> &a, vector<double> &b, vector<double> &c, vector<double> &d, vector<double> &e, vector<double> &r) {
    double h = (B - A) / (N);

    double lamPlusMu = lambda() + 2.0 * mu();
    
    for (int i = 0; i < a.size(); i++) {
        if (i == 0) {
            r[0] = (P_a * A);
            r[a.size() - 1] = -(P_b * B);
        }
            a[i] = 0;
            b[i] = 0;
            c[i] = 0;
            d[i] = 0;
            e[i] = 0;
        if (i % 2) {
            double Aa = riu(i-1), Ba = riu(i+1), l = lambda(), m = mu();

            vector<double> center = coefC(riu(i-1), riu(i), riu(i+1));
            double ai = center[0], bi = center[1], ci = center[2];

            vector<double> right = coefR(riu(i-1), riu(i), riu(i+1));
            double ria = right[0], rib = right[1], ric = right[2];

            vector<double> left = coefL(riu(i-1), riu(i), riu(i+1));
            double lia = left[0], lib = left[1], lic = left[2];

            b[i] += -0.16666666666666666*((l + 2*mu())*(6*Power(Aa,4)*ai*ria + 3*Power(Aa,2)*bi*rib + 4*Power(Aa,3)*(bi*ria + ai*rib) - 
        Power(Ba,2)*(6*ai*Power(Ba,2)*ria + 4*Ba*bi*ria + 4*ai*Ba*rib + 3*bi*rib))) - 
   (l*(3*Power(Aa,4)*ai*ria + 6*Aa*ci*rib + 2*Power(Aa,3)*(2*bi*ria + ai*rib) + 3*Power(Aa,2)*(2*ci*ria + bi*rib) - 
        Ba*(Ba*(3*ai*Power(Ba,2) + 4*Ba*bi + 6*ci)*ria + (2*ai*Power(Ba,2) + 3*Ba*bi + 6*ci)*rib)))/6. - 
   (l*(3*Power(Aa,4)*ai*ria + 2*Power(Aa,3)*(bi*ria + 2*ai*rib) + 6*Aa*bi*ric + 3*Power(Aa,2)*(bi*rib + 2*ai*ric) - 
        Ba*(Ba*(bi*(2*Ba*ria + 3*rib) + ai*Ba*(3*Ba*ria + 4*rib)) + 6*(ai*Ba + bi)*ric)))/6. + 
   (l + 2*mu())*(-0.25*(Power(Aa,4)*ai*ria) + (ai*Power(Ba,4)*ria)/4. - (Power(Aa,3)*(bi*ria + ai*rib))/3. + (Power(Ba,3)*(bi*ria + ai*rib))/3. - 
      (Power(Aa,2)*(ci*ria + bi*rib + ai*ric))/2. + (Power(Ba,2)*(ci*ria + bi*rib + ai*ric))/2. - Aa*(ci*rib + bi*ric) + Ba*(ci*rib + bi*ric) - 
      ci*ric*log(Aa) + ci*ric*log(Ba));

            c[i] += -((Aa - Ba)*(ai*(Aa + Ba) + bi)*(ai*(Power(Aa,2) + Power(Ba,2)) + (Aa + Ba)*bi + 2*ci)*l) - 
   ((6*Power(ai,2)*(Power(Aa,4) - Power(Ba,4)) + 8*ai*(Power(Aa,3) - Power(Ba,3))*bi + 3*(Aa - Ba)*(Aa + Ba)*Power(bi,2))*(l + 2*mu()))/6. + 
   (l + 2*mu())*(-0.08333333333333333*(Aa*(6*Aa*Power(bi,2) + 8*bi*(Power(Aa,2)*ai + 3*ci) + 3*Aa*ai*(Power(Aa,2)*ai + 4*ci))) + 
      (Ba*(6*Ba*Power(bi,2) + 8*bi*(ai*Power(Ba,2) + 3*ci) + 3*ai*Ba*(ai*Power(Ba,2) + 4*ci)))/12. - Power(ci,2)*log(Aa) + 
      Power(ci,2)*log(Ba));

            d[i] += -0.16666666666666666*(l*(3*Power(Aa,4)*ai*lia + 6*Aa*ci*lib + 2*Power(Aa,3)*(2*bi*lia + ai*lib) + 3*Power(Aa,2)*(2*ci*lia + bi*lib) - 
        Ba*(Ba*(3*ai*Power(Ba,2) + 4*Ba*bi + 6*ci)*lia + (2*ai*Power(Ba,2) + 3*Ba*bi + 6*ci)*lib))) - 
   (l*(3*Power(Aa,4)*ai*lia + 2*Power(Aa,3)*(bi*lia + 2*ai*lib) + 6*Aa*bi*lic + 3*Power(Aa,2)*(bi*lib + 2*ai*lic) - 
        Ba*(Ba*(bi*(2*Ba*lia + 3*lib) + ai*Ba*(3*Ba*lia + 4*lib)) + 6*(ai*Ba + bi)*lic)))/6. - 
   ((6*Power(Aa,4)*ai*lia + 3*Power(Aa,2)*bi*lib + 4*Power(Aa,3)*(bi*lia + ai*lib) - 
        Power(Ba,2)*(6*ai*Power(Ba,2)*lia + 4*Ba*bi*lia + 4*ai*Ba*lib + 3*bi*lib))*(l + 2*mu()))/6. + 
   (l + 2*mu())*(-0.25*(Power(Aa,4)*ai*lia) + (ai*Power(Ba,4)*lia)/4. - (Power(Aa,3)*(bi*lia + ai*lib))/3. + (Power(Ba,3)*(bi*lia + ai*lib))/3. - 
      (Power(Aa,2)*(ci*lia + bi*lib + ai*lic))/2. + (Power(Ba,2)*(ci*lia + bi*lib + ai*lic))/2. - Aa*(ci*lib + bi*lic) + Ba*(ci*lib + bi*lic) - 
      ci*lic*log(Aa) + ci*lic*log(Ba));
        } else {
            if (i > 0) {
                double Aa = riu(i-2), Bb = riu(i), l =lambda(), m = mu();

                vector<double> center = coefC(riu(i-2), riu(i-1), riu(i));
                double ccca = center[0], cccb = center[1], cccc = center[2];
                vector<double> right = coefR(riu(i-2), riu(i-1), riu(i));
                double rrra = right[0], rrrb = right[1], rrcc = right[2];

                vector<double> left = coefL(riu(i-2), riu(i-1), riu(i));
                double llla = left[0], lllb = left[1], lllc = left[2];

                a[i] += -0.16666666666666666*(l*(6*Aa*lllb*rrcc + 3*Power(Aa,4)*llla*rrra + 2*Power(Aa,3)*(lllb*rrra + 2*llla*rrrb) + 
        3*Power(Aa,2)*(2*llla*rrcc + lllb*rrrb) - Bb*(6*Bb*llla*rrcc + 6*lllb*rrcc + 3*Power(Bb,3)*llla*rrra + 2*Power(Bb,2)*lllb*rrra + 
           Bb*(4*Bb*llla + 3*lllb)*rrrb))) + (l*(-3*Power(Aa,4)*llla*rrra - 6*Aa*lllc*rrrb - 2*Power(Aa,3)*(2*lllb*rrra + llla*rrrb) - 
        3*Power(Aa,2)*(2*lllc*rrra + lllb*rrrb) + Bb*(Bb*(3*Power(Bb,2)*llla + 4*Bb*lllb + 6*lllc)*rrra + 
           (2*Power(Bb,2)*llla + 3*Bb*lllb + 6*lllc)*rrrb)))/6. - 
   ((l + 2*m)*(6*Power(Aa,4)*llla*rrra + 3*Power(Aa,2)*lllb*rrrb + 4*Power(Aa,3)*(lllb*rrra + llla*rrrb) - 
        Power(Bb,2)*(6*Power(Bb,2)*llla*rrra + 3*lllb*rrrb + 4*Bb*(lllb*rrra + llla*rrrb))))/6. + 
   (l + 2*m)*(-0.25*(Power(Aa,4)*llla*rrra) + (Power(Bb,4)*llla*rrra)/4. - (Power(Aa,3)*(lllb*rrra + llla*rrrb))/3. + 
      (Power(Bb,3)*(lllb*rrra + llla*rrrb))/3. - (Power(Aa,2)*(llla*rrcc + lllc*rrra + lllb*rrrb))/2. + 
      (Power(Bb,2)*(llla*rrcc + lllc*rrra + lllb*rrrb))/2. - Aa*(lllb*rrcc + lllc*rrrb) + Bb*(lllb*rrcc + lllc*rrrb) - lllc*rrcc*log(Aa) + 
      lllc*rrcc*log(Bb));

                b[i] += (l*(-3*Power(Aa,4)*ccca*llla - 6*Aa*cccc*lllb - 2*Power(Aa,3)*(2*cccb*llla + ccca*lllb) - 3*Power(Aa,2)*(2*cccc*llla + cccb*lllb) + 
        Bb*(Bb*(3*Power(Bb,2)*ccca + 4*Bb*cccb + 6*cccc)*llla + (2*Power(Bb,2)*ccca + 3*Bb*cccb + 6*cccc)*lllb)))/6. + 
   (l*(-3*Power(Aa,4)*ccca*llla - 2*Power(Aa,3)*(cccb*llla + 2*ccca*lllb) - 6*Aa*cccb*lllc - 3*Power(Aa,2)*(cccb*lllb + 2*ccca*lllc) + 
        Bb*(Bb*(3*Power(Bb,2)*ccca*llla + 2*Bb*cccb*llla + 4*Bb*ccca*lllb + 3*cccb*lllb) + 6*(Bb*ccca + cccb)*lllc)))/6. - 
   ((6*Power(Aa,4)*ccca*llla + 3*Power(Aa,2)*cccb*lllb + 4*Power(Aa,3)*(cccb*llla + ccca*lllb) - 
        Power(Bb,2)*(6*Power(Bb,2)*ccca*llla + 3*cccb*lllb + 4*Bb*(cccb*llla + ccca*lllb)))*(l + 2*m))/6. + 
   (l + 2*m)*(-0.25*(Power(Aa,4)*ccca*llla) + (Power(Bb,4)*ccca*llla)/4. - (Power(Aa,3)*(cccb*llla + ccca*lllb))/3. + 
      (Power(Bb,3)*(cccb*llla + ccca*lllb))/3. - (Power(Aa,2)*(cccc*llla + cccb*lllb + ccca*lllc))/2. + 
      (Power(Bb,2)*(cccc*llla + cccb*lllb + ccca*lllc))/2. - Aa*(cccc*lllb + cccb*lllc) + Bb*(cccc*lllb + cccb*lllc) - cccc*lllc*log(Aa) + 
      cccc*lllc*log(Bb));

                c[i] += -((Aa - Bb)*l*((Aa + Bb)*llla + lllb)*((Power(Aa,2) + Power(Bb,2))*llla + (Aa + Bb)*lllb + 2*lllc)) - 
   ((6*(Power(Aa,4) - Power(Bb,4))*Power(llla,2) + 8*(Power(Aa,3) - Power(Bb,3))*llla*lllb + 3*(Aa - Bb)*(Aa + Bb)*Power(lllb,2))*(l + 2*m))/
    6. + (l + 2*m)*(-0.08333333333333333*(Aa*(6*Aa*Power(lllb,2) + 8*lllb*(Power(Aa,2)*llla + 3*lllc) + 
           3*Aa*llla*(Power(Aa,2)*llla + 4*lllc))) + (Bb*(6*Bb*Power(lllb,2) + 8*lllb*(Power(Bb,2)*llla + 3*lllc) + 
           3*Bb*llla*(Power(Bb,2)*llla + 4*lllc)))/12. - Power(lllc,2)*log(Aa) + Power(lllc,2)*log(Bb));

            }
            if (i < a.size() - 1) {
                double Aa = riu(i), Bb = riu(i + 2), l = lambda(), m = mu();

                vector<double> center = coefC(riu(i), riu(i+1), riu(i+2));
                double cca = center[0], ccb = center[1], ccc = center[2];

                vector<double> right = coefR(riu(i), riu(i+1), riu(i+2));
                double rra = right[0], rrb = right[1], rrc = right[2];

                vector<double> left = coefL(riu(i), riu(i+1), riu(i+2));
                double lla = left[0], llb = left[1], llc = left[2];

                c[i] += -0.16666666666666666*((l + 2*m)*(6*(Power(Aa,4) - Power(Bb,4))*Power(rra,2) + 8*(Power(Aa,3) - Power(Bb,3))*rra*rrb + 
        3*(Aa - Bb)*(Aa + Bb)*Power(rrb,2))) - (Aa - Bb)*l*((Aa + Bb)*rra + rrb)*((Power(Aa,2) + Power(Bb,2))*rra + (Aa + Bb)*rrb + 2*rrc) + 
   (l + 2*m)*(-0.08333333333333333*(Aa*(3*Power(Aa,3)*Power(rra,2) + 8*Power(Aa,2)*rra*rrb + 24*rrb*rrc + 6*Aa*(Power(rrb,2) + 2*rra*rrc))) + 
      (Bb*(3*Power(Bb,3)*Power(rra,2) + 8*Power(Bb,2)*rra*rrb + 24*rrb*rrc + 6*Bb*(Power(rrb,2) + 2*rra*rrc)))/12. - 
      Power(rrc,2)*log(Aa) + Power(rrc,2)*log(Bb));

                d[i] +=(l*(-3*Power(Aa,4)*cca*rra - 6*Aa*ccc*rrb - 2*Power(Aa,3)*(2*ccb*rra + cca*rrb) - 3*Power(Aa,2)*(2*ccc*rra + ccb*rrb) + 
        Bb*(Bb*(3*Power(Bb,2)*cca + 4*Bb*ccb + 6*ccc)*rra + (2*Power(Bb,2)*cca + 3*Bb*ccb + 6*ccc)*rrb)))/6. - 
   ((l + 2*m)*(6*Power(Aa,4)*cca*rra + 3*Power(Aa,2)*ccb*rrb + 4*Power(Aa,3)*(ccb*rra + cca*rrb) - 
        Power(Bb,2)*(6*Power(Bb,2)*cca*rra + 3*ccb*rrb + 4*Bb*(ccb*rra + cca*rrb))))/6. + 
   (l*(-3*Power(Aa,4)*cca*rra - 2*Power(Aa,3)*(ccb*rra + 2*cca*rrb) - 6*Aa*ccb*rrc - 3*Power(Aa,2)*(ccb*rrb + 2*cca*rrc) + 
        Bb*(Bb*(3*Power(Bb,2)*cca*rra + 2*Bb*ccb*rra + 4*Bb*cca*rrb + 3*ccb*rrb) + 6*(Bb*cca + ccb)*rrc)))/6. + 
   (l + 2*m)*(-0.25*(Power(Aa,4)*cca*rra) + (Power(Bb,4)*cca*rra)/4. - (Power(Aa,3)*(ccb*rra + cca*rrb))/3. + 
      (Power(Bb,3)*(ccb*rra + cca*rrb))/3. - (Power(Aa,2)*(ccc*rra + ccb*rrb + cca*rrc))/2. + (Power(Bb,2)*(ccc*rra + ccb*rrb + cca*rrc))/2. - 
      Aa*(ccc*rrb + ccb*rrc) + Bb*(ccc*rrb + ccb*rrc) - ccc*rrc*log(Aa) + ccc*rrc*log(Bb)); 

                e[i] += (l*(-3*Power(Aa,4)*lla*rra - 6*Aa*llc*rrb - 2*Power(Aa,3)*(2*llb*rra + lla*rrb) - 3*Power(Aa,2)*(2*llc*rra + llb*rrb) + 
        Bb*(Bb*(3*Power(Bb,2)*lla + 4*Bb*llb + 6*llc)*rra + (2*Power(Bb,2)*lla + 3*Bb*llb + 6*llc)*rrb)))/6. - 
   ((l + 2*m)*(6*Power(Aa,4)*lla*rra + 3*Power(Aa,2)*llb*rrb + 4*Power(Aa,3)*(llb*rra + lla*rrb) - 
        Power(Bb,2)*(6*Power(Bb,2)*lla*rra + 3*llb*rrb + 4*Bb*(llb*rra + lla*rrb))))/6. + 
   (l*(-3*Power(Aa,4)*lla*rra - 2*Power(Aa,3)*(llb*rra + 2*lla*rrb) - 6*Aa*llb*rrc - 3*Power(Aa,2)*(llb*rrb + 2*lla*rrc) + 
        Bb*(Bb*(3*Power(Bb,2)*lla*rra + 2*Bb*llb*rra + 4*Bb*lla*rrb + 3*llb*rrb) + 6*(Bb*lla + llb)*rrc)))/6. + 
   (l + 2*m)*(-0.25*(Power(Aa,4)*lla*rra) + (Power(Bb,4)*lla*rra)/4. - (Power(Aa,3)*(llb*rra + lla*rrb))/3. + 
      (Power(Bb,3)*(llb*rra + lla*rrb))/3. - (Power(Aa,2)*(llc*rra + llb*rrb + lla*rrc))/2. + (Power(Bb,2)*(llc*rra + llb*rrb + lla*rrc))/2. - 
      Aa*(llc*rrb + llb*rrc) + Bb*(llc*rrb + llb*rrc) - llc*rrc*log(Aa) + llc*rrc*log(Bb));
            }

        }
    }
}