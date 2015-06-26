/*!
 * fluid_model_ppr.cpp
 * \brief Source of the Peng-Robinson model.
 * \author S. Vitale, G. Gori, M. Pini, A. Guardone, P. Colonna
 * \version 4.0.0 "Cardinal"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "./../include/fluid_model.hpp"

CPengRobinson::CPengRobinson() : CIdealGas() {
	a= 0.0;
	b =0.0;
	k = 0.0;
	TstarCrit = 0.0;
}

CPengRobinson::CPengRobinson(double gamma, double R, double Pstar, double Tstar, double w) : CIdealGas(gamma, R) {

	a = 0.45724*Gas_Constant*Gas_Constant*Tstar*Tstar/Pstar;
	b = 0.0778*Gas_Constant*Tstar/Pstar;
	TstarCrit = Tstar;
	Zed=1.0;

	if (w <= 0.49)
        k = 0.37464 + 1.54226 * w - 0.26992 * w*w;
        else
        k = 0.379642 + 1.48503 * w - 0.164423 * w*w + 0.016666 * w*w*w;


}

CPengRobinson::~CPengRobinson(void) { }


double CPengRobinson::alpha2(double T) {

	return ( 1 + k*(1 - sqrt(T/TstarCrit)))*( 1 + k*(1 - sqrt(T/TstarCrit)));
}

double CPengRobinson::T_v_h(double v, double h) {
	double fv, A, B, C, T, d;
	double sqrt2=sqrt(2.0);

	d = (v*v+2*b*v-b*b);
	fv = atanh( b*sqrt2 / (v + b));

	A = Gas_Constant*(1 / Gamma_Minus_One + v/(v-b)) - a*v*k*k / (TstarCrit * d);
	B = a*k*(k+1)/sqrt(TstarCrit) *( fv/(b*sqrt2) + 2*v/d );
	C = h + a*(1+k)*(1+k)*(fv/(b*sqrt2) + v/d);

	T = ( -B + sqrt(B*B + 4*A*C) ) / (2*A); /// Only positive root considered

	return T*T;
}

void CPengRobinson::SetTDState_rhoe (double rho, double e ) {

    double DpDd_T, DpDT_d, DeDd_T, Cv;
    double A, B, C, sqrt2, fv, a2T, rho2;

    Density = rho;
    StaticEnergy = e;

    rho2 = rho*rho;
    sqrt2=sqrt(2);

    fv = atanh( rho * b * sqrt2/(1 + rho*b));
    
    A = Gas_Constant / Gamma_Minus_One;
    B = a*k*(k+1)*fv/(b*sqrt2*sqrt(TstarCrit));
    C = a*(k+1)*(k+1)*fv/(b*sqrt2) + e;

    Temperature = ( -B + sqrt(B*B + 4*A*C) ) / (2*A); /// Only positive root considered
    Temperature *= Temperature;

    a2T = alpha2(Temperature);

    A = (1/rho2 + 2*b/rho - b*b);
    B = 1/rho-b;

    Pressure = Temperature*Gas_Constant / B - a*a2T / A;

    Entropy = Gas_Constant / Gamma_Minus_One*log(Temperature) + Gas_Constant*log(B) - a*sqrt(a2T) *k*fv/(b*sqrt2*sqrt(Temperature*TstarCrit));

    DpDd_T =  ( Temperature*Gas_Constant /(B*B    ) - 2*a*a2T*(1/rho + b) /( A*A ) ) /(rho2);

    DpDT_d = Gas_Constant /B + a*k / A * sqrt( a2T/(Temperature*TstarCrit) );

    Cv = Gas_Constant/Gamma_Minus_One + ( a*k*(k+1)*fv ) / ( 2*b*sqrt(2*Temperature*TstarCrit) );

    dPde_rho = DpDT_d/Cv;

    DeDd_T = - a*(1+k) * sqrt( a2T ) / A / (rho2);

    dPdrho_e = DpDd_T - dPde_rho*DeDd_T;

    SoundSpeed2 = dPdrho_e + Pressure/(rho2)*dPde_rho;

    dTde_rho = 1/Cv;

    Zed = Pressure/(Gas_Constant*Temperature*Density);


}

void CPengRobinson::SetTDState_PT (double P, double T ) {
	double toll= 1e-6;
	double A, B, Z, DZ=1.0, F, F1;
	double rho, fv, e;
	double sqrt2=sqrt(2);
	unsigned short nmax = 20, count=0;

	A= a*alpha2(T)*P/(T*Gas_Constant)/(T*Gas_Constant);
	B= b*P/(T*Gas_Constant);

  if (Zed > 0.1) Z = min(Zed, 0.99);
		else Z=0.99;
  
	do {
		F = Z*Z*Z + Z*Z*(B - 1.0) + Z*(A - 2*B - 3*B*B)  + (B*B*B + B*B - A*B);
		F1 = 3*Z*Z + 2*Z*(B - 1.0) + (A - 2*B - 3*B*B);
		DZ = F/F1;
		Z-= DZ;
	} while(abs(DZ)>toll && count < nmax);

	if (count == nmax) {
		cout << "Warning Newton-Raphson exceed number of max iteration in PT"<< endl;
		cout << "Compressibility factor  "<< Z << " would be substituted with "<< Zed<< endl;
	}
	// check if the solution is physical otherwise uses previous point  solution
	if (Z <= 1.0001 && Z >= 0.05 && count < nmax)
	    Zed = Z;


	rho= P/(Zed*Gas_Constant*T);
	fv = atanh( rho * b * sqrt2/(1 + rho*b));

        e = T*Gas_Constant/Gamma_Minus_One - a*(k+1)*sqrt( alpha2(T) )*fv / (b*sqrt2);

	SetTDState_rhoe(rho, e);
}

void CPengRobinson::SetTDState_Prho (double P, double rho ) {

    SetEnergy_Prho(P, rho);

	SetTDState_rhoe(rho, StaticEnergy);

}

void CPengRobinson::SetTDState_hs (double h, double s ) {

	double T, fv, sqrt2=sqrt(2.0), A;
	double f, v;
	double x1, x2, xmid, dx, fx1, fx2, fmid, rtb;
	double toll = 1e-9, FACTOR=0.2;
	double cons_s, cons_h;
	unsigned short countrtb=0, NTRY=10, ITMAX=100;

	A = Gas_Constant / Gamma_Minus_One;
	T = h*Gamma_Minus_One/Gas_Constant/Gamma;
	v = exp(-1/Gamma_Minus_One*log(T) + s/Gas_Constant);


	if (Zed<0.9999) {
		x1 = Zed*v;
		x2 = v;

	} else{
		x1 = 0.2*v;
		x2 = v;
	}


	T = T_v_h(x1, h);
	fv = atanh( b*sqrt2 / (x1 + b));
	fx1 = A*log(T) + Gas_Constant*log(x1 - b) - a*sqrt(alpha2(T)) *k*fv/(b*sqrt2*sqrt(T*TstarCrit)) - s;
	T = T_v_h(x2, h);
	fv = atanh( b*sqrt2 / (x2 + b));
	fx2 = A*log(T) + Gas_Constant*log(x2 - b) - a*sqrt(alpha2(T)) *k*fv/(b*sqrt2*sqrt(T*TstarCrit)) - s;

	// zbrac algorithm NR

	for (int j=1; j<=NTRY; j++) {
		if (fx1*fx2 > 0.0) {
			if (fabs(fx1) < fabs(fx2)) {
				x1 += FACTOR*(x1-x2);
				T = T_v_h(x1, h);
				fv = atanh( b*sqrt2/(x1 + b));
				fx1 = A*log(T) + Gas_Constant*log(x1 - b) - a*sqrt(alpha2(T)) *k*fv/(b*sqrt2*sqrt(T*TstarCrit)) - s;
			} else{
				x2 += FACTOR*(x2-x1);
				T = T_v_h(x2, h);
				fv = atanh( b*sqrt2/(x2 + b));
				fx2 = A*log(T) + Gas_Constant*log(x2 - b) - a*sqrt(alpha2(T)) *k*fv/(b*sqrt2*sqrt(T*TstarCrit)) - s;
			}
		}
	}

	// rtbis algorithm NR

	f=fx1;
	fmid=fx2;
	if (f*fmid >= 0.0) {
		cout<< "Root must be bracketed for bisection in rtbis"<< endl;
		SetTDState_rhoT(Density, Temperature);
	}
	rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
	do{
		xmid=rtb+(dx *= 0.5);
		T = T_v_h(xmid, h);
		fv = atanh( b* sqrt2/(xmid + b));
		fmid= A*log(T) + Gas_Constant*log(xmid - b) - a*sqrt(alpha2(T)) *k*fv/(b*sqrt2*sqrt(T*TstarCrit)) - s;

		if (fmid <= 0.0) rtb=xmid;
		countrtb++;
	}while(abs(fmid) > toll && countrtb<ITMAX);

	v = xmid;
	if (countrtb==ITMAX) {
		cout <<"Too many bisections in rtbis" << endl;
//			do{
//					fv = atanh( b/v* sqrt2/(1 + b/v));
//					T=T_v_h(v, h);
//					f = A*log(T) + Gas_Constant*log(v - b) - a*sqrt(alpha2(T)) *k*fv/(b*sqrt2*sqrt(T*TstarCrit)) - s;
//					f1= Gas_Constant/(v-b)+ a*sqrt(alpha2(T)) *k/(sqrt(T*TstarCrit)*(v*v - b*b - 2*v*b));
//					dv= f/f1;
//					v-= dv;
//					countnw++;
//			}while(abs(f/x2) > toll && countnw<ITMAXNW);
//
//		} else{
	}
	if (v!=v) {
		cout <<"not physical solution found, h and s input " << h << " "<< s << endl;
		SetTDState_rhoT(Density, Temperature);
	}

	T=T_v_h(v, h);
	SetTDState_rhoT(1/v, T);

	// consistency check
	cons_h= abs(((StaticEnergy + Pressure/Density) - h)/h);
	cons_s= abs((Entropy-s)/s);

	if (cons_h >1e-4 or cons_s >1e-4) {
		cout<< "TD consistency not verified in hs call"<< endl;
			 //cout <<"Before  "<< h <<" "<< s << endl;
			 //cout <<"After  "<< StaticEnergy + Pressure/Density <<" "<< Entropy << fmid <<" "<< f<< " "<< countrtb<<" "<< countnw<< endl;
			 //getchar();
	}
}


void CPengRobinson::SetEnergy_Prho (double P, double rho) {

    double ad;
    double A, B, C, T, vb1, vb2;
    vb1 = (1/rho -b);
    vb2 = (1/rho/rho + 2*b/rho - b*b);

    A =   Gas_Constant/vb1 - a*k*k/TstarCrit/vb2;

    B =   2*a*k*(k+1)/sqrt(TstarCrit)/vb2;

    C = - P - a*(1+k)*(1+k)/vb2;

    T = ( -B + sqrt(B*B - 4*A*C) ) / (2*A);
    T *= T;

    ad = a*(k+1)*sqrt( alpha2(T) ) / ( b*sqrt(2) ) * atanh( rho * b * sqrt(2)/(1 + rho*b) ) ;

    StaticEnergy = T * Gas_Constant / Gamma_Minus_One - ad;

}

void CPengRobinson::SetTDState_rhoT (double rho, double T) {
	double fv, e;

	fv = atanh( rho * b * sqrt(2)/(1 + rho*b));
	e = T*Gas_Constant/Gamma_Minus_One - a*(k+1)*sqrt( alpha2(T) ) / ( b*sqrt(2) ) * fv;
	SetTDState_rhoe(rho, e);
}




