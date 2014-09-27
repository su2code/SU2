/*!
 * fluid_model_ppr.cpp
 * \brief Source of the Peng-Robinson model.
 * \author: S.Vitale, G.Gori, M.Pini, A.Guardone, P.Colonna
 * \version 3.2.1 "eagle"
 *
 * SU2, Copyright (C) 2012-2014 Aerospace Design Laboratory (ADL).
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


double CPengRobinson::alpha2(double T){

	return ( 1 + k*(1 - sqrt(T/TstarCrit)))*( 1 + k*(1 - sqrt(T/TstarCrit)));
}

void CPengRobinson::SetTDState_rhoe (double rho, double e ) {

    double DpDd_T, DpDT_d,DeDd_T, Cv;
    double A, B, C, sqrt2=sqrt(2), fv, g, g1;

    Density = rho;
    StaticEnergy = e;

    fv = atanh( rho * b * sqrt2/(1 + rho*b));
    A = Gas_Constant / Gamma_Minus_One;
    B = a*k*(k+1)*fv/(b*sqrt2*sqrt(TstarCrit));
    C = - a*(k+1)*(k+1)*fv/(b*sqrt2) - e;
    g = rho * b * sqrt2/(1 + rho*b);
    g1 = b * sqrt2*((1 + rho*b) - rho*b)/((1 + rho*b)*(1 + rho*b));

    Temperature = ( -B + sqrt(B*B - 4*A*C) ) / (2*A); /// Only positive root considered
    Temperature *= Temperature;

    Pressure = rho*Temperature*Gas_Constant / (1 - rho*b) - a*alpha2(Temperature) / ( 1/rho/rho + 2*b/rho - b*b );

    Entropy = A*log(Temperature) + Gas_Constant*log(1/Density -b) - a*sqrt(alpha2(Temperature)) *k*fv/(b*sqrt2*sqrt(Temperature*TstarCrit));

    A = (1/rho/rho + 2*b/rho - b*b);
    B = - 0.5*k/sqrt(Temperature*TstarCrit); //(D alpha / DT)

    DpDd_T = - (- Temperature*Gas_Constant /((1/rho - b)*(1/rho - b)) + 2*a*alpha2(Temperature)*(1/rho + b) /( A*A))/(rho*rho);

    DpDT_d = Gas_Constant /(1/rho - b) - 2*a*sqrt(alpha2(Temperature))*B / A;

    Cv = Gas_Constant/Gamma_Minus_One - a/b/sqrt2 * ( 2*sqrt(alpha2(Temperature)) * B +k*B*sqrt(Temperature/TstarCrit) + 0.5*k*sqrt(alpha2(Temperature)/Temperature/TstarCrit)) * fv;

    dPde_rho = DpDT_d/Cv;

    DeDd_T = g1/((1-g)*(1-g));

    dPdrho_e = DpDd_T - dPde_rho*DeDd_T;

    SoundSpeed2 = dPdrho_e + Pressure/(Density*Density)*dPde_rho;

    dTde_rho = 1/Cv;

    Zed = Pressure/(Gas_Constant*Temperature*Density);


//DTDd_e = Gamma_Minus_One/Gas_Constant*a;
//
//DvDT_p = - DpDT_v / DpDv_T;
//
//    Cv = Gas_Constant/Gamma_Minus_One - a_c/Covolume/sqrt2 * ( 2* sqrt(alpha2) * B +k*B*sqrt(T) + 0.5*k*\sqrt(T)*Crit_Temperature ) * fv; /// R
//
//    Cp = Cv -T*Crit_Temperature * pow(DpDT_v, 2) / DpDv_T;
//
// 	SpeedSound = sqrt( - Cp/Cv * DpDv_T *rho*rho);

}

void CPengRobinson::SetTDState_PT (double P, double T ) {
	double toll= 1e-6;
	double A, B, Z, DZ, F, F1;
	double rho, fv, e;
	double sqrt2=sqrt(2);
	unsigned short nmax = 20, count=0;

//	cout <<"Before  "<< P <<" "<< T <<endl;

	A= a*alpha2(T)*P/(T*Gas_Constant)/(T*Gas_Constant);
	B= b*P/(T*Gas_Constant);

	if(Zed > 0.1)
			Z=min(Zed, 0.99);
		else
			Z=0.99;
	DZ= 1.0;
	do{
		F = Z*Z*Z + Z*Z*(B - 1.0) + Z*(A - 2*B - 3*B*B)  + (B*B*B + B*B - A*B);
		F1 = 3*Z*Z + 2*Z*(B - 1.0) + (A - 2*B - 3*B*B);
		DZ = F/F1;
		Z-= DZ;
	}while(abs(DZ)>toll && count < nmax);

	if (count == nmax){
		cout << "Warning Newton-Raphson exceed number of max iteration in PT"<<endl;
		cout << "Compressibility factor  "<< Z << " would be substituted with "<< Zed<<endl;
	}
	// check if the solution is physical otherwise uses previous point  solution
	if (Z <= 1.0001 && Z >= 0.05 && count < nmax)
	    Zed = Z;


	rho= P/(Zed*Gas_Constant*T);
	fv = atanh( rho * b * sqrt2/(1 + rho*b));

    e = T*Gas_Constant/Gamma_Minus_One + a*k*(k+1)*fv/(b*sqrt2*sqrt(TstarCrit))*sqrt(T) - a*(k+1)*(k+1)*fv/(b*sqrt2);
	SetTDState_rhoe(rho, e);

//	cout <<"After  "<< Pressure <<" "<< Temperature <<endl;

}

void CPengRobinson::SetTDState_Prho (double P, double rho ) {

    SetEnergy_Prho(P,rho);

	SetTDState_rhoe(rho, StaticEnergy);

}

void CPengRobinson::SetTDState_hs (double h, double s ){

	double fv, A, B, C, sqrt2=sqrt(2);
	double T, P, rho;
	double f, v;
	double dv = 1.0;
	double x1,x2, xmid, dx, fx1,fx2, fmid,rtb;
	double toll = 1e-6, FACTOR=0.2;
	unsigned short count=0, NTRY=10, ITMAX=40;

//	 cout <<"Before  "<< h <<" "<< s <<endl;

	A = Gas_Constant / Gamma_Minus_One;

	T = 1.0*h*Gamma_Minus_One/Gas_Constant/Gamma;
	v = exp(-1/Gamma_Minus_One*log(T) + s/Gas_Constant);
	if(Zed<0.9999){
		x1 = Zed*v;
		x2 = v;

	}else{
		x1 = 0.2*v;
		x2 = v;
	}

	rho =1/x1;
	P = rho*T*Gas_Constant / (1 - rho*b) - a*alpha2(T) / ( 1/rho/rho + 2*b/rho - b*b );
	fv = atanh( rho * b * sqrt2/(1 + rho*b));
	B = a*k*(k+1)*fv/(b*sqrt2*sqrt(TstarCrit));
	C = - a*(k+1)*(k+1)*fv/(b*sqrt2) - h + P/rho;
	T = ( -B + sqrt(B*B - 4*A*C) ) / (2*A); /// Only positive root considered
	T *= T;
	fx1 = A*log(T) + Gas_Constant*log(x1 - b) - a*sqrt(alpha2(T)) *k*fv/(b*sqrt2*sqrt(T*TstarCrit)) - s;
	rho =1/x2;
	P = rho*T*Gas_Constant / (1 - rho*b) - a*alpha2(T) / ( 1/rho/rho + 2*b/rho - b*b );
	fv = atanh( rho * b * sqrt2/(1 + rho*b));
	B = a*k*(k+1)*fv/(b*sqrt2*sqrt(TstarCrit));
	C = - a*(k+1)*(k+1)*fv/(b*sqrt2) - h + P/rho;
	T = ( -B + sqrt(B*B - 4*A*C) ) / (2*A); /// Only positive root considered
	T *= T;
	fx2 = A*log(T) + Gas_Constant*log(x2 - b) - a*sqrt(alpha2(T)) *k*fv/(b*sqrt2*sqrt(T*TstarCrit)) - s;

	// zbrac algorithm NR

	for (int j=1;j<=NTRY;j++) {
		if (fx1*fx2 > 0.0){
			if (fabs(fx1) < fabs(fx2)){
				x1 += FACTOR*(x1-x2);
				rho =1/x1;
				P = rho*T*Gas_Constant / (1 - rho*b) - a*alpha2(T) / ( 1/rho/rho + 2*b/rho - b*b );
				fv = atanh( rho * b * sqrt2/(1 + rho*b));
				B = a*k*(k+1)*fv/(b*sqrt2*sqrt(TstarCrit));
				C = - a*(k+1)*(k+1)*fv/(b*sqrt2) - h + P/rho;
				T = ( -B + sqrt(B*B - 4*A*C) ) / (2*A); /// Only positive root considered
				T *= T;
				fx1 = A*log(T) + Gas_Constant*log(x1 - b) - a*sqrt(alpha2(T)) *k*fv/(b*sqrt2*sqrt(T*TstarCrit)) - s;
			}else{
				x2 += FACTOR*(x2-x1);
				P = rho*T*Gas_Constant / (1 - rho*b) - a*alpha2(T) / ( 1/rho/rho + 2*b/rho - b*b );
				fv = atanh( rho * b * sqrt2/(1 + rho*b));
				B = a*k*(k+1)*fv/(b*sqrt2*sqrt(TstarCrit));
				C = - a*(k+1)*(k+1)*fv/(b*sqrt2) - h + P/rho;
				T = ( -B + sqrt(B*B - 4*A*C) ) / (2*A); /// Only positive root considered
				T *= T;
				fx2 = A*log(T) + Gas_Constant*log(x2 - b) - a*sqrt(alpha2(T)) *k*fv/(b*sqrt2*sqrt(T*TstarCrit)) - s;
			}
		}
	}

	  // rtbis algorithm NR

		f=fx1;
		fmid=fx2;
		if (f*fmid >= 0.0){
			cout<< "Root must be bracketed for bisection in rtbis"<<endl;
			SetTDState_rhoT(Density, Temperature);
		}
		rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
		do{
			xmid=rtb+(dx *= 0.5);
			rho =1/xmid;
			P = rho*T*Gas_Constant / (1 - rho*b) - a*alpha2(T) / ( 1/rho/rho + 2*b/rho - b*b );
			fv = atanh( rho * b * sqrt2/(1 + rho*b));
			B = a*k*(k+1)*fv/(b*sqrt2*sqrt(TstarCrit));
			C = - a*(k+1)*(k+1)*fv/(b*sqrt2) - h + P/rho;
			T = ( -B + sqrt(B*B - 4*A*C) ) / (2*A); /// Only positive root considered
			T *= T;
			fmid= A*log(T) + Gas_Constant*log(xmid - b) - a*sqrt(alpha2(T)) *k*fv/(b*sqrt2*sqrt(T*TstarCrit)) - s;
			if (fmid <= 0.0) rtb=xmid;
			count++;
		}while(abs(dx/x1) > toll && count<ITMAX);

		v = xmid;
		if(count==ITMAX){
			cout <<"Too many bisections in rtbis" <<endl;
		}


		rho = 1/v;
		rho =1/xmid;
		P = rho*T*Gas_Constant / (1 - rho*b) - a*alpha2(T) / ( 1/rho/rho + 2*b/rho - b*b );
		fv = atanh( rho * b * sqrt2/(1 + rho*b));
		B = a*k*(k+1)*fv/(b*sqrt2*sqrt(TstarCrit));
		C = - a*(k+1)*(k+1)*fv/(b*sqrt2) - h + P/rho;
		T = ( -B + sqrt(B*B - 4*A*C) ) / (2*A); /// Only positive root considered
		T *= T;

		SetTDState_rhoT(rho, T);

//		cout <<"After  "<< StaticEnergy + Pressure/Density <<" "<< Entropy << fmid <<" "<< f<< " "<< count<<endl;

//	A = Gas_Constant / Gamma_Minus_One;
//
//	T = 1.0*h*Gamma_Minus_One/Gas_Constant/Gamma;
//	v = exp(-1/Gamma_Minus_One*log(T) + s/Gas_Constant);
//	P = T*Gas_Constant / (v - b) - a*alpha2(T) / ( v*v + 2*b*v - b*b);
//	rho =1/v;
//
//	do{
//		fv = atanh( rho * b * sqrt2/(1 + rho*b));
//		B = a*k*(k+1)*fv/(b*sqrt2*sqrt(TstarCrit));
//		C = - a*(k+1)*(k+1)*fv/(b*sqrt2) - h + P/rho;
//		T = ( -B + sqrt(B*B - 4*A*C) ) / (2*A); /// Only positive root considered
//		T *= T;
//
//		f = A*log(T) + Gas_Constant*log(v - b) - a*sqrt(alpha2(T)) *k*fv/(b*sqrt2*sqrt(T*TstarCrit)) - s;
//		f1= Gas_Constant/(v-b)+ a*sqrt(alpha2(T)) *k/(sqrt(T*TstarCrit)*(v*v - b*b - 2*v*b));
//		dv= f/f1;
//		v-= dv;
//		rho = 1/v;
//		P = rho*T*Gas_Constant / (1 - rho*b) - a*alpha2(T) / ( 1/rho/rho + 2*b/rho - b*b );
//
//	}while(abs(dv) > toll);
//
//
//    Z = P/(Gas_Constant*T*rho);
//	// check if the solution is physical otherwise uses previous solution
//	if (Z <= 1.0001 && Z >= 0.05){
//		Zed = Z;
//        SetTDState_rhoT(rho, T);
//	}else{
//		SetTDState_rhoT(Density, Temperature);
//	}

}


void CPengRobinson::SetEnergy_Prho (double P, double rho ) {

	double ad;
    double A, B, C, T,vb1, vb2;
    vb1 = (1/rho -b);
    vb2 = (1/rho/rho + 2*b/rho - b*b);

    A =    Gas_Constant/vb1 - a*k*k/TstarCrit/vb2;

    B =   2*a*k*(k+1)/sqrt(TstarCrit)/vb2;

    C = - P - a*(1+k)*(1+k)/vb2;

    T = ( -B + sqrt(B*B - 4*A*C) ) / (2*A);
    T *= T;


    ad = a*sqrt(alpha2(T)/2)/b*(sqrt(alpha2(T)) + k*sqrt(T/TstarCrit))*atanh( rho * b * sqrt(2)/(1 + rho*b) ) ;
    StaticEnergy = T * Gas_Constant / Gamma_Minus_One - ad;

}

void CPengRobinson::SetTDState_rhoT (double rho, double T){
	double fv,e, sqrt2=sqrt(2);
	fv = atanh( rho * b * sqrt2/(1 + rho*b));
	e = T*Gas_Constant/Gamma_Minus_One + a*k*(k+1)*fv/(b*sqrt2*sqrt(TstarCrit))*sqrt(T) - a*(k+1)*(k+1)*fv/(b*sqrt2);
	SetTDState_rhoe(rho, e);
}




