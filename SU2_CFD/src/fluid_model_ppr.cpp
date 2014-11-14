/*!
 * fluid_model_ppr.cpp
 * \brief Source of the Peng-Robinson model.
 * \author: S.Vitale, G.Gori, M.Pini, A.Guardone, P.Colonna
 * \version 3.2.4 "eagle"
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
	double toll= 1e-9;
	double A, B, Z, DZ, F, F1;
	double rho, fv, e;
	double sqrt2=sqrt(2);
	A= a*alpha2(T)*P/(T*Gas_Constant)/(T*Gas_Constant);
	B= b*P/(T*Gas_Constant);

    Z= max(B, 1.01);
	DZ= 1.0;
	do{
		F = Z*Z*Z + Z*Z*(B - 1.0) + Z*(A - 2*B - 3*B*B)  + (B*B*B + B*B - A*B);
		F1 = 3*Z*Z + 2*Z*(B - 1.0) + (A - 2*B - 3*B*B);
		DZ = F/F1;
		Z-= DZ;
	}while(DZ>toll);

	rho= P/(Z*Gas_Constant*T);
	fv = atanh( rho * b * sqrt2/(1 + rho*b));

    e = T*Gas_Constant/Gamma_Minus_One + a*k*(k+1)*fv/(b*sqrt2*sqrt(TstarCrit))*sqrt(T) - a*(k+1)*(k+1)*fv/(b*sqrt2);
	SetTDState_rhoe(rho, e);
}

void CPengRobinson::SetTDState_Prho (double P, double rho ) {

    SetEnergy_Prho(P,rho);

	SetTDState_rhoe(rho, StaticEnergy);

}

void CPengRobinson::SetTDState_hs (double h, double s ){

	double fv, A, B, C, sqrt2=sqrt(2);
	double f, f1, v;
	double dv = 1.0;
	double toll =1e-9;

	Temperature = h*Gamma_Minus_One/Gas_Constant/Gamma;
	v = exp(-1/Gamma_Minus_One*log(Temperature) + s/Gas_Constant);
	Pressure = Temperature*Gas_Constant / (v - b) - a*alpha2(Temperature) / ( v*v + 2*b*v - b*b);
	Density =1/v;

	do{
		fv = atanh( Density * b * sqrt2/(1 + Density*b));
		A = Gas_Constant / Gamma_Minus_One;
		B = a*k*(k+1)*fv/(b*sqrt2*sqrt(TstarCrit));
		C = - a*(k+1)*(k+1)*fv/(b*sqrt2) - h + Pressure/Density;
		Temperature = ( -B + sqrt(B*B - 4*A*C) ) / (2*A); /// Only positive root considered
		Temperature *= Temperature;

		f = A*log(Temperature) + Gas_Constant*log(1/Density -b) - a*sqrt(alpha2(Temperature)) *k*fv/(b*sqrt2*sqrt(Temperature*TstarCrit)) - s;
		f1= Gas_Constant/(v-b)+ a*sqrt(alpha2(Temperature)) *k/(sqrt(Temperature*TstarCrit)*(v*v - b*b - 2*v*b));
		dv= f/f1;
		v-= dv;
		Density = 1/v;
		Pressure = Density*Temperature*Gas_Constant / (1 - Density*b) - a*alpha2(Temperature) / ( 1/Density/Density + 2*b/Density - b*b );

	}while(abs(dv) > toll);

	double e = Temperature*Gas_Constant/Gamma_Minus_One + a*k*(k+1)*fv/(b*sqrt2*sqrt(TstarCrit))*sqrt(Temperature) - a*(k+1)*(k+1)*fv/(b*sqrt2);
	SetTDState_rhoe(Density, e);


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




