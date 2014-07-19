/*!
 * fluid_model_ppr.cpp
 * \brief Source of the Peng-Robinson model.
 * \author: S.Vitale, G.Gori, M.Pini, A.Guardone, P.Colonna
 * \version 3.2.0 "eagle"
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

//	a= 0.0;
//	b =0.0;

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

    A = (1/rho/rho + 2*b/rho - b*b);
    B = - 0.5*k/sqrt(Temperature*TstarCrit); //(D alpha / DT)

    DpDd_T = - (- Temperature*Gas_Constant /((1/rho - b)*(1/rho - b)) + 2*a*alpha2(Temperature)*(1/rho + b) /( A*A))/(rho*rho);

    DpDT_d = Gas_Constant /(1/rho - b) - 2*a*sqrt(alpha2(Temperature))*B / A;

    Cv = Gas_Constant/Gamma_Minus_One - a/b/sqrt2 * ( 2*sqrt(alpha2(Temperature)) * B +k*B*sqrt(Temperature/TstarCrit) + 0.5*k*sqrt(alpha2(Temperature)/Temperature/TstarCrit)) * fv;

    DpDe_d = DpDT_d/Cv;

    DeDd_T = g1/((1-g)*(1-g));

    DpDd_e = DpDd_T - DpDe_d*DeDd_T;

    SoundSpeed2 = DpDd_e + Pressure/(Density*Density)*DpDe_d;

    DTDe_d = 1/Cv;
    //DTDd_e = Gamma_Minus_One/Gas_Constant*a;

    //DvDT_p = - DpDT_v / DpDv_T;

//    Cv = Gas_Constant/Gamma_Minus_One - a_c/Covolume/sqrt2 * ( 2* sqrt(alpha2) * B +k*B*sqrt(T) + 0.5*k*\sqrt(T)*Crit_Temperature ) * fv; /// R
//
//    Cp = Cv -T*Crit_Temperature * pow(DpDT_v, 2) / DpDv_T;
//
// 	SpeedSound = sqrt( - Cp/Cv * DpDv_T *rho*rho);

}

void CPengRobinson::SetTDState_PT (double P, double T ) {
//
//    /// T E P MUST BE Tr e Pr
//
//	double a1;
//	double a2;
//	double a3;
//	double A, B;
//	double alpha;
//	double* root  = new double [3];
//	double rho;
//	double ad;
//	double e;
//
//	alpha = pow(  ( 1 + k*(1 - sqrt(T)) ) , 2)
//    Mol_Attraction = a_c * alpha;
//
//	A = Mol_Attraction*P/pow(T * Gas_Constant, 2);
//	B = Covolume*P/(T* Gas_Constant);
//
//	a1 = - ( 1 - B );
//	a2 =   ( A - 3*B -2*B );
//	a3 = - ( A*B -B*B - B*B*B );
//
//	poly_solve_cubic(a1, a2, a3, root);
//
//    // check the root of the polynomial
//
//	rho = P / ( root[2]*T*Gas_Constant ) * Crit_Density;
//
//    ad = a_c/Covolume/sqrt(2) * k*k * atanh( rho * Covolume * sqrt(2)/(1 + rho*Covolume) ) ;
//    ad *= (1 + 2*k + k*k) + (3*k - 1) * k *sqrt(T);
//    e = T*Crit_Temperature*Gas_Constant / Gamma_Minus_One - ad;
//
//	SetTDState_de(rho, e);
//
//	delete [] root;
//
}

void CPengRobinson::SetTDState_Prho (double P, double rho ) {

    SetEnergy_Prho(P,rho);

	SetTDState_rhoe(rho, StaticEnergy);

}

void CPengRobinson::SetTDState_hs (double h, double s ){

/// WARNING TODO

}

//void CPengRobinson::SetTDState_Ps (double P, double s ) {
//	/// Working
//
//    double rho = 0.0;
//
//	SetTDState_Pd(P, rho);
//
//}

void CPengRobinson::SetEnergy_Prho (double P, double rho ) {

	double ad;
    double A, B, C, T,vb1, vb2;
    vb1 = (1/rho -b);
    vb2 = (1/rho/rho + 2*b/rho - b*b);

    A =    Gas_Constant/vb1 - a*k*k/TstarCrit/vb2;

    B =   2*a*k*(k+1)/sqrt(TstarCrit)/vb2;

    C = - P - a*(1+k)*(1+k)/vb2;

    T = ( -B + sqrt(B*B - 4*A*C) ) / (2*A); /// Only positive root considered
    T *= T;


    ad = a*sqrt(alpha2(T)/2)/b*(sqrt(alpha2(T)) + k*sqrt(T/TstarCrit))*atanh( rho * b * sqrt(2)/(1 + rho*b) ) ;
    StaticEnergy = T * Gas_Constant / Gamma_Minus_One - ad;

}





//// function for solving a cubic equation
//
//#define SWAP(a,b) do { double tmp = b ; b = a ; a = tmp ; } while(0)
//
//double poly_solve_cubic (double a, double b, double c,
//                      double *root) {
//
//  // the equation is in the form: x^3 + a x^2 + b x + c = 0
//
//  double q = (a * a - 3 * b);
//  double r = (2 * a * a * a - 9 * a * b + 27 * c);
//
//  double Q = q / 9;
//  double R = r / 54;
//
//  double Q3 = Q * Q * Q;
//  double R2 = R * R;
//
//  double CR2 = 729 * r * r;
//  double CQ3 = 2916 * q * q * q;
//
//  if (R == 0 && Q == 0)
//    {
//      root[0] = - a / 3 ;
//      root[1] = - a / 3 ;
//      root[2] = - a / 3 ;
//      return 3 ;
//    }
//  else if (CR2 == CQ3)
//    {
//      /* this test is actually R2 == Q3, written in a form suitable
//         for exact computation with integers */
//
//      /* Due to finite precision some double roots may be missed, and
//         considered to be a pair of complex roots z = x +/- epsilon i
//         close to the real axis. */
//
//      double sqrtQ = sqrt (Q);
//
//      if (R > 0)
//        {
//    	  root[0] = -2 * sqrtQ  - a / 3;
//    	  root[1] = sqrtQ - a / 3;
//    	  root[2] = sqrtQ - a / 3;
//        }
//      else
//        {
//    	  root[0] = - sqrtQ  - a / 3;
//    	  root[1] = - sqrtQ - a / 3;
//    	  root[2] = 2 * sqrtQ - a / 3;
//        }
//      return 3 ;
//    }
//  else if (R2 < Q3)
//    {
//      double sgnR = (R >= 0 ? 1 : -1);
//      double ratio = sgnR * sqrt (R2 / Q3);
//      double theta = acos (ratio);
//      double norm = -2 * sqrt (Q);
//      root[0] = norm * cos (theta / 3) - a / 3;
//      root[1] = norm * cos ((theta + 2.0 * M_PI) / 3) - a / 3;
//      root[2] = norm * cos ((theta - 2.0 * M_PI) / 3) - a / 3;
//
//      /* Sort root[0], root[1], root[2] into increasing order */
//
//      if (root[0] > root[1])
//        SWAP(root[0], root[1]) ;
//
//      if (root[1] > root[2])
//        {
//          SWAP(root[1], root[2]) ;
//
//          if (root[0] > root[1])
//            SWAP(root[0], root[1]) ;
//        }
//
//      return 3;
//    }
//  else
//    {
//      double sgnR = (R >= 0 ? 1 : -1);
//      double A = -sgnR * pow (fabs (R) + sqrt (R2 - Q3), 1.0/3.0);
//      double B = Q / A ;
//      root[0] = A + B - a / 3;
//      return 1;
//    }
//}
//
//





