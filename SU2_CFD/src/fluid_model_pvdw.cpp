/*!
 * fluid_model_pvdw.cpp
 * \brief Source of the Polytropic Van der Waals model.
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

#include "../include/fluid_model.hpp"

CVanDerWaalsGas::CVanDerWaalsGas() : CIdealGas() {
	a = 0.0;
	b = 0.0;
}


CVanDerWaalsGas::CVanDerWaalsGas(double gamma, double R, double Pstar, double Tstar): CIdealGas(gamma, R) {
    a = 27.0/64.0*Gas_Constant*Gas_Constant*Tstar*Tstar/Pstar;
	b = 1.0/8.0*Gas_Constant*Tstar/Pstar;
//	a = 0.0;
//    b = 0.0;

}


CVanDerWaalsGas::~CVanDerWaalsGas(void) { }


void CVanDerWaalsGas::SetTDState_rhoe (double rho, double e ) {
	Density = rho;
	StaticEnergy = e;

	Pressure = Gamma_Minus_One*Density/(1.0-Density*b)*(StaticEnergy + Density*a) - a*Density*Density;
    Temperature = (Pressure+Density*Density*a)*((1-Density*b)/(Density*Gas_Constant));
//	Temperature = Gamma_Minus_One/Gas_Constant*(StaticEnergy + Density*a);

    dPde_rho = Density*Gamma_Minus_One/(1.0 - Density*b);
    dPdrho_e = Gamma_Minus_One/(1.0 - Density*b)*((StaticEnergy + 2*Density*a) + Density*b*(StaticEnergy + Density*a)/(1.0 - Density*b)) - 2*Density*a;
    dTdrho_e = Gamma_Minus_One/Gas_Constant*a;
    dTde_rho = Gamma_Minus_One/Gas_Constant;

    SoundSpeed2 = dPdrho_e + Pressure/(Density*Density)*dPde_rho;

}

void CVanDerWaalsGas::SetTDState_PT (double P, double T ) {
//	double a1 = -(P*b+T)/P;
//	double a2 =  a/P;
//	double a3 = -a*b/P;
//	double* root  = new double [3];
//	poly_solve_cubic(a1,a2,a3,root);
//
//    // check the root of the polynomial
//	double rho = 1.0/root[2];
//    double ad = a*rho;
//    double e = T/Gamma_Minus_One -ad;
//
//	SetTDState_de(rho, e);

}

void CVanDerWaalsGas::SetTDState_Prho (double P, double rho ) {

	SetEnergy_Prho(P, rho);

	SetTDState_rhoe(rho, StaticEnergy);

}

void CVanDerWaalsGas::SetTDState_hs (double h, double s ){

/// WARNING TODO

}

//void CVanDerWaalsGas::SetTDState_Ps (double P, double s ) {
//	// an implicit equation must be solved for rho
//	// The equation is in the form:
//	// (P+a*rho^2)*(od-b)^Gamma = s, where od = 1/rho
//
//    double rho = 0.0;
//
//	SetTDState_Pd(P, rho);
//
//}

void CVanDerWaalsGas::SetEnergy_Prho (double P, double rho ) {
	double T = (P+rho*rho*a)*(1-rho*b)/(rho*Gas_Constant);
	StaticEnergy = T*Gas_Constant/Gamma_Minus_One - rho*a;

}

//void CVanDerWaalsGas::SetEntropy_Prho (double P, double rho ) {
//	double od = 1/rho;
//	double odmb = od -b;
//	Entropy = ( P +a/pow(od,2) ) *pow(odmb,Gamma);
//
//}


// function for solving a cubic equation

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







