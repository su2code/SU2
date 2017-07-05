/*!
 * fluid_model_ppr.cpp
 * \brief Source of the Peng-Robinson model.
 * \author S. Vitale, G. Gori, M. Pini, A. Guardone, P. Colonna
 * \version 5.0.0 "Raven"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2017 SU2, the open-source CFD code.
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

CPengRobinson::CPengRobinson(su2double gamma, su2double R, su2double Pstar, su2double Tstar, su2double w) : CIdealGas(gamma, R) {

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


su2double CPengRobinson::alpha2(su2double T) {

  return ( 1 + k*(1 - sqrt(T/TstarCrit)))*( 1 + k*(1 - sqrt(T/TstarCrit)));
}

su2double CPengRobinson::T_v_h(su2double v, su2double h) {
  su2double fv, A, B, C, T, d, atanh;
  su2double sqrt2=sqrt(2.0);

  d = (v*v+2*b*v-b*b);
  
  atanh = (log(1.0+( b*sqrt2 / (v + b))) - log(1.0-( b*sqrt2 / (v + b))))/2.0;
  
  fv = atanh;

  A = Gas_Constant*(1 / Gamma_Minus_One + v/(v-b)) - a*v*k*k / (TstarCrit * d);
  B = a*k*(k+1)/sqrt(TstarCrit) *( fv/(b*sqrt2) + 2*v/d );
  C = h + a*(1+k)*(1+k)*(fv/(b*sqrt2) + v/d);

  T = ( -B + sqrt(B*B + 4*A*C) ) / (2*A); /// Only positive root considered

  return T*T;
}

su2double CPengRobinson::T_P_rho(su2double P, su2double rho) {
  su2double A, B, C, T, vb1, vb2;
  vb1 = (1/rho -b);
  vb2 = (1/rho/rho + 2*b/rho - b*b);

  A =   Gas_Constant/vb1 - a*k*k/TstarCrit/vb2;

  B =   2*a*k*(k+1)/sqrt(TstarCrit)/vb2;

  C = - P - a*(1+k)*(1+k)/vb2;

  T = ( -B + sqrt(B*B - 4*A*C) ) / (2*A);
  T *= T;

  return T;
}

void CPengRobinson::SetTDState_rhoe (su2double rho, su2double e ) {

    su2double DpDd_T, DpDT_d, DeDd_T, Cv;
    su2double A, B, C, sqrt2, fv, a2T, rho2, atanh;

    Density = rho;
    StaticEnergy = e;

    rho2 = rho*rho;
    sqrt2=sqrt(2.0);

    atanh = (log(1.0+( rho * b * sqrt2/(1 + rho*b))) - log(1.0-( rho * b * sqrt2/(1 + rho*b))))/2.0;
  
    fv = atanh;
    
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

void CPengRobinson::SetTDState_PT (su2double P, su2double T ) {
  su2double toll= 1e-6;
  su2double A, B, Z, DZ=1.0, F, F1, atanh;
  su2double rho, fv, e;
  su2double sqrt2=sqrt(2.0);
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
  
  atanh = (log(1.0+( rho * b * sqrt2/(1 + rho*b))) - log(1.0-( rho * b * sqrt2/(1 + rho*b))))/2.0;
  
  fv = atanh;

  e = T*Gas_Constant/Gamma_Minus_One - a*(k+1)*sqrt( alpha2(T) )*fv / (b*sqrt2);

  SetTDState_rhoe(rho, e);
}

void CPengRobinson::SetTDState_Prho (su2double P, su2double rho ) {

  SetEnergy_Prho(P, rho);

  SetTDState_rhoe(rho, StaticEnergy);

}

void CPengRobinson::SetTDState_hs (su2double h, su2double s ) {

  su2double T, fv, sqrt2=sqrt(2.0), A;
  su2double f, v, atanh;
  su2double x1, x2, xmid, dx, fx1, fx2, fmid, rtb;
  su2double toll = 1e-9, FACTOR=0.2;
  su2double cons_s, cons_h;
  unsigned short countrtb=0, NTRY=10, ITMAX=100;

  A = Gas_Constant / Gamma_Minus_One;
  T = h*Gamma_Minus_One/Gas_Constant/Gamma;
  v = exp(-1/Gamma_Minus_One*log(T) + s/Gas_Constant);


  if (Zed<0.9999) {
    x1 = Zed*v;
    x2 = v;

  } else {
    x1 = 0.2*v;
    x2 = v;
  }


  T = T_v_h(x1, h);
  
  atanh = (log(1.0+( b*sqrt2 / (x1 + b))) - log(1.0-( b*sqrt2 / (x1 + b))))/2.0;
  fv = atanh;
  
  fx1 = A*log(T) + Gas_Constant*log(x1 - b) - a*sqrt(alpha2(T)) *k*fv/(b*sqrt2*sqrt(T*TstarCrit)) - s;
  T = T_v_h(x2, h);
  
  atanh = (log(1.0+( b*sqrt2 / (x2 + b))) - log(1.0-( b*sqrt2 / (x2 + b))))/2.0;
  fv = atanh;
  
  fx2 = A*log(T) + Gas_Constant*log(x2 - b) - a*sqrt(alpha2(T)) *k*fv/(b*sqrt2*sqrt(T*TstarCrit)) - s;

  // zbrac algorithm NR

  for (int j=1; j<=NTRY; j++) {
    if (fx1*fx2 > 0.0) {
      if (fabs(fx1) < fabs(fx2)) {
        x1 += FACTOR*(x1-x2);
        T = T_v_h(x1, h);
        atanh = (log(1.0+( b*sqrt2/(x1 + b))) - log(1.0-( b*sqrt2/(x1 + b))))/2.0;
        fv = atanh;
        fx1 = A*log(T) + Gas_Constant*log(x1 - b) - a*sqrt(alpha2(T)) *k*fv/(b*sqrt2*sqrt(T*TstarCrit)) - s;
      } else {
        x2 += FACTOR*(x2-x1);
        T = T_v_h(x2, h);
        atanh = (log(1.0+( b*sqrt2/(x2 + b))) - log(1.0-( b*sqrt2/(x2 + b))))/2.0;
        fv = atanh;
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
    atanh = (log(1.0+( b* sqrt2/(xmid + b))) - log(1.0-( b* sqrt2/(xmid + b))))/2.0;
    fv = atanh;
    fmid= A*log(T) + Gas_Constant*log(xmid - b) - a*sqrt(alpha2(T)) *k*fv/(b*sqrt2*sqrt(T*TstarCrit)) - s;

    if (fmid <= 0.0) rtb=xmid;
    countrtb++;
  }while(abs(fmid) > toll && countrtb<ITMAX);

  v = xmid;
  if (countrtb==ITMAX) {
    cout <<"Too many bisections in rtbis" << endl;
//      do{
//          atanh = (log(1.0+( b/v* sqrt2/(1 + b/v))) - log(1.0-( b/v* sqrt2/(1 + b/v))))/2.0;
//          fv = atanh;
//          T=T_v_h(v, h);
//          f = A*log(T) + Gas_Constant*log(v - b) - a*sqrt(alpha2(T)) *k*fv/(b*sqrt2*sqrt(T*TstarCrit)) - s;
//          f1= Gas_Constant/(v-b)+ a*sqrt(alpha2(T)) *k/(sqrt(T*TstarCrit)*(v*v - b*b - 2*v*b));
//          dv= f/f1;
//          v-= dv;
//          countnw++;
//      }while(abs(f/x2) > toll && countnw<ITMAXNW);
//
//    } else {
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

  if (cons_h >1e-4 || cons_s >1e-4) {
    cout<< "TD consistency not verified in hs call"<< endl;
       //cout <<"Before  "<< h <<" "<< s << endl;
       //cout <<"After  "<< StaticEnergy + Pressure/Density <<" "<< Entropy << fmid <<" "<< f<< " "<< countrtb<<" "<< countnw<< endl;
       //getchar();
  }
}


void CPengRobinson::SetEnergy_Prho (su2double P, su2double rho) {

    su2double ad;
    su2double A, B, C, T, vb1, vb2, atanh;
    vb1 = (1/rho -b);
    vb2 = (1/rho/rho + 2*b/rho - b*b);

    A =   Gas_Constant/vb1 - a*k*k/TstarCrit/vb2;

    B =   2*a*k*(k+1)/sqrt(TstarCrit)/vb2;

    C = - P - a*(1+k)*(1+k)/vb2;

    T = ( -B + sqrt(B*B - 4*A*C) ) / (2*A);
    T *= T;

    atanh = (log(1.0+( rho * b * sqrt(2.0)/(1 + rho*b) )) - log(1.0-( rho * b * sqrt(2.0)/(1 + rho*b) )))/2.0;
    ad = a*(k+1)*sqrt( alpha2(T) ) / ( b*sqrt(2.0) ) * atanh ;

    StaticEnergy = T * Gas_Constant / Gamma_Minus_One - ad;

}

void CPengRobinson::SetTDState_rhoT (su2double rho, su2double T) {
  su2double fv, e, atanh;

  atanh = (log(1.0+( rho * b * sqrt(2.0)/(1 + rho*b))) - log(1.0-( rho * b * sqrt(2.0)/(1 + rho*b))))/2.0;
  fv = atanh;
  e = T*Gas_Constant/Gamma_Minus_One - a*(k+1)*sqrt( alpha2(T) ) / ( b*sqrt(2.0) ) * fv;
  SetTDState_rhoe(rho, e);
}

void CPengRobinson::SetTDState_Ps (su2double P, su2double s) {

  su2double T, rho, v, cons_P, cons_s, fv, A, atanh;
  su2double x1,x2, fx1, fx2,f, fmid, rtb, dx, xmid, sqrt2=sqrt(2.0);
  su2double toll = 1e-5, FACTOR=0.2;
  unsigned short count=0, NTRY=10, ITMAX=100;

  A = Gas_Constant / Gamma_Minus_One;
  T   = exp(Gamma_Minus_One/Gamma* (s/Gas_Constant +log(P) -log(Gas_Constant)) );
  v = (T*Gas_Constant)/P;

  if(Zed<0.9999) {
    x1 = Zed*v;
    x2 = v;

  }else {
    x1 = 0.2*v;
    x2 = v;
  }
  T = T_P_rho(P,1.0/x1);
  
  atanh = (log(1.0 + ( b*sqrt2 / (x1 + b) )) - log(1.0-( b*sqrt2 / (x1 + b) )))/2.0;
  fv = atanh;
  
  fx1 = A*log(T) + Gas_Constant*log(x1 - b) - a*sqrt(alpha2(T)) *k*fv/(b*sqrt2*sqrt(T*TstarCrit)) - s;
  T = T_P_rho(P,1.0/x2);
  
  atanh = (log(1.0 + ( b*sqrt2 / (x2 + b) )) - log(1.0-( b*sqrt2 / (x2 + b) )))/2.0;
  fv = atanh;

  fx2 = A*log(T) + Gas_Constant*log(x2 - b) - a*sqrt(alpha2(T)) *k*fv/(b*sqrt2*sqrt(T*TstarCrit)) - s;

  // zbrac algorithm NR

  for (int j=1;j<=NTRY;j++) {
    if (fx1*fx2 > 0.0) {
      if (fabs(fx1) < fabs(fx2)) {
        x1 += FACTOR*(x1-x2);
        T = T_P_rho(P,1.0/x1);
        
        atanh = (log(1.0 + ( b*sqrt2 / (x1 + b) )) - log(1.0-( b*sqrt2 / (x1 + b) )))/2.0;
        fv = atanh;
        
        fx1 = A*log(T) + Gas_Constant*log(x1 - b) - a*sqrt(alpha2(T)) *k*fv/(b*sqrt2*sqrt(T*TstarCrit)) - s;
      }else {
        T = T_P_rho(P,1.0/x2);
        
        atanh = (log(1.0 + ( b*sqrt2 / (x2 + b) )) - log(1.0-( b*sqrt2 / (x2 + b) )))/2.0;
        fv = atanh;

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
    T = T_P_rho(P,1.0/xmid);
    
    atanh = (log(1.0 + ( b*sqrt2 / (xmid + b) )) - log(1.0-( b*sqrt2 / (xmid + b) )))/2.0;
    fv = atanh;
    
    fmid = A*log(T) + Gas_Constant*log(xmid - b) - a*sqrt(alpha2(T)) *k*fv/(b*sqrt2*sqrt(T*TstarCrit)) - s;
    if (fmid <= 0.0) rtb=xmid;
    count++;
    }while(abs(fmid) > toll && count<ITMAX);

    if(count==ITMAX) {
      cout <<"Too many bisections in rtbis" << endl;
    }

  rho = 1.0/xmid;
  T = T_P_rho(P, rho);
  SetTDState_rhoT(rho, T);
//  cout << xmid << " "<< T<< " "<< Pressure<< " "<< P << " "<< Entropy << " "<< s <<endl;

  cons_P= abs((Pressure -P)/P);
  cons_s= abs((Entropy-s)/s);

  if(cons_P >1e-3 || cons_s >1e-3) {
    cout<< "TD consistency not verified in hs call"<< endl;
  }

}


CPengRobinson_Generic::CPengRobinson_Generic() : CPengRobinson() {
  a= 0.0;
  b =0.0;
  k = 0.0;
  TstarCrit = 0.0;
}

CPengRobinson_Generic::CPengRobinson_Generic(su2double gamma, su2double R, su2double Pstar, su2double Tstar, su2double w) :
		CPengRobinson(gamma, R, Pstar, Tstar, w) {

  a = 0.45724*Gas_Constant*Gas_Constant*Tstar*Tstar/Pstar;
  b = 0.0778*Gas_Constant*Tstar/Pstar;
  TstarCrit = Tstar;
  Zed=1.0;

  if (w <= 0.49)
        k = 0.37464 + 1.54226 * w - 0.26992 * w*w;
        else
        k = 0.379642 + 1.48503 * w - 0.164423 * w*w + 0.016666 * w*w*w;


}

CPengRobinson_Generic::~CPengRobinson_Generic(void) { }

su2double CPengRobinson_Generic::alpha2(su2double T) {

	// alpha call corrected
  su2double alpha_2 = pow(1 + k - k*sqrt(T/TstarCrit), 2);

  return  ( alpha_2 );
}

su2double CPengRobinson_Generic::dalphadT(su2double T) {

	// alpha call corrected

  return ( -0.5*k / sqrt(T*TstarCrit));
}

su2double CPengRobinson_Generic::dalpha2dT2(su2double T) {

	// alpha call corrected

  return ( 0.25 * k * pow(T*TstarCrit, -1.5) * TstarCrit);
}


su2double CPengRobinson_Generic::T_v_h(su2double v, su2double h) {

	//corrected
  su2double fv, A, B, C, T, d, atanh, T_new;
  su2double toll = 1e-9, error = 1;
  unsigned int count_T = 0, ITMAX = 100;
  su2double sqrt2=sqrt(2.0);

  d = (v*v+2*b*v-b*b);

  HeatCapacity->Set_Cv0(TstarCrit);
  Cv0 = HeatCapacity->Get_Cv0();

  atanh = (log(1.0+( b*sqrt2 / (v + b))) - log(1.0-( b*sqrt2 / (v + b))))/2.0;
  fv = atanh;

  A = Cv0 + Gas_Constant*v/(v-b) - a*v*k*k / (TstarCrit * d);
  B = a*k*(k+1)/sqrt(TstarCrit) *( fv/(b*sqrt2) + 2*v/d );
  C = h + a*(1+k)*(1+k)*(fv/(b*sqrt2) + v/d);

  T_new = ( -B + sqrt(B*B + 4*A*C) ) / (2*A); /// Only positive root considered
  T_new *= T_new;

  do {
	  T = T_new;
	  HeatCapacity->Set_Cv0(T);
	  Cv0 = HeatCapacity->Get_Cv0();
	  Set_Cv(T, v);
//	  cout << T_new << " " << v << " " << Cv << endl;

      A = Cv + Gas_Constant*v/(v-b) - a*v*k*k / (TstarCrit * d);

      T_new = ( -B + sqrt(B*B + 4*A*C) ) / (2*A); /// Only positive root considered
      T_new = T_new * T_new;
      error = abs((T - T_new)/T);
      count_T++;

  } while (count_T < ITMAX && error > toll);

  if (count_T == ITMAX)
  	cout << "Too many iterations in T_v_h function" << endl;

  return T_new;

}


su2double CPengRobinson_Generic::T_P_rho(su2double P, su2double rho) {
  su2double A, B, C, T, vb1, vb2;

  vb1 = (1/rho -b);
  vb2 = (1/rho/rho + 2*b/rho - b*b);

  A =   Gas_Constant/vb1 - a*k*k/TstarCrit/vb2;

  B =   2*a*k*(k+1)/sqrt(TstarCrit)/vb2;

  C = - P - a*(1+k)*(1+k)/vb2;

  T = ( -B + sqrt(B*B - 4*A*C) ) / (2*A);
  T *= T;

  return T;
}

void CPengRobinson_Generic::SetTDState_rhoe (su2double rho, su2double e ) {

	// rhoe call corrected

    su2double DpDd_T, DpDT_d, DeDd_T, der, Temperature_new;
    su2double error = 1, toll = 1e-9;
    su2double A, B, C, sqrt2, fv, a2T, rho2, atanh;
    unsigned int count_T = 0, ITMAX = 100;

    Density = rho;
    StaticEnergy = e;

    rho2 = rho*rho;
    sqrt2=sqrt(2.0);
    atanh = (log(1.0+( rho * b * sqrt2/(1 + rho*b))) - log(1.0-( rho * b * sqrt2/(1 + rho*b))))/2.0;
	fv = atanh;

	HeatCapacity->Set_Cv0(TstarCrit);
	A = HeatCapacity->Get_Cv0();
	B = a*k*(k+1)*fv/(b*sqrt2*sqrt(TstarCrit));
	C = a*(k+1)*(k+1)*fv/(b*sqrt2) + e;

	Temperature_new = ( -B + sqrt(B*B + 4*A*C) ) / (2*A); /// Only positive root considered
	Temperature_new *= Temperature_new;

    do {
    	Temperature = Temperature_new;

    	HeatCapacity->Set_Cv0(Temperature);
    	Cv0 = HeatCapacity->Get_Cv0();
    	Set_Cv( Temperature, 1/rho);

		A = Cv;

		Temperature_new = ( -B + sqrt(B*B + 4*A*C) ) / (2*A); /// Only positive root considered
		Temperature_new *= Temperature_new;
		error = abs(Temperature - Temperature_new)/Temperature;
		count_T++;

    } while (count_T < ITMAX && error > toll);

    if (count_T == ITMAX) {
    	cout << "Too many iterations in rho_e call" << endl;
    	cout << "Error on Temperature: " << error << endl;
    }

    Temperature = Temperature_new;
	HeatCapacity->Set_Cv0(Temperature);
	Cv0 = HeatCapacity->Get_Cv0();
	Set_Cv( Temperature, 1/rho);




    a2T = alpha2(Temperature);

    A = (1/rho2 + 2*b/rho - b*b);
    B = 1/rho-b;

    Pressure = Temperature*Gas_Constant / B - a*a2T / A;

    SetGamma_Trho();


    if 	(Gamma > 6 ) {
    	cout << "Warning: Gamma value greater than 6, switch to CONSTANT_GAMMA" << endl;
     	cout << "Ideal gas correction implemented, Cp = Cv + R" << endl;
     	Cp = Cv + Gas_Constant;
    	Gamma = Cp/Cv;
    	Gamma_Minus_One = Gamma - 1;
//		getchar();
    }

    if 	( Gamma < 1) {
     	cout << "Warning: Gamma value lower than 1, switch to CONSTANT_GAMMA" << endl;
    	cout << "Ideal gas correction implemented, Cp = Cv + R" << endl;
    	Cp = Cv + Gas_Constant;
    	Gamma = Cp/Cv;
    	Gamma_Minus_One = Gamma - 1;

    }


    Entropy = Cv*log(Temperature) + Gas_Constant*log(B) - a*sqrt(a2T) *k*fv/(b*sqrt2*sqrt(Temperature*TstarCrit));

    DpDd_T =  ( Temperature*Gas_Constant /(B*B    ) - 2*a*a2T*(1/rho + b) /( A*A ) ) /(rho2);

    DpDT_d = Gas_Constant /B + a*k / A * sqrt( a2T/(Temperature*TstarCrit) );

    der = Cv + ( a *k*(k+1)*fv ) / ( 2*b*sqrt(2*Temperature*TstarCrit) );

    dPde_rho = DpDT_d/der;

    DeDd_T = - a*(1+k) * sqrt( a2T ) / A / (rho2);

    dPdrho_e = DpDd_T - dPde_rho*DeDd_T;

    SoundSpeed2 = dPdrho_e + Pressure/(rho2)*dPde_rho;

    dTde_rho = 1/der;

    Zed = Pressure/(Gas_Constant*Temperature*Density);

}

void CPengRobinson_Generic::SetTDState_PT (su2double P, su2double T ) {

// PT call corrected

  su2double toll= 1e-6;
  su2double A, B, Z, DZ=1.0, F, F1, atanh;
  su2double rho, fv, e;
  su2double sqrt2=sqrt(2.0);
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

  atanh = (log(1.0+( rho * b * sqrt2/(1 + rho*b))) - log(1.0-( rho * b * sqrt2/(1 + rho*b))))/2.0;

  fv = atanh;

  HeatCapacity->Set_Cv0(T);
  Cv0 = HeatCapacity->Get_Cv0();

  Set_Cv(T, 1/rho);

  e = T*Cv - a*(k+1)*sqrt( alpha2(T) )*fv / (b*sqrt2);

  SetTDState_rhoe(rho, e);
}

void CPengRobinson_Generic::SetTDState_Prho (su2double P, su2double rho ) {
	// Prho call corrected

  SetEnergy_Prho(P, rho);


  SetTDState_rhoe(rho, StaticEnergy);

}

void CPengRobinson_Generic::SetTDState_hs (su2double h, su2double s ) {

  su2double T, fv, fv_new, sqrt2=sqrt(2.0), A, rho, rho_new;
  su2double f, v, atanh, P;
  su2double x1, x2, xmid, dx, fx1, fx2, fmid, rtb;
  su2double toll = 1e-9, FACTOR=0.2;
  su2double cons_s, cons_h;

  su2double error=1, error_v = 1, T_new, v_new;
  unsigned short countrtb=0, NTRY=50, ITMAX=100, count_T=0, count_v = 0;

/*
  HeatCapacity->Set_Cv0(TstarCrit);
  Cp = HeatCapacity->Get_Cv0 ();
  T_new = abs(h)/Cp;


	  do{
		T = T_new;
		HeatCapacity->Set_Cv0 (T);
		Cv0 = HeatCapacity->Get_Cv0 ();
		rho_new = 1/(exp((s - Cv0*log(T) + a*sqrt(alpha2(T)) *k/(b*sqrt(2*T*TstarCrit)))/ Gas_Constant)+ b) ;

		do {
			rho = rho_new;
			Set_Cv(T,1/rho);
			fv = (log(1.0+( b*sqrt2 / (1/rho + b))) - log(1.0-( b*sqrt2 / (1/rho + b))))/2.0;

				rho_new = s - Cv*log(T) + a*sqrt(alpha2(T)) *k*fv/(b*sqrt(2*T*TstarCrit));
			    rho_new = exp(rho_new/Gas_Constant) + b;
			    rho_new = 1/rho_new;

			error_v = abs(rho_new - rho)/rho_new;
//			cout << "2iter " << rho_new << " " << Cv << endl;
			count_v++;

		}while(error_v >1e-5 && count_v < ITMAX);

		error_v = 1;
		count_v = 0;

		T_new = T_v_h(1/rho_new, h);
		error = abs(T - T_new)/T;
		count_T++;

	  }while(error >toll && count_T < ITMAX);

//  cout << T << " " << rho_new << endl;
//  getchar();
   SetTDState_rhoT(1/v_new, T_new);

   cons_h= abs((StaticEnergy + Pressure/Density -h)/h);
   cons_s= abs((Entropy-s)/s);



   if(cons_h >1e-3 || cons_s >1e-3) {
     cout<< "TD consistency not verified in hs call"<< endl;
   }

}*/

    HeatCapacity->Set_Cv0(TstarCrit);
    A = HeatCapacity->Get_Cv0 ();
    T = abs(h)/A;
    v = exp(-Cv*log(T)/Gas_Constant + s/Gas_Constant);


    if (Zed<0.9999) {
      x1 = Zed*v;
      x2 = v;

    } else {
      x1 = 0.2*v;
      x2 = v;
    }


    T = T_v_h(x1, h);

    atanh = (log(1.0+( b*sqrt2 / (x1 + b))) - log(1.0-( b*sqrt2 / (x1 + b))))/2.0;
    fv = atanh;

    fx1 = Cv*log(T) + Gas_Constant*log(x1 - b) - a*sqrt(alpha2(T)) *k*fv/(b*sqrt2*sqrt(T*TstarCrit)) - s;
    T = T_v_h(x2, h);

    atanh = (log(1.0+( b*sqrt2 / (x2 + b))) - log(1.0-( b*sqrt2 / (x2 + b))))/2.0;
    fv = atanh;

    fx2 = Cv*log(T) + Gas_Constant*log(x2 - b) - a*sqrt(alpha2(T)) *k*fv/(b*sqrt2*sqrt(T*TstarCrit)) - s;

    // zbrac algorithm NR

    for (int j=1; j<=NTRY; j++) {
      if (fx1*fx2 > 0.0) {
        if (fabs(fx1) < fabs(fx2)) {
          x1 += FACTOR*(x1-x2);
          T = T_v_h(x1, h);
          atanh = (log(1.0+( b*sqrt2/(x1 + b))) - log(1.0-( b*sqrt2/(x1 + b))))/2.0;
          fv = atanh;
          fx1 = Cv*log(T) + Gas_Constant*log(x1 - b) - a*sqrt(alpha2(T)) *k*fv/(b*sqrt2*sqrt(T*TstarCrit)) - s;
        } else {
          x2 += FACTOR*(x2-x1);
          T = T_v_h(x2, h);
          atanh = (log(1.0+( b*sqrt2/(x2 + b))) - log(1.0-( b*sqrt2/(x2 + b))))/2.0;
          fv = atanh;
          fx2 = Cv*log(T) + Gas_Constant*log(x2 - b) - a*sqrt(alpha2(T)) *k*fv/(b*sqrt2*sqrt(T*TstarCrit)) - s;
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
      atanh = (log(1.0+( b* sqrt2/(xmid + b))) - log(1.0-( b* sqrt2/(xmid + b))))/2.0;
      fv = atanh;
      fmid= Cv*log(T) + Gas_Constant*log(xmid - b) - a*sqrt(alpha2(T)) *k*fv/(b*sqrt2*sqrt(T*TstarCrit)) - s;

      if (fmid <= 0.0) rtb=xmid;
      countrtb++;
    }while(abs(fmid) > toll && countrtb<ITMAX);

    v = xmid;
    if (countrtb==ITMAX) {
      cout <<"Too many bisections in rtbis" << endl;
  //      do{
  //          atanh = (log(1.0+( b/v* sqrt2/(1 + b/v))) - log(1.0-( b/v* sqrt2/(1 + b/v))))/2.0;
  //          fv = atanh;
  //          T=T_v_h(v, h);
  //          f = A*log(T) + Gas_Constant*log(v - b) - a*sqrt(alpha2(T)) *k*fv/(b*sqrt2*sqrt(T*TstarCrit)) - s;
  //          f1= Gas_Constant/(v-b)+ a*sqrt(alpha2(T)) *k/(sqrt(T*TstarCrit)*(v*v - b*b - 2*v*b));
  //          dv= f/f1;
  //          v-= dv;
  //          countnw++;
  //      }while(abs(f/x2) > toll && countnw<ITMAXNW);
  //
  //    } else {
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

    if (cons_h >1e-4 || cons_s >1e-4) {
      cout<< "TD consistency not verified in hs call"<< endl;
         //cout <<"Before  "<< h <<" "<< s << endl;
         //cout <<"After  "<< StaticEnergy + Pressure/Density <<" "<< Entropy << fmid <<" "<< f<< " "<< countrtb<<" "<< countnw<< endl;
         //getchar();
    }
  }



void CPengRobinson_Generic::SetEnergy_Prho (su2double P, su2double rho) {

	// Energy _ Prho call corrected
    su2double ad;
    su2double A, B, C, T, vb1, vb2, atanh;
    vb1 = (1/rho -b);
    vb2 = (1/rho/rho + 2*b/rho - b*b);

    A =   Gas_Constant/vb1 - a*k*k/TstarCrit/vb2;

    B =   2*a*k*(k+1)/sqrt(TstarCrit)/vb2;

    C = - P - a*(1+k)*(1+k)/vb2;

    T = ( -B + sqrt(B*B - 4*A*C) ) / (2*A);
    T *= T;

    atanh = (log(1.0+( rho * b * sqrt(2.0)/(1 + rho*b) )) - log(1.0-( rho * b * sqrt(2.0)/(1 + rho*b) )))/2.0;
    ad = a*(k+1)*sqrt( alpha2(T) ) / ( b*sqrt(2.0) ) * atanh ;

    HeatCapacity->Set_Cv0(T);
    Cv0 = HeatCapacity->Get_Cv0();

    Set_Cv(T, 1/rho);

    StaticEnergy = T * Cv - ad;


}

void CPengRobinson_Generic::SetTDState_rhoT (su2double rho, su2double T) {

	// rhoT call corrected
  su2double fv, e, atanh, a2T;

  a2T = alpha2(T);

  atanh = (log(1.0+( rho * b * sqrt(2.0)/(1 + rho*b))) - log(1.0-( rho * b * sqrt(2.0)/(1 + rho*b))))/2.0;
  fv = atanh;

  HeatCapacity->Set_Cv0(T);
  Cv0 = HeatCapacity->Get_Cv0();
  Set_Cv(T, 1/rho);

  e = T * Cv - a*(k+1)*sqrt( a2T ) / ( b*sqrt(2.0) ) * fv;

  SetTDState_rhoe(rho, e);


}

void CPengRobinson_Generic::SetTDState_Ps (su2double P, su2double s) {

  su2double T, rho, v, cons_P, cons_s, fv, A, atanh, T_new, v_new;
  su2double error = 0, error_v = 0;
  su2double x1,x2, fx1, fx2,f, fmid, rtb, dx, xmid, sqrt2=sqrt(2.0);

  su2double toll = 1e-9, FACTOR=0.2;
  unsigned short count=0, NTRY=10, ITMAX=100, count_T = 0, count_v = 0;

  T_new = TstarCrit;// exp(Gamma_Minus_One/Gamma* (s/Gas_Constant +log(P) -log(Gas_Constant)) );
//getchar();

	  do{
		T = T_new;
		HeatCapacity->Set_Cv0 (T);
		Cv0 = HeatCapacity->Get_Cv0 ();
		v_new = exp((s - Cv0*log(T)) / Gas_Constant)+ b ;

		do {
			v = v_new;
			Set_Cv(T,v);

			fv = (log(1.0+( b*sqrt2 / (v + b))) - log(1.0-( b*sqrt2 / (v + b))))/2.0;
			v_new = (s - Cv*log(T) + a*sqrt(alpha2(T)) *k*fv/(b*sqrt2*sqrt(T*TstarCrit)))/Gas_Constant;
			v_new = exp(v_new) + b;
			error_v = abs(v - v_new)/v;
			count_v++;

		}while(error_v >toll && count_v < ITMAX);

		error_v = 1;
		count_v = 0;

		T_new = T_P_rho(P, 1/v);


		error = abs(T - T_new)/T;
		count_T++;
	  }while(error >toll && count_T < ITMAX);

   SetTDState_rhoT(1/v, T);


  cons_P= abs((Pressure -P)/P);
  cons_s= abs((Entropy-s)/s);

  if(cons_P >1e-3 || cons_s >1e-3) {
    cout<< "TD consistency not verified in hs call"<< endl;
  }

}

void CPengRobinson_Generic::SetGamma_Trho () {

// Gamma call corrected

  su2double dPodT, dPodv, CpmCv, daT, a2T;

  daT = dalphadT(Temperature);
  a2T = alpha2(Temperature);

  dPodT = 2 * a*sqrt(a2T) * daT;
  dPodT = dPodT / (1/Density/Density + 2*b/Density - b*b);
  dPodT = Gas_Constant/ (1/Density - b) - dPodT;

  dPodv = -b*b + 2* b / Density + 1/Density/Density;
  dPodv = + 2* a * a2T * (1/ Density + b) / pow(dPodv, 2);
  dPodv = dPodv -Gas_Constant * Temperature / (1/Density - b)/ (1/Density - b);

  CpmCv = -Temperature * pow(dPodT, 2)/dPodv;

  Cp = Cv + CpmCv;

  Gamma = Cp/Cv;

}


void CPengRobinson_Generic::Set_Cv (su2double T, su2double v) {
	// cv call corrected
  su2double CvmCv0, a2T, daT, da2T;

  a2T = alpha2(T);
  da2T = dalpha2dT2(T);
  daT = dalphadT(T);

  CvmCv0 = sqrt(a2T) * da2T;
  CvmCv0 = CvmCv0 + pow( daT, 2);
  CvmCv0 = -2 * a * CvmCv0 * T;

  CvmCv0 = CvmCv0 * 0.5 * (log( v + b - sqrt(2)*b) - log( v + b + sqrt(2)*b)) / sqrt(2) / b;

  Cv = Cv0 + CvmCv0;


}


