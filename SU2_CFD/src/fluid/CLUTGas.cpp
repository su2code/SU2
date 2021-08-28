/*!
 * \file CLUTGas.cpp
 * \brief Source of the interpolation table fluid model.
 * \author F. Dittmann
 * \version 7.0.7 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/fluid/CLUTGas.hpp"

CLUTGas::CLUTGas(bool CompEntropy) : CFluidModel() {
  ComputeEntropy = CompEntropy;
}

//saturation energy 
su2double CLUTGas::esat_rho(su2double rho) const
{   
    return coefesat[0] + coefesat[1]*rho + coefesat[2]*pow(rho,ONE2) + coefesat[3]*pow(rho,ONE3);
}

//interpolation
su2double CLUTGas::interpolateTable(su2double xi, su2double yi, const su2double x[2], const su2double y[2], const su2double z[Nx*Ny]) const
{   
    //fractional indices assuming equal spacing
    su2double ix = (xi-x[0])/(x[1]-x[0]) *(Nx-1); su2double iy = (yi-y[0])/(y[1]-y[0]) *(Ny-1);
    //bounded left and right indices (out of bound indices will result in extrapolation)
    int ixl = min(max(0,(int) SU2_TYPE::GetValue(ix)), Nx-2); int ixr = ixl+1;
    int iyl = min(max(0,(int) SU2_TYPE::GetValue(iy)), Ny-2); int iyr = iyl+1;      
    //linear interpolation
    su2double zil = z[ixl*Ny+iyl] + (iy - (su2double) iyl) *(z[ixl*Ny+iyr] - z[ixl*Ny+iyl]);
    su2double zir = z[ixr*Ny+iyl] + (iy - (su2double) iyl) *(z[ixr*Ny+iyr] - z[ixr*Ny+iyl]);
    return zil + (ix - (su2double) ixl) *(zir - zil);
}

//root finding
su2double CLUTGas::rootFunc(su2double x0, function<su2double(su2double)>func) const
{
    //initial guess
    int n = 0;
    su2double x = x0;        
    su2double y = func(x);
    x0  = 1.01*x0;
    su2double y0, dy, dx;
    //secant method
    while (abs(y)>1e-9*x0 && n<20) 
    {
        y0 = func(x0);
	    dy = y - y0;
	    dx = x - x0;
	    y = y0;
   	 	x = x0; 
	    x0 = x0 - y0*dx/dy;
        n++;
    }
    return x;
}

//---state variables from interpolation of rho-e table--- 
su2double CLUTGas::z_rhoe(su2double rhoi, su2double ei, const su2double z_TABLE[Nx*Ny]) const
{
    su2double Dei = ei - esat_rho(rhoi);
    return interpolateTable(rhoi,Dei,Rho,De,z_TABLE);
}

//---state variables from root finding of interpolation of rho-e table-- 
su2double CLUTGas::e_rhoz(su2double rhoi, su2double zi, const su2double z_TABLE[Nx*Ny]) const
{
    su2double ei0 = esat_rho(rhoi) + De[0];
    return rootFunc(ei0,[this,rhoi,zi,z_TABLE](su2double x){return z_rhoe(rhoi,x,z_TABLE)-zi;});
}

//---state variables from root finding of root finding of interpolation of rho-e table--- 
su2double CLUTGas::rho_zz(su2double z1i, su2double z2i, const su2double z1_TABLE[Nx*Ny], const su2double z2_TABLE[Nx*Ny]) const
{
    su2double rhoi0 = Rho[0];
    return rootFunc(rhoi0,[this,z1i,z2i,z1_TABLE,z2_TABLE](su2double x){return z_rhoe(x,e_rhoz(x,z1i,z1_TABLE),z2_TABLE)-z2i;});
}


//---cheap set state call---

void CLUTGas::SetTDState_rhoe(su2double rho, su2double e)
{
    Density = rho;
    StaticEnergy = e;
    Pressure = z_rhoe(rho,e,P_TABLE);
    Temperature = z_rhoe(rho,e,T_TABLE);
    //StaticEnthalpy = h_rhoe(rho,e);

    SoundSpeed2 = z_rhoe(rho,e,a2_TABLE);
    dPdrho_e = z_rhoe(rho,e,dPdrho_e_TABLE);
    dPde_rho = z_rhoe(rho,e,dPde_rho_TABLE);
    dTdrho_e = z_rhoe(rho,e,dTdrho_e_TABLE);
    dTde_rho = z_rhoe(rho,e,dTde_rho_TABLE);

    Cv = z_rhoe(rho,e,cv_TABLE);
    Cp = z_rhoe(rho,e,cp_TABLE);

    if (ComputeEntropy) Entropy = z_rhoe(rho,e,s_TABLE);
}


//---not so cheap set state calls---

void CLUTGas::SetEnergy_Prho(su2double P, su2double rho)
{
    StaticEnergy = e_rhoz(rho,P,P_TABLE);
}

void CLUTGas::SetTDState_Prho(su2double P, su2double rho)
{
    SetTDState_rhoe(rho,e_rhoz(rho,P,P_TABLE));
}

void CLUTGas::SetTDState_rhoT(su2double rho, su2double T)
{   
    SetTDState_rhoe(rho,e_rhoz(rho,T,T_TABLE));
}

void CLUTGas::SetTDState_rhoh(su2double rho, su2double h)
{   
    SetTDState_rhoe(rho,e_rhoz(rho,h,h_TABLE));
}


//---expensive set state calls---

void CLUTGas::SetTDState_PT(su2double P, su2double T)
{
    SetTDState_Prho(P,rho_zz(P,T,P_TABLE,T_TABLE));   
}

void CLUTGas::SetTDState_Ps(su2double P, su2double s)
{
    SetTDState_Prho(P,rho_zz(P,s,P_TABLE,s_TABLE));    
}

void CLUTGas::SetTDState_hs(su2double h, su2double s)
{
    SetTDState_rhoh(rho_zz(h,s,h_TABLE,s_TABLE),h);    
}
