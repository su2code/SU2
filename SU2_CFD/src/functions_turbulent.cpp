#include "../include/functions_turbulent.hpp"

SpalartAllmarasConstants::SpalartAllmarasConstants(){
  /*--- Spalart-Allmaras closure constants ---*/
  cv1_3 = pow(7.1,3.0);
  k2 = pow(0.41,2.0);
  cb1 = 0.1355;
  cw2 = 0.3;
  cw3_6 = pow(2.0,6.0);
  sigma = 2./3.;
  cb2 = 0.622;
  cb2_sigma = cb2/sigma;
  cw1 = cb1/k2+(1+cb2)/sigma;
}

/* Computes the spalart-allmaras source term.
 the outputs are 
 (Production, Destruction, CrossProduction) */
void SpalartAllmarasSourceTerm(void){
  /*--- Computation of vorticity ---*/
  Vorticity = (PrimVar_Grad_i[2][0]-PrimVar_Grad_i[1][1])*(PrimVar_Grad_i[2][0]-PrimVar_Grad_i[1][1]);
  if (nDim == 3) Vorticity += ( (PrimVar_Grad_i[3][1]-PrimVar_Grad_i[2][2])*(PrimVar_Grad_i[3][1]-PrimVar_Grad_i[2][2]) + (PrimVar_Grad_i[1][2]-PrimVar_Grad_i[3][0])*(PrimVar_Grad_i[1][2]-PrimVar_Grad_i[3][0]) );
  Omega = sqrt(Vorticity);
  
  /*--- Rotational correction term ---*/
  if (rotating_frame) {
    div = PrimVar_Grad_i[1][0] + PrimVar_Grad_i[2][1];
    if (nDim == 3) div += PrimVar_Grad_i[3][2];
    StrainMag = 0.0;
    // add diagonals
    StrainMag += pow(PrimVar_Grad_i[1][0] - 1.0/3.0*div,2.0);
    StrainMag += pow(PrimVar_Grad_i[2][1] - 1.0/3.0*div,2.0);
    if (nDim == 3) StrainMag += pow(PrimVar_Grad_i[3][2] - 1.0/3.0*div,2.0);
    // add off diagonals
    StrainMag += 2.0*pow(0.5*(PrimVar_Grad_i[1][1]+PrimVar_Grad_i[2][0]),2.0);
    if (nDim == 3) {
      StrainMag += 2.0*pow(0.5*(PrimVar_Grad_i[1][2]+PrimVar_Grad_i[3][0]),2.0);
      StrainMag += 2.0*pow(0.5*(PrimVar_Grad_i[2][2]+PrimVar_Grad_i[3][1]),2.0);
    }
    StrainMag = sqrt(2.0*StrainMag);
    Omega += 2.0*min(0.0,StrainMag-Omega);
  }
  
  if (dist_i > 1e-10) {
    
    /*--- Production term ---*/
    dist_i_2 = dist_i*dist_i;
    nu = Laminar_Viscosity_i/Density_i;
    Ji = TurbVar_i[0]/nu;
    Ji_2 = Ji*Ji;
    Ji_3 = Ji_2*Ji;
    fv1 = Ji_3/(Ji_3+cv1_3);
    fv2 = 1.0 - Ji/(1.0+Ji*fv1);
    S = Omega;
    inv_k2_d2 = 1.0/(k2*dist_i_2);
    
    Shat = S + TurbVar_i[0]*fv2*inv_k2_d2;
    inv_Shat = 1.0/max(Shat, 1.0e-10);
    
    /*--- Production term ---*/
    
    if (!transition) Production = cb1*Shat*TurbVar_i[0]*Volume;
    else Production = cb1*Shat*TurbVar_i[0]*Volume*intermittency;
    
    /*--- Destruction term ---*/
    
    r = min(TurbVar_i[0]*inv_Shat*inv_k2_d2,10.0);
    g = r + cw2*(pow(r,6.0)-r);
    g_6 =	pow(g,6.0);
    glim = pow((1.0+cw3_6)/(g_6+cw3_6),1.0/6.0);
    fw = g*glim;
    
    if (!transition) Destruction = cw1*fw*TurbVar_i[0]*TurbVar_i[0]/dist_i_2*Volume;
    else Destruction = cw1*fw*TurbVar_i[0]*TurbVar_i[0]/dist_i_2*Volume*min(max(intermittency,0.1),1.0);
    
    /*--- Diffusion term ---*/
    
    norm2_Grad = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      norm2_Grad += TurbVar_Grad_i[0][iDim]*TurbVar_Grad_i[0][iDim];
    CrossProduction = cb2_sigma*norm2_Grad*Volume;
    
    val_residual[0] = Production - Destruction + CrossProduction;
    
    /*--- Implicit part ---*/
    
    /*--- Production term ---*/
    dfv1 = 3.0*Ji_2*cv1_3/(nu*pow(Ji_3+cv1_3,2.));
    dfv2 = -(1/nu-Ji_2*dfv1)/pow(1.+Ji*fv1,2.);
    if ( Shat <= 1.0e-10 ) dShat = 0.0;
    else dShat = (fv2+TurbVar_i[0]*dfv2)*inv_k2_d2;
    val_Jacobian_i[0][0] += cb1*(TurbVar_i[0]*dShat+Shat)*Volume;
    
    /*--- Destruction term ---*/
    dr = (Shat-TurbVar_i[0]*dShat)*inv_Shat*inv_Shat*inv_k2_d2;
    if (r == 10.0) dr = 0.0;
    dg = dr*(1.+cw2*(6.*pow(r,5.)-1.));
    dfw = dg*glim*(1.-g_6/(g_6+cw3_6));
    val_Jacobian_i[0][0] -= cw1*(dfw*TurbVar_i[0] +	2.*fw)*TurbVar_i[0]/dist_i_2*Volume;
};