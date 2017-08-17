void GetPMatrix_R(su2double *val_density, su2double *val_velocity,
	su2double *val_soundspeed, su2double *val_normal, 
	su2double *val_p_tensor, su2double sqvel, su2double Gamma_Minus_One, 
	int nDim, int nVar) {
  
  su2double rhooc, rhoxc;
  su2double DN0, DN1, DN2;
  
  rhooc = *val_density / *val_soundspeed;
  rhoxc = (*val_density * *val_soundspeed) / Gamma_Minus_One;
  
  if (nDim == 2) {
    
    //sqvel = val_velocity[0]*val_velocity[0]+val_velocity[1]*val_velocity[1];
    DN0 = *val_density*val_normal[0];
    DN1 = *val_density*val_normal[1];
    
    val_p_tensor[0] = 1.0;
    val_p_tensor[1]=0.0;
    val_p_tensor[2]=0.5*rhooc;
    val_p_tensor[3]=0.5*rhooc;
    
    val_p_tensor[nVar+0]=val_velocity[0];
    val_p_tensor[nVar+1]=DN1;
    val_p_tensor[nVar+2]=0.5*(val_velocity[0]*rhooc+DN0);
    val_p_tensor[nVar+3]=0.5*(val_velocity[0]*rhooc-DN0);
    
    val_p_tensor[2*nVar+0]=val_velocity[1];
    val_p_tensor[2*nVar+1]=-DN0;
    val_p_tensor[2*nVar+2]=0.5*(val_velocity[1]*rhooc+DN1);
    val_p_tensor[2*nVar+3]=0.5*(val_velocity[1]*rhooc-DN1);
    
    val_p_tensor[3*nVar+0]=0.5*sqvel;
    val_p_tensor[3*nVar+1]=val_velocity[0]*DN1-val_velocity[1]*DN0;
    val_p_tensor[3*nVar+2]=0.5*(0.5*sqvel*rhooc+val_velocity[0]*DN0+val_velocity[1]*DN1+rhoxc);
    val_p_tensor[3*nVar+3]=0.5*(0.5*sqvel*rhooc-val_velocity[0]*DN0-val_velocity[1]*DN1+rhoxc);
    
  }
  else {
    
    //sqvel = val_velocity[0]*val_velocity[0]+val_velocity[1]*val_velocity[1]+val_velocity[2]*val_velocity[2];
    DN0 = *val_density*val_normal[0];
    DN1 = *val_density*val_normal[1];
    DN2 = *val_density*val_normal[2];
    
    val_p_tensor[0]=val_normal[0];
    val_p_tensor[1]=val_normal[1];
    val_p_tensor[2]=val_normal[2];
    val_p_tensor[3]=0.5*rhooc;
    val_p_tensor[4]=0.5*rhooc;
    
    val_p_tensor[nVar+0]=val_velocity[0]*val_normal[0];
    val_p_tensor[nVar+1]=val_velocity[0]*val_normal[1]-DN2;
    val_p_tensor[nVar+2]=val_velocity[0]*val_normal[2]+DN1;
    val_p_tensor[nVar+3]=0.5*(val_velocity[0]*rhooc+DN0);
    val_p_tensor[nVar+4]=0.5*(val_velocity[0]*rhooc-DN0);
    
    val_p_tensor[2*nVar+0]=val_velocity[1]*val_normal[0]+DN2;
    val_p_tensor[2*nVar+1]=val_velocity[1]*val_normal[1];
    val_p_tensor[2*nVar+2]=val_velocity[1]*val_normal[2]-DN0;
    val_p_tensor[2*nVar+3]=0.5*(val_velocity[1]*rhooc+DN1);
    val_p_tensor[2*nVar+4]=0.5*(val_velocity[1]*rhooc-DN1);
    
    val_p_tensor[3*nVar+0]=val_velocity[2]*val_normal[0]-DN1;
    val_p_tensor[3*nVar+1]=val_velocity[2]*val_normal[1]+DN0;
    val_p_tensor[3*nVar+2]=val_velocity[2]*val_normal[2];
    val_p_tensor[3*nVar+3]=0.5*(val_velocity[2]*rhooc+DN2);
    val_p_tensor[3*nVar+4]=0.5*(val_velocity[2]*rhooc-DN2);
    
    val_p_tensor[4*nVar+0]=0.5*sqvel*val_normal[0]+val_velocity[1]*DN2-val_velocity[2]*DN1;
    val_p_tensor[4*nVar+1]=0.5*sqvel*val_normal[1]-val_velocity[0]*DN2+val_velocity[2]*DN0;
    val_p_tensor[4*nVar+2]=0.5*sqvel*val_normal[2]+val_velocity[0]*DN1-val_velocity[1]*DN0;
    val_p_tensor[4*nVar+3]=0.5*(0.5*sqvel*rhooc+(val_velocity[0]*DN0+val_velocity[1]*DN1+val_velocity[2]*DN2)+rhoxc);
    val_p_tensor[4*nVar+4]=0.5*(0.5*sqvel*rhooc-(val_velocity[0]*DN0+val_velocity[1]*DN1+val_velocity[2]*DN2)+rhoxc);
    
  }
  
}

void GetPMatrix_inv_R(su2double *val_density, su2double *val_velocity,
    su2double *val_soundspeed, su2double *val_normal, 
	su2double *val_invp_tensor, su2double sqvel,
	su2double Gamma_Minus_One, int nDim, int nVar) {
  
  su2double rhoxc, c2, gm1, DN0, DN1, DN2, gm1_o_c2, gm1_o_rhoxc;

  rhoxc = *val_density * *val_soundspeed;
  c2 = *val_soundspeed * *val_soundspeed;
  gm1 = Gamma_Minus_One;
  DN0= val_normal[0] / *val_density;
  DN1= val_normal[1] / *val_density;
  DN2 = val_normal[2] / *val_density;
  gm1_o_c2 = gm1/c2;
  gm1_o_rhoxc = gm1/rhoxc;

  if (nDim == 3) {
    
    val_invp_tensor[0]=val_normal[0]-DN2*val_velocity[1] + DN1*val_velocity[2] -val_normal[0]*0.5*gm1_o_c2*sqvel;
    val_invp_tensor[1]=val_normal[0]*gm1_o_c2*val_velocity[0];
    val_invp_tensor[2]=DN2 +val_normal[0]*gm1_o_c2*val_velocity[1];
    val_invp_tensor[3]=-DN1 +val_normal[0]*gm1_o_c2*val_velocity[2];
    val_invp_tensor[4]=-val_normal[0]*gm1_o_c2;

    val_invp_tensor[1*nVar+0]=val_normal[1]+DN2*val_velocity[0] - DN0*val_velocity[2] - val_normal[1]*0.5*gm1_o_c2*sqvel;
    val_invp_tensor[1*nVar+1]=-DN2 + val_normal[1]*gm1_o_c2*val_velocity[0];
    val_invp_tensor[1*nVar+2]=val_normal[1]*gm1_o_c2*val_velocity[1];
    val_invp_tensor[1*nVar+3]=DN0 + val_normal[1]*gm1_o_c2*val_velocity[2];
    val_invp_tensor[1*nVar+4]=-val_normal[1]*gm1_o_c2;

    val_invp_tensor[2*nVar+0]=val_normal[2]-DN1*val_velocity[0] +DN0*val_velocity[1] -val_normal[2]*0.5*gm1_o_c2*sqvel;
    val_invp_tensor[2*nVar+1]=DN1 + val_normal[2]*gm1_o_c2*val_velocity[0];
    val_invp_tensor[2*nVar+2]=-DN0 +val_normal[2]*gm1_o_c2*val_velocity[1];
    val_invp_tensor[2*nVar+3]=val_normal[2]*gm1_o_c2*val_velocity[2];
    val_invp_tensor[2*nVar+4]=-val_normal[2]*gm1_o_c2;

    val_invp_tensor[3*nVar+0]=-(DN0*val_velocity[0]+DN1*val_velocity[1]+DN2*val_velocity[2]) +0.5*gm1_o_rhoxc*sqvel;
    val_invp_tensor[3*nVar+1]=DN0 -gm1_o_rhoxc*val_velocity[0];
    val_invp_tensor[3*nVar+2]=DN1 -gm1_o_rhoxc*val_velocity[1];
    val_invp_tensor[3*nVar+3]=DN2 -gm1_o_rhoxc*val_velocity[2];
    val_invp_tensor[3*nVar+4]=gm1_o_rhoxc;

    val_invp_tensor[4*nVar+0]=(DN0*val_velocity[0]+DN1*val_velocity[1]+DN2*val_velocity[2]) +0.5*gm1_o_rhoxc*sqvel;
    val_invp_tensor[4*nVar+1]=-DN0 -gm1_o_rhoxc*val_velocity[0];
    val_invp_tensor[4*nVar+2]=-DN1 -gm1_o_rhoxc*val_velocity[1];
    val_invp_tensor[4*nVar+3]=-DN2 -gm1_o_rhoxc*val_velocity[2];
    val_invp_tensor[4*nVar+4]=gm1_o_rhoxc;
    
  }
  if (nDim == 2) {
    
    val_invp_tensor[0] = 1.0-0.5*gm1_o_c2*sqvel;
    val_invp_tensor[1]=gm1_o_c2*val_velocity[0];
    val_invp_tensor[2]=gm1_o_c2*val_velocity[1];
    val_invp_tensor[3]=-gm1_o_c2;

    val_invp_tensor[1*nVar+0]=-DN1*val_velocity[0]+DN0*val_velocity[1];
    val_invp_tensor[1*nVar+1]=DN1;
    val_invp_tensor[1*nVar+2]=-DN0;
    val_invp_tensor[1*nVar+3]=0.0;

    val_invp_tensor[2*nVar+0]=-DN0*val_velocity[0]-DN1*val_velocity[1]+0.5*gm1_o_rhoxc*sqvel;
    val_invp_tensor[2*nVar+1]=DN0-gm1_o_rhoxc*val_velocity[0];
    val_invp_tensor[2*nVar+2]=DN1-gm1_o_rhoxc*val_velocity[1];
    val_invp_tensor[2*nVar+3]=gm1_o_rhoxc;

    val_invp_tensor[3*nVar+0]=DN0*val_velocity[0]+DN1*val_velocity[1]+0.5*gm1_o_rhoxc*sqvel;
    val_invp_tensor[3*nVar+1]=-DN0-gm1_o_rhoxc*val_velocity[0];
    val_invp_tensor[3*nVar+2]=-DN1-gm1_o_rhoxc*val_velocity[1];
    val_invp_tensor[3*nVar+3]=gm1_o_rhoxc;
    
  }
}

void GetInviscidProjFlux_R(su2double *val_density, su2double *val_velocity,
     su2double *val_pressure, su2double *val_enthalpy, su2double *val_normal,
     su2double *val_Proj_Flux, int nDim) {

    su2double rhou, rhov, rhow;

  if (nDim == 2) {

    rhou = (*val_density)*val_velocity[0]*val_normal[0];
    rhov = (*val_density)*val_velocity[1]*val_normal[1];

    val_Proj_Flux[0] = rhou;
    val_Proj_Flux[1] = rhou*val_velocity[0]+(*val_pressure)*val_normal[0];
    val_Proj_Flux[2] = rhou*val_velocity[1];
    val_Proj_Flux[3] = rhou*(*val_enthalpy);

    val_Proj_Flux[0] += rhov;
    val_Proj_Flux[1] += rhov*val_velocity[0];
    val_Proj_Flux[2] += rhov*val_velocity[1]+(*val_pressure)*val_normal[1];
    val_Proj_Flux[3] += rhov*(*val_enthalpy);
    
  }
  else {
    
    rhou = (*val_density)*val_velocity[0]*val_normal[0];
    rhov = (*val_density)*val_velocity[1]*val_normal[1];
    rhow = (*val_density)*val_velocity[2]*val_normal[2];

    val_Proj_Flux[0] = rhou;
    val_Proj_Flux[1] = rhou*val_velocity[0]+(*val_pressure)*val_normal[0];
    val_Proj_Flux[2] = rhou*val_velocity[1];
    val_Proj_Flux[3] = rhou*val_velocity[2];
    val_Proj_Flux[4] = rhou*(*val_enthalpy);

    val_Proj_Flux[0] += rhov;
    val_Proj_Flux[1] += rhov*val_velocity[0];
    val_Proj_Flux[2] += rhov*val_velocity[1]+(*val_pressure)*val_normal[1];
    val_Proj_Flux[3] += rhov*val_velocity[2];
    val_Proj_Flux[4] += rhov*(*val_enthalpy);

    val_Proj_Flux[0] += rhow;
    val_Proj_Flux[1] += rhow*val_velocity[0];
    val_Proj_Flux[2] += rhow*val_velocity[1];
    val_Proj_Flux[3] += rhow*val_velocity[2]+(*val_pressure)*val_normal[2];
    val_Proj_Flux[4] += rhow*(*val_enthalpy);

  }

}

void GetInviscidProjJac_R(su2double *val_velocity, su2double *val_energy,
     su2double *val_normal, su2double val_scale,
     su2double **val_Proj_Jac_Tensor, su2double Gamma, 
     int nDim) {
  unsigned short iDim, jDim;
  su2double proj_vel, sqvel, phi, a1, a2;

  a2 = Gamma-1.0;

  switch (nDim) {
  case 2:
    sqvel    = val_velocity[0]*val_velocity[0]+val_velocity[1]*val_velocity[1];
    proj_vel = val_velocity[0]*val_normal[0]+val_velocity[1]*val_normal[1];
    phi = 0.5*a2*sqvel;
    a1 = Gamma*(*val_energy)-phi;

    val_Proj_Jac_Tensor[0][0] = 0.0;
    val_Proj_Jac_Tensor[0][1] = val_scale*val_normal[0];
    val_Proj_Jac_Tensor[0][2] = val_scale*val_normal[1];
    val_Proj_Jac_Tensor[0][3] = 0.0;
    val_Proj_Jac_Tensor[1][0] = val_scale*(val_normal[0]*phi - val_velocity[0]*proj_vel);
    val_Proj_Jac_Tensor[1][1] = val_scale*(val_normal[0]*val_velocity[0]-a2*val_normal[0]*val_velocity[0]);
    val_Proj_Jac_Tensor[1][2] = val_scale*(val_normal[1]*val_velocity[0]-a2*val_normal[0]*val_velocity[1]);
    val_Proj_Jac_Tensor[1][1]+= val_scale*proj_vel;
    val_Proj_Jac_Tensor[1][3] = val_scale*a2*val_normal[0];
    val_Proj_Jac_Tensor[2][0] = val_scale*(val_normal[1]*phi - val_velocity[1]*proj_vel);
    val_Proj_Jac_Tensor[2][1] = val_scale*(val_normal[0]*val_velocity[1]-a2*val_normal[1]*val_velocity[0]);
    val_Proj_Jac_Tensor[2][2] = val_scale*(val_normal[1]*val_velocity[1]-a2*val_normal[1]*val_velocity[1]);
    val_Proj_Jac_Tensor[2][2]+= val_scale*proj_vel;
    val_Proj_Jac_Tensor[2][3] = val_scale*a2*val_normal[1];

    val_Proj_Jac_Tensor[3][0] = val_scale*proj_vel*(phi-a1);
    val_Proj_Jac_Tensor[3][1] = val_scale*(val_normal[0]*a1-a2*val_velocity[0]*proj_vel);
    val_Proj_Jac_Tensor[3][2] = val_scale*(val_normal[1]*a1-a2*val_velocity[1]*proj_vel);
    val_Proj_Jac_Tensor[3][3] = val_scale*Gamma*proj_vel;

  break; // end 2

  case 3: // nDim==3
    sqvel    = val_velocity[0]*val_velocity[0] 
             + val_velocity[1]*val_velocity[1]
             + val_velocity[2]*val_velocity[2];
    proj_vel = val_velocity[0]*val_normal[0] 
             + val_velocity[1]*val_normal[1]
             + val_velocity[2]*val_normal[2];
    phi = 0.5*a2*sqvel;
    a1 = Gamma*(*val_energy)-phi;

    val_Proj_Jac_Tensor[0][0] = 0.0;
    val_Proj_Jac_Tensor[0][1] = val_scale*val_normal[0];
    val_Proj_Jac_Tensor[0][2] = val_scale*val_normal[1];
    val_Proj_Jac_Tensor[0][3] = val_scale*val_normal[2];
    val_Proj_Jac_Tensor[0][4] = 0.0;

    val_Proj_Jac_Tensor[1][0] = val_scale*(val_normal[0]*phi - val_velocity[0]*proj_vel);
    val_Proj_Jac_Tensor[1][1] = val_scale*(val_normal[0]*val_velocity[0]-a2*val_normal[0]*val_velocity[0]);
    val_Proj_Jac_Tensor[1][2] = val_scale*(val_normal[1]*val_velocity[0]-a2*val_normal[0]*val_velocity[1]);
    val_Proj_Jac_Tensor[1][3] = val_scale*(val_normal[2]*val_velocity[0]-a2*val_normal[0]*val_velocity[2]);
    val_Proj_Jac_Tensor[1][1]+= val_scale*proj_vel;
    val_Proj_Jac_Tensor[1][4] = val_scale*a2*val_normal[0];

    val_Proj_Jac_Tensor[2][0] = val_scale*(val_normal[1]*phi - val_velocity[1]*proj_vel);
    val_Proj_Jac_Tensor[2][1] = val_scale*(val_normal[0]*val_velocity[1]-a2*val_normal[1]*val_velocity[0]);
    val_Proj_Jac_Tensor[2][2] = val_scale*(val_normal[1]*val_velocity[1]-a2*val_normal[1]*val_velocity[1]);
    val_Proj_Jac_Tensor[2][3] = val_scale*(val_normal[2]*val_velocity[1]-a2*val_normal[1]*val_velocity[2]);
    val_Proj_Jac_Tensor[2][2]+= val_scale*proj_vel;
    val_Proj_Jac_Tensor[2][4] = val_scale*a2*val_normal[1];

    val_Proj_Jac_Tensor[3][0] = val_scale*(val_normal[2]*phi - val_velocity[2]*proj_vel);
    val_Proj_Jac_Tensor[3][1] = val_scale*(val_normal[0]*val_velocity[2]-a2*val_normal[2]*val_velocity[0]);
    val_Proj_Jac_Tensor[3][2] = val_scale*(val_normal[1]*val_velocity[2]-a2*val_normal[2]*val_velocity[1]);
    val_Proj_Jac_Tensor[3][3] = val_scale*(val_normal[2]*val_velocity[2]-a2*val_normal[2]*val_velocity[2]);
    val_Proj_Jac_Tensor[3][3]+= val_scale*proj_vel;
    val_Proj_Jac_Tensor[3][4] = val_scale*a2*val_normal[2];

    val_Proj_Jac_Tensor[4][0] = val_scale*proj_vel*(phi-a1);
    val_Proj_Jac_Tensor[4][1] = val_scale*(val_normal[0]*a1-a2*val_velocity[0]*proj_vel);
    val_Proj_Jac_Tensor[4][2] = val_scale*(val_normal[1]*a1-a2*val_velocity[1]*proj_vel);
    val_Proj_Jac_Tensor[4][3] = val_scale*(val_normal[2]*a1-a2*val_velocity[2]*proj_vel);
    val_Proj_Jac_Tensor[4][4] = val_scale*Gamma*proj_vel;
 break; // end 3

 default: // general code
  sqvel=0.0, proj_vel=0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    sqvel    += val_velocity[iDim]*val_velocity[iDim];
    proj_vel += val_velocity[iDim]*val_normal[iDim];
  }

  phi = 0.5*a2*sqvel;
  a1 = Gamma*(*val_energy)-phi;

   val_Proj_Jac_Tensor[0][0] = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    val_Proj_Jac_Tensor[0][iDim+1] = val_scale*val_normal[iDim];
  val_Proj_Jac_Tensor[0][nDim+1] = 0.0;

  for (iDim = 0; iDim < nDim; iDim++) {
    val_Proj_Jac_Tensor[iDim+1][0] = val_scale*(val_normal[iDim]*phi - val_velocity[iDim]*proj_vel);
    for (jDim = 0; jDim < nDim; jDim++)
      val_Proj_Jac_Tensor[iDim+1][jDim+1] = val_scale*(val_normal[jDim]*val_velocity[iDim]-a2*val_normal[iDim]*val_velocity[jDim]);
    val_Proj_Jac_Tensor[iDim+1][iDim+1] += val_scale*proj_vel;
    val_Proj_Jac_Tensor[iDim+1][nDim+1] = val_scale*a2*val_normal[iDim];
  }

  val_Proj_Jac_Tensor[nDim+1][0] = val_scale*proj_vel*(phi-a1);
  for (iDim = 0; iDim < nDim; iDim++)
    val_Proj_Jac_Tensor[nDim+1][iDim+1] = val_scale*(val_normal[iDim]*a1-a2*val_velocity[iDim]*proj_vel);
  val_Proj_Jac_Tensor[nDim+1][nDim+1] = val_scale*Gamma*proj_vel;
 }
}
