#include "../include/solver_boom.hpp"
#include <cmath>
#include <algorithm>
#include <complex>
#include <valarray>
#include <iostream>
#include <fstream>

typedef std::complex<su2double> Complex;
typedef std::valarray<Complex> CArray;

CBoom_AugBurgers::CBoom_AugBurgers(){

}

CBoom_AugBurgers::CBoom_AugBurgers(CSolver *solver, CConfig *config, CGeometry *geometry){

  int rank, nProcessor = 1;

  rank = SU2_MPI::GetRank();
  nProcessor = SU2_MPI::GetSize();

  Kind_Boom_Cost = config->GetKind_ObjFunc();
  CFL_reduce = config->GetBoom_cfl_reduce();
  Kind_Step = config->GetBoom_step_type();
  Step_size = config->GetBoom_step_size();
  Step_growth = config->GetBoom_step_growth();
  AD_Mode = config->GetAD_Mode();

  /*---Make sure to read in hard-coded values in the future!---*/

  nDim = geometry->GetnDim();

  /*---Flight variables---*/
  flt_h = config->GetBoom_flt_h(); // altitude [m]
  flt_M = config->GetMach();
  flt_psi = M_PI/2.;  // heading angle [rad]
  flt_gamma = 0.; // flight path angle [rad]
  flt_mu = asin(1./flt_M); // Mach angle [rad]

  /*---Atmosphere variables---*/
  atm_g = config->GetGamma();

  /*---Values from config file---*/
  const su2double deg2rad = M_PI/180.;
  n_prop = config->GetBoom_N_prof();
  if(nDim == 2){
    ray_N_phi = 1;
    ray_phi = new su2double[1];
    ray_phi[0] = 0.0;
  }
  else{
    ray_N_phi = config->GetBoom_N_phi();
    ray_phi = config->GetBoom_phi();
    if(ray_phi == NULL){
      ray_N_phi = 1;
      ray_phi = new su2double[1];
      ray_phi[0] = 0.0;
    }
    else{
      for(unsigned short iPhi = 0; iPhi < ray_N_phi; iPhi++){
        ray_phi[iPhi] *= deg2rad;
      }
    }
  }
  ray_r0 = config->GetBoom_r0();
  scale_L = config->GetRefLength();

  /*---Cost function---*/
  PLdB = new su2double[ray_N_phi];

  /*---Set reference pressure, make sure signal is dimensional---*/
  su2double Pressure_FreeStream=config->GetPressure_FreeStream();
  su2double Pressure_Ref;
  if (config->GetRef_NonDim() == DIMENSIONAL) {
    Pressure_Ref      = 1.0;
  }
  else if (config->GetRef_NonDim() == FREESTREAM_PRESS_EQ_ONE) {
    Pressure_Ref      = Pressure_FreeStream;     // Pressure_FreeStream = 1.0
  }
  else if (config->GetRef_NonDim() == FREESTREAM_VEL_EQ_MACH) {
    Pressure_Ref      = atm_g*Pressure_FreeStream; // Pressure_FreeStream = 1.0/Gamma
  }
  else if (config->GetRef_NonDim() == FREESTREAM_VEL_EQ_ONE) {
    Pressure_Ref      = flt_M*flt_M*atm_g*Pressure_FreeStream; // Pressure_FreeStream = 1.0/(Gamma*(M_inf)^2)
  }
  
  if(rank == MASTER_NODE && !AD_Mode)
    cout << "Pressure_Ref = " << Pressure_Ref << ", Pressure_FreeStream = " << Pressure_FreeStream << endl;

  /*---Perform search on domain to determine where line intersects boundary---*/
  if(rank == MASTER_NODE)
    cout << "Search for start of line." << endl;
  SearchLinear(config, geometry, ray_r0, ray_phi);

  /*---Walk through neighbors to determine all points containing line---*/
  if(rank == MASTER_NODE)
    cout << "Extract line." << endl;
  for(unsigned short iPhi = 0; iPhi < ray_N_phi; iPhi++){
    if(nPanel[iPhi] > 0){
      ExtractLine(geometry, ray_r0, iPhi);
    }
  }

  /*---Initialize some signal parameters---*/
  nPointID = new unsigned long[ray_N_phi];
  PointID = new unsigned long*[ray_N_phi];
  signal.len = new unsigned long[ray_N_phi];
  signal.x = new su2double*[ray_N_phi];
  signal.p_prime = new su2double*[ray_N_phi];
  Coord_original = new su2double**[ray_N_phi];
    
  for( unsigned short iPhi = 0; iPhi < ray_N_phi; iPhi++){
    PointID[iPhi] = NULL;
    signal.x[iPhi] = NULL;
    signal.p_prime[iPhi] = NULL;
  }
    
  if(AD_Mode) AD::StartRecording();

  /*---Interpolate pressures along line---*/
  if(rank == MASTER_NODE)
    cout << "Extract pressure signature." << endl;
    
  nPointID_loc = 0;
  for(unsigned short iPhi = 0; iPhi < ray_N_phi; iPhi++){
    nPointID[iPhi] = 0;
    if(nPanel[iPhi] > 0){      
      ExtractPressure(solver, config, geometry, iPhi);
    }
  }

  nPointID_proc = new unsigned long[nProcessor];
  unsigned long iPanel, panelCount, totSig, maxSig;
  unsigned long *nPanel_loc = new unsigned long[nProcessor];
  nPointID_tot = 0;
  for(unsigned short iPhi = 0; iPhi < ray_N_phi; iPhi++){
    for(iPanel = 0; iPanel < nPanel[iPhi]; iPanel++){
      signal.p_prime[iPhi][iPanel] = signal.p_prime[iPhi][iPanel]*Pressure_Ref - Pressure_FreeStream;
    }

    totSig = 0;
    maxSig = 0;

    SU2_MPI::Allreduce(&nPanel[iPhi], &totSig, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&nPanel[iPhi], &maxSig, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
    SU2_MPI::Gather(&nPanel[iPhi], 1, MPI_UNSIGNED_LONG, nPanel_loc, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Gather(&nPointID[iPhi], 1, MPI_UNSIGNED_LONG, nPointID_proc, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);

    su2double* Buffer_Recv_Press = NULL;
    su2double* Buffer_Recv_x = NULL;
    su2double* Buffer_Send_Press = new su2double[maxSig];
    su2double* Buffer_Send_x = new su2double[maxSig];
    if(rank == MASTER_NODE){
      Buffer_Recv_Press = new su2double[nProcessor*maxSig];
      Buffer_Recv_x = new su2double[nProcessor*maxSig];
    }
    
    for(iPanel = 0; iPanel < nPanel[iPhi]; iPanel++){
      Buffer_Send_x[iPanel] = signal.x[iPhi][iPanel];
      Buffer_Send_Press[iPanel] = signal.p_prime[iPhi][iPanel];
    }

    SU2_MPI::Gather(Buffer_Send_Press, maxSig, MPI_DOUBLE, Buffer_Recv_Press,  maxSig , MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Gather(Buffer_Send_x, maxSig, MPI_DOUBLE, Buffer_Recv_x,  maxSig , MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);

    if (rank == MASTER_NODE)
      cout << "Gathered signal data to MASTER_NODE." << endl;

    if (rank == MASTER_NODE){
      panelCount = 0;
      nPanel[iPhi] = totSig;
      signal.x[iPhi] = new su2double[nPanel[iPhi]];
      signal.p_prime[iPhi] = new su2double[nPanel[iPhi]];
      for(unsigned int iProcessor = 0; iProcessor < nProcessor; iProcessor++){
        for(iPanel = 0; iPanel < nPanel_loc[iProcessor]; iPanel++){
          signal.x[iPhi][panelCount] = Buffer_Recv_x[iProcessor*maxSig+iPanel];
          signal.p_prime[iPhi][panelCount] = Buffer_Recv_Press[iProcessor*maxSig+iPanel];
          panelCount++;
        }
        nPointID_tot += nPointID_proc[iProcessor];
      }

      /*---Sort signal in order of x-coordinate---*/
      cout << "Sorting signal data. " << nPanel[iPhi] << " points to sort." << endl;
      MergeSort(signal.x[iPhi], signal.p_prime[iPhi], 0, totSig-1);

      /*---Check for duplicate points---*/
      signal.len[iPhi] = nPanel[iPhi];
      for(iPanel = 1; iPanel < nPanel[iPhi]; iPanel++){
        if(abs(signal.x[iPhi][iPanel-1]-signal.x[iPhi][iPanel]) < 1.0E-6){
          for(unsigned long jPanel = iPanel; jPanel < nPanel[iPhi]; jPanel++){
            signal.x[iPhi][jPanel-1] = signal.x[iPhi][jPanel];
            signal.p_prime[iPhi][jPanel-1] = signal.p_prime[iPhi][jPanel];
          }
          iPanel--;
          nPanel[iPhi]--;
        }
      }

      if(nPanel[iPhi] != totSig){
        cout << "Eliminating duplicate points." << endl;
        su2double *xtmp = new su2double[nPanel[iPhi]], *ptmp = new su2double[nPanel[iPhi]];
        for(iPanel = 0; iPanel < nPanel[iPhi]; iPanel++){
          xtmp[iPanel] = signal.x[iPhi][iPanel];
          ptmp[iPanel] = signal.p_prime[iPhi][iPanel];
        }
        signal.len[iPhi] = nPanel[iPhi];
        signal.x[iPhi] = new su2double[nPanel[iPhi]];
        signal.p_prime[iPhi] = new su2double[nPanel[iPhi]];
        for(iPanel = 0; iPanel < nPanel[iPhi]; iPanel++){
          signal.x[iPhi][iPanel] = xtmp[iPanel];
          signal.p_prime[iPhi][iPanel] = ptmp[iPanel];
        }
        delete [] xtmp;
        delete [] ptmp;
        totSig = nPanel[iPhi];
      }
    }
  }

  /*---Initialize sensitivities---*/
  dJdU = NULL;
  dJdX = NULL;
  if(AD_Mode){
    
    dJdU = new su2double**[ray_N_phi];
    dJdX = new su2double**[ray_N_phi];
    
    for(unsigned short iPhi = 0; iPhi < ray_N_phi; iPhi++){
      dJdU[iPhi] = new su2double* [nDim+3];
      dJdX[iPhi] = new su2double* [nDim];
      
      for(int iDim = 0; iDim < nDim+3; iDim++){
        dJdU[iPhi][iDim] = NULL;
        if(nPointID[iPhi] > 0){
          dJdU[iPhi][iDim] = new su2double[nPointID[iPhi]];
          for(iPanel = 0;  iPanel < nPointID[iPhi]; iPanel++){
            dJdU[iPhi][iDim][iPanel] = 0.0;
          }
        }
      }
      
      for(int iDim = 0; iDim < nDim; iDim++){
        dJdX[iPhi][iDim] = NULL;
        if(nPointID[iPhi] > 0){
          dJdX[iPhi][iDim] = new su2double[nPointID[iPhi]];
          for(iPanel = 0;  iPanel < nPointID[iPhi]; iPanel++){
            dJdX[iPhi][iDim][iPanel] = 0.0;
          }
        }
      }
        
    }

    if (rank==MASTER_NODE) cout << "Sensitivities initialized." << endl;
    
  }
    
  if(rank == MASTER_NODE) cout << "ABE initialized." << endl;

}

CBoom_AugBurgers::~CBoom_AugBurgers(void){
    
  /*---Clear up memory from dJ/dU and dJ/dX---*/
  if(dJdU != NULL){
    for(unsigned short iPhi = 0; iPhi < ray_N_phi; iPhi++){
      if(dJdU[iPhi] != NULL){
        for (unsigned short i = 0; i < nDim+3; i++){
          if(dJdU[iPhi][i] != NULL) delete [] dJdU[iPhi][i];
        }
        for (unsigned short i = 0; i < nDim; i++){
          if(dJdX[iPhi][i] != NULL) delete [] dJdX[iPhi][i];
        }
        delete [] dJdU[iPhi];
        delete [] dJdX[iPhi];
      }
    }
    delete [] dJdU;
    delete [] dJdX;
  }
  if(nPointID != NULL) delete [] nPointID;

}

void MergeSort(su2double x[], su2double p[], int l, int r){
  if(l < r){
    int m = l + (r-l)/2;
    MergeSort(x, p, l, m);
    MergeSort(x, p, m+1, r);
    merge(x, p, l, m, r);
  }
}

void MergeSort(su2double y[], unsigned long k[], int l, int r){
  if(l < r){
    int m = l + (r-l)/2;
    MergeSort(y, k, l, m);
    MergeSort(y, k, m+1, r);
    merge(y, k, l, m, r);
  }
}

void merge(su2double x[], su2double p[], int l, int m, int r){
  int i, j, k;
  int n1 = m-l+1;
  int n2 = r-m;

  su2double L[n1], R[n2], Lp[n1], Rp[n2];

  /*--- Temporary arrays ---*/
  for(i = 0; i < n1; i++){
    L[i]  = x[l+i];
    Lp[i] = p[l+i];
  }
  for(j = 0; j < n2; j++){
    R[j]  = x[m+1+j];
    Rp[j] = p[m+1+j];
  }

  i = 0;
  j = 0;
  k = l;

  /*--- Begin sorting ---*/
  while(i < n1 && j < n2){
    if(L[i] <= R[j]){
      x[k] = L[i];
      p[k] = Lp[i];
      i++;
    }
    else{
      x[k] = R[j];
      p[k] = Rp[j];
      j++;
    }
    k++;
  }

  /*--- Fill rest of arrays ---*/
  while(i < n1){
    x[k] = L[i];
    p[k] = Lp[i];
    i++;
    k++;
  }

  while(j < n2){
    x[k] = R[j];
    p[k] = Rp[j];
    j++;
    k++;
  }
}

void merge(su2double x[], unsigned long p[], int l, int m, int r){
  int i, j, k;
  int n1 = m-l+1;
  int n2 = r-m;

  su2double L[n1], R[n2];
  unsigned long Lp[n1], Rp[n2];

  /*--- Temporary arrays ---*/
  for(i = 0; i < n1; i++){
    L[i]  = x[l+i];
    Lp[i] = p[l+i];
  }
  for(j = 0; j < n2; j++){
    R[j]  = x[m+1+j];
    Rp[j] = p[m+1+j];
  }

  i = 0;
  j = 0;
  k = l;

  /*--- Begin sorting ---*/
  while(i < n1 && j < n2){
    if(L[i] <= R[j]){
      x[k] = L[i];
      p[k] = Lp[i];
      i++;
    }
    else{
      x[k] = R[j];
      p[k] = Rp[j];
      j++;
    }
    k++;
  }

  /*--- Fill rest of arrays ---*/
  while(i < n1){
    x[k] = L[i];
    p[k] = Lp[i];
    i++;
    k++;
  }

  while(j < n2){
    x[k] = R[j];
    p[k] = Rp[j];
    j++;
    k++;
  }
}

void Isoparameters(unsigned short nDim, unsigned short nDonor, su2double *X, su2double *xj, su2double *isoparams) {
  short iDonor,iDim,k; // indices
  su2double tmp, tmp2;

  su2double x[nDim+1];
  su2double x_tmp[nDim+1];
  su2double Q[nDonor*nDonor];
  su2double R[nDonor*nDonor];
  su2double A[(nDim+1)*nDonor];
  su2double *A2    = NULL;
  su2double x2[nDim+1];
  
  bool test[nDim+1];
  bool testi[nDim+1];
  
  su2double eps = 1E-10;
  
  short n = nDim+1;

  if (nDonor>2) {
    /*--- Create Matrix A: 1st row all 1's, 2nd row x coordinates, 3rd row y coordinates, etc ---*/
    /*--- Right hand side is [1, \vec{x}']'---*/
    for (iDonor=0; iDonor<nDonor; iDonor++) {
      isoparams[iDonor]=0;
      A[iDonor] = 1.0;
      for (iDim=0; iDim<nDim; iDim++)
        A[(iDim+1)*nDonor+iDonor]=X[iDim*nDonor+iDonor];
    }

    x[0] = 1.0;
    for (iDim=0; iDim<nDim; iDim++)
      x[iDim+1]=xj[iDim];

    /*--- Eliminate degenerate rows:
     * for example, if z constant including the z values will make the system degenerate
     * TODO: improve efficiency of this loop---*/
    test[0]=true; // always keep the 1st row
    for (iDim=1; iDim<n; iDim++) {
      // Test this row against all previous
      test[iDim]=true; // Assume that it is not degenerate
      for (k=0; k<iDim; k++) {
        tmp=0; tmp2=0;
        for (iDonor=0;iDonor<nDonor;iDonor++) {
          tmp+= A[iDim*nDonor+iDonor]*A[iDim*nDonor+iDonor];
          tmp2+=A[k*nDonor+iDonor]*A[k*nDonor+iDonor];
        }
        tmp  = pow(tmp,0.5);
        tmp2 = pow(tmp2,0.5);
        testi[k]=false;
        for (iDonor=0; iDonor<nDonor; iDonor++) {
          // If at least one ratio is non-matching row iDim is not degenerate w/ row k
          if (abs(A[iDim*nDonor+iDonor]/tmp-A[k*nDonor+iDonor]/tmp2) > eps)
            testi[k]=true;
        }
        // If any of testi (k<iDim) are false, row iDim is degenerate
        test[iDim]=(test[iDim] && testi[k]);
      }
      if (!test[iDim]) n--;
    }

    /*--- Initialize A2 now that we might have a smaller system --*/
    A2 = new su2double[n*nDonor];
    iDim=0;
    /*--- Copy only the rows that are non-degenerate ---*/
    for (k=0; k<n; k++) {
      if (test[k]) {
        for (iDonor=0;iDonor<nDonor;iDonor++ ) {
          A2[nDonor*iDim+iDonor]=A[nDonor*k+iDonor];
        }
        x2[iDim]=x[k];
        iDim++;
      }
    }

    /*--- Initialize Q,R to 0 --*/
    for (k=0; k<nDonor*nDonor; k++) {
      Q[k]=0;
      R[k]=0;
    }
    /*--- TODO: make this loop more efficient ---*/
    /*--- Solve for rectangular Q1 R1 ---*/
    for (iDonor=0; iDonor<nDonor; iDonor++) {
      tmp=0;
      for (iDim=0; iDim<n; iDim++)
        tmp += (A2[iDim*nDonor+iDonor])*(A2[iDim*nDonor+iDonor]);

      R[iDonor*nDonor+iDonor]= pow(tmp,0.5);
      if (tmp>eps && iDonor<n) {
        for (iDim=0; iDim<n; iDim++)
          Q[iDim*nDonor+iDonor]=A2[iDim*nDonor+iDonor]/R[iDonor*nDonor+iDonor];
      }
      else if (tmp!=0) {
        for (iDim=0; iDim<n; iDim++)
          Q[iDim*nDonor+iDonor]=A2[iDim*nDonor+iDonor]/tmp;
      }
      for (iDim=iDonor+1; iDim<nDonor; iDim++) {
        tmp=0;
        for (k=0; k<n; k++)
          tmp+=A2[k*nDonor+iDim]*Q[k*nDonor+iDonor];

        R[iDonor*nDonor+iDim]=tmp;

        for (k=0; k<n; k++)
          A2[k*nDonor+iDim]=A2[k*nDonor+iDim]-Q[k*nDonor+iDonor]*R[iDonor*nDonor+iDim];
      }
    }
    /*--- x_tmp = Q^T * x2 ---*/
    for (iDonor=0; iDonor<nDonor; iDonor++)
      x_tmp[iDonor]=0.0;
    for (iDonor=0; iDonor<nDonor; iDonor++) {
      for (iDim=0; iDim<n; iDim++)
        x_tmp[iDonor]+=Q[iDim*nDonor+iDonor]*x2[iDim];
    }

    /*--- solve x_tmp = R*isoparams for isoparams: upper triangular system ---*/
    for (iDonor = n-1; iDonor>=0; iDonor--) {
      if (R[iDonor*nDonor+iDonor]>eps)
        isoparams[iDonor]=x_tmp[iDonor]/R[iDonor*nDonor+iDonor];
      else
        isoparams[iDonor]=0;
      for (k=0; k<iDonor; k++)
        x_tmp[k]=x_tmp[k]-R[k*nDonor+iDonor]*isoparams[iDonor];
    }
  }
  else {
    /*-- For 2-donors (lines) it is simpler: */
    tmp =  pow(X[0*nDonor+0]- X[0*nDonor+1],2.0);
    tmp += pow(X[1*nDonor+0]- X[1*nDonor+1],2.0);
    tmp = sqrt(tmp);

    tmp2 = pow(X[0*nDonor+0] - xj[0],2.0);
    tmp2 += pow(X[1*nDonor+0] - xj[1],2.0);
    tmp2 = sqrt(tmp2);
    isoparams[1] = tmp2/tmp;

    tmp2 = pow(X[0*nDonor+1] - xj[0],2.0);
    tmp2 += pow(X[1*nDonor+1] - xj[1],2.0);
    tmp2 = sqrt(tmp2);
    isoparams[0] = tmp2/tmp;
  }

  /*--- Isoparametric coefficients have been calculated. Run checks to eliminate outside-element issues ---*/
  if (nDonor==4 && nDim == 2) {
    //-- Bilinear coordinates, bounded by [-1,1] ---
    su2double xi, eta;
    xi = -isoparams[0]+isoparams[1]+isoparams[2]-isoparams[3];
    eta = 2.0*(isoparams[2]-isoparams[0]) - xi;
    if (xi>1.0) xi=1.0;
    if (xi<-1.0) xi=-1.0;
    if (eta>1.0) eta=1.0;
    if (eta<-1.0) eta=-1.0;
    isoparams[0]=0.25*(1-xi)*(1-eta);
    isoparams[1]=0.25*(1+xi)*(1-eta);
    isoparams[2]=0.25*(1+xi)*(1+eta);
    isoparams[3]=0.25*(1-xi)*(1+eta);
  }
  if (nDonor<4) {
    tmp = 0.0; // value for normalization
    tmp2=0; // check for maximum value, to be used to id nearest neighbor if necessary
    k=0; // index for maximum value
    for (iDonor=0; iDonor< nDonor; iDonor++) {
      if (isoparams[iDonor]>tmp2) {
        k=iDonor;
        tmp2=isoparams[iDonor];
      }
      // [0,1]
      if (isoparams[iDonor]<0) isoparams[iDonor]=0;
      if (isoparams[iDonor]>1) isoparams[iDonor] = 1;
      tmp +=isoparams[iDonor];
    }
    if (tmp>0)
      for (iDonor=0; iDonor< nDonor; iDonor++)
        isoparams[iDonor]=isoparams[iDonor]/tmp;
    else {
      isoparams[k] = 1.0;
    }
  }
  
  if (A2 != NULL) delete [] A2;

}

void CBoom_AugBurgers::SearchLinear(CConfig *config, CGeometry *geometry, 
               const su2double r0, const su2double *phi){
  
  /*--- Loop over boundary markers ---*/
  unsigned long nMarker = config->GetnMarker_All();
  unsigned long iMarker, iVertex, iPoint;
  unsigned long jElem;
  unsigned short iElem, nElem;

  bool inside=false;

  su2double x = 0.0, y = 0.0, z = 0.0;
  su2double *Coord = NULL;
  su2double *pp0 = NULL, *pp1 = NULL;

  nPanel = new unsigned long[ray_N_phi];
  for(unsigned short i = 0; i < ray_N_phi; i++){
    nPanel[i] = 0;
  }

  pp0 = new su2double[nDim];
  pp1 = new su2double[nDim];

  /*--- Search on boundaries ---*/  
  for(iMarker = 0; iMarker < nMarker; iMarker++){
    /*--- Only look at farfield boundary (or send/recv for parallel computations) ---*/
    if(config->GetMarker_All_KindBC(iMarker) == FAR_FIELD   || 
      config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE ||
      config->GetMarker_All_KindBC(iMarker) == SYMMETRY_PLANE){
      for(iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++){
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Coord = geometry->node[iPoint]->GetCoord();
        y = SU2_TYPE::GetValue(Coord[1]);
        if(nDim == 3) z = SU2_TYPE::GetValue(Coord[2]);

        /*--- Only look at points below aircraft ---*/
        if((nDim == 2 && y < 0.0) || (nDim == 3 && z < 0.0)){
          x = SU2_TYPE::GetValue(Coord[0]);
          nElem = geometry->node[iPoint]->GetnElem();
          for(iElem = 0; iElem < nElem; iElem++){
            jElem = geometry->node[iPoint]->GetElem(iElem);
            if(jElem < geometry->GetnElem()){
              for(unsigned short iPhi = 0; iPhi < ray_N_phi; iPhi++){
                inside = InsideElem(geometry, r0, phi[iPhi], jElem, pp0, pp1);
                if(inside){
                  if(nPanel[iPhi] == 0){
                    nPanel[iPhi] = 1;
                    pointID_original = new unsigned long*[ray_N_phi];
                    pointID_original[iPhi] = new unsigned long[1];
                    pointID_original[iPhi][0] = jElem;
                  }

                  else{
                    nPanel[iPhi]++;

                    unsigned long *pointID_tmp = new unsigned long[nPanel[iPhi]-1];
                    for(unsigned long i = 0; i < nPanel[iPhi]-1; i++){
                      pointID_tmp[i] = pointID_original[iPhi][i];
                    }
                    delete [] pointID_original[iPhi];

                    pointID_original[iPhi] = new unsigned long[nPanel[iPhi]];
                    for(unsigned long i = 0; i < nPanel[iPhi]-1; i++){
                      pointID_original[iPhi][i] = pointID_tmp[i];
                    }
                    delete [] pointID_tmp;

                    pointID_original[iPhi][nPanel[iPhi]-1] = jElem;
                  }

                  break;
                }
              }
            }
          }
        }
      }
    }
  }

  if(pp0 != NULL) delete [] pp0;
  if(pp1 != NULL) delete [] pp1;

}

void CBoom_AugBurgers::ExtractLine(CGeometry *geometry, const su2double r0, unsigned short iPhi){
  bool inside, inside_iPanel, addPanel, end = false;
  unsigned short iElem, nElem;
  unsigned long jElem, jElem_m1, nElem_tot = geometry->GetnElem();

  unsigned long *pointID_tmp;
  su2double *pp0 = new su2double[nDim], *pp1 = new su2double[nDim];

  while(!end){
    inside_iPanel = false;
    for(unsigned long jPanel = 0; jPanel < nPanel[iPhi]; jPanel++){
      jElem_m1 = pointID_original[iPhi][jPanel];
      nElem = geometry->elem[jElem_m1]->GetnNeighbor_Elements();

      for(iElem = 0; iElem < nElem; iElem++){
        addPanel = true;
        inside = false;
        jElem = geometry->elem[jElem_m1]->GetNeighbor_Elements(iElem);
        /*--- Check element ---*/
        if(jElem < nElem_tot){
          inside = InsideElem(geometry, r0, ray_phi[iPhi], jElem, pp0, pp1);
          if(inside){
            for(unsigned long iPanel = 0; iPanel < nPanel[iPhi]; iPanel++){
              if(jElem == pointID_original[iPhi][iPanel]){
                addPanel = false;
                break;
              }
            }
            if(addPanel){ // If not a repeat panel
              nPanel[iPhi]++;

              pointID_tmp = new unsigned long[nPanel[iPhi]-1];
              for(unsigned long i = 0; i < nPanel[iPhi]-1; i++){
                pointID_tmp[i] = pointID_original[iPhi][i];
              }
              delete [] pointID_original[iPhi];

              pointID_original[iPhi] = new unsigned long[nPanel[iPhi]];
              for(unsigned long i = 0; i < nPanel[iPhi]-1; i++){
                pointID_original[iPhi][i] = pointID_tmp[i];
              }
              delete [] pointID_tmp;
              pointID_original[iPhi][nPanel[iPhi]-1] = jElem;
              inside_iPanel = true;
            }

          }
        }
      }
    }
    if(!inside_iPanel){
      end = true;
    }
  }

}

void CBoom_AugBurgers::ExtractPressure(CSolver *solver, CConfig *config, CGeometry *geometry, unsigned short iPhi){
  unsigned short iDim, iNode, nNode;
  unsigned long iElem, jElem, jNode, jjNode, iPoint;
  unsigned long nNode_list, *jNode_list;
  su2double rho, rho_ux, rho_uy, rho_uz, rho_E, TKE;
  su2double ux, uy, uz, StaticEnergy;
  bool addNode;

  su2double **isoparams = NULL;
  su2double *X_donor = NULL;
  su2double *pp0 = new su2double[nDim], *pp1 = new su2double[nDim];

  isoparams = new su2double*[nPanel[iPhi]];

  for(unsigned long i = 0; i < nPanel[iPhi]; i++){
    jElem = pointID_original[iPhi][i];
    nNode = geometry->elem[jElem]->GetnNodes();
    nPointID[iPhi] += nNode;
    isoparams[i] = NULL;
  }
  PointID[iPhi] = new unsigned long[nPointID[iPhi]];
  jNode_list = new unsigned long[nPointID[iPhi]];
  nNode_list = 0;

  /*--- Compile list of all nodes ---*/
  for(unsigned long i = 0; i < nPanel[iPhi]; i++){
    jElem = pointID_original[iPhi][i];
    nNode = geometry->elem[jElem]->GetnNodes();
    for(iNode = 0; iNode < nNode; iNode++){
      jNode = geometry->elem[jElem]->GetNode(iNode);

      if(nNode_list == 0){
        jNode_list[nNode_list] = jNode;
        nNode_list++;
      }
      else{
        addNode = true;
        for(unsigned long ii = 0; ii < nNode_list; ii++){
          if(jNode == jNode_list[ii]){
            addNode = false;
            break;
          }
        }
        if(addNode){
          jNode_list[nNode_list] = jNode;
          nNode_list++;
        }
      }
    }
  }

  /*---Register coordinates as input for adjoint computation---*/
  if (AD_Mode){
    for(iPoint = 0; iPoint < nNode_list; iPoint++){
      jNode = jNode_list[iPoint];
      for(iDim = 0; iDim < nDim; iDim++){
        AD::RegisterInput(geometry->node[jNode]->GetCoord()[iDim] );
      }
    }
  }
    
  /*---Extract coordinates---*/
  pp0 = new su2double[nDim];
  pp1 = new su2double[nDim];
  signal.x[iPhi] = new su2double[nPanel[iPhi]];
  signal.p_prime[iPhi] = new su2double[nPanel[iPhi]];
  Coord_original[iPhi] = new su2double*[nPanel[iPhi]];
  for(unsigned long i = 0; i < nPanel[iPhi]; i++){
    Coord_original[iPhi][i] = new su2double[nDim];
    jElem = pointID_original[iPhi][i];
    addNode = InsideElem(geometry, ray_r0, ray_phi[iPhi], jElem, pp0, pp1);
    Coord_original[iPhi][i][0] = (pp0[0] + pp1[0])/2.0;
    if(nDim == 2){
      Coord_original[iPhi][i][1] = -ray_r0;
    }
      
    else{
      Coord_original[iPhi][i][1] = ray_r0*sin(ray_phi[iPhi]);
      Coord_original[iPhi][i][2] = -ray_r0*cos(ray_phi[iPhi]);
    }
      
    /*--- Get info needed for isoparameter computation ---*/
    nNode = geometry->elem[jElem]->GetnNodes();
    X_donor = new su2double[nDim*nNode];
    for(iNode = 0; iNode < nNode; iNode++){
      jNode = geometry->elem[jElem]->GetNode(iNode);
      for(iDim = 0; iDim < nDim; iDim++){
        X_donor[iDim*nNode + iNode] = geometry->node[jNode]->GetCoord(iDim);
      }
    }
      
    /*--- Compute isoparameters ---*/
    isoparams[i] = new su2double[nNode];
    Isoparameters(nDim, nNode, X_donor, Coord_original[iPhi][i], isoparams[i]);
      
    /*--- x-locations of nearfield signal ---*/
    signal.x[iPhi][i] = Coord_original[iPhi][i][0];
    signal.p_prime[iPhi][i] = 0.0;
  }

    
  /*--- Now interpolate pressure ---*/
  for(iPoint = 0; iPoint < nNode_list; iPoint++){
    jNode = jNode_list[iPoint];
    PointID[iPhi][iPoint] = geometry->node[jNode]->GetGlobalIndex();

    /*---Extract conservative flow data---*/
    rho = solver->node[jNode]->GetSolution(nDim);
    rho_ux = solver->node[jNode]->GetSolution(nDim+1);
    rho_uy = solver->node[jNode]->GetSolution(nDim+2);
    if(nDim == 3) rho_uz = solver->node[jNode]->GetSolution(nDim+3);
    rho_E = solver->node[jNode]->GetSolution(2*nDim+1);
    TKE = 0.0;

    /*---Register conservative variables as input for adjoint computation---*/
    if (AD_Mode){
      AD::RegisterInput(rho );
      AD::RegisterInput(rho_ux );
      AD::RegisterInput(rho_uy );
      if (nDim==3) AD::RegisterInput(rho_uz );
      AD::RegisterInput(rho_E );
      AD::RegisterInput(TKE );
    }

    ux = rho_ux/rho;
    uy = rho_uy/rho;
    uz = 0.0;
    if (nDim == 3) uz = rho_uz/rho;
    StaticEnergy = rho_E/rho-0.5*(ux*ux+uy*uy+uz*uz)-TKE;

    /*--- Check if node is part of any elements ---*/
    for(iElem = 0; iElem < nPanel[iPhi]; iElem++){
      jElem = pointID_original[iPhi][iElem];
      nNode = geometry->elem[jElem]->GetnNodes();
      for(iNode = 0; iNode < nNode; iNode++){
        jjNode = geometry->elem[jElem]->GetNode(iNode);
        /*--- If node surrounds element, add contribution of node pressure ---*/
        if(jNode == jjNode){
          signal.p_prime[iPhi][iElem] += (config->GetGamma()-1)*rho*StaticEnergy*isoparams[iElem][iNode];
        }
      }
    }
  }

  nPointID[iPhi] = nNode_list;

  /*--- Clean up interpolation variables ---*/
  if(isoparams != NULL){
    for(iElem = 0; iElem < nPanel[iPhi]; iElem++){
      if(isoparams[iElem] != NULL){
        delete [] isoparams[iElem];
      }
    }
    delete [] isoparams;
  }

  if(X_donor != NULL) delete [] X_donor;
  
}

bool CBoom_AugBurgers::InsideElem(CGeometry *geometry, su2double r0, su2double phi, unsigned long jElem, su2double *pp0, su2double *pp1){
  bool inside = false;
  unsigned long iPoint, jNode;
  unsigned short iNode, nNode, count, intersect;

  su2double *ppp0 = new su2double[nDim];
  su2double *ppp1 = new su2double[nDim];

  /*--- Now determine if line intersects element ---*/
  if(nDim == 2){ // Check if line intersects any edge

    /*--- Store node coordinates ---*/
    nNode = geometry->elem[jElem]->GetnNodes();
    su2double **Coord_elem = new su2double*[nNode];
    for(iNode = 0; iNode < nNode; iNode++){
      jNode = geometry->elem[jElem]->GetNode(iNode);

      Coord_elem[iNode] = new su2double[nDim];
      for(unsigned short iDim = 0; iDim < nDim; iDim++){
        Coord_elem[iNode][iDim] = geometry->node[jNode]->GetCoord(iDim);
      }
    }

    count = 0;
    for(unsigned short iEdge = 0; iEdge < nNode; iEdge++){
      unsigned short iEdge_p1 = iEdge + 1;
      if(iEdge == nNode-1) iEdge_p1 = 0;
        intersect = Intersect2D(r0, Coord_elem[iEdge], Coord_elem[iEdge_p1], ppp0, ppp1);
        if(intersect == 1){
          if(count == 0){
            pp0[0] = ppp0[0];
            pp0[1] = ppp0[1];
          }
          else{
            pp1[0] = ppp0[0];
            pp1[1] = ppp0[1];
          }
        }
        else if(intersect == 2){
          pp0[0] = ppp0[0];
          pp0[1] = ppp0[1];
          pp1[0] = ppp1[0];
          pp1[1] = ppp1[1];
          inside = true;
          break;
        }
        count += intersect;
        if(count > 1){
          inside = true;
          break;
        }
    }

    if(Coord_elem != NULL){
      for(iNode = 0; iNode < nNode; iNode++){
        if(Coord_elem[iNode] != NULL){
          delete [] Coord_elem[iNode];
        }
      }
      delete [] Coord_elem;
    }
  }

  else{ // Check if line intersects face

    su2double **Coord_face;
    count = 0;
    unsigned short nNodeFace = 0; 
    for(unsigned short iFace = 0; iFace < geometry->elem[jElem]->GetnFaces(); iFace++){
      nNodeFace = geometry->elem[jElem]->GetnNodesFace(iFace);
      Coord_face = new su2double*[nNodeFace];
      for(iNode = 0; iNode < nNodeFace; iNode++){
        Coord_face[iNode] = new su2double[3];
        iPoint = geometry->elem[jElem]->GetNode(geometry->elem[jElem]->GetFaces(iFace, iNode));
        for(unsigned short iDim = 0; iDim < 3; iDim++){
          Coord_face[iNode][iDim] = geometry->node[iPoint]->GetCoord(iDim);
        }
      }
      intersect = Intersect3D(r0, phi, nNodeFace, Coord_face, ppp0);
      if(intersect == 1){
        if(count == 0){
          pp0[0] = ppp0[0];
          pp0[1] = ppp0[1];
          pp0[2] = ppp0[2];
        }
        else{
          pp1[0] = pp0[0];
          pp1[1] = pp0[1];
          pp1[2] = pp0[2];
        }
        count++;
        if(count > 1){
          inside = true;
          break;
        }
      }
    }

    if(Coord_face != NULL){
      for(iNode = 0; iNode < nNodeFace; iNode++){
        if(Coord_face[iNode] != NULL){
          delete [] Coord_face[iNode];
        }
      }
      delete [] Coord_face;
    }
      
  }
    
  delete [] ppp0;
  delete [] ppp1;

  return inside;
}

int CBoom_AugBurgers::Intersect2D(su2double r0, su2double *Coord_i, su2double *Coord_ip1, su2double *pp0, su2double *pp1){

  /*--- Interpolate if point between yi and yi+1 ---*/
  if((Coord_i[1] > -r0 && Coord_ip1[1] < -r0) || (Coord_i[1] < -r0 && Coord_ip1[1] > -r0)){
    su2double t = (-r0-Coord_i[1])/(Coord_ip1[1]-Coord_i[1]);
    pp0[0] = Coord_i[0] + t*(Coord_ip1[0]-Coord_i[0]);
    return 1;
  }
  /*--- Colinear segments at r0 ---*/
  else if(abs(Coord_i[1] + r0) < 1.0E-8 && abs(Coord_ip1[1] + r0) < 1.0E-8){
    pp0[0] = Coord_i[0];
    pp0[1] = -r0;
    pp1[0] = Coord_ip1[0];
    pp1[1] = -r0;
    return 2;
  }
  else{
    return 0;
  }

}

int CBoom_AugBurgers::Intersect3D(su2double r0, su2double phi, int nCoord, su2double **Coord_i, su2double *pp1){
  
  su2double y0 = r0*sin(phi), z0 = -r0*cos(phi);
  su2double ymin = 1.0E9, ymax = -1.0E9, zmin = 1.0E9, zmax = -1.0E9;

  /*--- First check simple bounding box ---*/
  for(int iCoord = 0; iCoord < nCoord; iCoord++){
    if(Coord_i[iCoord][1] < ymin) ymin = Coord_i[iCoord][1];
    if(Coord_i[iCoord][1] > ymax) ymax = Coord_i[iCoord][1];
    if(Coord_i[iCoord][2] < zmin) zmin = Coord_i[iCoord][2];
    if(Coord_i[iCoord][2] > zmax) zmax = Coord_i[iCoord][2];
  }

  if(y0 < ymin || y0 > ymax || z0 < zmin || z0 > zmax){
    return 0;
  }

  /*--- If inside bounding box, check sum of angles ---*/
  su2double d0, d1, d2;
  su2double a_x, a_y, b_x, b_y;
  su2double c;
  su2double deg = 0.0;
  bool cw;

  for(int iCoord = 0; iCoord < nCoord; iCoord++){
    int i = iCoord, ip = iCoord+1;
    if(i == nCoord-1) ip = 0;
    /*--- Vector magnitudes ---*/
    d0 = sqrt((Coord_i[i][1]-Coord_i[ip][1])*(Coord_i[i][1]-Coord_i[ip][1]) + (Coord_i[i][2]-Coord_i[ip][2])*(Coord_i[i][2]-Coord_i[ip][2]));
    d1 = sqrt((y0-Coord_i[ip][1])*(y0-Coord_i[ip][1]) + (z0-Coord_i[ip][2])*(z0-Coord_i[ip][2]));
    d2 = sqrt((Coord_i[i][1]-y0)*(Coord_i[i][1]-y0) + (Coord_i[i][2]-z0)*(Coord_i[i][2]-z0));
    /*--- Vector directions ---*/
    a_x = Coord_i[i][1] - y0;
    a_y = Coord_i[i][2] - z0;
    b_x = Coord_i[ip][1] - y0;
    b_y = Coord_i[ip][2] - z0;
    /*--- Clockwise or counterclockwise ---*/
    c = b_y*a_x - b_x*a_y;
    cw = (c < 0);
    deg += acos((d1*d1+d2*d2-d0*d0)/(2.0*d1*d2))*180./M_PI;
  }

  if(abs(deg - 360.) <= 3.){
    /*--- Get info needed for isoparameter computation ---*/
    su2double *Coord = new su2double[2];
    Coord[0] = y0; Coord[1] = z0;
    su2double *X_donor = new su2double[2*nCoord];
    for(int iCoord = 0; iCoord < nCoord; iCoord++){
      for(int iDim = 0; iDim < 2; iDim++){  
        X_donor[iDim*nCoord + iCoord] = Coord_i[iCoord][iDim+1];
      }
    }

    /*--- Compute isoparameters ---*/
    su2double *isoparams = new su2double[nCoord];
    Isoparameters(2, nCoord, X_donor, Coord, isoparams);

    /*--- Interpolate x-coord ---*/
    pp1[0] = 0.0;
    pp1[1] = y0;
    pp1[2] = z0;
    for(int iCoord = 0; iCoord < nCoord; iCoord++){
      pp1[0] += isoparams[iCoord]*Coord_i[iCoord][0];
    }

    if(isoparams != NULL) delete [] isoparams;
    if(Coord != NULL) delete [] Coord;
    if(X_donor != NULL) delete [] X_donor;
      
    return 1;
  }
  else{
    return 0;
  }

}

void CBoom_AugBurgers::AtmosISA(su2double& h0, su2double& T, su2double& a, su2double& p,
                      su2double& rho, su2double& g){
  /*---Calculate temperature, speed of sound, pressure, and density at a given
       altitude using standard atmosphere---*/
  const su2double GMR = 34.163195;    // hydrostatic constant
  const su2double R = 287.058;    // gas constant of air [m^2/s^2*K]
  const su2double Re_earth = 6369.;    // equatorial radius [km]

  su2double h = (h0/1000.)*Re_earth/(h0/1000.+Re_earth);    // geometric to geopotential altitude

  su2double htab[8] = {0.0, 11.0, 20.0, 32.0, 47.0, 51.0, 71.0, 84.852};
  su2double ttab[8] = {288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.946};
  su2double ptab[8] = {1.0, 2.233611E-1, 5.403295E-2, 8.5666784E-3, 1.0945601E-3,
                    6.6063531E-4, 3.9046834E-5, 3.68501E-6};
  su2double gtab[8] = {-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0, 0.0};
  su2double tgrad, tbase, tlocal, deltah, theta, pdelta, sigma;

  int i = 1, j=8, k;

  // Binary search
  while(j > i+1){
    k = (i+j)/2;
    if(h < htab[k-1]){
      j = k;
    }
    else{
      i = k;
    }
  }

  tgrad = gtab[i-1];
  tbase = ttab[i-1];
  deltah = h - htab[i-1];
  tlocal = tbase + tgrad*deltah;
  theta = tlocal/ttab[0];    // temperature ratio wrt sea-level

  if(abs(tgrad) < 1.0E-6){
    pdelta = ptab[i-1]*exp(-GMR*deltah/tbase);    // pressure ratio wrt sea-level
  }
  else{
    pdelta = ptab[i-1]*pow((tbase/tlocal), GMR/tgrad);
  }

  sigma = pdelta/theta;    // density ratio wrt sea-level

  T = theta*288.15;
  a = sqrt(g*R*T);
  p = pdelta*101325.;
  rho = sigma*1.225;

}

void CBoom_AugBurgers::HumidityISO(su2double& z0, su2double& p_inf, su2double& T_inf, su2double& h){

  su2double T01 = 273.16,  // Reference temperature (Humidity)
            Pr  = 101325.; // Reference pressure (Absorption)
  su2double logPsat_Pr = -6.8346*pow(T01/T_inf,1.261) + 4.6151;  // ISO standard for saturation vapor pressure

  su2double hr_prof[2][61] = {{0, 304.8, 609.6, 914.4, 1219.2,  1524,  1828.8,  2133.6,  2438.4,  2743.2,  3048,
                               3352.8,  3657.6,  3962.4,  4267.2,  4572,  4876.8,  5181.6,  5486.4,  5791.2,  6096,  
                               6400.8,  6705.6,  7010.4,  7315.2,  7620,  7924.8,  8229.6,  8534.4,  8839.2,  9144,  
                               9448.8,  9753.6,  10058.4, 10363.2, 10668, 10972.8, 11277.6, 11582.4, 11887.2, 12192, 
                               12496.8, 12801.6, 13106.4, 13411.2, 13716, 14020.8, 14325.6, 14630.4, 14935.2, 15240, 
                               15544.8, 16154.4, 16459.2, 16764, 17068.8, 17373.6, 17678.4, 17983.2, 18288, 18592.8},
                              {59.62, 60.48, 62.03, 63.83, 65.57, 67.06, 68.2,  68.97, 69.41, 69.62, 69.72, 69.83, 
                               70.05, 70.46, 71.12, 72.04, 73.19, 74.48, 75.77, 76.9,  77.61, 77.66, 76.77, 74.75, 
                               71.47, 66.96, 61.38, 55.07, 48.44, 41.95, 36.01, 30.95, 27.01, 24.38, 23.31, 24.29, 
                               28.6,  25.61, 22.1,  19.04, 16.42, 14.2,  12.34, 10.8,  9.53,  8.49,  7.64,  6.95,  
                               6.4, 5.97,  5.63,  5.38,  5.08,  5.02,  5.01,  5.03,  5.09,  5.18,  5.28,  5.4, 5.5}};
  su2double hr = 5.5;  // Interpolate relative humidity from profile (hr_prof is the default profile for sBOOM)
  if(z0 < 1.0E-8){
    hr = 59.62;
  }
  else{
    for(unsigned short i = 0; i < 60; i++){
      if(hr_prof[0][i] <= z0 && hr_prof[0][i+1] >= z0){
        hr = hr_prof[1][i] + (z0-hr_prof[0][i])*(hr_prof[1][i+1]-hr_prof[1][i])/(hr_prof[0][i+1]-hr_prof[0][i]);
        break;
      }
    }
  }

  h = hr * Pr/p_inf * pow(10.,logPsat_Pr);             // Absolute humidity (percent)

}

void CBoom_AugBurgers::SetKindSens(unsigned short kind_sensitivity){
  kind_sens = kind_sensitivity;
  int rank = SU2_MPI::GetRank();
  if(rank == MASTER_NODE){
    if(kind_sens == mesh_sens){      cout << endl; cout << "Computing mesh sensitivity." << endl;}
    else if(kind_sens == flow_sens){ cout << endl; cout << "Computing flow sensitivity." << endl;}
  }
}

void CBoom_AugBurgers::Run(CConfig *config){
    
  int rank = SU2_MPI::GetRank();
    
  if(rank == MASTER_NODE){
    for(unsigned short iPhi = 0; iPhi < ray_N_phi; iPhi++){
      cout << "Propagating signal for phi = " << ray_phi[iPhi] << ". " << endl;
      PropagateSignal(iPhi);
    }
  }
    
  if(rank == MASTER_NODE) cout << "Propagation complete." << endl;
    
  su2double Objective_Function = 0.0;
  if (rank == MASTER_NODE){
    for(unsigned short iPhi = 0; iPhi < ray_N_phi; iPhi++){
      Objective_Function += PLdB[iPhi]/su2double(ray_N_phi); // Normalize by number of propagated signals
    }
  }
    
  /*---Write boom strength metrics to file---*/
  if (rank == MASTER_NODE){
    ofstream sigFile;
    sigFile.precision(15);
    sigFile.open("boomSU2", ios::out);
    sigFile << "Objective_Function= " << Objective_Function << endl;
    sigFile.close();
  }
    
  /*---Extract sensitivities for discrete adjoint---*/
  if(AD_Mode){
    if (rank==MASTER_NODE){
      SU2_TYPE::SetDerivative(Objective_Function,1.0);
    }else{
      SU2_TYPE::SetDerivative(Objective_Function,0.0);
    }
    AD::StopRecording();
    AD::ComputeAdjoint();
        
    if (rank==MASTER_NODE) cout<<"Finished computing boom adjoint."<<endl;
        
    su2double extracted_derivative;
      
    if(kind_sens == mesh_sens){
      /*---Mesh sensitivities---*/
      for(unsigned short iPhi = 0; iPhi < ray_N_phi; iPhi++){
        for(unsigned int iSig=0; iSig<nPointID[iPhi]; iSig++){
          for(unsigned short i =0; i< nDim; i++){
            dJdX[iPhi][i][iSig]=SU2_TYPE::GetDerivative(extracted_derivative);
          }
        }
      }
    }
      
    if(kind_sens == flow_sens){
      /*---Flow sensitivities---*/
      for(unsigned short iPhi = 0; iPhi < ray_N_phi; iPhi++){
        for(unsigned int iSig=0; iSig<nPointID[iPhi]; iSig++){
          for(unsigned short i =0; i< nDim+3; i++){
            dJdU[iPhi][i][iSig]=SU2_TYPE::GetDerivative(extracted_derivative);
          }
        }
      }
    }
      
    AD::ClearAdjoints();
        
    if(rank==MASTER_NODE) cout<<"Finished extracting derivatives."<<endl;
        
    /*---Write sensitivities to file---*/
    WriteSensitivities();
    
  }
}

void CBoom_AugBurgers::PropagateSignal(unsigned short iPhi){

  unsigned long iIter = 0;
  ground_flag = false;

  while(!ground_flag){

  	Preprocessing(iPhi, iIter);
    Scaling(iPhi);
    Relaxation(iPhi, iIter);
    Attenuation(iPhi);
   	Nonlinearity(iPhi);

    ray_z -= dz;

    iIter++;
  }


  if(!AD_Mode){
    cout.width(5); cout << iIter;
    cout.width(12); cout.precision(6); cout << ray_z;
    cout.width(12); cout.precision(6); cout << p_peak << endl;
    cout << "Signal propagated in " << iIter << " iterations." << endl;
      
    WriteGroundPressure(iPhi);
  }

  if(Kind_Boom_Cost==BOOM_LOUD){
    PerceivedLoudness(iPhi);
  }
  else{
    AcousticEnergy(iPhi);
  }

}

void CBoom_AugBurgers::Preprocessing(unsigned short iPhi, unsigned long iIter){
  
  su2double R = 287.058,         // Gas constant
            mu0 = 1.846E-5,      // Reference viscosity
            kappa0 = 2.624E-2,   // Reference thermal conduction coefficient
            T0  = 300.,          // Reference temperature (Sutherland)
            Ts  = 110.4,         // Reference temperature (Sutherland)
            Ta  = 245.4,         // Reference temperature (Sutherland)
            Tb  = 27.6,          // Reference temperature (Sutherland)
            Tr  = 293.15,        // Reference temperature (Absorption)
            Pr  = 101325.;       // Reference pressure (Absorption)

  su2double h;

  /*---Preprocess signal for first iteration---*/
  if(iIter == 0){

    AtmosISA(flt_h, T_inf, c0, p_inf, rho0, atm_g);
    p0 = p_inf;
    flt_U = flt_M*c0;
    ray_z    = flt_h - ray_r0*cos(ray_phi[iPhi]);
    AtmosISA(ray_z, T_inf, c0, p_inf, rho0, atm_g);
    beta = (atm_g + 1.)/2.;
    HumidityISO(ray_z, p_inf, T_inf, h);
    c0 *= (1+0.0016*h); // Humidity correction to speed of sound

    /*---First create uniform grid since algorithm requires it---*/
    CreateUniformGridSignal(iPhi);
    ground_flag = false;

  	signal.P       = new su2double[signal.len[iPhi]];
  	signal.t       = new su2double[signal.len[iPhi]];
  	signal.t_prime = new su2double[signal.len[iPhi]];
  	signal.tau     = new su2double[signal.len[iPhi]];
    signal.taud    = new su2double[signal.len[iPhi]];
    p_peak = 0.;
    for(unsigned long i = 0; i < signal.len[iPhi]; i++){
      p_peak = max(p_peak, signal.p_prime[iPhi][i]);
    }

    M_a = p_peak/(rho0*pow(c0,2));
    if(!AD_Mode) cout << "Acoustic Mach number: " << M_a << endl;

  	for(unsigned long i = 0; i < signal.len[iPhi]; i++){
      signal.t[i]       = (signal.x[iPhi][i]-signal.x[iPhi][0])/flt_U;
  	}
    f0 = flt_U/(signal.x[iPhi][signal.len[iPhi]-1]-signal.x[iPhi][0]); // Characteristic time governed by length of signal
    w0 = 2.*M_PI*f0;

    for(unsigned long i = 0; i < signal.len[iPhi]; i++){
      signal.P[i] = signal.p_prime[iPhi][i]/p0;
      signal.tau[i] = w0*signal.t[i];
    }
    dt = signal.t[1] - signal.t[0];
    dtau = signal.tau[1] - signal.tau[0];

    CreateInitialRayTube(iPhi);

    /*---Initial signal---*/
    if(!AD_Mode){
      ofstream sigFile;
      char filename [64];
      SPRINTF (filename, "nearfield_%d.dat", SU2_TYPE::Int(iPhi));
      sigFile.precision(15);
      sigFile.open(filename, ios::out);
      for(int j = 0; j < signal.len[iPhi]; j++){
        sigFile << signal.tau[j]/w0 << "\t" << signal.P[j] << endl;
      }
      sigFile.close();

      cout << " Iter" << "        z[m]" << "   p_max[Pa]" << endl;
    }

  }

  else{
    AtmosISA(ray_z, T_inf, c0, p_inf, rho0, atm_g);
    HumidityISO(ray_z, p_inf, T_inf, h);
    c0 *= (1+0.0016*h); // Humidity correction to speed of sound
  }

  /*---Compute other coefficients needed for solution of ABE---*/
  xbar      = rho0*pow(c0,3)/(beta*w0*p0);
  mu        = mu0*pow(T_inf/T0,1.5)*(T0+Ts)/(T_inf+Ts);
  kappa     = kappa0*pow(T_inf/T0,1.5)*(T0+Ta*exp(-Tb/T0))/(T_inf+Ta*exp(-Tb/T_inf));
  delta     = mu/rho0*(4./3. + 0.6 + pow((atm_g-1.),2)*kappa/(atm_g*R*mu));
  alpha0_tv = delta*pow(w0,2)/(2.*pow(c0,3));
  Gamma     = 1./(alpha0_tv*xbar);

  su2double A_nu_O2 = 0.01278 * pow(T_inf/Tr,-2.5) * exp(-2239.1/T_inf);
  su2double A_nu_N2 = 0.1068  * pow(T_inf/Tr,-2.5) * exp(-3352./T_inf);

  su2double f_nu_O2 = p_inf/Pr * (24. + 4.04E4*h*(0.02+h)/(0.391+h));
  su2double f_nu_N2 = p_inf/Pr * sqrt(Tr/T_inf) * (9. + 280.*h*exp(-4.170*(pow(Tr/T_inf, 1./3.) - 1.)));

  theta_nu_O2 = f0/f_nu_O2;
  theta_nu_N2 = f0/f_nu_N2;

  C_nu_O2 = A_nu_O2 * f0 * xbar * theta_nu_O2;
  C_nu_N2 = A_nu_N2 * f0 * xbar * theta_nu_N2;

  /*--- Stop recording during step size computation to avoid implicit dependence of controller on state ---*/
  AD::StopRecording();

  /*--- Step size controller ---*/
  DetermineStepSize(iPhi, iIter);

  /*--- Continue recording ---*/
  AD::StartRecording();

  /*--- Output some information about the propagation ---*/
  if((!AD_Mode) && (iIter%50 == 0)){
    cout.width(5); cout << iIter;
    cout.width(12); cout.precision(6); cout << ray_z;
    cout.width(12); cout.precision(6); cout << p_peak << endl;
  }

}

void CBoom_AugBurgers::CreateUniformGridSignal(unsigned short iPhi){

  /*---Find start and end of significant part of signal---*/
  unsigned long istart = 0, iend = signal.len[iPhi]-1;
  for(unsigned long i = 0; i < signal.len[iPhi]; i++){
    if(abs(signal.p_prime[iPhi][i]) > 1.0E-3){
      istart = i;
      break;
    }
  }

  for(unsigned long i = signal.len[iPhi]-1; i >= istart+1; i--){
    if(abs(signal.p_prime[iPhi][i]) > 1.0E-3){
      iend = i;
      break;
    }
  }
  scale_L = signal.x[iPhi][iend] - signal.x[iPhi][istart];
  unsigned long len_new = iend + 1 - istart;

  /*---Average spacing and sampling frequency---*/
  dx_avg = scale_L/su2double(len_new);
  fs = flt_U/dx_avg;
  if(fs < 7250. && len_new < 1000){ // Make sure signal is adequately (at least somewhat) refined
    dx_avg = scale_L/1000.;
    fs = flt_U/dx_avg;
  }

  /*---Create new temp signal---*/
  unsigned long j = 0;
  len_new = ceil((signal.x[iPhi][signal.len[iPhi]-1]-signal.x[iPhi][0])/dx_avg);
  unsigned long len1 = len_new;
  su2double *xtmp = new su2double[len_new];
  su2double *ptmp = new su2double[len_new];
  su2double p_max = 0.;
  xtmp[0] = signal.x[iPhi][0];
  ptmp[0] = signal.p_prime[iPhi][0];
  unsigned long j0 = 0;
  for(unsigned long i = 1; i < len_new; i++){
    xtmp[i] = xtmp[i-1] + dx_avg;

    /*---Interpolate signal---*/
    if(xtmp[i] > signal.x[iPhi][signal.len[iPhi]-1]){
      j = signal.len[iPhi]-2;
      j0 = j;
      break;
    }
    else{
      for(unsigned long j = j0; j < signal.len[iPhi]-1; j++){
        if(xtmp[i] > signal.x[iPhi][j] && xtmp[i] < signal.x[iPhi][j+1]){
          j0 = j;
          break;
        }
      }
    }
    ptmp[i] = signal.p_prime[iPhi][j0] + (xtmp[i] - signal.x[iPhi][j0]) * (signal.p_prime[iPhi][j0+1]-signal.p_prime[iPhi][j0])/(signal.x[iPhi][j0+1]-signal.x[iPhi][j0]);
    p_max = max(p_max, abs(ptmp[i]));
    
  }

  /*---Store new signal---*/
  unsigned long len_pad = 5*ceil(scale_L/dx_avg);
  su2double dp_dx_end = abs(ptmp[len_new-1]-ptmp[len_new-2])/dx_avg;
  if(ptmp[len_new-1] > 0) dp_dx_end = - dp_dx_end;
  len_new = len_new+ceil(len_pad*1.75);
  signal.len[iPhi] = len_new;
  signal.x[iPhi] = new su2double[signal.len[iPhi]];
  signal.p_prime[iPhi] = new su2double[signal.len[iPhi]];
  unsigned long i0 = len_pad, i1 = i0+len1, i2 = signal.len[iPhi];
  /*---Zero-pad front of signal---*/
  for(unsigned long i = 0; i < i0; i++){
    signal.x[iPhi][i] = xtmp[0]-dx_avg*su2double(i0-i);
    signal.p_prime[iPhi][i] = 0.;
  }
  /*---Interpolated signal---*/
  for(unsigned long i = i0; i < i1; i++){
    signal.x[iPhi][i] = xtmp[i-i0];
    signal.p_prime[iPhi][i] = ptmp[i-i0];
  }
  /*---Recompress aft of signal---*/
  for(unsigned long i = i1; i < i2; i++){
    signal.x[iPhi][i] = signal.x[iPhi][i1-1]+dx_avg*su2double(i+1-i1);
    if((signal.p_prime[iPhi][i-1]+dp_dx_end*dx_avg) / signal.p_prime[iPhi][i-1] > 0.){ // If no sign change
      signal.p_prime[iPhi][i] = signal.p_prime[iPhi][i-1]+dp_dx_end*dx_avg;
    }
    else{
      signal.p_prime[iPhi][i] = 0.;
    }
  }

  if(!AD_Mode){
    cout << "Signal interpolated and padded, now contains " << signal.len[iPhi] << " points." << endl;
    cout << "Length scale of waveform = " << scale_L << " m." << endl;
    cout << "Sample frequency = " << fs << " Hz." << endl;
  }

  delete [] xtmp;
  delete [] ptmp;

}

void CBoom_AugBurgers::CreateInitialRayTube(unsigned short iPhi){

  ray_x = new su2double[4];
  ray_y = new su2double[4];
  ray_gamma = new su2double[2];
  ray_theta = new su2double[2];
  ray_c0 = new su2double[2];
  ray_nu = new su2double[2];
  ray_dt = 1.0E-1;
  ray_dphi = 1.0E-1;

  ray_theta[0] = asin(sin(flt_mu)*sin(flt_gamma) - cos(flt_mu)*cos(flt_gamma)*cos(ray_phi[iPhi]));
  ray_theta[1] = asin(sin(flt_mu)*sin(flt_gamma) - cos(flt_mu)*cos(flt_gamma)*cos(ray_phi[iPhi]+ray_dphi));

  ray_nu[0] = flt_psi - atan2(cos(flt_mu)*sin(ray_phi[iPhi]), sin(flt_mu)*cos(flt_gamma) + cos(flt_mu)*sin(flt_gamma)*cos(ray_phi[iPhi]));
  ray_nu[1] = flt_psi - atan2(cos(flt_mu)*sin(ray_phi[iPhi]+ray_dphi), sin(flt_mu)*cos(flt_gamma) + cos(flt_mu)*sin(flt_gamma)*cos(ray_phi[iPhi]+ray_dphi));

  ray_c0[0] = c0/cos(ray_theta[0]);  // TODO: Add wind contribution
  ray_c0[1] = c0/cos(ray_theta[1]);  // TODO: Add wind contribution
  ray_theta[0] = acos(c0/ray_c0[0]);
  ray_theta[1] = acos(c0/ray_c0[1]);

  su2double dx_dz, dy_dz, uwind = 0., vwind = 0.;
  dx_dz  = (c0*cos(ray_theta[0])*sin(ray_nu[0]) - uwind)/(c0*sin(ray_theta[0]));
  dy_dz  = (c0*cos(ray_theta[0])*cos(ray_nu[0]) - vwind)/(c0*sin(ray_theta[0]));
  ray_x[0] = ray_x[3] = dx_dz*ray_r0*cos(ray_phi[iPhi]);
  ray_y[0] = ray_y[1] = dy_dz*ray_r0*cos(ray_phi[iPhi]);
  dy_dz  = (c0*cos(ray_theta[1])*cos(ray_nu[1]) - vwind)/(c0*sin(ray_theta[1]));
  ray_x[1] = ray_x[2] = ray_x[0] + ray_dt*flt_M*c0;
  ray_y[2] = ray_y[3] = dy_dz*ray_r0*cos(ray_phi[iPhi]+ray_dphi);

  su2double u[2] = {ray_x[2]-ray_x[0], ray_y[2]-ray_y[0]};
  su2double v[2] = {ray_x[3]-ray_x[1], ray_y[3]-ray_y[1]};
  su2double c    = (u[0]*v[1] - u[1]*v[0]);
  su2double A_h   = 0.5*sqrt(pow(c,2));

  ray_A = c0*A_h*tan(ray_theta[0])/(ray_c0[0]);  // TODO: Add wind contribution
 
}

void CBoom_AugBurgers::DetermineStepSize(unsigned short iPhi, unsigned long iIter){

  su2double dsigma_non, dsigma_tv, dsigma_relO, dsigma_relN, dsigma_A, dsigma_rc, dsigma_c;
  su2double dx_dz, dy_dz, ds_dz, ds, z_new;
  su2double uwind = 0., vwind = 0.;

  if(Kind_Step == FIXED){
    dz = Step_size;
    dx_dz  = (c0*cos(ray_theta[0])*sin(ray_nu[0]) - uwind)/(c0*sin(ray_theta[0]));
    dy_dz  = (c0*cos(ray_theta[0])*cos(ray_nu[0]) - vwind)/(c0*sin(ray_theta[0]));
    ds_dz   = sqrt(pow(dx_dz,2)+pow(dy_dz,2)+1);
    ds = ds_dz*dz;
    dsigma = ds/xbar;
    /*---Get peak pressure for output---*/
    p_peak = 0.;
    for(unsigned long i = 0; i < signal.len[iPhi]-1; i++){
      p_peak = max(signal.P[i]*p0, p_peak);
    }
  }

  else if(Kind_Step == EXPONENTIAL){
    dz = Step_size*exp(iIter*Step_growth);
    dx_dz  = (c0*cos(ray_theta[0])*sin(ray_nu[0]) - uwind)/(c0*sin(ray_theta[0]));
    dy_dz  = (c0*cos(ray_theta[0])*cos(ray_nu[0]) - vwind)/(c0*sin(ray_theta[0]));
    ds_dz   = sqrt(pow(dx_dz,2)+pow(dy_dz,2)+1);
    ds = ds_dz*dz;
    dsigma = ds/xbar;
    /*---Get peak pressure for output---*/
    p_peak = 0.;
    for(unsigned long i = 0; i < signal.len[iPhi]-1; i++){
      p_peak = max(signal.P[i]*p0, p_peak);
    }
  }

  else{
    /*---Restrict dsigma to avoid multivalued waveforms---*/
    su2double dp, max_dp = 1.0E-9;
    p_peak = 0.;
    for(unsigned long i = 0; i < signal.len[iPhi]-1; i++){
      dp = signal.P[i+1]-signal.P[i];
      max_dp = max(dp, max_dp);
      p_peak = max(signal.P[i]*p0, p_peak);
    }

    dsigma_non = 0.5*pow(dtau,2.)/(max_dp*dtau + 2./Gamma); // dsigma < dtau^2/max|dp|*dtau + 2/Gamma

    /*---Restrict dsigma based on thermoviscous effects---*/
    dsigma_tv = 0.1*Gamma/signal.len[iPhi];

    /*---Restrict dsigma based on relaxation---*/
    dsigma_relO = 0.1*(1.+signal.len[iPhi]*pow(theta_nu_O2,2))/(C_nu_O2*signal.len[iPhi])*min(1.,2.*M_PI/(sqrt(signal.len[iPhi])*theta_nu_O2));
    dsigma_relN = 0.1*(1.+signal.len[iPhi]*pow(theta_nu_N2,2))/(C_nu_N2*signal.len[iPhi])*min(1.,2.*M_PI/(sqrt(signal.len[iPhi])*theta_nu_N2));

    /*---Restrict dsigma based on spreading---*/
    ds = xbar*1.0E-5;
    dx_dz  = (c0*cos(ray_theta[0])*sin(ray_nu[0]) - uwind)/(c0*sin(ray_theta[0]));
    dy_dz  = (c0*cos(ray_theta[0])*cos(ray_nu[0]) - vwind)/(c0*sin(ray_theta[0]));
    ds_dz   = sqrt(pow(dx_dz,2)+pow(dy_dz,2)+1);
    dz = ds/ds_dz;
    z_new = ray_z - dz;

    su2double x_new[4], y_new[4], A_new;

    x_new[0] = ray_x[0] + dx_dz*dz;
    y_new[0] = ray_y[0] + dy_dz*dz;
    x_new[1] = ray_x[1] + dx_dz*dz;
    y_new[1] = ray_y[1] + dy_dz*dz;

    dx_dz  = (c0*cos(ray_theta[1])*sin(ray_nu[1]) - uwind)/(c0*sin(ray_theta[1]));
    dy_dz  = (c0*cos(ray_theta[1])*cos(ray_nu[1]) - vwind)/(c0*sin(ray_theta[1]));

    x_new[2] = ray_x[2] + dx_dz*dz;
    y_new[2] = ray_y[2] + dy_dz*dz;
    x_new[3] = ray_x[3] + dx_dz*dz;
    y_new[3] = ray_y[3] + dy_dz*dz;

    su2double u[2] = {x_new[2]-x_new[0], y_new[2]-y_new[0]};
    su2double v[2] = {x_new[3]-x_new[1], y_new[3]-y_new[1]};
    su2double c    = (u[0]*v[1] - u[1]*v[0]);
    su2double A_h   = 0.5*sqrt(pow(c,2));

    A_new = c0*A_h*tan(ray_theta[0])/(ray_c0[0]);  // TODO: Add wind contribution
    su2double dA_dsigma = (A_new-ray_A)/(1.0E-5);

    dsigma_A = abs(0.05*ray_A/dA_dsigma);

    /*---Restrict dsigma based on stratification---*/
    su2double Tp1, cp1, pp1, rhop1, hp1;
    AtmosISA(z_new, Tp1, cp1, pp1, rhop1, atm_g);
    HumidityISO(z_new, pp1, Tp1, hp1);
    cp1 *= (1+0.0016*hp1);

    su2double drhoc_dsigma = (rhop1*cp1-rho0*c0)/(1.0E-5);
    dsigma_rc = abs(0.05*rho0*c0/drhoc_dsigma);

    su2double dc_dsigma = (cp1-c0)/(1.0E-5);
    dsigma_c = abs(0.05*c0/dc_dsigma);

    /*---Pick minimum dsigma---*/

    dsigma = dsigma_non;
    dsigma = min(dsigma, dsigma_tv);
    dsigma = min(dsigma, dsigma_relO);
    dsigma = min(dsigma, dsigma_relN);
    dsigma = min(dsigma, dsigma_A);
    dsigma = min(dsigma, dsigma_rc);
    dsigma = min(dsigma, dsigma_c);
    dsigma *= CFL_reduce;

    ds = xbar*dsigma;
    dx_dz  = (c0*cos(ray_theta[0])*sin(ray_nu[0]) - uwind)/(c0*sin(ray_theta[0]));
    dy_dz  = (c0*cos(ray_theta[0])*cos(ray_nu[0]) - vwind)/(c0*sin(ray_theta[0]));
    ds_dz   = sqrt(pow(dx_dz,2)+pow(dy_dz,2)+1);
    dz = ds/ds_dz;
  }

  /*---Check for intersection with ground plane---*/
  if(ray_z-dz < 0.0){
    dz = ray_z;
    ds = ds_dz*dz;
    dsigma = ds/xbar;
    ground_flag = true;
  }

}

void CBoom_AugBurgers::Nonlinearity(unsigned short iPhi){

  /*---Roe scheme for nonlinear term---*/
  su2double *Ptmp = new su2double[signal.len[iPhi]];
  su2double p_i, p_ip, p_im, f_i, f_ip, f_im, A_ip, A_im;
  /*---Predictor---*/
  for(unsigned long i = 1; i < signal.len[iPhi]-1; i++){

    p_i = signal.P[i];
    p_ip = signal.P[i+1];
    p_im = signal.P[i-1];

    f_i = -0.5*pow(p_i,2);
    f_ip = -0.5*pow(p_ip,2);
    f_im = -0.5*pow(p_im,2);

    if(abs(p_ip-p_i) < 1.0E-6) A_ip = -0.5*(p_ip + p_i);
    else                       A_ip = (f_ip - f_i)/(p_ip - p_i);
    if(abs(p_i-p_im) < 1.0E-6) A_im = -0.5*(p_im + p_i);
    else                       A_im = (f_i - f_im)/(p_i - p_im);
    
    f_ip = 0.5*(f_ip + f_i - abs(A_ip)*(p_ip - p_i));
    f_im = 0.5*(f_im + f_i - abs(A_im)*(p_i - p_im));

    Ptmp[i] = p_i - dsigma/(2.*dtau)*(f_ip - f_im);

  }
  Ptmp[0] = signal.P[0];
  Ptmp[signal.len[iPhi]-1] = signal.P[signal.len[iPhi]-1];

  /*---Corrector---*/
  su2double p_iml, p_imr, p_ipl, p_ipr, f_l, f_r;
  for(unsigned long i = 2; i < signal.len[iPhi]-2; i++){

    p_iml = Ptmp[i-1] + 0.5*(Ptmp[i-1]-Ptmp[i-2]);
    p_imr = Ptmp[i] - 0.5*(Ptmp[i]-Ptmp[i-1]);
    p_ipl = Ptmp[i] + 0.5*(Ptmp[i]-Ptmp[i-1]);
    p_ipr = Ptmp[i+1] - 0.5*(Ptmp[i+1]-Ptmp[i]);

    f_l = -0.5*pow(p_iml,2);
    f_r = -0.5*pow(p_imr,2);
    if(abs(p_imr - p_iml) < 1.0E-6) A_im = -0.5*(p_iml + p_imr);
    else                            A_im = (f_r - f_l)/(p_imr - p_iml);

    f_im = 0.5*(f_l + f_r - abs(A_im)*(p_imr - p_iml));
    f_l = -0.5*pow(p_ipl,2);
    f_r = -0.5*pow(p_ipr,2);
    if(abs(p_ipr - p_ipl) < 1.0E-6) A_ip = -0.5*(p_ipl + p_ipr);
    else                            A_ip = (f_r - f_l)/(p_ipr - p_ipl);
    f_ip = 0.5*(f_l + f_r - abs(A_ip)*(p_ipr - p_ipl));

    signal.P[i] = signal.P[i] - dsigma/dtau*(f_ip - f_im);

  }

  delete [] Ptmp;

}

void CBoom_AugBurgers::Attenuation(unsigned short iPhi){

  su2double lambda = dsigma/(2.*Gamma*pow(dtau,2)),
            alpha, alpha_p;
  su2double *y = new su2double[signal.len[iPhi]];

  alpha = 0.5;
  alpha_p = 1.-alpha;

  /*---Tridiagonal matrix-vector multiplication B_tv * Pk (see Cleveland thesis)---*/
  for(unsigned long i = 1; i < signal.len[iPhi]-1; i++){

    y[i] = alpha_p*lambda*(signal.P[i+1]+signal.P[i-1]) + (1.-2.*alpha_p*lambda)*signal.P[i];

  }

  /*---Solve for Pk+1 via Thomas algorithm for tridiagonal matrix---*/
  su2double *ci = new su2double[signal.len[iPhi]-1],
            *di = new su2double[signal.len[iPhi]-1];

  ci[0] = 0.;
  di[0] = signal.P[0];
  for(unsigned long i = 1; i < signal.len[iPhi]-1; i++){

    ci[i] = -alpha*lambda/((1.+2.*alpha*lambda) + alpha*lambda*ci[i-1]);
    di[i] = (y[i] + alpha*lambda*di[i-1])/((1.+2.*alpha*lambda) + alpha*lambda*ci[i-1]);

  }

  for(int i = signal.len[iPhi]-2; i >= 1; i--){

    signal.P[i] = di[i] - ci[i]*signal.P[i+1];
  }

  delete [] y;
  delete [] ci;
  delete [] di;
}

void CBoom_AugBurgers::Relaxation(unsigned short iPhi, unsigned long iIter){

  su2double lambda[2] = {dsigma*C_nu_O2/pow(dtau,2), dsigma*C_nu_N2/pow(dtau,2)},
            mur[2] = {theta_nu_O2/(2.*dtau), theta_nu_N2/(2.*dtau)},
            alpha, alpha_p, a, b, c;
  su2double *y = new su2double[signal.len[iPhi]], *ci, *di;

  alpha = 0.5;
  alpha_p = 1.-alpha;

  /*---Compute effect of different relaxation modes---*/
  for(unsigned short j = 0; j < 2; j++){
    /*---Tridiagonal matrix-vector multiplication B_nu * Pk (see Cleveland thesis)---*/
    for(unsigned long i = 1; i < signal.len[iPhi]-1; i++){

      y[i] = (alpha_p*lambda[j]-mur[j])*signal.P[i-1] + (1.-2.*alpha_p*lambda[j])*signal.P[i] + (alpha_p*lambda[j]+mur[j])*signal.P[i+1];

    }
    
    a = -(alpha*lambda[j] + mur[j]);
    b = 1. + 2.*alpha*lambda[j];
    c = -(alpha*lambda[j] - mur[j]);
    
    /*---Solve for Pk+1 via Thomas algorithm for tridiagonal matrix---*/
    ci  = new su2double[signal.len[iPhi]-1];
    di  = new su2double[signal.len[iPhi]-1];
    ci[0] = 0.;
    di[0] = signal.P[0];
    for(unsigned long i = 1; i < signal.len[iPhi]-1; i++){

      ci[i] = c/(b - a*ci[i-1]);
      di[i] = (y[i] - a*di[i-1])/(b - a*ci[i-1]);

    }

    for(int i = signal.len[iPhi]-2; i >= 1; i--){

      signal.P[i] = di[i] - ci[i]*signal.P[i+1];

    }

    delete [] ci;
    delete [] di;

  }

  delete [] y;
}

void CBoom_AugBurgers::Scaling(unsigned short iPhi){

  /*---Compute ray tube area at sigma+dsigma---*/
  su2double dx_dz, dy_dz, ds_dz, ds, z_new;
  su2double uwind = 0., vwind = 0.;
  ds = xbar*dsigma;
  dx_dz  = (c0*cos(ray_theta[0])*sin(ray_nu[0]) - uwind)/(c0*sin(ray_theta[0]));
  dy_dz  = (c0*cos(ray_theta[0])*cos(ray_nu[0]) - vwind)/(c0*sin(ray_theta[0]));
  ds_dz   = sqrt(pow(dx_dz,2)+pow(dy_dz,2)+1);
  dz = ds/ds_dz;
  z_new = ray_z - dz;

  su2double x_new[4], y_new[4], A_new;

  x_new[0] = ray_x[0] + dx_dz*dz;
  y_new[0] = ray_y[0] + dy_dz*dz;
  x_new[1] = ray_x[1] + dx_dz*dz;
  y_new[1] = ray_y[1] + dy_dz*dz;

  dx_dz  = (c0*cos(ray_theta[1])*sin(ray_nu[1]) - uwind)/(c0*sin(ray_theta[1]));
  dy_dz  = (c0*cos(ray_theta[1])*cos(ray_nu[1]) - vwind)/(c0*sin(ray_theta[1]));

  x_new[2] = ray_x[2] + dx_dz*dz;
  y_new[2] = ray_y[2] + dy_dz*dz;
  x_new[3] = ray_x[3] + dx_dz*dz;
  y_new[3] = ray_y[3] + dy_dz*dz;

  su2double u[2] = {x_new[2]-x_new[0], y_new[2]-y_new[0]};
  su2double v[2] = {x_new[3]-x_new[1], y_new[3]-y_new[1]};
  su2double c    = (u[0]*v[1] - u[1]*v[0]);
  su2double A_h   = 0.5*sqrt(pow(c,2));

  A_new = c0*A_h*tan(ray_theta[0])/(ray_c0[0]);  // TODO: Add wind contribution

  /*---Compute change in pressure from exact solution to general ray tube area eqn---*/
  /*---and due to atmospheric stratification---*/
  su2double Tp1, cp1, pp1, rhop1, hp1;
  AtmosISA(z_new, Tp1, cp1, pp1, rhop1, atm_g);
  HumidityISO(z_new, pp1, Tp1, hp1);
  cp1 *= (1+0.0016*hp1);

  for(unsigned long i = 0; i < signal.len[iPhi]; i++){

    signal.P[i] = sqrt((rhop1*cp1*ray_A)/(rho0*c0*A_new))*signal.P[i];

  }

  /*---Set new ray properties (except z, we'll do that later)---*/
  for(unsigned short i = 0; i < 4; i++){
    ray_x[i] = x_new[i];
    ray_y[i] = y_new[i];
  }
  ray_A = A_new;

  /*---Snell's law for wave normals---*/
  ray_theta[0] = acos(cp1*cos(ray_theta[0])/c0);
  ray_theta[1] = acos(cp1*cos(ray_theta[1])/c0);

}

void CBoom_AugBurgers::PerceivedLoudness(unsigned short iPhi){

  su2double p_ref = 20.E-6;             // [Pa]

  unsigned long n_sample, len_new; //n_sample = ceil(14500*signal.len[iPhi]*dx_avg/(flt_U)), // fmax*N/Fs
  unsigned short n_band = 41, N;

  su2double *w, *p_of_w, *p_of_t;

  su2double *wtmp, *ptmp;

  su2double fc[41]    = {1.25, 1.6, 2.0, 2.5, 3.15,
                         4., 5., 6.3, 8., 10.,
                         12.5, 16., 20., 25., 31.5,
                         40., 50., 63., 80., 100.,
                         125., 160., 200., 250., 315.,
                         400., 500., 630., 800., 1000.,
                         1250., 1600., 2000., 2500., 3150.,
                         4000., 5000., 6300., 8000., 10000.,
                         12500.}, 
            f_min[41] = {1.12, 1.41, 1.78, 2.24, 2.82,
                         3.55, 4.47, 5.62, 7.08, 8.91,
                         11.2, 14.1, 17.8, 22.4, 28.2,
                         35.5, 44.7, 56.2, 70.8, 89.1,
                         112., 141., 178., 224., 282.,
                         355., 447., 562., 708., 891.,
                         1120., 1410., 1780., 2240., 2820.,
                         3550., 4470., 5620., 7080., 8910.,
                         11200.}, 
            f_max[41] = {1.41, 1.78, 2.24, 2.82,
                         3.55, 4.47, 5.62, 7.08, 8.91,
                         11.2, 14.1, 17.8, 22.4, 28.2,
                         35.5, 44.7, 56.2, 70.8, 89.1,
                         112., 141., 178., 224., 282.,
                         355., 447., 562., 708., 891.,
                         1120., 1410., 1780., 2240., 2820.,
                         3550., 4470., 5620., 7080., 8910.,
                         11200., 14100.},
            E_band[41], SPL_band[41];

  /*--- Upsample if necessary ---*/
  if(fs < 29000.){ // Nyquist criterion for Mark VII
    if(!AD_Mode) cout << "Upsampling signal." << endl;

    /*--- Zero-fill signal ---*/
    N = ceil(29000./fs);
    len_new = N*(signal.len[iPhi]-1)+1;
    ptmp = new su2double[len_new];
    su2double t;
    for(unsigned long i = 0; i < signal.len[iPhi]-1; i++){
      ptmp[i*N] = signal.P[i];
      for(unsigned short j = 1; j < N; j++){
        t = signal.tau[i]+su2double(j)/su2double(N)*dtau;
        ptmp[i*N+j] = 0.;
        /*--- Sinc interpolation ---*/
        for(unsigned long k = 0; k < signal.len[iPhi]; k++){
          ptmp[i*N+j] += signal.P[k]*sin(M_PI*(t-signal.tau[k])*fs/w0)/(M_PI*(t-signal.tau[k])*fs/w0);
        }
      }
    }
    ptmp[len_new-1] = signal.P[signal.len[iPhi]-1];

    fs *= N;
    dtau /= N;

    /*---Write upsampled signal---*/
    if(!AD_Mode){
      ofstream sigFile;
      char filename [64];
      SPRINTF (filename, "ground_sinc_%d.dat", SU2_TYPE::Int(iPhi));
      sigFile.precision(15);
      sigFile.open(filename, ios::out);
      if(iPhi == 0) sigFile << "# phi, T, p" << endl;
      for(int j = 0; j < len_new; j++){
        sigFile << su2double(j)*dtau/w0 << "\t" << ptmp[j]*p0 << endl;
      }
      sigFile.close();
    }

  }
  else{
    len_new = signal.len[iPhi];
    ptmp = new su2double[len_new];
    for(unsigned long i = 0; i < signal.len[iPhi]; i++){
      ptmp[i] = signal.P[i];
    }
  }

  /*-- Zero-pad signal ---*/
  if(!AD_Mode) cout << "Zero-padding signal." << endl;
  m_pow_2  = ceil(log(su2double(len_new))/log(2.)); // Next power of 2 (for FFT)
  n_sample = pow(2,m_pow_2);
  w        = new su2double[n_sample/2];
  p_of_w   = new su2double[n_sample];
  p_of_t   = new su2double[n_sample];
  for(unsigned long i = 0; i < n_sample; i++){
    if(i < n_sample/2){w[i] = su2double(i)*fs/(su2double(n_sample));}
    p_of_w[i] = 0.;

    if(i < len_new){p_of_t[i] = ptmp[i];}
    else{p_of_t[i] = 0.;}
  }
  delete [] ptmp;

  /*--- Compute frequency domain signal ---*/
  if(!AD_Mode) cout << "Performing Fourier Transform." << endl;
  FFT(m_pow_2, p_of_t, p_of_w);

  if(!AD_Mode) cout << "Interpolating signal at band edges." << endl;
  n_sample = n_sample/2;
  wtmp = new su2double[n_sample+42];
  ptmp = new su2double[n_sample+42];
  unsigned short k = 0;
  for(unsigned short i = 0; i < n_sample+42; i++){
    if(i < 41){
      wtmp[i] = f_min[i];
      if(f_min[i] < w[0]){ // Interpolate between DC and f0
        ptmp[i] = p_of_w[0] + f_min[i] * (p_of_w[1]-p_of_w[0])/w[1];
      }
      else{
        while(f_min[i] > w[k]){
          k++;
        }
        ptmp[i] = p_of_w[k-1] + (f_min[i]-w[k-1]) * (p_of_w[k]-p_of_w[k-1])/(w[k]-w[k-1]);
      }
    }
    else if(i == 41){
      wtmp[i] = f_max[40];
      while(f_max[40] > w[k]){
        k++;
      }
      ptmp[i] = p_of_w[k-1] + (f_max[40]-w[k-1]) * (p_of_w[k]-p_of_w[k-1])/(w[k]-w[k-1]);
    }
    else{
      wtmp[i] = w[i-42];
      ptmp[i] = p_of_w[i-42];
    }
  }
  n_sample = n_sample + 42;
  MergeSort(wtmp, ptmp, 0, n_sample-1);

  /*--- Compute 1/3-oct band energies and pressure levels ---*/
  if(!AD_Mode) cout << "Computing 1/3-oct band pressure levels." << endl;
  if(!AD_Mode) cout << " Band" << "      E[Pa^2*s]" << endl;

  k = 0;
  for(unsigned short j = 1; j < n_sample; j++){
    if(wtmp[j] > f_min[0]){
      k = j-1;
      break;
    }
  }
  for(unsigned short j = 0; j < n_band; j++){
    E_band[j] = 0.0;
    while(wtmp[k] < f_max[j]){
      E_band[j] += 0.5*(ptmp[k]+ptmp[k+1])*(wtmp[k+1]-wtmp[k]);
      k++;
    }
      
    if(!AD_Mode){
      cout.width(5); cout << j+1;
      cout.width(15); cout.precision(6); cout << E_band[j] << endl;
    }
      
    E_band[j] /= 0.07;  // Critical time of human ear is 0.07 s (Shepherd, 1991)
    SPL_band[j] = 10.*log10(E_band[j]/pow(p_ref,2.)) - 3.0;
  }

  /*--- Use band pressure levels to compute perceived loudness ---*/
  if(!AD_Mode) cout << "Computing perceived loudness (MarkVII)." << endl;
  MarkVII(iPhi, SPL_band, fc, n_band);

  /*--- Clean up ---*/
  delete [] w;
  delete [] p_of_w;
  delete [] p_of_t;
  delete [] wtmp;
  delete [] ptmp;

}

void CBoom_AugBurgers::FFT(unsigned long m, su2double *x, su2double *y){

  /*---Source: http://paulbourke.net/miscellaneous/dft---*/

  long n,i,i1,j,k,i2,l,l1,l2;
  su2double c1,c2,tx,ty,t1,t2,u1,u2,z;

  /* Calculate the number of points */
  n = 1;
  for (i=0;i<m;i++){
    n *= 2;
  }

  /* Do the bit reversal */
  i2 = n >> 1;
  j = 0;
  for (i=0;i<n-1;i++) {
    if (i < j) {
       tx = x[i];
       ty = y[i];
       x[i] = x[j];
       y[i] = y[j];
       x[j] = tx;
       y[j] = ty;
    }
    k = i2;
    while (k <= j) {
       j -= k;
       k >>= 1;
    }
    j += k;
  }

  /* Compute the FFT */
  c1 = -1.0; 
  c2 = 0.0;
  l2 = 1;
  for (l=0;l<m;l++) {
    l1 = l2;
    l2 <<= 1;
    u1 = 1.0; 
    u2 = 0.0;
    for (j=0;j<l1;j++) {
       for (i=j;i<n;i+=l2) {
          i1 = i + l1;
          t1 = u1 * x[i1] - u2 * y[i1];
          t2 = u1 * y[i1] + u2 * x[i1];
          x[i1] = x[i] - t1; 
          y[i1] = y[i] - t2;
          x[i] += t1;
          y[i] += t2;
       }
       z =  u1 * c1 - u2 * c2;
       u2 = u1 * c2 + u2 * c1;
       u1 = z;
    }
    c2 = -sqrt((1.0 - c1) / 2.0);
    c1 = sqrt((1.0 + c1) / 2.0);
  }

  /* Scaling for forward transform */
  for (i=0;i<n/2+1;i++) {
    y[i] = 2.*(x[i]*x[i] + y[i]*y[i])*pow(p0*dtau/w0,2.);
  }
  y[0] /= 2;

}

void CBoom_AugBurgers::FourierTransform(unsigned short iPhi, su2double *w, su2double *p_of_w, su2double& p_dc, unsigned short n_sample){

  su2double t1, t2, y1, y2;
  su2double p_real, p_imag;
  unsigned long N = signal.len[iPhi];

  /*---Initialize frequency domain signal ---*/
  for(unsigned short i = 0; i < n_sample; i++){
    w[i] = su2double(i)*flt_U/(su2double(N)*dx_avg);  // f = n*Fs/N
  }

  /*--- Transform ---*/

  for(unsigned long j = 0; j < n_sample; j++){

    p_real = 0.;
    p_imag = 0.;

    for(unsigned long i = 0; i < N-1; i++){
      t1 = signal.tau[i]/w0; t2 = signal.tau[i+1]/w0;
      y1 = signal.P[i]*p0;   y2 = signal.P[i+1]*p0;

      p_real += 0.5*(y2*cos(2*M_PI*w[j]*t2) + y1*cos(2*M_PI*w[j]*t1))*(t2-t1);
      p_imag -= 0.5*(y2*sin(2*M_PI*w[j]*t2) + y1*sin(2*M_PI*w[j]*t1))*(t2-t1);

    }

    p_of_w[j] = 2.*(p_real*p_real + p_imag*p_imag);

  }

  /*--- DC pressure ---*/
  p_real = 0.;
  for(unsigned long i = 0; i < N-1; i++){
    t1 = signal.tau[i]/w0; t2 = signal.tau[i+1]/w0;
    y1 = signal.P[i]*p0;   y2 = signal.P[i+1]*p0;

    p_real += 0.5*(y2 + y1)*(t2-t1);
  }
  p_dc = p_real*p_real;

}

void CBoom_AugBurgers::MarkVII(unsigned short iPhi, su2double *SPL_band, su2double *fc, unsigned short n_band){

  short band;
  su2double L_eq[41], sonband;
  su2double A, B, llb, ulb, xb, sonmax, sonsum, F;

  /*--- Sone -> F factor table from Stevens, 1972 ---*/
  su2double ffactr[186] = {0.181, 0.100, 0.196, 0.122, 0.212, 0.140, 0.230, 0.158,
                           0.248, 0.174, 0.269, 0.187, 0.290, 0.200, 0.314, 0.212,
                           0.339, 0.222, 0.367, 0.232, 0.396, 0.241, 0.428, 0.250,
                           0.463, 0.259, 0.500, 0.267, 0.540, 0.274, 0.583, 0.281,
                           0.630, 0.287, 0.680, 0.293, 0.735, 0.298, 0.794, 0.303,
                           0.857, 0.308, 0.926, 0.312, 1.000, 0.316, 1.080, 0.319,
                           1.170, 0.320, 1.260, 0.322, 1.360, 0.322, 1.470, 0.320,
                           1.590, 0.319, 1.710, 0.317, 1.850, 0.314, 2.000, 0.311,
                           2.160, 0.308, 2.330, 0.304, 2.520, 0.300, 2.720, 0.296,
                           2.940, 0.292, 3.180, 0.288, 3.430, 0.284, 3.700, 0.279,
                           4.000, 0.275, 4.320, 0.270, 4.670, 0.266, 5.040, 0.262,
                           5.440, 0.258, 5.880, 0.253, 6.350, 0.248, 6.860, 0.244,
                           7.410, 0.240, 8.000, 0.235, 8.640, 0.230, 9.330, 0.226,
                           10.10, 0.222, 10.90, 0.217, 11.80, 0.212, 12.70, 0.208,
                           13.70, 0.204, 14.80, 0.200, 16.00, 0.197, 17.30, 0.195,
                           18.70, 0.194, 20.20, 0.193, 21.80, 0.192, 23.50, 0.191,
                           25.40, 0.190, 27.40, 0.190, 29.60, 0.190, 32.00, 0.190,
                           34.60, 0.190, 37.30, 0.190, 40.30, 0.191, 43.50, 0.191,
                           47.00, 0.192, 50.80, 0.193, 54.90, 0.194, 59.30, 0.195,
                           64.00, 0.197, 69.10, 0.199, 74.70, 0.201, 80.60, 0.203,
                           87.10, 0.205, 94.10, 0.208, 102.0, 0.210, 110.0, 0.212,
                           119.0, 0.215, 128.0, 0.217, 138.0, 0.219, 149.0, 0.221,
                           161.0, 0.223, 174.0, 0.224, 188.0, 0.225, 203.0, 0.226,
                           219.0, 0.227},

  /*--- Leq -> sone table from Jackson, 1973 ---*/
            son[140]    = {0.079, 0.087, 0.097, 0.107, 0.118,
                           0.129, 0.141, 0.153, 0.166, 0.181,
                           0.196, 0.212, 0.230, 0.248, 0.269,
                           0.290, 0.314, 0.339, 0.367, 0.396,
                           0.428, 0.463, 0.500, 0.540, 0.583,
                           0.630, 0.680, 0.735, 0.794, 0.857,
                           0.926, 1.000, 1.080, 1.170, 1.260,
                           1.360, 1.470, 1.590, 1.710, 1.850,
                           2.000, 2.160, 2.330, 2.520, 2.720,
                           2.940, 3.180, 3.430, 3.700, 4.000,
                           4.320, 4.670, 5.040, 5.440, 5.880,
                           6.350, 6.860, 7.410, 8.000, 8.640,
                           9.330, 10.10, 10.90, 11.80, 12.70,
                           13.70, 14.80, 16.00, 17.30, 18.70,
                           20.20, 21.80, 23.50, 25.40, 27.40,
                           29.60, 32.00, 34.60, 37.30, 40.30,
                           43.50, 47.00, 50.80, 54.90, 59.30,
                           64.00, 69.10, 74.70, 80.60, 87.10,
                           94.10, 102.0, 110.0, 119.0, 128.0,
                           138.0, 149.0, 161.0, 174.0, 188.0,
                           203.0, 219.0, 237.0, 256.0, 276.0,
                           299.0, 323.0, 348.0, 376.0, 406.0,
                           439.0, 474.0, 512.0, 553.0, 597.0,
                           645.0, 697.0, 752.0, 813.0, 878.0,
                           948.0, 1024., 1106., 1194., 1290.,
                           1393., 1505., 1625., 1756., 1896.,
                           2048., 2212., 2389., 2580., 2787.,
                           3010., 3251., 3511., 3792., 4096.},

  /*--- Band limits and coefficients from Jackson, 1973 ---*/
            ll[8]       = {86.5,  85.0,  83.5,  82.0,  80.5,  79.0,  77.5,  76.0},
            ul[8]       = {131.5, 130.0, 128.5, 127.0, 125.5, 124.0, 122.5, 121.0},
            X[8]        = {10.5,  9.0,   7.5,   6.0,   4.5,   3.0,   1.5,   0.0};

  for(unsigned short i = 0; i < n_band; i++){
    band = i+1;
    if(band < 27){
      /*--- Get band limits ---*/
      if(band < 20){
        llb = ll[0];
        ulb = ul[0];
        xb = X[0];
        B = 160. - (160. - SPL_band[i])*log10(80.)/log10(fc[i]);
        SPL_band[i] = B;
        fc[i] = 80.0;
      }
      else{
        llb = ll[i-18];
        ulb = ul[i-18];
        xb = X[i-18];
      }

      /*--- Compute equivalent loudness ---*/
      if(SPL_band[i] <= llb){
        A = 115. - (115. - SPL_band[i])*log10(400.)/log10(fc[i]);
      }
      else if(SPL_band[i] > ulb){
        A = 160. - (160. - SPL_band[i])*log10(400.)/log10(fc[i]);
      }
      else{
        A = SPL_band[i] - xb;
      }
      L_eq[i] = A - 8.0;
    }

    /*--- Compute equivalent loudness ---*/
    else if(band >= 27 && band <= 31) L_eq[i] = SPL_band[i] - 8.0;
    else if(band >= 32 && band <= 34) L_eq[i] = SPL_band[i] - 2.*(35. - su2double(band));
    else if(band >= 35 && band <= 39) L_eq[i] = SPL_band[i];
    else                              L_eq[i] = SPL_band[i] + 4.*(39. - su2double(band));
  }

  /*--- Interpolate sone based on Leq ---*/
  sonmax = 0.0;
  sonsum = 0.0;
  for(unsigned short i = 0; i < n_band; i++){
    band = floor(L_eq[i]);
    if(band > 139) sonband = son[139];
    else if(band < 1) sonband = 0;
    else sonband = son[band-1] + (L_eq[i] - su2double(band))*(son[band] - son[band-1]); // Leq in sone table increment by 1, so no division
    
    sonsum += sonband;
    sonmax = max(sonmax, sonband);
  } 

  /*--- Interpolate F factor at max(sone) ---*/
  F = 0.;
  if(sonmax >= 219.) F = 0.227;
  else if(sonmax < ffactr[0]) F = 0.;
  else{
    for(unsigned short i = 0; i < 92; i++){
      if(sonmax >= ffactr[2*i] && sonmax <= ffactr[2*i+2]){
        F = ffactr[2*i+1] + (sonmax - ffactr[2*i])*(ffactr[2*i+3] - ffactr[2*i+1])/(ffactr[2*i+2] - ffactr[2*i]);
        break;
      }
    }
  }

  /*-- Calculate loudness ---*/
  PLdB[iPhi] = 32. + 9.*log(sonmax + F*(sonsum - sonmax))/log(2.);

}

void CBoom_AugBurgers::AcousticEnergy(unsigned short iPhi){

  if(!AD_Mode) cout << "Computing acoustic energy." << endl;

  PLdB[iPhi] = 0.;
  for(int j = 1; j < signal.len[iPhi]; j++){
    if(signal.P[j]*signal.P[j-1] < 0.0){ // if sign change, do double triangular integration
      /*--- Find root of line segment ---*/
      su2double tau0 = signal.tau[j-1] + (-signal.P[j-1])
                    *(signal.tau[j]-signal.tau[j-1])/(signal.P[j]-signal.P[j-1]);
      PLdB[iPhi] += (0.5*(signal.P[j-1]*signal.P[j-1])*(tau0-signal.tau[j-1])
                      + 0.5*(signal.P[j]*signal.P[j])*(signal.tau[j]-tau0))*p0*p0/w0;
    }
    else{ // otherwise, do trapezoidal integration
      PLdB[iPhi] += 0.5*(signal.P[j]*signal.P[j]+signal.P[j-1]*signal.P[j-1])
                    *(signal.tau[j]-signal.tau[j-1])*p0*p0/w0;
    }
  }

}

void CBoom_AugBurgers::WriteGroundPressure(unsigned short iPhi){
  /*---Final signal file---*/
  ofstream sigFile;
  char filename [64];
  SPRINTF (filename, "ground_%d.dat", SU2_TYPE::Int(iPhi));
  sigFile.precision(15);
  sigFile.open(filename, ios::out);
  if(iPhi == 0) sigFile << "# phi, T, p" << endl;
  for(int j = 0; j < signal.len[iPhi]; j++){
    sigFile << signal.tau[j]/w0 << "\t" << signal.P[j]*p0 << endl;
  }
  sigFile.close();
}

void CBoom_AugBurgers::WriteSensitivities(){
  unsigned long iVar, iSig, Max_nPointID, nVar, Global_Index, Total_Index;
  ofstream Boom_AdjointFile;

  int rank = 0, iProcessor, nProcessor = 1;
  rank = SU2_MPI::GetRank();
  nProcessor = SU2_MPI::GetSize();
    
  unsigned long Buffer_Send_nPointID[1], *Buffer_Recv_nPointID = NULL;

  if(rank == MASTER_NODE) Buffer_Recv_nPointID = new unsigned long [nProcessor];
    
  /*--- Determine number of variables based on sensitivity mode ---*/
  if(kind_sens == mesh_sens) nVar = nDim;
  if(kind_sens == flow_sens) nVar = nDim+3;

  /*--- Open proper file ---*/
  if(rank == MASTER_NODE){
    Boom_AdjointFile.precision(15);
    if(kind_sens == mesh_sens)      Boom_AdjointFile.open("Adj_Boom_dJdX.dat", ios::out);
    else if(kind_sens == flow_sens) Boom_AdjointFile.open("Adj_Boom_dJdU.dat", ios::out);
  }

  for(unsigned short iPhi = 0; iPhi < ray_N_phi; iPhi++){

    Buffer_Send_nPointID[0] = nPointID[iPhi]; 

    SU2_MPI::Gather(&Buffer_Send_nPointID, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nPointID, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD); // send the number of vertices at each process to the master
    SU2_MPI::Allreduce(&nPointID[iPhi], &Max_nPointID, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD); // find the max num of vertices over all processes

    /*--- Pack sensitivity values in each processor and send to root ---*/
    su2double *Buffer_Send_dJdU = new su2double [Max_nPointID*nVar];
    unsigned long *Buffer_Send_GlobalIndex = new unsigned long [Max_nPointID];
    /*--- Zero send buffers ---*/
    for(unsigned int i =  0; i < Max_nPointID*nVar; i++){
      Buffer_Send_dJdU[i] = 0.0;
    }
    for(unsigned int i = 0; i < Max_nPointID; i++){
      Buffer_Send_GlobalIndex[i] = 0;
    }
    su2double *Buffer_Recv_dJdU = NULL;
    unsigned long *Buffer_Recv_GlobalIndex = NULL;

    if(rank == MASTER_NODE) {
      Buffer_Recv_dJdU = new su2double [nProcessor*Max_nPointID*nVar];
      Buffer_Recv_GlobalIndex = new unsigned long[nProcessor*Max_nPointID];
    }

    /*--- Fill send buffers with dJ/dX or dJ/dU ---*/
    for(iVar = 0; iVar < nVar; iVar++){
      for(iSig = 0; iSig < nPointID[iPhi]; iSig++){
        if(kind_sens == mesh_sens)      Buffer_Send_dJdU[iVar*nPointID[iPhi]+iSig] = dJdX[iPhi][iVar][iSig];
        else if(kind_sens == flow_sens) Buffer_Send_dJdU[iVar*nPointID[iPhi]+iSig] = dJdU[iPhi][iVar][iSig];
      }
    }

    for(iSig = 0; iSig < nPointID[iPhi]; iSig++){
       Buffer_Send_GlobalIndex[iSig] = PointID[iPhi][iSig];
    }

    SU2_MPI::Gather(Buffer_Send_dJdU, Max_nPointID*nVar, MPI_DOUBLE, Buffer_Recv_dJdU,  Max_nPointID*nVar , MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Gather(Buffer_Send_GlobalIndex, Max_nPointID, MPI_UNSIGNED_LONG, Buffer_Recv_GlobalIndex, Max_nPointID , MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);

    if(rank == MASTER_NODE){

      /*--- Loop through all of the collected data and write each node's values ---*/
      for(iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        for(iSig = 0; iSig < Buffer_Recv_nPointID[iProcessor]; iSig++) {
          Global_Index = Buffer_Recv_GlobalIndex[iProcessor*Max_nPointID+iSig];
          Boom_AdjointFile  << scientific << Global_Index << "\t";

          for(iVar = 0; iVar < nVar; iVar++){
            /*--- Current index position and global index ---*/
            Total_Index  = iProcessor*Max_nPointID*nVar + iVar*Buffer_Recv_nPointID[iProcessor]  + iSig;

            /*--- Write to file---*/
            Boom_AdjointFile << scientific <<  Buffer_Recv_dJdU[Total_Index]   << "\t";
          }
          Boom_AdjointFile  << endl;

        }
      }

      delete [] Buffer_Recv_nPointID;
      delete [] Buffer_Recv_dJdU;
      delete [] Buffer_Recv_GlobalIndex;
    }

    delete [] Buffer_Send_dJdU;
    delete [] Buffer_Send_GlobalIndex;

  }

  if (rank == MASTER_NODE) cout << "\nFinished writing boom adjoint file." << endl;

}

