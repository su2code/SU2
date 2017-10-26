#include "../include/SUBoom.hpp"
#include <cmath>
#include <algorithm>
#include <complex>
#include <valarray>
#include <iostream>
#include <fstream>

typedef std::complex<su2double> Complex;
typedef std::valarray<Complex> CArray;

SUBoom::SUBoom(){

}

SUBoom::SUBoom(CSolver *solver, CConfig *config, CGeometry *geometry){

  int rank, nProcessor = 1;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
#endif

  /*---Make sure to read in hard-coded values in the future!---*/

  nDim = geometry->GetnDim();

  /*---Flight variables---*/
  flt_h = 15240; // altitude [m]
  flt_M = config->GetMach();
  flt_psi = 0.;  // heading angle [deg]
  flt_gamma = 0.; // flight path angle [deg]

  /*---Atmosphere variables---*/
  atm_g = config->GetGamma();
  atm_noise_flag = 0;

  /*---Scale factors---*/
  scale_L = config->GetRefLength();

  /*---Values from config file---*/
  const su2double deg2rad = M_PI/180.;
  n_prof = config->GetBoom_N_prof();
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
  tol_dphi = config->GetBoom_Tol_dphi();
  tol_dr = config->GetBoom_Tol_dr();
  tol_m = config->GetBoom_Tol_m();
  tol_dp = config->GetBoom_Tol_dp();
  tol_l = config->GetBoom_Tol_l();

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
  
  if(rank == MASTER_NODE)
    cout << "Pressure_Ref = " << Pressure_Ref << ", Pressure_FreeStream = " << Pressure_FreeStream << endl;

  /*---Perform search on domain to determine where line intersects boundary---*/
  if(rank == MASTER_NODE)
    cout << "Search for start of line." << endl;
  SearchLinear(config, geometry, ray_r0, ray_phi, ray_N_phi);

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
  signal.original_len = new unsigned long[ray_N_phi];
  signal.final_len = new unsigned long[ray_N_phi];
  signal.x = new su2double*[ray_N_phi];
  signal.original_p = new su2double*[ray_N_phi];
  signal.original_T = new su2double*[ray_N_phi];
  signal.final_p = new su2double*[ray_N_phi];
  signal.final_T = new su2double*[ray_N_phi];

  /*---Interpolate pressures along line---*/
  if(rank == MASTER_NODE)
    cout << "Extract pressure signature." << endl;
  nPointID_loc = 0;
  for(unsigned short iPhi = 0; iPhi < ray_N_phi; iPhi++){
    nPointID[iPhi] = 0;
    if(nPanel[iPhi] > 0){      
      ExtractPressure(solver, config, geometry, iPhi);
    }
    else{
      PointID[iPhi] = new unsigned long[1];
      PointID[iPhi][0] = -1;
    }
  }

  nPointID_proc = new unsigned long[nProcessor];
  unsigned long iPanel, panelCount, totSig, maxSig;
  unsigned long *nPanel_loc = new unsigned long[nProcessor];
  nPointID_tot = 0;
  for(unsigned short iPhi = 0; iPhi < ray_N_phi; iPhi++){
    for(iPanel = 0; iPanel < nPanel[iPhi]; iPanel++){
      signal.original_p[iPhi][iPanel] = signal.original_p[iPhi][iPanel]*Pressure_Ref - Pressure_FreeStream;
    }

    totSig = 0;
    maxSig = 0;
#ifdef HAVE_MPI
    SU2_MPI::Allreduce(&nPanel[iPhi], &totSig, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&nPanel[iPhi], &maxSig, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
    SU2_MPI::Gather(&nPanel[iPhi], 1, MPI_UNSIGNED_LONG, nPanel_loc, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Gather(&nPointID[iPhi], 1, MPI_UNSIGNED_LONG, nPointID_proc, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
#endif

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
      Buffer_Send_Press[iPanel] = signal.original_p[iPhi][iPanel];
    }

#ifdef HAVE_MPI
    SU2_MPI::Gather(Buffer_Send_Press, maxSig, MPI_DOUBLE, Buffer_Recv_Press,  maxSig , MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Gather(Buffer_Send_x, maxSig, MPI_DOUBLE, Buffer_Recv_x,  maxSig , MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#endif

    if (rank == MASTER_NODE)
      cout << "Gathered signal data to MASTER_NODE." << endl;

    if (rank == MASTER_NODE){
      panelCount = 0;
      nPanel[iPhi] = totSig;
      signal.x[iPhi] = new su2double[nPanel[iPhi]];
      signal.original_p[iPhi] = new su2double[nPanel[iPhi]];
      signal.original_T[iPhi] = new su2double[nPanel[iPhi]];
      for(unsigned int iProcessor = 0; iProcessor < nProcessor; iProcessor++){
        for(iPanel = 0; iPanel < nPanel_loc[iProcessor]; iPanel++){
          signal.x[iPhi][panelCount] = Buffer_Recv_x[iProcessor*maxSig+iPanel];
          signal.original_p[iPhi][panelCount] = Buffer_Recv_Press[iProcessor*maxSig+iPanel];
          panelCount++;
        }
        nPointID_tot += nPointID_proc[iProcessor];
      }

      /*---Sort signal in order of x-coordinate---*/
      cout << "Sorting signal data. " << nPanel[iPhi] << " points to sort." << endl;
      MergeSort(signal.x[iPhi], signal.original_p[iPhi], 0, totSig-1);

      /*---Check for duplicate points---*/
      signal.original_len[iPhi] = nPanel[iPhi];
      for(iPanel = 1; iPanel < nPanel[iPhi]; iPanel++){
        if(abs(signal.x[iPhi][iPanel-1]-signal.x[iPhi][iPanel]) < 1.0E-12){
          for(unsigned long jPanel = iPanel; jPanel < nPanel[iPhi]; jPanel++){
            signal.x[iPhi][jPanel-1] = signal.x[iPhi][jPanel];
            signal.original_p[iPhi][jPanel-1] = signal.original_p[iPhi][jPanel];
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
          ptmp[iPanel] = signal.original_p[iPhi][iPanel];
        }
        signal.original_len[iPhi] = nPanel[iPhi];
        signal.x[iPhi] = new su2double[nPanel[iPhi]];
        signal.original_T[iPhi] = new su2double[nPanel[iPhi]];
        signal.original_p[iPhi] = new su2double[nPanel[iPhi]];
        for(iPanel = 0; iPanel < nPanel[iPhi]; iPanel++){
          signal.x[iPhi][iPanel] = xtmp[iPanel];
          signal.original_p[iPhi][iPanel] = ptmp[iPanel];
        }
        delete [] xtmp;
        delete [] ptmp;
        totSig = nPanel[iPhi];
      }
    }
  }

  /*---Now write to file---*/
  if(rank == MASTER_NODE){
    nPanel_tot = 0;
    ofstream sigFile;
    sigFile.precision(15);
    sigFile.open("signal_original.dat");
    sigFile << "# phi, x, p" << endl;
    for(unsigned short iPhi = 0; iPhi < ray_N_phi; iPhi++){
      for(iPanel = 0; iPanel < nPanel[iPhi]; iPanel++){
        sigFile << scientific << ray_phi[iPhi] << "\t";
        sigFile << scientific << signal.x[iPhi][iPanel] << "\t";
        sigFile << scientific << signal.original_p[iPhi][iPanel]   << "\t";
        sigFile << endl;
      }
      nPanel_tot += nPanel[iPhi];
    }
    sigFile.close();
    cout << "Signal written. Total nPanel = " << nPanel_tot << "." << endl;

  }

  /*---Initialize sensitivities---*/
  if(config->GetAD_Mode()){
    dJdU = new su2double**[ray_N_phi];
    for(unsigned short iPhi = 0; iPhi < ray_N_phi; iPhi++){
      dJdU[iPhi] = new su2double* [nDim+3];
      for(int iDim = 0; iDim < nDim+3 ; iDim++){
        if(nPointID[iPhi] > 0){
          dJdU[iPhi][iDim] = new su2double[nPointID[iPhi]];
          for(iPanel = 0;  iPanel< nPointID[iPhi]; iPanel++){
            dJdU[iPhi][iDim][iPanel] = 0.0;
          }
        }
        else{
          dJdU[iPhi][iDim] = new su2double[1];
          dJdU[iPhi][iDim][0] = 0.0;
        }
      }
    }

    if (rank==MASTER_NODE)
      cout << "Sensitivities initialized." << endl;
  }

}

SUBoom::~SUBoom(void){

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

void QuickSort(su2double x[], su2double p[], int l, int r){
  int i = l, j = r;
  su2double tmp, pivot=x[(l+r)/2];

  while(i <= j){
    while(x[i] < pivot) i++;
    while(x[j] > pivot) j--;
    if(i <= j){
        tmp = x[i];
        x[i] = x[j];
        x[j] = tmp;
        tmp = p[i];
        p[i] = p[j];
        p[j] = tmp;
        i++;
        j--;
    }
  }
  if(l < j) QuickSort(x, p, l, j);
  if(i < r) QuickSort(x, p, i, r);
}

void SUBoom::SearchLinear(CConfig *config, CGeometry *geometry, 
               const su2double r0, const su2double *phi, unsigned short nPhi){
  
  /*--- Loop over boundary markers ---*/
  unsigned long nMarker = config->GetnMarker_All();
  unsigned long iMarker, iVertex, iPoint;
  unsigned long jElem;
  unsigned short iElem, nElem;

  bool inside=false;

  su2double Minf = config->GetMach();
  su2double x = 0.0, y = 0.0, z = 0.0;
  su2double *Coord;
  su2double *p0, *p1;

  nPanel = new unsigned long[nPhi];
  for(unsigned short i = 0; i < nPhi; i++){
    nPanel[i] = 0;
  }

  p0 = new su2double[nDim];
  p1 = new su2double[nDim];

  /*--- Search on boundaries ---*/  
  for(iMarker = 0; iMarker < nMarker; iMarker++){
    /*--- Only look at farfield boundary (or send/recv for parallel computations) ---*/
    if(config->GetMarker_All_KindBC(iMarker) == FAR_FIELD || config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE){
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
              for(unsigned short iPhi = 0; iPhi < nPhi; iPhi++){
                inside = InsideElem(geometry, r0, phi[iPhi], jElem, p0, p1);
                if(inside){
                  if(nPanel[iPhi] == 0){
                    nPanel[iPhi] = 1;
                    pointID_original = new unsigned long*[nPhi];
                    pointID_original[iPhi] = new unsigned long[1];
                    Coord_original = new su2double**[nPhi];
                    Coord_original[iPhi] = new su2double*[1];
                    Coord_original[iPhi][0] = new su2double[nDim];

                    pointID_original[iPhi][0] = jElem;
                    Coord_original[iPhi][0][0] = (p0[0] + p1[0])/2.0;

                    if(nDim == 2){
                      Coord_original[iPhi][0][1] = -r0;
                    }

                    else{
                      Coord_original[iPhi][0][1] = r0*sin(ray_phi[iPhi]);
                      Coord_original[iPhi][0][2] = -r0*cos(ray_phi[iPhi]);
                    }
                  }

                  else{
                    nPanel[iPhi]++;

                    unsigned long *pointID_tmp = new unsigned long[nPanel[iPhi]-1];
                    su2double **Coord_tmp = new su2double*[nPanel[iPhi]-1];
                    for(unsigned long i = 0; i < nPanel[iPhi]-1; i++){
                      Coord_tmp[i] = new su2double[nDim];
                      pointID_tmp[i] = pointID_original[iPhi][i];
                      Coord_tmp[i][0] = Coord_original[iPhi][i][0];
                      Coord_tmp[i][1] = Coord_original[iPhi][i][1];

                      delete [] Coord_original[iPhi][i];
                    }
                    delete [] pointID_original[iPhi];
                    delete [] Coord_original[iPhi];

                    pointID_original[iPhi] = new unsigned long[nPanel[iPhi]];
                    Coord_original[iPhi] = new su2double*[nPanel[iPhi]];
                    for(unsigned long i = 0; i < nPanel[iPhi]-1; i++){
                      Coord_original[iPhi][i] = new su2double[nDim];
                      pointID_original[iPhi][i] = pointID_tmp[i];
                      Coord_original[iPhi][i][0] = Coord_tmp[i][0];
                      Coord_original[iPhi][i][1] = Coord_tmp[i][1];

                      delete [] Coord_tmp[i];
                    }
                    delete [] pointID_tmp;
                    delete [] Coord_tmp;

                    Coord_original[iPhi][nPanel[iPhi]-1] = new su2double[nDim];
                    pointID_original[iPhi][nPanel[iPhi]-1] = jElem;
                    Coord_original[iPhi][nPanel[iPhi]-1][0] = (p0[0] + p1[0])/2.0;

                    if(nDim == 2){
                      Coord_original[iPhi][nPanel[iPhi]-1][1] = -r0;
                    }

                    else{
                      Coord_original[iPhi][nPanel[iPhi]-1][1] = -r0*sin(ray_phi[iPhi]);
                      Coord_original[iPhi][nPanel[iPhi]-1][2] = -r0*cos(ray_phi[iPhi]);
                    }
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

  delete [] Coord;
  delete [] p0;
  delete [] p1;

}

void SUBoom::ExtractLine(CGeometry *geometry, const su2double r0, unsigned short iPhi){
  bool inside, inside_iPanel, addPanel, end = false;
  unsigned short iElem, nElem;
  unsigned long jElem, jElem_m1, nElem_tot = geometry->GetnElem();
  su2double x_i, x_m1;

  unsigned long *pointID_tmp;
  su2double **Coord_tmp;
  su2double *p0 = new su2double[nDim], *p1 = new su2double[nDim];

  
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
          inside = InsideElem(geometry, r0, ray_phi[iPhi], jElem, p0, p1);
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
              Coord_tmp = new su2double*[nPanel[iPhi]-1];
              for(unsigned long i = 0; i < nPanel[iPhi]-1; i++){
                Coord_tmp[i] = new su2double[nDim];
                pointID_tmp[i] = pointID_original[iPhi][i];
                Coord_tmp[i][0] = Coord_original[iPhi][i][0];
                Coord_tmp[i][1] = Coord_original[iPhi][i][1];

                delete [] Coord_original[iPhi][i];
              }
              delete [] pointID_original[iPhi];
              delete [] Coord_original[iPhi];

              pointID_original[iPhi] = new unsigned long[nPanel[iPhi]];
              Coord_original[iPhi] = new su2double*[nPanel[iPhi]];
              for(unsigned long i = 0; i < nPanel[iPhi]-1; i++){
                Coord_original[iPhi][i] = new su2double[nDim];
                pointID_original[iPhi][i] = pointID_tmp[i];
                Coord_original[iPhi][i][0] = Coord_tmp[i][0];
                Coord_original[iPhi][i][1] = Coord_tmp[i][1];

                delete [] Coord_tmp[i];
              }
              delete [] pointID_tmp;
              delete [] Coord_tmp;

              Coord_original[iPhi][nPanel[iPhi]-1] = new su2double[nDim];
              pointID_original[iPhi][nPanel[iPhi]-1] = jElem;
              Coord_original[iPhi][nPanel[iPhi]-1][0] = (p0[0] + p1[0])/2.0;
              if(nDim == 2){
                Coord_original[iPhi][nPanel[iPhi]-1][1] = -r0;
              }
              else{
                Coord_original[iPhi][nPanel[iPhi]-1][1] = -r0*sin(ray_phi[iPhi]);
                Coord_original[iPhi][nPanel[iPhi]-1][2] = -r0*cos(ray_phi[iPhi]);
              }

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

void SUBoom::ExtractPressure(CSolver *solver, CConfig *config, CGeometry *geometry, unsigned short iPhi){
  unsigned short iDim, iNode, nNode;
  unsigned long iElem, jElem, jNode, jjNode, iPoint;
  unsigned long nNode_list, *jNode_list;
  unsigned long pointCount = 0;
  su2double rho, rho_ux, rho_uy, rho_uz, rho_E, TKE;
  su2double ux, uy, uz, StaticEnergy;
  bool addNode;

  su2double **isoparams;
  su2double *X_donor;
  su2double *Coord = new su2double[nDim];

  isoparams = new su2double*[nPanel[iPhi]];

  for(unsigned long i = 0; i < nPanel[iPhi]; i++){
    jElem = pointID_original[iPhi][i];
    nNode = geometry->elem[jElem]->GetnNodes();
    nPointID[iPhi] += nNode;
  }
  PointID[iPhi] = new unsigned long[nPointID[iPhi]];
  jNode_list = new unsigned long[nPointID[iPhi]];
  nNode_list = 0;

  signal.x[iPhi] = new su2double[nPanel[iPhi]];
  signal.original_p[iPhi] = new su2double[nPanel[iPhi]];
  for(unsigned long i = 0; i < nPanel[iPhi]; i++){
    /*--- Get info needed for isoparameter computation ---*/
    jElem = pointID_original[iPhi][i];
    nNode = geometry->elem[jElem]->GetnNodes();
    for(unsigned short j = 0; j < nDim; j++){
      Coord[j] = Coord_original[iPhi][i][j];
    }
    X_donor = new su2double[nDim*nNode];
    for(iNode = 0; iNode < nNode; iNode++){
      jNode = geometry->elem[jElem]->GetNode(iNode);
      for(iDim = 0; iDim < nDim; iDim++){  
        X_donor[iDim*nNode + iNode] = geometry->node[jNode]->GetCoord(iDim);
      }

      /*--- Compile list of all nodes ---*/
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

    /*--- Compute isoparameters ---*/
    isoparams[i] = new su2double[nNode];
    Isoparameters(nDim, nNode, X_donor, Coord, isoparams[i]);

    /*--- x-locations of nearfield signal ---*/
    signal.x[iPhi][i] = Coord_original[iPhi][i][0];
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

    //Register conservative variables as input for adjoint computation
    if (config->GetAD_Mode()){
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
          signal.original_p[iPhi][iElem] += (config->GetGamma()-1)*rho*StaticEnergy*isoparams[iElem][iNode];
        }
      }
    }
  }

  nPointID[iPhi] = nNode_list;

  /*--- Clean up interpolation variables ---*/
  for(iElem = 0; iElem < nPanel[iPhi]; iElem++){
    delete [] isoparams[iElem];
  }
  delete [] isoparams;
  
}

bool SUBoom::InsideElem(CGeometry *geometry, su2double r0, su2double phi, unsigned long jElem, su2double *p0, su2double *p1){
  bool inside = false;
  unsigned long iPoint, jNode;
  unsigned short iNode, nNode, count, intersect;

  su2double *pp0 = new su2double[nDim];
  su2double *pp1 = new su2double[nDim];

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
        intersect = Intersect2D(r0, Coord_elem[iEdge], Coord_elem[iEdge_p1], pp0, pp1);
        if(intersect == 1){
          if(count == 0){
            p0[0] = pp0[0];
            p0[1] = pp0[1];
          }
          else{
            p1[0] = pp0[0];
            p1[1] = pp0[1];
          }
        }
        else if(intersect == 2){
          p0[0] = pp0[0];
          p0[1] = pp0[1];
          p1[0] = pp1[0];
          p1[1] = pp1[1];
          inside = true;
          break;
        }
        count += intersect;
        if(count > 1){
          inside = true;
          break;
        }
    }

    for(iNode = 0; iNode < nNode; iNode++){
      delete [] Coord_elem[iNode];
    }
    delete [] Coord_elem;
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
      intersect = Intersect3D(r0, phi, nNodeFace, Coord_face, pp0);
      if(intersect == 1){
        if(count == 0){
          p0[0] = pp0[0];
          p0[1] = pp0[1];
          p0[2] = pp0[2];
        }
        else{
          p1[0] = pp0[0];
          p1[1] = pp0[1];
          p1[2] = pp0[2];
        }
        count++;
        if(count > 1){
          inside = true;
          break;
        }
      }
    }

    for(iNode = 0; iNode < nNodeFace; iNode++){
      delete [] Coord_face[iNode];
    }
    delete [] Coord_face;
  }

  return inside;
}

int SUBoom::Intersect2D(su2double r0, su2double *Coord_i, su2double *Coord_ip1, su2double *p0, su2double *p1){

  /*--- Interpolate if point between yi and yi+1 ---*/
  if((Coord_i[1] > -r0 && Coord_ip1[1] < -r0) || (Coord_i[1] < -r0 && Coord_ip1[1] > -r0)){
    su2double t = (-r0-Coord_i[1])/(Coord_ip1[1]-Coord_i[1]);
    p0[0] = Coord_i[0] + t*(Coord_ip1[0]-Coord_i[0]);
    return 1;
  }
  /*--- Colinear segments at r0 ---*/
  else if(abs(Coord_i[1] + r0) < 1.0E-8 && abs(Coord_ip1[1] + r0) < 1.0E-8){
    p0[0] = Coord_i[0];
    p0[1] = -r0;
    p1[0] = Coord_ip1[0];
    p1[1] = -r0;
    return 2;
  }
  else{
    return 0;
  }

}

int SUBoom::Intersect3D(su2double r0, su2double phi, int nCoord, su2double **Coord_i, su2double *p1){
  
  /*--- t = -(P0.N + d)/(V.N) where P0 is origin of ray, N is plane normal, d is offset vector of plane, V is direction of ray ---*/
  su2double normal[3] = {Coord_i[0][1]*Coord_i[1][2] - Coord_i[0][2]*Coord_i[1][1],
                      Coord_i[0][2]*Coord_i[1][0] - Coord_i[0][0]*Coord_i[1][2],
                      Coord_i[0][0]*Coord_i[1][1] - Coord_i[0][1]*Coord_i[1][0]};

  /*--- We only care about rays in +x direction, so only use x components in denominator---*/
  if(abs(normal[0]) < 1.0E-8){ // Ray and plane are parallel
    return 0;
  }
  else{
    unsigned short i0, i1, i2;
    su2double p0[3] = {-1.0, -r0*sin(phi), -r0*cos(phi)};
    su2double d = -Coord_i[0][0]*normal[0] - Coord_i[0][1]*normal[1] - Coord_i[0][2]*normal[2];
    su2double t = -(p0[0]*normal[0] + p0[1]*normal[1] + p0[2]*normal[2] + d)/(normal[0]);

    p1[0] = p0[0] + t;
    p1[1] = p0[1];
    p1[2] = p0[2];

    /*--- Project onto yz plane ---*/
    for(unsigned short i = 0; i < nCoord-2; i++){
      if(nCoord == 3){
        i0 = 0;
        i1 = 1;
        i2 = 2;
      }
      else{ // if 4 points, check if inside triangle 0-1-2 and 2-3-0
        i0 = 2*i;
        i1 = i0+1;
        i2 = (i0+2)%nCoord;
      }
      su2double u0 = p1[1] - Coord_i[i0][1];
      su2double v0 = p1[2] - Coord_i[i0][2];
      su2double u1 = Coord_i[i1][1] - Coord_i[i0][1];
      su2double u2 = Coord_i[i2][1] - Coord_i[i0][1];
      su2double v1 = Coord_i[i1][2] - Coord_i[i0][2];
      su2double v2 = Coord_i[i2][2] - Coord_i[i0][2];
      su2double alpha = -1.0, beta;

      if(u1 < 1.0E-8){
        beta = u0 / u2;
        if(0 <= beta && beta <= 1){
          alpha = (v0 - beta*v2)/v1;
        }
      }
      else{
        beta = (v0*u1 - u0*v1)/(v2*u1 - u2*v1);
        if(0 <= beta && beta <= 1){
          alpha = (u0 - beta*u2)/u1;
        }

      }

      if(alpha >= 0.0 && beta >= 0.0 && (alpha+beta) <= 1.0){
        return 1;
      }
    }



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
      for (iDim=0; iDim<n; iDim++)
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
  if (nDonor==4) {
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

void SUBoom::AtmosISA(su2double& h0, su2double& T, su2double& a, su2double& p,
                      su2double& rho, su2double& g){
  /*---Calculate temperature, speed of sound, pressure, and density at a given
       altitude using standard atmosphere---*/
  const su2double GMR = 34.163195;    // hydrostatic constant
  const su2double R = 287.058;    // gas constant of air [m^2/s^2*K]
  const su2double Re_earth = 6371.;    // equatorial radius [km]

  su2double h = (h0/1000.)*Re_earth/(h0/1000.+Re_earth);    // geometric to geopotential altitude

  su2double htab[8] = {0.0, 11.0, 20.0, 32.0, 47.0, 51.0, 71.0, 84.852};
  su2double ttab[8] = {288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.946};
  su2double ptab[8] = {1.0, 2.233611E-1, 5.403295E-2, 8.5666784E-3, 1.0945601E-3,
                    6.6063531E-4, 3.9046834E-5, 3.68501E-6};
  su2double gtab[8] = {-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0, 0.0};
  su2double tgrad, tbase, tlocal, deltah, theta, delta, sigma;

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

  if(tgrad == 0.0){
    delta = ptab[i-1]*exp(-GMR*deltah/tbase);    // pressure ratio wrt sea-level
  }
  else{
    delta = ptab[i-1]*pow((tbase/tlocal), GMR/tgrad);
  }

  sigma = delta/theta;    // density ratio wrt sea-level

  T = theta*288.15;
  a = sqrt(g*R*T);
  p = delta*101325.;
  rho = sigma*1.225;

}

void SUBoom::ConditionAtmosphericData(){
  su2double g = atm_g;    // specific heat ratio
  su2double h0 = flt_h;    // flight altitude [m]
  su2double T, p;

  /*---Conditions at flight altitude---*/
  AtmosISA(h0, T_inf, a_inf, p_inf, rho_inf, g);

  /*---Atmospheric conditions profile---*/
  z = new su2double[n_prof];
  a_of_z = new su2double[n_prof];
  p_of_z = new su2double[n_prof];
  rho_of_z = new su2double[n_prof];
  for(unsigned int i = 0; i < n_prof; i++){
    z[i] = h0*su2double(i)/(su2double(n_prof)-1.);
    AtmosISA(z[i], T, a_of_z[i], p_of_z[i], rho_of_z[i], g);
  }
}

void SUBoom::ScaleFactors(){

  scale_T = scale_L/(flt_M*a_inf);    // flow over aircraft [s]
  scale_p = p_inf;    // ambient [Pa]
  scale_m = scale_p/scale_T;    // slope of boom signal [Pa/s]
  scale_z = flt_h;    // altitude [m]

  scale_C1 = scale_p;    // [Pa]
  scale_C2 = scale_T;    // [s]

}

void SUBoom::DistanceToTime(){
  int len;
  for(unsigned short iPhi = 0; iPhi < ray_N_phi; iPhi++){
    len = signal.original_len[iPhi];
    for(int i = 0; i < len; i++){
      signal.original_T[iPhi][i] = signal.x[iPhi][i]/(a_inf*flt_M);
    }
  }
}

void SUBoom:: Sph2Cart(su2double& nx, su2double& ny, su2double& nz, su2double az, su2double elev,
              su2double r){
  /*---Compute spherical coordinates from elevation, azimuth, and radius---*/
  nx = r*cos(elev)*cos(az);
  ny = r*cos(elev)*sin(az);
  nz = r*sin(elev);
}

void SUBoom::InitialWaveNormals(){

  su2double a0;
  su2double nx, ny, nz;
  const su2double deg2rad = M_PI/180.;
  su2double M = flt_M;    // Mach
  su2double psi = flt_psi*deg2rad;    // heading angle [rad]
  su2double g = flt_gamma*deg2rad;    // climb angle [rad]
  su2double phi[ray_N_phi];    // azimuth angles of wave normal from vertical plane [rad]
  su2double phi_tube[ray_N_phi];
  su2double mu = asin(1./M);    // Mach angle [rad]
  flt_mu = mu;

  /*---Create arrays---*/
  ray_nu = new su2double*[ray_N_phi];
  ray_c0 = new su2double*[ray_N_phi];
  ray_theta0 = new su2double*[ray_N_phi];
  for(int i = 0; i < ray_N_phi; i++){
    ray_nu[i] = new su2double[2];
    ray_c0[i] = new su2double[2];
    ray_theta0[i] = new su2double[2];
  }

  /*---Initialize ray parameters---*/
  for(int i = 0; i < ray_N_phi; i++){
    ray_nu[i][0] = ray_nu[i][1] = 0.;
    ray_theta0[i][0] = ray_theta0[i][1] = 0.;
    ray_c0[i][0] = ray_c0[i][1] = 0.;

    phi[i] = ray_phi[i];
    phi_tube[i] = phi[i] + tol_dphi*deg2rad;
  }

  /*---Compute initial ray parameters---*/
  a0 = a_of_z[n_prof-1];

  for(int i = 0; i < ray_N_phi; i++){
    /*---Primary rays---*/
    ray_nu[i][0] = psi - atan2(cos(mu)*sin(phi[i]),
                               sin(mu)*cos(g) + cos(mu)*sin(g)*cos(phi[i]));
    ray_theta0[i][0] = asin(sin(mu)*sin(g) - cos(mu)*cos(g)*cos(phi[i]));
    ray_c0[i][0] = a0/cos(ray_theta0[i][0]);

    /*---Ray tube corners---*/
    ray_nu[i][1] = psi - atan2(cos(mu)*sin(phi_tube[i]),
                               sin(mu)*cos(g) + cos(mu)*sin(g)*cos(phi_tube[i]));
    ray_theta0[i][1] = asin(sin(mu)*sin(g) - cos(mu)*cos(g)*cos(phi_tube[i]));
    ray_c0[i][1] = a0/cos(ray_theta0[i][1]);
  }

  /*---Compute heading unit vector---*/
  Sph2Cart(nx, ny, nz, M_PI/2.-psi, g, 1.0);
  flt_heading[0] = nx;
  flt_heading[1] = ny;
  flt_heading[2] = nz;

}

su2double *derivs(su2double x, int m, su2double y[], SUBoom::RayData data){
  su2double *dydx;
  su2double a, rho, p, T;
  su2double g = 1.4;
  su2double xdim = x*data.L;

  SUBoom boom;
  boom.AtmosISA(xdim, T, a, p, rho, g);

  su2double theta = acos(a/data.c0);
  su2double num = cos(theta);
  su2double denom = sin(theta);

  dydx = new su2double[3];
  dydx[0] = (num*sin(data.nu))/denom;
  dydx[1] = (num*cos(data.nu))/denom;
  dydx[2] = (data.L/(data.T*a))/denom;

  return dydx;
}

su2double *SUBoom::rk4(su2double x0, int m, su2double y0[], su2double dx, RayData data,
                    su2double *f(su2double x, int m, su2double y[], RayData data)){
  su2double *f0, *f1, *f2, *f3;
  su2double x1, x2, x3;
  su2double *y, *y1, *y2, *y3;

  SUBoom boom;

  // Intermediate steps
  f0 = f(x0, m, y0, data);
  x1 = x0 + dx/2.0;
  y1 = new su2double[m];
  for(int i = 0; i < m; i++){
    y1[i] = y0[i] + dx*f0[i]/2.0;
  }

  f1 = f(x1, m, y1, data);
  x2 = x0 + dx/2.0;
  y2 = new su2double[m];
  for(int i = 0; i < m; i++){
    y2[i] = y0[i] + dx*f1[i]/2.0;
  }

  f2 = f(x2, m, y2, data);
  x3 = x0 + dx;
  y3 = new su2double[m];
  for(int i = 0; i < m; i++){
    y3[i] = y0[i] + dx*f2[i];
  }

  f3 = f(x3, m, y3, data);

  // Now estimate solution
  y = new su2double[m];
  for(int i = 0; i < m; i++){
    y[i] = y0[i] + dx*(f0[i] + 2.0*f1[i] + 2.0*f2[i] + f3[i])/6.0;
  }

  delete [] f0;
  delete [] f1;
  delete [] f2;
  delete [] f3;
  delete [] y1;
  delete [] y2;
  delete [] y3;

  return y;
}

su2double *SUBoom::SplineGetDerivs(su2double x[], su2double y[], int N){
  su2double *ks;    // Derivatives vector, using Wikipedia notation q'(x)=k
  su2double *a, *b, *c, *d;
  su2double *cp, *dp;

  /*---Create a, b, c, d vectors---*/
  a = new su2double[N];
  b = new su2double[N];
  c = new su2double[N];
  d = new su2double[N];
  cp = new su2double[N];
  dp = new su2double[N];

  /*---Create tridiagonal system---*/
  for(int i = 1; i < N-1; i++){
    a[i] = 1/(x[i]-x[i-1]);
    b[i] = 2*(1/(x[i]-x[i-1]) + 1/(x[i+1]-x[i]));
    c[i] = 1/(x[i+1]-x[i]);
    d[i] = 3*((y[i]-y[i-1])/((x[i]-x[i-1])*(x[i]-x[i-1]))
                + (y[i+1]-y[i])/((x[i+1]-x[i])*(x[i+1]-x[i])));
  }
  b[0] = 2/(x[1]-x[0]);
  c[0] = 1/(x[1]-x[0]);
  d[0] = 3*(y[1]-y[0])/((x[1]-x[0])*(x[1]-x[0]));

  a[N-1] = 1/(x[N-1]-x[N-2]);
  b[N-1] = 2/(x[N-1]-x[N-2]);
  d[N-1] = 3*(y[N-1]-y[N-2])/((x[N-1]-x[N-2])*(x[N-1]-x[N-2]));

  /*---Forward sweep---*/
  cp[0] = c[0]/b[0];
  dp[0] = d[0]/b[0];
  for(int i = 1; i < N; i++){
    cp[i] = c[i]/(b[i]-a[i]*cp[i-1]);
    dp[i] = (d[i]-a[i]*dp[i-1])/(b[i]-a[i]*cp[i-1]);
  }

  /*---Back substitution---*/
  ks = new su2double[N];
  ks[N-1] =dp[N-1];
  for(int i = N-2; i > -1; i--){
    ks[i] = dp[i]-cp[i]*ks[i+1];
  }

  /*---Clear up memory---*/
  delete [] a;
  delete [] b;
  delete [] c;
  delete [] d;
  delete [] cp;
  delete [] dp;

  return ks;
}

void SUBoom::RayTracer(unsigned short iPhi){

  /*---Scale factors---*/
  su2double L = flt_h;
  su2double T = scale_T;
  su2double a0, rho0, p0;
  su2double a, rho, p;
  su2double r0[3];
  su2double *f;
  su2double dz = (z[0] - z[1]);

  //GetAtmosphericData(a0, rho0, p0, L);
  a0 = a_inf;
  rho0 = rho_inf;
  p0 = p_inf;

  /*---Class for passing data to RK4 solver---*/
  RayData data;
  data.L = L;
  data.T = T;

  /*---Create arrays---*/
  f = new su2double[3];

  x_of_z = new su2double*[4];
  y_of_z = new su2double*[4];
  t_of_z = new su2double*[4];
  ray_theta = new su2double*[4];
  for(unsigned short j = 0; j < 4; j++){
    x_of_z[j] = new su2double[n_prof];
    y_of_z[j] = new su2double[n_prof];
    t_of_z[j] = new su2double[n_prof];
    ray_theta[j] = new su2double[n_prof];
  }

  /*---Primary ray---*/
  data.c0 = ray_c0[iPhi][0];
  data.nu = ray_nu[iPhi][0];
  r0[0] = r0[1] = r0[2] = 0.;
  x_of_z[0][n_prof-1] = r0[0];
  y_of_z[0][n_prof-1] = r0[1];
  t_of_z[0][n_prof-1] = r0[2];
  ray_theta[0][n_prof-1] = acos(a_of_z[n_prof-1]/data.c0);
  for(int j = n_prof-1; j > 0; j--){
    f = rk4(z[j]/data.L, 3, r0, dz/data.L, data, derivs);
    x_of_z[0][j-1] = -f[0];
    y_of_z[0][j-1] = -f[1];
    t_of_z[0][j-1] = -f[2];

    r0[0] = -x_of_z[0][j-1];
    r0[1] = -y_of_z[0][j-1];
    r0[2] = -t_of_z[0][j-1];

    /*---Derived data---*/
    ray_theta[0][j-1] = acos(a_of_z[j-1]/data.c0);
  }

  /*---Ray tube corners: {0, +dheading}---*/
  r0[0] = flt_heading[0]*tol_dr;
  r0[1] = flt_heading[1]*tol_dr;
  r0[2] = flt_heading[2]*tol_dr;
  x_of_z[1][n_prof-1] = r0[0];
  y_of_z[1][n_prof-1] = r0[1];
  ray_theta[1][n_prof-1] = acos(a_of_z[n_prof-1]/data.c0);
  for(int j = n_prof-1; j > 0; j--){
    f = rk4(z[j]/data.L, 3, r0, dz/data.L, data, derivs);
    x_of_z[1][j-1] = -f[0];
    y_of_z[1][j-1] = -f[1];
    t_of_z[1][j-1] = -f[2];

    r0[0] = -x_of_z[1][j-1];
    r0[1] = -y_of_z[1][j-1];
    r0[2] = -t_of_z[1][j-1];

    /*---Derived data---*/
    ray_theta[1][j-1] = acos(a_of_z[j-1]/data.c0);
  }

  /*---Ray tube corners: {+dphi, 0}---*/
  data.c0 = ray_c0[iPhi][1];
  data.nu = ray_nu[iPhi][1];
  r0[0] = sin(data.nu)*tol_dr;
  r0[1] = cos(data.nu)*tol_dr;
  r0[2] = 0.;
  x_of_z[2][n_prof-1] = r0[0];
  y_of_z[2][n_prof-1] = r0[1];
  t_of_z[2][n_prof-1] = r0[2];
  ray_theta[2][n_prof-1] = acos(a_of_z[n_prof-1]/data.c0);
  for(int j = n_prof-1; j > 0; j--){
    f = rk4(z[j]/data.L, 3, r0, dz/data.L, data, derivs);
    x_of_z[2][j-1] = -f[0];
    y_of_z[2][j-1] = -f[1];
    t_of_z[2][j-1] = -f[2];

    r0[0] = -x_of_z[2][j-1];
    r0[1] = -y_of_z[2][j-1];
    r0[2] = -t_of_z[2][j-1];

    /*---Derived data---*/
    ray_theta[2][j-1] = acos(a_of_z[j-1]/data.c0);
  }

  /*---Ray tube corners: {+dphi, +dheading}---*/
  r0[0] = flt_heading[0]*tol_dr + sin(data.nu)*tol_dr;
  r0[1] = flt_heading[1]*tol_dr + cos(data.nu)*tol_dr;
  r0[2] = flt_heading[2]*tol_dr;
  x_of_z[3][n_prof-1] = r0[0];
  y_of_z[3][n_prof-1] = r0[1];
  ray_theta[3][n_prof-1] = acos(a_of_z[n_prof-1]/data.c0);
  for(int j = n_prof-1; j > 0; j--){
    f = rk4(z[j]/data.L, 3, r0, dz/data.L, data, derivs);
    x_of_z[3][j-1] = -f[0];
    y_of_z[3][j-1] = -f[1];
    t_of_z[3][j-1] = -f[2];

    r0[0] = -x_of_z[3][j-1];
    r0[1] = -y_of_z[3][j-1];
    r0[2] = -t_of_z[3][j-1];

    /*---Derived data---*/
    ray_theta[3][j-1] = acos(a_of_z[j-1]/data.c0);
  }

  /*---Clear up memory---*/
  delete [] f;

}

void SUBoom::RayTubeArea(unsigned short iPhi){

  su2double Ah, x_int, y_int, z_int;
  su2double corners[4][3];
  int M;

  ray_A = new su2double[n_prof];

  /*---Loop over rays---*/
  M = n_prof;
  Ah = 0;
  for(int j = 0; j < M; j++){
    for(int k = 0; k < 4; k++){
      for(int kk = 0; kk < 3; kk++){
        corners[k][kk] = 0.0;
      }
    }
    z_int = z[j]/scale_z;
    for(int k = 0; k < 4; k++){
  	  x_int = x_of_z[k][j];
      y_int = y_of_z[k][j];
      corners[k][0] = x_int;
      corners[k][1] = y_int;
      corners[k][2] = z_int;
    }
    su2double u[3] = {corners[3][0]-corners[0][0], corners[3][1]-corners[0][1], corners[3][2]-corners[0][2]};
    su2double v[3] = {corners[2][0]-corners[1][0], corners[2][1]-corners[1][1], corners[2][2]-corners[1][2]};
    /*---Cross product---*/
    su2double c[3] = {u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0]};
    Ah = 0.5*sqrt(pow(c[0],2)+pow(c[1],2)+pow(c[2],2));
    ray_A[j] = Ah*a_of_z[j]*tan(ray_theta[0][j])/ray_c0[iPhi][0];
  }

}

su2double SUBoom::matchr(int j, su2double h_L, su2double r0){
  su2double f;
  f = sqrt(pow(x_of_z[0][j],2) + pow(y_of_z[0][j],2) + pow(z[j]/scale_z-1.0,2))*h_L - r0/scale_L;
  return f;
}

void SUBoom::FindInitialRayTime(){

  su2double h_L = scale_z/scale_L;
  su2double eps;

  su2double f[n_prof], t[n_prof];
  su2double f0, f1, t0, t1, tmp, trat;
  //int flag;

  su2double a, b;
  su2double *ks;

  ks = new su2double[n_prof];

  for(unsigned int j = 0; j < n_prof; j++){
    t[j] = t_of_z[0][j];
    f[j] = matchr(j, h_L, ray_r0);
  }

  /*---Get spline---*/
  ks = SplineGetDerivs(t, f, n_prof);

  /*---Secant method---*/
  f0 = f[(n_prof-1)/3];
  f1 = f[2*(n_prof-1)/3];
  t0 = t[(n_prof-1)/3];
  t1 = t[2*(n_prof-1)/3];
  eps = 1.E6;
  while(eps > 1.E-8){
    tmp = t1 - f1*(t1-t0)/(f1-f0);
    t0 = t1;
    t1 = tmp;
    f0 = f1;
    // Find interval which contains t1
    unsigned int j_sp = n_prof-1;
    for(int j = n_prof-1; j > -1; j--){
    if(t1 < t[j]){
      j_sp = j+1;
      break;
    }
    }
    if(j_sp == 0){j_sp = 1;}

    trat = (t1-t[j_sp-1])/(t[j_sp]-t[j_sp-1]);
    a = ks[j_sp-1]*(t[j_sp]-t[j_sp-1]) - (f[j_sp]-f[j_sp-1]);
    b = -ks[j_sp]*(t[j_sp]-t[j_sp-1]) + (f[j_sp]-f[j_sp-1]);

    f1 = (1.0-trat)*f[j_sp-1] + trat*f[j_sp] + trat*(1-trat)*(a*(1-trat)+b*trat);
    eps = abs(t1-t0);
  }

  ray_t0 = t1;

  /*---Clear up memory---*/
  delete [] ks;

}

void SUBoom::ODETerms(unsigned short iPhi){

  su2double *cn;
  su2double *dadt, *drhodt, *dAdt;
  su2double *dcndt;
  su2double g = atm_g;

  cn = new su2double[n_prof];

  dadt = new su2double[n_prof];
  drhodt = new su2double[n_prof];
  dAdt = new su2double[n_prof];
  dcndt = new su2double[n_prof];

  for (unsigned int j = 0; j < n_prof; j++){
    cn[j] = ray_c0[iPhi][0]*cos(ray_theta[0][j]);
  }

  /*--- Spline interpolation of a, rho, A---*/
  dadt = SplineGetDerivs(t_of_z[0], a_of_z, n_prof);
  drhodt = SplineGetDerivs(t_of_z[0], rho_of_z, n_prof);
  dAdt = SplineGetDerivs(t_of_z[0], ray_A, n_prof);
  dcndt = SplineGetDerivs(t_of_z[0], cn, n_prof);

  ray_C1 = new su2double[n_prof];
  ray_C2 = new su2double[n_prof];
  for(unsigned int j = 0; j < n_prof; j++){
    ray_C1[j] = ((g+1.)/(2.*g))*a_of_z[j]/(p_of_z[j]*cn[j])*scale_p; // TESTING scaling by pressure (See Thomas paper)
    ray_C2[j] = 0.5*((3./a_of_z[j])*dadt[j] + drhodt[j]/rho_of_z[j] - (2./cn[j])*dcndt[j] - (1./ray_A[j])*dAdt[j]);
  }
  ray_dC1 = SplineGetDerivs(t_of_z[0], ray_C1, n_prof);
  ray_dC2 = SplineGetDerivs(t_of_z[0], ray_C2, n_prof);

  delete [] ray_c0[iPhi];
  delete [] ray_nu[iPhi];
  delete [] ray_theta0[iPhi];
  
  delete [] cn;
  delete [] dadt;
  delete [] drhodt;
  delete [] dAdt;
  delete [] dcndt;
}

void SUBoom::CreateSignature(unsigned short iPhi){
  int len = signal.original_len[iPhi];
  su2double pp[2][len-1];
  su2double ll[len-1];
  su2double mm[len-1];

  /*---Pressure signal to waveform---*/
  /*---Collect data into segments---*/
  for(int i = 0; i < len-1; i++){
    pp[0][i] = signal.original_p[iPhi][i]/scale_p;  // [Pa] (relative)
    pp[1][i] = signal.original_p[iPhi][i+1]/scale_p;
    ll[i] = (signal.original_T[iPhi][i+1] - signal.original_T[iPhi][i])/scale_T;  // [s]
    mm[i] = (pp[1][i] - pp[0][i])/ll[i]; // [Pa/s]
  }

  /*---Remove shocks---*/
  int M = len-1;
  int i = 0;
  while(i <= M-1){
    if(mm[i] > tol_m/scale_m){// || mm[i] < -tol_m/scale_m){  // shock present
      /*---Remove segment i---*/
      for(int j = i; j < M; j++){
        pp[0][j] = pp[0][j+1];
        pp[1][j] = pp[1][j+1];
	      ll[j] = ll[j+1];
	      mm[j] = mm[j+1];
        //ll[i] = tol_l;
        //mm[i] = (pp[1][i] - pp[0][i])/ll[i];
      }
      i--;
      M--;
    }
    else if(mm[i] < -tol_m/scale_m && ll[i] < tol_l/scale_T){  // "expansion shock" present
      /*---Remove segment i---*/
      ll[i] = tol_l/scale_T;
      mm[i] = (pp[1][i] - pp[0][i])/ll[i];
    }
    i++;
  }

  /*---Record signal---*/
  signal.dp = new su2double[M];
  signal.m = new su2double[M];
  signal.l = new su2double[M];
  signal.M = M;

  signal.dp[0] = pp[0][0];
  signal.m[0] = mm[0];
  signal.l[0] = ll[0];
  for(int i = 1; i < M; i++){
    signal.dp[i] = pp[0][i] - pp[1][i-1];
    signal.m[i] = mm[i];
    signal.l[i] = ll[i];
  }

}

su2double **WaveformToPressureSignal(su2double fvec[], unsigned int M, int &Msig){
  su2double TT[2][M], pp[2][M];
  su2double m[M], dp[M], l[M];
  su2double **sig;

  /*---Get m, dp, l from fvec---*/
  for(unsigned int i = 0; i < M; i++){
    m[i] = fvec[i];
    dp[i] = fvec[i+M];
    l[i] = fvec[i+2*M];
  }

  /*---First segment---*/
  TT[0][0] = 0;
  TT[1][0] = l[0];
  pp[0][0] = dp[0];
  pp[1][0] = pp[0][0] + m[0]*l[0];

  /*---Build signal---*/
  for(unsigned int seg = 1; seg < M; seg++){
    /*---First node---*/
    TT[0][seg] = TT[1][seg-1];
    pp[0][seg] = pp[1][seg-1] + dp[seg];

    /*---Second node---*/
    TT[1][seg] = TT[0][seg] + l[seg];
    pp[1][seg] = pp[0][seg] + m[seg]*l[seg];
  }

  if(Msig < 0){
    /*---Just return pressure value at segments---*/
    sig = new su2double*[2];
    for(unsigned short j = 0; j < 2; j++){sig[j] = new su2double[M];}

    for(unsigned int j = 0; j < M; j++){
        sig[0][j] = pp[0][j];
        sig[1][j] = pp[1][j];

    }
  }
  else{
    /*---Build 2-D vector---*/
    sig = new su2double*[2];
    for(unsigned int j = 0; j < 2; j++){sig[j] = new su2double[2*M];}

    /*---First segment---*/
    sig[0][0] = TT[0][0];
    sig[0][1] = TT[1][0];
    sig[1][0] = pp[0][0];
    sig[1][1] = pp[1][0];

    int i = 2;
    for(unsigned int seg = 1; seg < M; seg++){
      /*---First node---*/
      if(pp[0][seg] != pp[1][seg-1]){  // pressure jump at this juncture
        sig[0][i] = TT[0][seg];
        sig[1][i] = pp[0][seg];
        i++;
      }

      sig[0][i] = TT[1][seg];
      sig[1][i] = pp[1][seg];
      i++;
      Msig = i;
    }
  }

  return sig;

}

su2double *derivsProp(su2double t, int m, su2double y[], SUBoom::RayData data){
  su2double *dydt;

  int M = data.M;
  su2double diff_m[M], diff_dp[M];
  su2double **current_signal;
  su2double C1, C2;
  int Msig = -1;

  dydt = new su2double[3*M];

  /*---Get final pressure value for diff_dp[end]---*/
  current_signal = new su2double*[2];
  for(int i = 0; i < 2; i++){current_signal[i] = new su2double[M];}
  current_signal = WaveformToPressureSignal(y, M, Msig);

  for(int i = 0; i < M; i++){
    if(i == 0){
      diff_m[i] = y[i];
      diff_dp[i] = y[i+M] + y[i+M+1];
    }
    else if(i == M-1){
      diff_m[i] = y[i] + y[i-1];
      diff_dp[i] =  y[i+M] - current_signal[1][M-1];
    }
    else{
      diff_m[i] = y[i] + y[i-1];
      diff_dp[i] = y[i+M] + y[i+M+1];
    }
  }

  /*---Get coefficients for derivatives---*/
  C1 = EvaluateSpline(t,5,data.t,data.C1,data.dC1);
  C2 = EvaluateSpline(t,5,data.t,data.C2,data.dC2);

  for(int i = 0; i < 3*M; i++){
    if(i < M){
      dydt[i] = C1*y[i]*y[i] + C2*y[i];
    }
    else if(i  < 2*M){
      dydt[i] = 0.5*C1*y[i]*diff_m[i-M] + C2*y[i];
    }
    else{
      dydt[i] = -0.5*C1*diff_dp[i-2*M] - C1*y[i-2*M]*y[i];
    }
  }

  /*---Free memory---*/
  for(int i = 0; i < 2; i++){
    delete [] current_signal[i];
  }
  delete [] current_signal;

  return dydt;
}

su2double *SUBoom::ClipLambdaZeroSegment(su2double fvec[], int &M){
  su2double m[M], dp[M], l[M];
  su2double dp_seg;
  su2double *fvec_new, **current_signal;
  int N = M;
  int Msig = -1;

  /*---Get pressure values to get pressure gap---*/
  current_signal = new su2double*[2];
  for(int j = 0; j < 2; j++){current_signal[j] = new su2double[M];}

  /*---Decompose f vector---*/
  //fvec_new = new su2double[3*M];
  for(int j = 0; j < 3*M; j++){
    //fvec_new[j] = fvec[j];
    if(j < M){
	    m[j] = fvec[j];
    }
    else if(j < 2*M){
	    dp[j-M] = fvec[j];
    }
    else{
	    l[j-2*M] = fvec[j];
    }
  }
  /*---Remove segments with l = 0---*/
  int i = 0;
  while(i <= N-1){
    if(l[i] <= tol_l/scale_T || m[i] >= tol_m/scale_m){// || m[i] <= -tol_m){
      /*---Record pressure gap---*/
      current_signal = WaveformToPressureSignal(fvec, N, Msig);
      dp_seg = dp[i] + (current_signal[1][i] - current_signal[0][i]);
      /*---Add to next segment if needed---*/
      if(dp_seg > tol_dp/scale_p){
          if(i < N-1){
              dp[i+1] = dp[i+1] + dp_seg;
          }
      }

      N--;
      /*---Remove segment---*/
      for(int j = i; j < N; j++){
          m[j] = m[j+1];
          dp[j] = dp[j+1];
          l[j] = l[j+1];
      }
      i--;
      //fvec_new = new su2double[3*N];
      for(int j = 0; j < N; j++){
          fvec[j] = m[j];
          fvec[j+N] = dp[j];
          fvec[j+2*N] = l[j];
      }

    }
    i++;

  }

  fvec_new = new su2double[3*N];
  for(int j = 0; j < N; j++){
    fvec_new[j] = m[j];
    fvec_new[j+N] = dp[j];
    fvec_new[j+2*N] = l[j];
  }

  M = N;

  /*---Free memory---*/
  for(int j = 0; j < 2; j++){
    delete [] current_signal[j];
  }
  delete [] current_signal;

  return fvec_new;

}

su2double EvaluateSpline(su2double x, int N, su2double t[], su2double fit[], su2double coeffs[]){
  su2double x_min, x_max;
  su2double x1, x2;
  su2double g1, g2, y1, y2;
  su2double a, b, rat;
  int j;

  su2double C;

  x_min = t[N-1]; x_max = t[0];

  if(x < x_min){j = N-2;}
  else if(x > x_max){j = 0;}
  else{
    for(int k = 0; k < N-1; k++){
      if(x >= t[k+1] && x <= t[k]){
	    j = k;
        break;
      }
    }
  }

  x1 = t[j+1]; x2 = t[j];
  g1 = coeffs[j+1]; g2 = coeffs[j];
  y1 = fit[j+1];    y2 = fit[j];
  rat = (x-x1)/(x2-x1);
  a = g1*(x2-x1) - (y2-y1);
  b = -g2*(x2-x1) + (y2-y1);
  C = (1.-rat)*y1 + rat*y2 + rat*(1-rat)*(a*(1.-rat) + b*rat);

  return C;

}

void SUBoom::PropagateSignal(unsigned short iPhi){
  su2double t0, tf, dt;
  unsigned int j0, j0_dat;
  RayData data;
  su2double *fvec;
  su2double *f;
  su2double **ground_signal;
  int M = signal.M;


    t0 = ray_t0;
    /*---Assemble f vector and ray data for integration---*/
    fvec = new su2double[3*signal.M];
    for(unsigned int j = 0; j < 3*signal.M; j++){
      if(j < signal.M){
	      fvec[j] = signal.m[j];
      }
      else if(j < 2*signal.M){
	      fvec[j] = signal.dp[j-signal.M];
      }
      else{
	      fvec[j] = signal.l[j-2*signal.M];
      }
    }
    data.t = new su2double[5];
    data.C1 = new su2double[5];
    data.C2 = new su2double[5];
    data.dC1 = new su2double[5];
    data.dC2 = new su2double[5];

    data.M = signal.M;
    data.scale_C1 = scale_C1;
    data.scale_C2 = scale_C2;

    //fvec = signal.fvec;
    for(unsigned int j = n_prof-2; j >= 0; j--){
        if(t_of_z[0][j] >= ray_t0){
            j0 = j+1;
            break;
        }
    }
    for(unsigned int j = j0; j > 0; j--){
      /*---Integrate---*/
      if(j >= n_prof-3) j0_dat = n_prof-3;
      else if(j < 2) j0_dat = 2;
      else j0_dat = j;
      for(unsigned short jj = 0; jj < 5; jj++){
          data.t[jj] = t_of_z[0][j0_dat-2+jj];
          data.C1[jj] = ray_C1[j0_dat-2+jj];
          data.C2[jj] = ray_C2[j0_dat-2+jj];
          data.dC1[jj] = ray_dC1[j0_dat-2+jj];
          data.dC2[jj] = ray_dC2[j0_dat-2+jj];
      }
      tf = t_of_z[0][j-1];
      dt = (tf - t0);
      f = rk4(t0, 3*M, fvec, dt, data, derivsProp);
      fvec = ClipLambdaZeroSegment(f, M);
      data.M = M;
      t0 = tf;
    }

    /*---Waveform to pressure signal---*/
    ground_signal = new su2double*[2];
    for(unsigned short ii = 0; ii < 2; ii++){
      ground_signal[ii] = new su2double[2*M];
    }
    int Msig = M;
    ground_signal = WaveformToPressureSignal(fvec, M, Msig);

    signal.final_T[iPhi] = new su2double[Msig];
    signal.final_p[iPhi] = new su2double[Msig];
    signal.final_len[iPhi] = Msig;

    /*---Final signal and boom strength---*/
    ofstream sigFile;
    sigFile.precision(15);
    sigFile.open("signal_final.dat", ios::out);
    if(iPhi == 0) sigFile << "# phi, T, p" << endl;
    p_max[iPhi] = -1.0E10;
    p_int2[iPhi] = 0.0;
    for(int j = 0; j < Msig; j++){
      signal.final_T[iPhi][j] = ground_signal[0][j]*scale_T;
      signal.final_p[iPhi][j] = ground_signal[1][j]*scale_p;
      sigFile << signal.final_T[iPhi][j] << "\t" << signal.final_p[iPhi][j] << endl;
      if(signal.final_p[iPhi][j] > p_max[iPhi]) p_max[iPhi] = signal.final_p[iPhi][j];
      if(j > 0) p_int2[iPhi] = p_int2[iPhi] + 0.5*(signal.final_p[iPhi][j]*signal.final_p[iPhi][j]+signal.final_p[iPhi][j-1]*signal.final_p[iPhi][j-1])
                        *(signal.final_T[iPhi][j]-signal.final_T[iPhi][j-1]);
    }
    sigFile.close();
    p_rise[iPhi] = signal.final_p[iPhi][0];
    p_rise2[iPhi] = -signal.final_p[iPhi][Msig-1];
    cout << "p_rise = " << p_rise[iPhi] << ", p_max = " << p_max[iPhi] << ", p_int2 = " << p_int2[iPhi] << "." << endl;

  /*---Clean up---*/
  for(unsigned short j = 0; j < 4; j++){
      delete [] x_of_z[j];
      delete [] y_of_z[j];
      delete [] t_of_z[j];
      delete [] ray_theta[j];
  }

  delete [] x_of_z;
  delete [] y_of_z;
  delete [] t_of_z;
  delete [] ray_theta;
  delete [] ray_A;
  delete [] ray_C1;
  delete [] ray_C2;
  delete [] ray_dC1;
  delete [] ray_dC2;

  /*---Clear up memory from RayData class---*/
  delete [] data.t;
  delete [] data.C1;
  delete [] data.C2;
  delete [] data.dC1;
  delete [] data.dC2;

}

void SUBoom::WriteSensitivities(){
  unsigned long iVar, iSig, Max_nPointID, nVar, Global_Index, Total_Index;
  ofstream Boom_AdjointFile;

  int rank, iProcessor, nProcessor;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
  unsigned long Buffer_Send_nPointID[1], *Buffer_Recv_nPointID = NULL;

  if (rank == MASTER_NODE) Buffer_Recv_nPointID= new unsigned long [nProcessor];

    if(rank == MASTER_NODE){
      Boom_AdjointFile.precision(15);
      Boom_AdjointFile.open("Adj_Boom.dat", ios::out);
    }

  for(unsigned short iPhi = 0; iPhi < ray_N_phi; iPhi++){

    Buffer_Send_nPointID[0]=nPointID[iPhi]; 
#ifdef HAVE_MPI
    SU2_MPI::Gather(&Buffer_Send_nPointID, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nPointID, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD); //send the number of vertices at each process to the master
    SU2_MPI::Allreduce(&nPointID[iPhi],&Max_nPointID,1,MPI_UNSIGNED_LONG,MPI_MAX,MPI_COMM_WORLD); //find the max num of vertices over all processes
#endif

    nVar = nDim+3;

    /* pack sensitivity values in each processor and send to root */
    su2double *Buffer_Send_dJdU = new su2double [Max_nPointID*nVar];
    unsigned long *Buffer_Send_GlobalIndex = new unsigned long [Max_nPointID];
    /*---Zero send buffers---*/
    for (unsigned int i=0; i <Max_nPointID*nVar; i++){
      Buffer_Send_dJdU[i]=0.0;
    }
    for (unsigned int i=0; i <Max_nPointID; i++){
      Buffer_Send_GlobalIndex[i]=0;
    }
    su2double *Buffer_Recv_dJdU = NULL;
    unsigned long *Buffer_Recv_GlobalIndex = NULL;

    if (rank == MASTER_NODE) {
      Buffer_Recv_dJdU = new su2double [nProcessor*Max_nPointID*nVar];
      Buffer_Recv_GlobalIndex = new unsigned long[nProcessor*Max_nPointID];
    }

    /*---Fill send buffers with dJ/dU---*/
    for (iVar=0; iVar<nVar; iVar++){
        for(iSig=0; iSig<nPointID[iPhi]; iSig++){
            Buffer_Send_dJdU[iVar*nPointID[iPhi]+iSig] = dJdU[iPhi][iVar][iSig];
          }
      }

    for (iSig=0; iSig<nPointID[iPhi]; iSig++){
       Buffer_Send_GlobalIndex[iSig] = PointID[iPhi][iSig];
    }

#ifdef HAVE_MPI
    SU2_MPI::Gather(Buffer_Send_dJdU, Max_nPointID*nVar, MPI_DOUBLE, Buffer_Recv_dJdU,  Max_nPointID*nVar , MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Gather(Buffer_Send_GlobalIndex,Max_nPointID, MPI_UNSIGNED_LONG, Buffer_Recv_GlobalIndex, Max_nPointID , MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
#endif

    if (rank == MASTER_NODE){

      /*--- Loop through all of the collected data and write each node's values ---*/
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        for (iSig = 0; iSig < Buffer_Recv_nPointID[iProcessor]; iSig++) {
          Global_Index = Buffer_Recv_GlobalIndex[iProcessor*Max_nPointID+iSig];
          /*--- Check dJdU[][iVar==0] and only write if greater than some tolerance ---*/
          Boom_AdjointFile  << scientific << Global_Index << "\t";

          for (iVar = 0; iVar < nVar; iVar++){
            /*--- Current index position and global index ---*/
            Total_Index  = iProcessor*Max_nPointID*nVar + iVar*Buffer_Recv_nPointID[iProcessor]  + iSig;

            /*--- Write to file---*/
            Boom_AdjointFile << scientific <<  Buffer_Recv_dJdU[Total_Index]   << "\t";
          }
          Boom_AdjointFile  << endl;

        }
      }

      delete [] Buffer_Recv_dJdU;
      delete [] Buffer_Recv_GlobalIndex;
    }

  delete [] Buffer_Send_dJdU;
  delete [] Buffer_Send_GlobalIndex;

  }

  if(rank == MASTER_NODE)
    Boom_AdjointFile.close();

  /*---Clear up  memory from dJdU---*/
  for(unsigned short iPhi = 0; iPhi < ray_N_phi; iPhi++){
    for (unsigned short i=0; i<nDim+3; i++){
      delete [] dJdU[iPhi][i];
    }
    delete [] dJdU[iPhi];
    delete [] PointID[iPhi];

  }
  delete [] dJdU;
  delete [] PointID;
  delete [] nPointID;

  if (rank == MASTER_NODE)
    cout << "\nFinished writing boom adjoint file." << endl;

}

