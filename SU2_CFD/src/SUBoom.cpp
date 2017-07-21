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

////  if (rank == MASTER_NODE){
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
  scale_L = config->GetRefLengthMoment();

  /*---Ray variables---*/
  if(nDim == 2) ray_N_phi = 1;
  else ray_N_phi = 1; // TODO: read in for 3D
  ray_phi = new su2double[ray_N_phi];
  for(unsigned int i = 0; i < ray_N_phi; i++){
    //ray_phi[i] = phi[i];
    ray_phi[i] = 0.0;
  }

  /*---Tolerance variables---*/
//  char cstr [200];
  string str;
  ifstream tolfile;

  tolfile.open("boom.in", ios::in);
  if (tolfile.fail()) {
    if(rank == MASTER_NODE)
      cout << "There is no boom.in file. Using default tolerances for boom propagation. " << endl;
    ray_r0 = 1.0;
    n_prof = 100001;
    tol_dphi = 1.0E-3;
    tol_dr = 1.0E-3;
    tol_m = 1.0E6;

    tol_dp = 1.0E-6;
    tol_l = 1.0E-4;
  }
  else{
    tolfile >> str >> ray_r0;
    tolfile >> str >> n_prof;
    tolfile >> str >> tol_dphi;
    tolfile >> str >> tol_dr;
    tolfile >> str >> tol_m;
    tolfile >> str >> tol_dp;
    tolfile >> str >> tol_l;
  }
////  }

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
  for(int i = 0; i < ray_N_phi; i++){
    if(startline[i]){
      ExtractLine(geometry, ray_r0, ray_phi, ray_N_phi);
    }
    else{
      nPanel = 0;
    }
  }

  /*---Interpolate pressures along line---*/
  if(rank == MASTER_NODE)
    cout << "Extract pressure signature." << endl;
  for(int i = 0; i < ray_N_phi; i++){
    if(startline[i]){      ExtractPressure(solver, config, geometry);
    }
  }

  unsigned long iPanel;
  for(iPanel = 0; iPanel < nPanel; iPanel++){
    signal.original_p[iPanel] = signal.original_p[iPanel]*Pressure_Ref - Pressure_FreeStream;
  }

  unsigned long totSig = 0;
  unsigned long maxSig = 0;
  unsigned long *nPanel_loc = new unsigned long[nProcessor];
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nPanel, &totSig, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&nPanel, &maxSig, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Gather(&nPanel, 1, MPI_UNSIGNED_LONG, nPanel_loc, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
#endif

  su2double* Buffer_Recv_Press = NULL;
  su2double* Buffer_Recv_x = NULL;
  su2double* Buffer_Send_Press = new su2double[maxSig];
  su2double* Buffer_Send_x = new su2double[maxSig];
  if(rank == MASTER_NODE){
    Buffer_Recv_Press = new su2double[nProcessor*maxSig];
    Buffer_Recv_x = new su2double[nProcessor*maxSig];
  }
  
  for(iPanel = 0; iPanel < nPanel; iPanel++){
    Buffer_Send_x[iPanel] = signal.x[iPanel];
    Buffer_Send_Press[iPanel] = signal.original_p[iPanel];
  }

#ifdef HAVE_MPI
  SU2_MPI::Gather(Buffer_Send_Press, maxSig, MPI_DOUBLE, Buffer_Recv_Press,  maxSig , MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_x, maxSig, MPI_DOUBLE, Buffer_Recv_x,  maxSig , MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#endif

  if (rank == MASTER_NODE)
    cout << "Gathered signal data to MASTER_NODE." << endl;

  if (rank == MASTER_NODE){
    ofstream sigFile;
    sigFile.precision(15);
    sigFile.open("signal_original.dat");
    sigFile << "# x, p" << endl;

  ////  unsigned long Total_Index;

  ////  for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
  ////    for (iPanel = 0; iPanel < Buffer_Recv_sigCount[iProcessor]; iPanel++) {
  ////      /*--- Current index position and global index ---*/
  ////      Total_Index  = iProcessor*totSig + iPanel;

  ////      signal.x[panelCount] = Buffer_Recv_x[Total_Index];
  ////      signal.original_p[panelCount] = Buffer_Recv_Press[Total_Index];

  ////      panelCount++;
  ////    }
  ////  }

    unsigned long panelCount = 0;
    nPanel = totSig;
    signal.original_len = nPanel;
    signal.x = new su2double[nPanel];
    signal.original_p = new su2double[nPanel];
    signal.original_T = new su2double[nPanel];
    for(unsigned int iProcessor = 0; iProcessor < nProcessor; iProcessor++){
      for(iPanel = 0; iPanel < nPanel_loc[iProcessor]; iPanel++){
        signal.x[panelCount] = Buffer_Recv_x[iProcessor*maxSig+iPanel];
        signal.original_p[panelCount] = Buffer_Recv_Press[iProcessor*maxSig+iPanel];
        panelCount++;
      }
    }

    /*---Sort signal in order of x-coordinate---*/
    cout << "Sorting signal data." << endl;
    MergeSort(signal.x, signal.original_p, 0, totSig-1);

    /*---Check for duplicate points---*/
    for(iPanel = 1; iPanel < nPanel; iPanel++){
      if(abs(signal.x[iPanel-1]-signal.x[iPanel]) < 1.0E-8){
        for(unsigned long jPanel = iPanel; jPanel < nPanel; jPanel++){
          signal.x[jPanel-1] = signal.x[jPanel];
          signal.original_p[jPanel-1] = signal.original_p[jPanel];
        }
        iPanel--;
        nPanel--;
      }
    }

    if(nPanel != totSig){
      cout << "Eliminating duplicate points." << endl;
      su2double *xtmp = new su2double[nPanel], *ptmp = new su2double[nPanel];
      for(iPanel = 0; iPanel < nPanel; iPanel++){
        xtmp[iPanel] = signal.x[iPanel];
        ptmp[iPanel] = signal.original_p[iPanel];
      }
      signal.x = new su2double[nPanel];
      signal.original_p = new su2double[nPanel];
      for(iPanel = 0; iPanel < nPanel; iPanel++){
        signal.x[iPanel] = xtmp[iPanel];
        signal.original_p[iPanel] = ptmp[iPanel];
      }
      delete [] xtmp;
      delete [] ptmp;
      totSig = nPanel;
    }

    /*---Now write to file---*/
    for(iPanel = 0; iPanel < nPanel; iPanel++){
      sigFile << scientific << signal.x[iPanel] << "\t";
      sigFile << scientific << signal.original_p[iPanel]   << "\t";
      sigFile << endl;
    }
    sigFile.close();
    cout << "Signal written. nPanel = " << nPanel << "." << endl;

  }

  /*---Initialize sensitivities---*/
  if(config->GetAD_Mode()){
    dJdU = new su2double* [nDim+3];
    for(int iDim = 0; iDim < nDim+3 ; iDim++){
      dJdU[iDim] = new su2double[nPointID];
      for(iPanel = 0;  iPanel< nPointID; iPanel++){
        dJdU[iDim][iPanel] = 0.0;
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
  unsigned long iMarker, iVertex, iPoint, *iPointmin, itmp;
  unsigned long jElem;
  unsigned short iElem, nElem, nNearest, *iPointmax, *ixmin;

  bool inside=false;

  su2double Minf = config->GetMach();
  su2double x = 0.0, y = 0.0, z = 0.0;
  su2double *Coord;
  su2double *y0, *z0;
  su2double yy, zz;
  su2double r2, r02 = r0*r0, *r2min, *xmin, mintmp;
  su2double sq2 = sqrt(2.0), sq2_2 = sq2/2.0;
  su2double mu = asin(1.0/Minf), cotmu = 1.0/tan(mu);
  su2double *p0, *p1;


  if(nDim == 2){
    nNearest = 4;
    y0 = new su2double[1];
    r2min = new su2double[nNearest];
    xmin = new su2double[nNearest];
    iPointmin = new unsigned long[nNearest];
    iPointmax = new unsigned short[1];
    ixmin = new unsigned short[1];
    startline = new bool[1];
    endline = new bool[1];
    p0 = new su2double[2];
    p1 = new su2double[2];
    y0[0] = r0;
    for(unsigned short j = 0; j < nNearest; j++){
      r2min[j] = 1.0E6;
      xmin[j] = 1.0E6;
    }
    iPointmax[0] = 0;
    ixmin[0] = 0;
    startline[0] = false;
    endline[0] = false;
  }
  else{
    nNearest = 8;
    y0 = new su2double[nPhi];
    z0 = new su2double[nPhi];
    r2min = new su2double[nNearest*nPhi];
    xmin = new su2double[nNearest*nPhi];
    iPointmin = new unsigned long[nNearest*nPhi];
    iPointmax = new unsigned short[nPhi];
    ixmin = new unsigned short[nPhi];
    startline = new bool[nPhi];
    endline = new bool[nPhi];
    p0 = new su2double[3];
    p1 = new su2double[3];
    for(int i = 0; i < nPhi; i++){
      y0[i] = r0*sin(phi[i]);
      z0[i] = r0*cos(phi[i]);
      for(unsigned short j = 0; j < nNearest; j++){
        r2min[nNearest*i+j] = 1.0E6;
        xmin[nNearest*i+j] = 1.0E6;
      }
      iPointmax[i] = 0;
      ixmin[i] = 0;
      startline[i] = false;
      endline[i] = false;
    }
  }

  for(iMarker = 0; iMarker < nMarker; iMarker++){
    /*--- Only look at farfield boundary (or send/recv for parallel computations) ---*/
    if(config->GetMarker_All_KindBC(iMarker) == FAR_FIELD || config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE){
      for(iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++){
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        /*--- Make sure point is in domain ---*/
        if(geometry->node[iPoint]->GetDomain()){
          Coord = geometry->node[iPoint]->GetCoord();
          y = SU2_TYPE::GetValue(Coord[1]);
          if(nDim == 3) z = SU2_TYPE::GetValue(Coord[2]);
          /*--- Only look at points below aircraft ---*/
          if((nDim == 2 && y < 0.0) || (nDim == 3 && z < 0.0)){
            r2 = y*y;
            if(nDim == 3) r2 += z*z;
            /*--- Limit search to points within strip ---*/
            if(r2 > r02*sq2_2 && r2 < r02*sq2){
              x = SU2_TYPE::GetValue(Coord[0]);

              if(nDim == 2){
                yy = abs(y);
                r2 = (yy-y0[0])*(yy-y0[0]);
                if(r2 < r2min[iPointmax[0]]){
                  iPointmin[iPointmax[0]] = iPoint;
                  r2min[iPointmax[0]] = r2;
                  xmin[iPointmax[0]] = x;
                  startline[0] = true;
                }
                for(unsigned short j = 0; j < nNearest; j++){
                  if(r2min[j] > r2min[iPointmax[0]]) iPointmax[0] = j;
                }
              }

              else{
                yy = abs(y);
                zz = abs(z);
                for(int i = 0; i < nPhi; i++){
                  r2 = (yy-y0[i])*(yy-y0[i]) + (zz-z0[i])*(zz-z0[i]);
                  if(r2 < r2min[2*i+1]){
                    iPointmin[2*i+1] = iPoint;
                    r2min[2*i+1] = r2;
                    xmin[2*i+1] = x;
                    startline[i] = true;
                    /*--- Now sort min values based on r2 ---*/
                    if(r2min[2*i+1] < r2min[2*i]){
                      itmp = iPointmin[2*i+1];
                      iPointmin[2*i+1] = iPointmin[2*i];
                      iPointmin[2*i] = itmp;

                      mintmp = r2min[2*i+1];
                      r2min[2*i+1] = r2min[2*i];
                      r2min[2*i] = mintmp;

                      mintmp = xmin[2*i+1];
                      xmin[2*i+1] = xmin[2*i];
                      xmin[2*i] = mintmp;
                    }
                  }
                }
              }

            }
          }
        }
      }
    }
  }

  /*--- Reorder iPointmin by x location ---*/
  for(int i = 0; i < nPhi; i++){
    MergeSort(xmin, iPointmin, i*nNearest, i*nNearest+nNearest-1);
  }

  if(nDim == 2){
    if(startline[0]){
      for(unsigned short iNearest = 0; iNearest < nNearest; iNearest++){
        Coord = geometry->node[iPointmin[iNearest]]->GetCoord();
        nElem = geometry->node[iPointmin[iNearest]]->GetnElem();
        for(iElem = 0; iElem < nElem; iElem++){
          jElem = geometry->node[iPointmin[iNearest]]->GetElem(iElem);
          inside = InsideElem(geometry, r0, 0.0, jElem, p0, p1);
          if(inside){
            nPanel = 1;
            pointID_original = new unsigned long[nPanel];
            Coord_original = new su2double*[nPanel];
            Coord_original[0] = new su2double[nDim];

            pointID_original[0] = jElem;
            Coord_original[0][0] = (p0[0] + p1[0])/2.0;
            Coord_original[0][1] = -r0;

            break;
          }
        }
        if(inside) break;
      }
      if(!inside) startline[0] = false;
    }
  }

  else if(nDim == 3){

  }
 

  delete [] iPointmin;
  delete [] Coord;
  delete [] y0;
  if(nDim == 3) delete [] z0;
  delete [] r2min;
  delete [] p0;
  delete [] p1;

}

void SUBoom::ExtractLine(CGeometry *geometry, const su2double r0, const su2double *phi, unsigned short nPhi){
  bool inside, boundary, end = false;
  unsigned short iElem, nElem;
  unsigned long jElem, jElem_m1, nElem_tot = geometry->GetnElem(), nPoint_tot = geometry->GetnPointDomain();
  su2double x_i, x_m1;

  unsigned long *pointID_tmp;
  su2double **Coord_tmp;
  su2double *p0 = new su2double[nDim], *p1 = new su2double[nDim];

  while(!end){
    if(nDim == 2){

      jElem_m1 = pointID_original[nPanel-1];
      x_m1 = geometry->elem[jElem_m1]->GetCG(0);
      nElem = geometry->elem[jElem_m1]->GetnNeighbor_Elements();
      inside = false;

      for(iElem = 0; iElem < nElem; iElem++){
        jElem = geometry->elem[jElem_m1]->GetNeighbor_Elements(iElem);
        /*--- Don't extract boundary elements ---*/
        boundary = false;
        /*for(unsigned short iPoint = 0; iPoint < geometry->elem[jElem]->GetnNodes(); iPoint++){
          unsigned long jPoint = geometry->elem[jElem]->GetNode(iPoint);
          //cout << "jPoint = " << jPoint << endl;//", nPoint_tot = " << nPoint_tot << endl;
          //if(jPoint < nPoint_tot) boundary = true;
          if(!geometry->node[jPoint]->GetDomain()){
            boundary = true;
            cout << "jPoint = " << jPoint << ", jElem = " << jElem << endl;
          }
        }*/
        //if(!boundary){
          if(jElem < nElem_tot){
            x_i = geometry->elem[jElem]->GetCG(0);

            if(x_i > x_m1){
              inside = InsideElem(geometry, r0, 0.0, jElem, p0, p1);
              if(inside){
                nPanel++;

                pointID_tmp = new unsigned long[nPanel-1];
                Coord_tmp = new su2double*[nPanel-1];
                for(unsigned long i = 0; i < nPanel-1; i++){
                  Coord_tmp[i] = new su2double[nDim];
                  pointID_tmp[i] = pointID_original[i];
                  Coord_tmp[i][0] = Coord_original[i][0];
                  Coord_tmp[i][1] = Coord_original[i][1];

                  delete [] Coord_original[i];
                }
                delete [] pointID_original;
                delete [] Coord_original;

                pointID_original = new unsigned long[nPanel];
                Coord_original = new su2double*[nPanel];
                for(unsigned long i = 0; i < nPanel-1; i++){
                  Coord_original[i] = new su2double[nDim];
                  pointID_original[i] = pointID_tmp[i];
                  Coord_original[i][0] = Coord_tmp[i][0];
                  Coord_original[i][1] = Coord_tmp[i][1];

                  delete [] Coord_tmp[i];
                }
                delete [] pointID_tmp;
                delete [] Coord_tmp;

                Coord_original[nPanel-1] = new su2double[nDim];
                pointID_original[nPanel-1] = jElem;
                Coord_original[nPanel-1][0] = (p0[0] + p1[0])/2.0;
                Coord_original[nPanel-1][1] = -r0;

                break;
              }
            }
          }
        //}
      }
      if(!inside){
        end = true;
      }
    }
  }

  cout << "nPanel extracted = " << nPanel << endl;

}

void SUBoom::ExtractPressure(CSolver *solver, CConfig *config, CGeometry *geometry){
  unsigned short iDim, iNode, nNode;
  unsigned long jElem, jNode;
  unsigned long pointCount = 0;
  su2double rho, rho_ux, rho_uy, rho_uz, rho_E, TKE;
  su2double rho_i, rho_ux_i, rho_uy_i, rho_uz_i, rho_E_i, TKE_i;
  su2double ux, uy, uz, StaticEnergy, p;

  su2double *isoparams;
  su2double *X_donor;
  su2double *Coord = new su2double[nDim];

  nPointID = 0;
  for(unsigned long i = 0; i < nPanel; i++){
    jElem = pointID_original[i];
    nPointID += geometry->elem[jElem]->GetnNodes();
  }
  PointID = new unsigned long[nPointID];

  signal.original_len = nPanel;
  signal.x = new su2double[nPanel];
  signal.original_p = new su2double[nPanel];
  for(unsigned long i = 0; i < nPanel; i++){
    /*--- Get info needed for isoparameter computation ---*/
    jElem = pointID_original[i];
    nNode = geometry->elem[jElem]->GetnNodes();
    for(unsigned short j = 0; j < nDim; j++){
      Coord[j] = Coord_original[i][j];
    }
    X_donor = new su2double[nDim*nNode];
    for(iNode = 0; iNode < nNode; iNode++){
      jNode = geometry->elem[jElem]->GetNode(iNode);
      for(iDim = 0; iDim < nDim; iDim++){  
        X_donor[iDim*nNode + iNode] = geometry->node[jNode]->GetCoord(iDim);
      }
    }

    /*--- Compute isoparameters ---*/
    isoparams = new su2double[nNode];
    Isoparameters(nDim, nNode, X_donor, Coord, isoparams);

    /*--- Now interpolate pressure ---*/
    p = 0.0;
    rho_i = 0.0;
    rho_ux_i = 0.0; rho_uy_i = 0.0; rho_uz_i = 0.0;
    rho_E_i = 0.0; TKE_i = 0.0;
    for(iNode = 0; iNode < nNode; iNode++){
      if(isoparams[iNode]*isoparams[iNode] > 0.0){
        jNode = geometry->elem[jElem]->GetNode(iNode);

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

        /*---Compute pressure---*/
        rho_i += rho*isoparams[iNode];
        rho_ux_i += rho_ux*isoparams[iNode];
        rho_uy_i += rho_uy*isoparams[iNode];
        if(nDim == 3) rho_uz_i += rho_uz*isoparams[iNode];
        rho_E_i += rho_E*isoparams[iNode];
        TKE_i += TKE*isoparams[iNode];

        /*ux = rho_ux/rho;
        uy = rho_uy/rho;
        uz = 0.0;
        if(nDim == 3) uz= rho_uz/rho;
        StaticEnergy =  rho_E/rho-0.5*(ux*ux+uy*uy+uz*uz)-TKE;
        p += (config->GetGamma()-1)*rho*StaticEnergy*isoparams[iNode];*/

        PointID[pointCount] = geometry->node[jNode]->GetGlobalIndex();
        pointCount++;
      }
    }
    
    ux = rho_ux_i/rho_i;
    uy = rho_uy_i/rho_i;
    uz = 0.0;
    if(nDim == 3) uz= rho_uz_i/rho_i;
    StaticEnergy =  rho_E_i/rho_i-0.5*(ux*ux+uy*uy+uz*uz)-TKE_i;
    p = (config->GetGamma()-1)*rho_i*StaticEnergy;

    signal.x[i] = Coord[0];
    signal.original_p[i] = p;

  }
  
}

bool SUBoom::InsideElem(CGeometry *geometry, su2double r0, su2double phi, unsigned long jElem, su2double *p0, su2double *p1){
  bool inside = false;
  unsigned long iPoint;
  unsigned short iNode, nNode, count, intersect;

  nNode = geometry->elem[jElem]->GetnNodes();

  su2double **Coord_elem = new su2double*[nNode];
  su2double *pp0 = new su2double[nDim];
  su2double *pp1 = new su2double[nDim];

  /*--- Store node coordinates ---*/
  for(iNode = 0; iNode < nNode; iNode++){
    iPoint = geometry->elem[jElem]->GetNode(iNode);
    Coord_elem[iNode] = new su2double[nDim];
    for(unsigned short iDim = 0; iDim < nDim; iDim++){
      Coord_elem[iNode][iDim] = geometry->node[iPoint]->GetCoord(iDim);
    }
  }

  /*--- Now determine if line intersects element ---*/
  if(nDim == 2){
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
      }
      count += intersect;
      if(count > 1){
        inside = true;
        break;
      }
    }
  }

  for(iNode = 0; iNode < nNode; iNode++){
    delete [] Coord_elem[iNode];
  }
  delete [] Coord_elem;

  return inside;
}

int SUBoom::Intersect2D(su2double r0, su2double *Coord_i, su2double *Coord_ip1, su2double *p0, su2double *p1){

  su2double line[2][2] = {{-1.0,-r0},{1.0E3,-r0}};
  su2double u[2] = {line[1][0]-line[0][0], line[1][1]-line[0][1]};
  su2double v[2] = {Coord_ip1[0]-Coord_i[0], Coord_ip1[1]-Coord_i[1]};
  su2double w[2] = {line[0][0]-Coord_i[0], line[0][1]-Coord_i[1]};
  su2double upv = u[0]*v[1] - u[1]*v[0], 
            upw = u[0]*w[1] - u[1]*w[0], 
            vpw = v[0]*w[1] - v[1]*w[0]; // Perp dot product

  if(abs(upv) < 1E-8){ // Segments are parallel
    if(abs(upw) > 1.0E-8 || abs(vpw) > 1.0E-8){ // Not colinear
      return 0;
    }
    /*--- Get overlap of collinear segments ---*/
    su2double t0, t1;
    su2double w2[2] = {line[1][0]-Coord_i[0], line[1][1]-Coord_i[1]};
    if(abs(v[0]) > 1.0E-8){
      t0 = w[0]/v[0];
      t1 = w2[0]/v[0];
    }
    else{
      t0 = w[1]/v[1];
      t1 = w2[1]/v[1];
    }
    if(t0 > t1){
      su2double t = t0;
      t0 = t1;
      t1 = t;
    }

    t0 = t0<0.0? 0.0 : t0; // Clip min to 0
    t1 = t1>1.0? 1.0 : t1; // Clip min to 1

    p0[0] = Coord_i[0] + t0*v[0];
    p0[1] = Coord_i[1] + t0*v[1];
    p1[0] = Coord_i[0] + t1*v[0];
    p1[1] = Coord_i[1] + t1*v[1];

    return 2;
  }

  /*--- Segments not parallel and may intersect at a point ---*/
  su2double s = vpw/upv, t = upw/upv;
  if(s < 0 || s > 1 || t < 0 || t > 1){
    return 0;
  }

  p0[0] = line[0][0] + s*u[0];
  p0[1] = line[0][1] + s*u[1];

  return 1;

}

int SUBoom::Intersect3D(){
  
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
    if (xi>1.0) xi=1.0;
    if (xi<-1.0) xi=-1.0;
    if (eta>1.0) eta=1.0;
    if (eta<-1.0) eta=-1.0;
    isoparams[0]=0.25*(1-xi)*(1-eta);
    isoparams[1]=0.25*(1+xi)*(1-eta);
    isoparams[2]=0.25*(1+xi)*(1+eta);
    isoparams[3]=0.25*(1-xi)*(1+eta);
  }
/*  if (nDonor<4) {
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
  }*/
  
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
  rho_of_z = new su2double[n_prof];
  for(unsigned int i = 0; i < n_prof; i++){
    z[i] = h0*su2double(i)/(su2double(n_prof)-1.);
    AtmosISA(z[i], T, a_of_z[i], p, rho_of_z[i], g);
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

void SUBoom::GetAtmosphericData(su2double& a, su2double& rho, su2double& p, su2double h){
  /*---Get speed of sound, density, and pressure at an altitude---*/
  int i_h = -1;
  for(int i = 0; i < n_prof; i++){
    if(h == z[i]){
      i_h = i;
      break;
    }
  }
  a = a_of_z[i_h];
  p = p_inf;
  rho = rho_of_z[i_h];
  //p = p_of_z[i_h];
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

    phi[i] = ray_phi[i]*deg2rad;
    phi_tube[i] = phi[i] + tol_dphi*deg2rad;
  }

  /*---Compute initial ray parameters---*/
  //GetAtmosphericData(a0, rho0, p0, h);
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

void SUBoom::RayTracer(){

  /*---Scale factors---*/
  su2double L = flt_h;
  su2double T = scale_T;
  su2double a0, rho0, p0;
  su2double a, rho, p;
  //  su2double theta[n_prof];
  su2double r0[3];
  su2double *f, *x, *y, *t;
  su2double *kx, *ky, *kz;
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
  x = new su2double[n_prof];
  y = new su2double[n_prof];
  t = new su2double[n_prof];
  kx = new su2double[n_prof];
  ky = new su2double[n_prof];
  kz = new su2double[n_prof];

  x_of_z = new su2double**[ray_N_phi];
  y_of_z = new su2double**[ray_N_phi];
  t_of_z = new su2double**[ray_N_phi];
  dxdt = new su2double**[ray_N_phi];
  dydt = new su2double**[ray_N_phi];
  dzdt = new su2double**[ray_N_phi];
  theta = new su2double**[ray_N_phi];
  for(unsigned int i = 0; i < ray_N_phi; i++){
    x_of_z[i] = new su2double*[4];
    y_of_z[i] = new su2double*[4];
    t_of_z[i] = new su2double*[4];
    dxdt[i] = new su2double*[4];
    dydt[i] = new su2double*[4];
    dzdt[i] = new su2double*[4];
    theta[i] = new su2double*[4];
    for(unsigned short j = 0; j < 4; j++){
      x_of_z[i][j] = new su2double[n_prof];
      y_of_z[i][j] = new su2double[n_prof];
      t_of_z[i][j] = new su2double[n_prof];
      dxdt[i][j] = new su2double[n_prof];
      dydt[i][j] = new su2double[n_prof];
      dzdt[i][j] = new su2double[n_prof];
      theta[i][j] = new su2double[n_prof];
    }
  }

   for(unsigned int i = 0; i < ray_N_phi; i++){

    /*---Primary ray---*/
    data.c0 = ray_c0[i][0];
    data.nu = ray_nu[i][0];
    r0[0] = r0[1] = r0[2] = 0.;

    for(int j = n_prof-1; j > -1; j--){
      f = rk4(z[j]/data.L, 3, r0, dz/data.L, data, derivs);
      x[j] = x_of_z[i][0][j] = -f[0];
      y[j] = y_of_z[i][0][j] = -f[1];
      t[j] = t_of_z[i][0][j] = -f[2];

      r0[0] = -x[j];
      r0[1] = -y[j];
      r0[2] = -t[j];

      /*---Derived data---*/
      GetAtmosphericData(a, rho, p, z[j]);
      theta[i][0][j] = acos(a/data.c0);
    }

    kx = SplineGetDerivs(t, x, n_prof);
    ky = SplineGetDerivs(t, y, n_prof);
    kz = SplineGetDerivs(t, z, n_prof);
    for(unsigned int ik = 0; ik < n_prof; ik++){
      dxdt[i][0][ik] = kx[ik];
      dydt[i][0][ik] = ky[ik];
      dzdt[i][0][ik] = kz[ik];
    }

    /*---Ray tube corners: {0, +dheading}---*/
    r0[0] = flt_heading[0]*tol_dr;
    r0[1] = flt_heading[1]*tol_dr;
    r0[2] = flt_heading[2]*tol_dr;
    for(int j = n_prof-1; j > -1; j--){
      f = rk4(z[j]/data.L, 3, r0, dz/data.L, data, derivs);
      x[j] = x_of_z[i][1][j] = -f[0];
      y[j] = y_of_z[i][1][j] = -f[1];
      t[j] = -f[2];

      r0[0] = -x[j];
      r0[1] = -y[j];
      r0[2] = -t[j];

      /*---Derived data---*/
      GetAtmosphericData(a, rho, p, z[j]);
      theta[i][1][j] = acos(a/data.c0);
    }

    kx = SplineGetDerivs(t, x, n_prof);
    ky = SplineGetDerivs(t, y, n_prof);
    for(unsigned int ik = 0; ik < n_prof; ik++){
      dxdt[i][1][ik] = kx[ik];
      dydt[i][1][ik] = ky[ik];
    }

    /*---Ray tube corners: {+dphi, 0}---*/
    data.c0 = ray_c0[i][1];
    data.nu = ray_nu[i][1];
    r0[0] = r0[1] = r0[2] = 0.;

    for(int j = n_prof-1; j > -1; j--){
      f = rk4(z[j]/data.L, 3, r0, dz/data.L, data, derivs);
      x[j] = x_of_z[i][2][j] = -f[0];
      y[j] = y_of_z[i][2][j] = -f[1];
      t[j] = t_of_z[i][2][j] = -f[2];

      r0[0] = -x[j];
      r0[1] = -y[j];
      r0[2] = -t[j];

      /*---Derived data---*/
      GetAtmosphericData(a, rho, p, z[j]);
      theta[i][2][j] = acos(a/data.c0);
    }

    kx = SplineGetDerivs(t, x, n_prof);
    ky = SplineGetDerivs(t, y, n_prof);
    kz = SplineGetDerivs(t, z, n_prof);
    for(unsigned int ik = 0; ik < n_prof; ik++){
      dxdt[i][2][ik] = kx[ik];
      dydt[i][2][ik] = ky[ik];
      dzdt[i][2][ik] = kz[ik];
    }

    /*---Ray tube corners: {+dphi, +dheading}---*/
    r0[0] = flt_heading[0]*tol_dr;
    r0[1] = flt_heading[1]*tol_dr;
    r0[2] = flt_heading[2]*tol_dr;
    for(int j = n_prof-1; j > -1; j--){
      f = rk4(z[j]/data.L, 3, r0, dz/data.L, data, derivs);
      x[j] = x_of_z[i][3][j] = -f[0];
      y[j] = y_of_z[i][3][j] = -f[1];
      t[j] = -f[2];

      r0[0] = -x[j];
      r0[1] = -y[j];
      r0[2] = -t[j];

      /*---Derived data---*/
      GetAtmosphericData(a, rho, p, z[j]);
      theta[i][3][j] = acos(a/data.c0);
    }

    kx = SplineGetDerivs(t, x, n_prof);
    ky = SplineGetDerivs(t, y, n_prof);
    for(unsigned int ik = 0; ik < n_prof; ik++){
      dxdt[i][3][ik] = kx[ik];
      dydt[i][3][ik] = ky[ik];
    }

  }

  cout << "Clearing up ray tracer memory" << endl;

  /*---Clear up memory---*/
  delete [] f;
  delete [] x;
  delete [] y;
  delete [] t;
  delete [] kx;
  delete [] ky;
  delete [] kz;

}

void SUBoom::RayTubeArea(){

  su2double Ah, x_int, y_int, z_int;
  su2double corners[4][3];
  int M;

  ray_A = new su2double*[ray_N_phi];
  for(int i = 0; i < ray_N_phi; i++){
    ray_A[i] = new su2double[n_prof];
  }

  /*---Loop over rays---*/
  for(int i = 0; i < ray_N_phi; i++){
    M = n_prof;
    Ah = 0;
    //ray_A[i][M-1] = 0.0;
    for(int j = 0; j < M; j++){
      for(int k = 0; k < 4; k++){
	    for(int kk = 0; kk < 3; kk++){
	    corners[k][kk] = 0.0;
	    }
      }
      z_int = z[j]/scale_z;
      for(int k = 0; k < 4; k++){
    	  x_int = x_of_z[i][k][j];
	      y_int = y_of_z[i][k][j];
	      corners[k][0] = x_int;
	      corners[k][1] = y_int;
	      corners[k][2] = z_int;
      }
      su2double u[3] = {corners[3][0]-corners[0][0], corners[3][1]-corners[0][1], corners[3][2]-corners[0][2]};
      su2double v[3] = {corners[2][0]-corners[1][0], corners[2][1]-corners[1][1], corners[2][2]-corners[1][2]};
      /*---Cross product---*/
      su2double c[3] = {u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0]};
      Ah = 0.5*sqrt(pow(c[0],2)+pow(c[1],2)+pow(c[2],2));
      ray_A[i][j] = Ah*a_of_z[j]*tan(theta[i][0][j])/ray_c0[i][0];
    }
  }

}

su2double SUBoom::matchr(int i, int j, su2double h_L, su2double r0){
  su2double f;
  f = sqrt(pow(x_of_z[i][0][j],2) + pow(y_of_z[i][0][j],2) + pow(z[j]/scale_z-1.0,2))*h_L - r0;
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

  ray_t0 = new su2double[ray_N_phi];

  ks = new su2double[n_prof];

  for(unsigned int i = 0; i < ray_N_phi; i++){
    for(unsigned int j = 0; j < n_prof; j++){
      t[j] = t_of_z[i][0][j];
      f[j] = matchr(i, j, h_L, ray_r0);
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

    ray_t0[i] = t1;

  }

  /*---Clear up memory---*/
  delete [] ks;

}

void SUBoom::ODETerms(){

  su2double t[n_prof], A[n_prof];
  su2double *cn;
  su2double *dadt, *drhodt, *dAdt;
  su2double *dcndt;
  su2double g = atm_g;

  cn = new su2double[n_prof];

  dadt = new su2double[n_prof];
  drhodt = new su2double[n_prof];
  dAdt = new su2double[n_prof];
  dcndt = new su2double[n_prof];

  ray_C1 = new su2double*[ray_N_phi];
  ray_C2 = new su2double*[ray_N_phi];
  ray_dC1 = new su2double*[ray_N_phi];
  ray_dC2 = new su2double*[ray_N_phi];
  for(unsigned int i = 0; i < ray_N_phi; i++){
    for (unsigned int j = 0; j < n_prof; j++){
      t[j] = t_of_z[i][0][j];
      A[j] = ray_A[i][j];
      cn[j] = ray_c0[i][0]*cos(theta[i][0][j]);
    }

    /*--- Spline interpolation of a, rho, A---*/
    dadt = SplineGetDerivs(t, a_of_z, n_prof);
    drhodt = SplineGetDerivs(t, rho_of_z, n_prof);
    dAdt = SplineGetDerivs(t, A, n_prof);
    dcndt = SplineGetDerivs(t, cn, n_prof);

    ray_C1[i] = new su2double[n_prof];
    ray_C2[i] = new su2double[n_prof];
    for(unsigned int j = 0; j < n_prof; j++){
      ray_C1[i][j] = ((g+1.)/(2.*g))*a_of_z[j]/cn[j];
      //if(A[j] > 1E-16){
          ray_C2[i][j] = 0.5*((3./a_of_z[j])*dadt[j] + drhodt[j]/rho_of_z[j] - (2./cn[j])*dcndt[j] - (1./A[j])*dAdt[j]);
      //}
      //else{
        //ray_C2[i][j] = 0.5*((3./a_of_z[j])*dadt[j] + drhodt[j]/rho_of_z[j] - (2./cn[j])*dcndt[j]);
      //}
      //cout << "A = " << A[j] << ", dAdt = " << dAdt[j] << ", C2 = " << ray_C2[i][j] << endl;
    }
    ray_dC1[i] = SplineGetDerivs(t, ray_C1[i], n_prof);
    ray_dC2[i] = SplineGetDerivs(t, ray_C2[i], n_prof);

    delete [] ray_c0[i];
    delete [] ray_nu[i];
    delete [] ray_theta0[i];
  }

  delete [] cn;
  delete [] dadt;
  delete [] drhodt;
  delete [] dAdt;
  delete [] dcndt;
  delete [] a_of_z;
  delete [] rho_of_z;
  delete [] ray_c0;
  delete [] ray_nu;
  delete [] ray_theta0;

}

void SUBoom::DistanceToTime(){
  int len = signal.original_len;
  for(int i = 0; i < len; i++){
    signal.original_T[i] = signal.x[i]/(a_inf*flt_M);
  }
}

void SUBoom::CreateSignature(){
  int len = signal.original_len;
  su2double pp[2][len-1];
  su2double ll[len-1];
  su2double mm[len-1];

  /*---Pressure signal to waveform---*/
  /*---Collect data into segments---*/
  for(int i = 0; i < len-1; i++){
    pp[0][i] = signal.original_p[i]/scale_p;  // [Pa] (relative)
    pp[1][i] = signal.original_p[i+1]/scale_p;
    ll[i] = (signal.original_T[i+1] - signal.original_T[i])/scale_T;  // [s]
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
      i -= 1;
      M -= 1;
    }
    else if(mm[i] < -tol_m/scale_m){  // "expansion shock" present
      /*---Remove segment i---*/
      ll[i] = tol_l/scale_T;
      mm[i] = (pp[1][i] - pp[0][i])/ll[i];
    }
    i += 1;
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
        i = i+1;
      }

      sig[0][i] = TT[1][seg];
      sig[1][i] = pp[1][seg];
      i = i+1;
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
  C1 = EvaluateSpline(t,9,data.t,data.C1,data.dC1);
  C2 = EvaluateSpline(t,9,data.t,data.C2,data.dC2);

  for(int i = 0; i < 3*M; i++){
    if(i < M){
      dydt[i] = C1*pow(y[i],2) + C2*y[i];
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
  for(int i = 0; i < 2; i++){current_signal[i] = new su2double[M+1];}

  /*---Decompose f vector---*/
  for(int j = 0; j < 3*M; j++){
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
    if(l[i] <= tol_l || m[i] >= tol_m || m[i] <= -tol_m){
      /*---Record pressure gap---*/
      current_signal = WaveformToPressureSignal(fvec, N, Msig);
      dp_seg = dp[i] + (current_signal[1][i] - current_signal[0][i]);
      /*---Add to next segment if needed---*/
      if(dp_seg > tol_dp){
          if(i < N-1){
              dp[i+1] = dp[i+1] + dp_seg;
          }
      }

      N -= 1;
      /*---Remove segment---*/
      for(int j = i; j < N; j++){
          m[j] = m[j+1];
          dp[j] = dp[j+1];
          l[j] = l[j+1];
      }
      //i -= 1;
      fvec_new = new su2double[3*N];
      for(int j = 0; j < N; j++){
          fvec[j] = m[j];
          fvec[j+N] = dp[j];
          fvec[j+2*N] = l[j];
      }

    }
    i += 1;

  }

  fvec_new = new su2double[3*N];
  for(int j = 0; j < N; j++){
    fvec_new[j] = m[j];
    fvec_new[j+N] = dp[j];
    fvec_new[j+2*N] = l[j];
  }

  M = N;

  /*---Free memory---*/
  for(int i = 0; i < 2; i++){
    delete [] current_signal[i];
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

void SUBoom::PropagateSignal(){
  su2double t0, tf, dt;
  unsigned int j0;
  RayData data;
  su2double *fvec;
  su2double *f;
  su2double **ground_signal;
  int M = signal.M;

  for(unsigned int i = 0; i < ray_N_phi; i++){
    t0 = ray_t0[i];

    /*---Assemble f vector and ray data for integration---*/
    signal.fvec = new su2double[3*signal.M];
    for(unsigned int j = 0; j < 3*signal.M; j++){
      if(j < signal.M){
	      signal.fvec[j] = signal.m[j];
      }
      else if(j < 2*signal.M){
	      signal.fvec[j] = signal.dp[j-signal.M];
      }
      else{
	      signal.fvec[j] = signal.l[j-2*signal.M];
      }
    }
    data.t = new su2double[9];
    data.C1 = new su2double[9];
    data.C2 = new su2double[9];
    data.dC1 = new su2double[9];
    data.dC2 = new su2double[9];

    data.M = signal.M;
    data.scale_C1 = scale_C1;
    data.scale_C2 = scale_C2;

    fvec = signal.fvec;
    for(int j = n_prof-1; j > -1; j--){
        if(t_of_z[i][0][j] >= ray_t0[i]){
            j0 = j+1;
            break;
        }
    }
    for(unsigned int j = j0; j > 0; j--){
      /*---Integrate---*/
      int j0_dat;
      if(j >= n_prof-4) j0_dat = n_prof-4;
      else if(j < 4) j0_dat = 4;
      else j0_dat = j;
      for(unsigned short jj = 0; jj < 9; jj++){
          data.t[jj] = t_of_z[i][0][j0_dat-4+jj];
          data.C1[jj] = ray_C1[i][j0_dat-4+jj];
          data.C2[jj] = ray_C2[i][j0_dat-4+jj];
          data.dC1[jj] = ray_dC1[i][j0_dat-4+jj];
          data.dC2[jj] = ray_dC2[i][j0_dat-4+jj];
      }
      tf = t_of_z[i][0][j-1];
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

    signal.final_T = new su2double[Msig];
    signal.final_p = new su2double[Msig];
    signal.final_M = Msig;

    /*---Final signal and boom strength---*/
    ofstream sigFile;
    sigFile.precision(15);
    sigFile.open("signal_final.dat", ios::out);
    sigFile << "# T, p" << endl;
    p_max = -1E10;
    p_int2 = 0;
    for(int j = 0; j < Msig; j++){
      signal.final_T[j] = ground_signal[0][j]*scale_T;
      signal.final_p[j] = ground_signal[1][j]*scale_p;
      sigFile << signal.final_T[j] << "\t" << signal.final_p[j] << endl;
      if(signal.final_p[j] > p_max) p_max = signal.final_p[j];
      if(j > 0) p_int2 = p_int2 + 0.5*(signal.final_p[j]*signal.final_p[j]+signal.final_p[j-1]*signal.final_p[j-1])
                        *(signal.final_T[j]-signal.final_T[j-1]);
    }
    sigFile.close();
    p_rise = signal.final_p[0];
    if(signal.final_p[0] > -signal.final_p[Msig-1]) p_rise2 = signal.final_p[0];
    else p_rise2 = signal.final_p[Msig-1];
    cout << "p_rise = " << p_rise << ", p_max = " << p_max << ", p_int2 = " << p_int2 << endl;

    /*---Write boom strength metrics to file---*/
    sigFile.open("pboomSU2", ios::out);
    sigFile << p_max << "," << p_rise << "," << p_rise2 << "," << p_int2 << endl;
    sigFile.close();
  }

  /*---Clean up---*/
  for(unsigned int i = 0; i < ray_N_phi; i++){
  for(unsigned short j = 0; j < 4; j++){
      delete [] x_of_z[i][j];
      delete [] y_of_z[i][j];
      delete [] t_of_z[i][j];
      delete [] dxdt[i][j];
      delete [] dydt[i][j];
      delete [] dzdt[i][j];
      delete [] theta[i][j];
  }
  delete [] x_of_z[i];
  delete [] y_of_z[i];
  delete [] t_of_z[i];
  delete [] dxdt[i];
  delete [] dydt[i];
  delete [] dzdt[i];
  delete [] theta[i];
  delete [] ray_A[i];
  delete [] ray_C1[i];
  delete [] ray_C2[i];
  delete [] ray_dC1[i];
  delete [] ray_dC2[i];
  }

  delete [] x_of_z;
  delete [] y_of_z;
  delete [] t_of_z;
  delete [] dxdt;
  delete [] dydt;
  delete [] dzdt;
  delete [] theta;
  delete [] ray_A;
  delete [] ray_C1;
  delete [] ray_C2;
  delete [] ray_dC1;
  delete [] ray_dC2;
  delete [] ray_t0;

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

  Buffer_Send_nPointID[0]=nPointID; 
#ifdef HAVE_MPI
  SU2_MPI::Gather(&Buffer_Send_nPointID, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nPointID, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD); //send the number of vertices at each process to the master
  SU2_MPI::Allreduce(&nPointID,&Max_nPointID,1,MPI_UNSIGNED_LONG,MPI_MAX,MPI_COMM_WORLD); //find the max num of vertices over all processes
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
      for(iSig=0; iSig<nPointID; iSig++){
          Buffer_Send_dJdU[iVar*nPointID+iSig] = dJdU[iVar][iSig];
        }
    }

  for (iSig=0; iSig<nPointID; iSig++){
     Buffer_Send_GlobalIndex[iSig] = PointID[iSig];
  }

#ifdef HAVE_MPI
  SU2_MPI::Gather(Buffer_Send_dJdU, Max_nPointID*nVar, MPI_DOUBLE, Buffer_Recv_dJdU,  Max_nPointID*nVar , MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_GlobalIndex,Max_nPointID, MPI_UNSIGNED_LONG, Buffer_Recv_GlobalIndex, Max_nPointID , MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
#endif

  if (rank == MASTER_NODE){
  Boom_AdjointFile.precision(15);
  Boom_AdjointFile.open("Adj_Boom.dat", ios::out);

  /*--- Loop through all of the collected data and write each node's values ---*/
  for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
    for (iSig = 0; iSig < Buffer_Recv_nPointID[iProcessor]; iSig++) {
        Global_Index = Buffer_Recv_GlobalIndex[iProcessor*Max_nPointID+iSig];
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

  Boom_AdjointFile.close();
  delete [] Buffer_Recv_dJdU;
  delete [] Buffer_Recv_GlobalIndex;
  }

  delete [] Buffer_Send_dJdU;
  delete [] Buffer_Send_GlobalIndex;

  /*---Clear up  memory from dJdU---*/
  for (unsigned short i=0; i<nDim+3; i++){
    delete [] dJdU[i];
  }
  delete [] dJdU;
  delete [] PointID;

  if (rank == MASTER_NODE)
    cout << "\nFinished writing boom adjoint file." << endl;

}

searchTree::searchTree(){

}

searchTree::searchTree(unsigned short nDim, unsigned long nPoints, const su2double **coord,
             const su2double Minf, const su2double r0, const su2double *phi, unsigned short nPhi){

  su2double x = 0.0, y = 0.0, z = 0.0;
  su2double r2 = 0.0, r02 = r0*r0;
  su2double sq2 = sqrt(2.0), sq2_2 = sq2/2.0;
  su2double mu = asin(1.0/Minf), cotmu = 1.0/tan(mu);

  unsigned long *ind_tmp;
  su2double *x_tmp, *y_tmp, *z_tmp;

  nDimTree = nDim - 1;
  nPointsTree = 0;

  treeNode nodes;
  nodes.ind = new unsigned long[1];
  nodes.x   = new su2double[1];
  nodes.y   = new su2double[1];
  if(nDimTree == 2) nodes.z = new su2double[1];

  /*--- Extract points and determine number of nodes in tree ---*/
  for(int i = 0; i < nPoints; i++){
    x = coord[i][0];
    y = coord[i][1];
    if(nDimTree == 2) z = coord[i][2];
    r2 = y*y + z*z;
    /*--- Check if point is within radius range, below aircraft, ahead of Mach angle ---*/
    if(r2 > sq2_2*r02 && r2 < sq2*r02 && ((nDimTree == 1 && y < 0 && x > -sq2_2*y*cotmu) || (nDimTree == 2 && z < 0 && x > -sq2_2*z*cotmu))){
      /*--- Modify arrays ---*/
      nPointsTree++;
      if(nPointsTree == 1){
        nodes.ind[0] = i;
        nodes.x[0] = x;
        nodes.y[0] = y;
        if(nDimTree == 2) nodes.z[0] = z;
      }
      else{
        ind_tmp = new unsigned long[nPointsTree-1];
        x_tmp   = new su2double[nPointsTree-1];
        y_tmp   = new su2double[nPointsTree-1];
        if(nDimTree == 2) z_tmp = new su2double[nPointsTree-1];
        
        for(int j = 0; j < nPointsTree-1; j++){
          ind_tmp[j] = nodes.ind[j];
          x_tmp[j] = nodes.x[j];
          y_tmp[j] = nodes.y[j];
          if(nDimTree == 2) z_tmp[j] = nodes.z[j];
        }

        nodes.ind = new unsigned long[nPointsTree];
        nodes.x   = new su2double[nPointsTree];
        nodes.y   = new su2double[nPointsTree];
        if(nDimTree == 2) nodes.z = new su2double[nPointsTree];
        
        for(int j = 0; j < nPointsTree-1; j++){
          nodes.ind[j] = ind_tmp[j];
          nodes.x[j] = x_tmp[j];
          nodes.y[j] = y_tmp[j];
          if(nDimTree == 2) nodes.z[j] = z_tmp[j];
        }

        nodes.ind[nPointsTree-1] = i;
        nodes.x[nPointsTree-1] = x;
        nodes.y[nPointsTree-1] = y;
        if(nDimTree == 2) nodes.z[nPointsTree-1] = z;

      }

    }
  }

  /*--- Sort nodes by location ---*/
  if(nPointsTree > 0){
    /*--- First create arrays for parents, children, and (r,phi) in 2D search ---*/
    nodes.root = -1;
    nodes.parent = new long[nPointsTree];
    nodes.left   = new long[nPointsTree];
    nodes.right  = new long[nPointsTree];
    if(nDimTree == 2){
      nodes.r   = new su2double[nPointsTree];
      nodes.phi = new su2double[nPointsTree];
    }
    for(int i = 0; i < nPointsTree; i++){
      nodes.parent[i] = -1;
      nodes.left[i]   = -1;
      nodes.right[i]  = -1;
      if(nDimTree == 2){
        nodes.r[i]   = sqrt(nodes.y[i]*nodes.y[i] + nodes.z[i]*nodes.z[i]);
        nodes.phi[i] = abs(atan(nodes.y[i]/nodes.z[i]));  // Only care about magnitude of angle
      }
    }

    /*--- Now sort arrays ---*/
    MergeSort(nodes, 0, nPointsTree-1, nDimTree);
  }

  /*--- Now build tree ---*/
  if(nDimTree == 1 && nPointsTree > 0){
    buildRBTree();
  }
  else if(nDimTree == 2 && nPointsTree > 0){
    buildQuadTree();
  }

  /*--- Search the tree ---*/
  if(nDimTree == 1){
    int k = 1000;
    count = new int[1];
    count[0] = 0;
    kNN = new int*[1];
    yNN = new su2double*[1];
    kNN[0] = new int[k];
    yNN[0] = new su2double[k];
    for(int i = 0; i < k; i++){
      kNN[0][i] = -1;
    }
    kNNRB(nodes.root, k, r0);
  }
  else if(nDimTree == 2){
    kNNQuad();
  }

  // TODO: Extract points along line (determine intersection)

}

void MergeSort(searchTree::treeNode nodes, int l, int r, unsigned short nDimTree){
  if(l < r){
    int m = l + (r-l)/2;
    MergeSort(nodes, l, m, nDimTree);
    MergeSort(nodes, m+1, r, nDimTree);
    merge(nodes, l, m, r, nDimTree);
  }
}

void merge(searchTree::treeNode nodes, int l, int m, int r, unsigned short nDimTree){
  int i, j, k;
  int n1 = m-l+1;
  int n2 = r-m;
  int n3;
  if(nDimTree == 1) n3 = 2;
  else n3 = 5;

  int Li[n1], Ri[n2];
  su2double **L, **R;
  L = new su2double*[n3];
  R = new su2double*[n3];
  for(int i = 0; i < n3; i++){
    L[i] = new su2double[n1];
    R[i] = new su2double[n2];
  }

  /*--- Temporary arrays ---*/
  for(i = 0; i < n1; i++){
    Li[i] = nodes.ind[l+i];
    L[0][i] = nodes.x[l+i];
    L[1][i] = nodes.y[l+i];
    if(nDimTree == 2){
      L[2][i] = nodes.z[l+i];
      L[3][i] = nodes.r[l+i];
      L[4][i] = nodes.phi[l+i];
    }
  }
  for(j = 0; j < n2; j++){
    Ri[j] = nodes.ind[m+1+j];
    R[0][j] = nodes.x[m+1+j];
    R[1][j] = nodes.y[m+1+j];
    if(nDimTree == 2){
      R[2][j] = nodes.z[m+1+j];
      R[3][j] = nodes.r[m+1+j];
      R[4][j] = nodes.phi[m+1+j];
    }
  }

  i = 0;
  j = 0;
  k = l;

  /*--- Begin sorting ---*/
  while(i < n1 && j < n2){
    if((nDimTree == 1 && L[1][i] <= R[1][j]) || (nDimTree == 2 && L[4][i] < R[4][j])){
      nodes.ind[k] = Li[i];
      nodes.x[k]   = L[0][i];
      nodes.y[k]   = L[1][i];
      if(nDimTree == 2){
        nodes.z[k]   = L[2][i];
        nodes.r[k]   = L[3][i];
        nodes.phi[k] = L[4][i];
      }
      i++;
    }
    else{
      nodes.ind[k] = Ri[j];
      nodes.x[k]   = R[0][j];
      nodes.y[k]   = R[1][j];
      if(nDimTree == 2){
        nodes.z[k]   = R[2][j];
        nodes.r[k]   = R[3][j];
        nodes.phi[k] = R[4][j];
      }
      j++;
    }
    k++;
  }

  /*--- Fill rest of arrays ---*/
  while(i < n1){
    nodes.ind[k] = Li[i];
    nodes.x[k]   = L[0][i];
    nodes.y[k]   = L[1][i];
    if(nDimTree == 2){
      nodes.z[k]   = L[2][i];
      nodes.r[k]   = L[3][i];
      nodes.phi[k] = L[4][i];
    }
    i++;
    k++;
  }

  while(j < n2){
    nodes.ind[k] = Ri[j];
    nodes.x[k]   = R[0][j];
    nodes.y[k]   = R[1][j];
    if(nDimTree == 2){
      nodes.z[k]   = R[2][j];
      nodes.r[k]   = R[3][j];
      nodes.phi[k] = R[4][j];
    }
    j++;
    k++;
  }
}

void searchTree::buildRBTree(){
  // TODO: Recursive build
  for(int i = 0; i < nPointsTree; i++){
    insertRB(i, nodes.root);
    insertRB_1(i);
  }
}

void searchTree::insertRB(long i, long j){
  if(nodes.root == -1)
    nodes.root = i;

  if(j == -1)
    return;
  else if(nodes.y[i] <= nodes.y[j]){
    nodes.parent[i] = j;
    insertRB(i,nodes.left[j]);
    if(nodes.left[j] == -1) nodes.left[j] = i;
  }
  else{
    nodes.parent[i] = j;
    insertRB(i,nodes.right[j]);
    if(nodes.right[j] == -1) nodes.right[j] = i;
  }
}

void searchTree::insertRB_1(long i){
  if(nodes.parent[i] == -1)
    nodes.color[i] = 0;
  else
    insertRB_2(i);
}

void searchTree::insertRB_2(long i){
  if(nodes.parent[i] == 0)
    return;
  else
    insertRB_3(i);
}

void searchTree::insertRB_3(long i){
  long u, g; // Uncle and grandparent
  if(nodes.parent[i] != -1)
    g = nodes.parent[nodes.parent[i]];
  else
    g = -1;

  if(g == -1)
    u = -1;
  else if(nodes.parent[i] == nodes.left[g])
    u = nodes.right[g];
  else
    u = nodes.left[g];

  if(u != -1 && nodes.color[u] == 1){
    nodes.color[nodes.parent[i]] = 0;
    nodes.color[u] = 0;
    nodes.color[g] = 1;
    insertRB_1(g);
  }
  else
    insertRB_4(i);
}

void searchTree::insertRB_4(long i){
  long g; // Grandparent
  if(nodes.parent[i] != -1)
    g = nodes.parent[nodes.parent[i]];
  else
    g = -1;

  if(g != -1){
    if(i == nodes.right[nodes.parent[i]] && nodes.parent[i] == nodes.left[g]){
      rotateLeft(i);
      i = nodes.left[i];
    }

    else if(i == nodes.left[nodes.parent[i]] && nodes.parent[i] == nodes.right[g]){
      rotateRight(i);
      i = nodes.right[i];
    }
  }

  insertRB_5(i);
}

void searchTree::insertRB_5(long i){
  long g; // Grandparent
  if(nodes.parent[i] != -1)
    g = nodes.parent[nodes.parent[i]];
  else
    g = -1;

  nodes.color[nodes.parent[i]] = 0;
  if(g != -1) nodes.color[g] = 1;
  if(i == nodes.left[nodes.parent[i]])
    rotateRight(nodes.parent[i]);
  else
    rotateLeft(nodes.parent[i]);
}

void searchTree::rotateLeft(long i){
  long saved_p, saved_left;
  saved_p = nodes.parent[i];
  saved_left = nodes.left[i];
  nodes.right[saved_p] = saved_left;
  
  nodes.parent[i] = nodes.parent[saved_p];
  
  if(nodes.parent[saved_p] == -1)
    nodes.root = i;
  else if(saved_p == nodes.left[nodes.parent[saved_p]])
    nodes.left[nodes.parent[saved_p]] = i;
  else
    nodes.right[nodes.parent[saved_p]] = i;

  nodes.left[i] = saved_p;
  nodes.parent[saved_p] = i;
}


void searchTree::rotateRight(long i){
  long saved_p, saved_right;
  saved_p = nodes.parent[i];
  saved_right = nodes.right[i];
  nodes.left[saved_p] = saved_right;
  
  nodes.parent[i] = nodes.parent[saved_p];
  
  if(nodes.parent[saved_p] == -1)
    nodes.root = i;
  else if(saved_p == nodes.right[nodes.parent[saved_p]])
    nodes.right[nodes.parent[saved_p]] = i;
  else 
    nodes.left[nodes.parent[saved_p]] = i;

  nodes.right[i] = saved_p;
  nodes.parent[saved_p] = i;
}

void searchTree::kNNRB(long i, int k, const su2double r0){
  if(count[0] < k){
    kNN[0][count[0]] = i;
    yNN[0][count[0]] = abs(r0-nodes.y[i]);
    count++;
//    MergeSort(yNN[0], kNN[0], 0, count[0]-1);
  }
  else{
    if(abs(r0-nodes.y[i]) < abs(r0-yNN[0][count[0]-1])){
      kNN[0][count[0]] = i;
      yNN[0][count[0]] = abs(r0-nodes.y[i]);
//      MergeSort(yNN[0], kNN[0], 0, count[0]-1);
    }
  }

  if(r0 <= nodes.y[i] && nodes.left[i] != -1)
    kNNRB(nodes.left[i], k, r0);
  else if(nodes.right[i] != -1)
    kNNRB(nodes.right[i], k, r0);
}

void searchTree::buildQuadTree(){

}

void searchTree::kNNQuad(){

}