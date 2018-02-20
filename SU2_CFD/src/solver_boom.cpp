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

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
#endif

  /*---Make sure to read in hard-coded values in the future!---*/

  nDim = geometry->GetnDim();

  /*---Flight variables---*/
  flt_h = config->GetBoom_flt_h(); // altitude [m]
  flt_M = config->GetMach();
  flt_psi = 0.;  // heading angle [deg]
  flt_gamma = 0.; // flight path angle [deg]

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
  signal.len = new unsigned long[ray_N_phi];
  signal.x = new su2double*[ray_N_phi];
  signal.p_prime = new su2double*[ray_N_phi];

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
      signal.p_prime[iPhi][iPanel] = signal.p_prime[iPhi][iPanel]*Pressure_Ref - Pressure_FreeStream;
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
      Buffer_Send_Press[iPanel] = signal.p_prime[iPhi][iPanel];
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
        if(abs(signal.x[iPhi][iPanel-1]-signal.x[iPhi][iPanel]) < 1.0E-12){
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
        sigFile << scientific << signal.p_prime[iPhi][iPanel]   << "\t";
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

CBoom_AugBurgers::~CBoom_AugBurgers(void){

}

void CBoom_AugBurgers::SearchLinear(CConfig *config, CGeometry *geometry, 
               const su2double r0, const su2double *phi){
  
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

  nPanel = new unsigned long[ray_N_phi];
  for(unsigned short i = 0; i < ray_N_phi; i++){
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
              for(unsigned short iPhi = 0; iPhi < ray_N_phi; iPhi++){
                inside = InsideElem(geometry, r0, phi[iPhi], jElem, p0, p1);
                if(inside){
                  if(nPanel[iPhi] == 0){
                    nPanel[iPhi] = 1;
                    pointID_original = new unsigned long*[ray_N_phi];
                    pointID_original[iPhi] = new unsigned long[1];
                    Coord_original = new su2double**[ray_N_phi];
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
                      if(nDim == 3) Coord_tmp[i][2] = Coord_original[iPhi][i][2];

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
                      if(nDim == 3) Coord_original[iPhi][i][2] = Coord_tmp[i][2];

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
                      Coord_original[iPhi][nPanel[iPhi]-1][1] = r0*sin(ray_phi[iPhi]);
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

void CBoom_AugBurgers::ExtractLine(CGeometry *geometry, const su2double r0, unsigned short iPhi){
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
                if(nDim == 3) Coord_tmp[i][2] = Coord_original[iPhi][i][2];

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
                if(nDim == 3) Coord_original[iPhi][i][2] = Coord_tmp[i][2];

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
                Coord_original[iPhi][nPanel[iPhi]-1][1] = r0*sin(ray_phi[iPhi]);
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

void CBoom_AugBurgers::ExtractPressure(CSolver *solver, CConfig *config, CGeometry *geometry, unsigned short iPhi){
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
  signal.p_prime[iPhi] = new su2double[nPanel[iPhi]];
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
          signal.p_prime[iPhi][iElem] += (config->GetGamma()-1)*rho*StaticEnergy*isoparams[iElem][iNode];
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

bool CBoom_AugBurgers::InsideElem(CGeometry *geometry, su2double r0, su2double phi, unsigned long jElem, su2double *p0, su2double *p1){
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

int CBoom_AugBurgers::Intersect2D(su2double r0, su2double *Coord_i, su2double *Coord_ip1, su2double *p0, su2double *p1){

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

int CBoom_AugBurgers::Intersect3D(su2double r0, su2double phi, int nCoord, su2double **Coord_i, su2double *p1){
  
  su2double y0 = r0*sin(phi), z0 = -r0*cos(phi);
  su2double ymin = 1.0E9, ymax = -1.0E9, zmin = 1.0E9, zmax = -1.0E9;

  /*--- First check simple bounding box ---*/
  //cout << "Check bounding box";
  for(int iCoord = 0; iCoord < nCoord; iCoord++){
    if(Coord_i[iCoord][1] < ymin) ymin = Coord_i[iCoord][1];
    if(Coord_i[iCoord][1] > ymax) ymax = Coord_i[iCoord][1];
    if(Coord_i[iCoord][2] < zmin) zmin = Coord_i[iCoord][2];
    if(Coord_i[iCoord][2] > zmax) zmax = Coord_i[iCoord][2];
  }

  if(y0 < ymin || y0 > ymax || z0 < zmin || z0 > zmax){
    //cout << "y0 = " << y0 << ", z0 = " << z0 << ", ymin = " << ymin << ", ymax = " << ymax << ", zmin = " << zmin << ", zmax = " << zmax << endl;
    return 0;
  }

  /*--- If inside bounding box, check sum of angles ---*/
  //cout << "Check sum of angles" << endl;
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
    p1[0] = 0.0;
    p1[1] = y0;
    p1[2] = z0;
    for(int iCoord = 0; iCoord < nCoord; iCoord++){
      p1[0] += isoparams[iCoord]*Coord_i[iCoord][0];
    }

    delete [] isoparams;
    delete [] Coord;
    delete [] X_donor;
    return 1;
  }
  else{
    return 0;
  }

}

void AtmosISA(su2double& h0, su2double& T, su2double& a, su2double& p,
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

void CBoom_AugBurgers::ConditionAtmosphericData(){
  su2double g = atm_g;    // specific heat ratio
  su2double h0 = flt_h;    // flight altitude [m]
  su2double T, p;

  /*---Conditions at flight altitude---*/
  AtmosISA(h0, T_inf, a_inf, p_inf, rho_inf, g);

  /*---Atmospheric conditions profile---*/
  z = new su2double[n_prop];
  a_of_z = new su2double[n_prop];
  p_of_z = new su2double[n_prop];
  rho_of_z = new su2double[n_prop];
  for(unsigned int i = 0; i < n_prop; i++){
    z[i] = h0*su2double(i)/(su2double(n_prop)-1.);
    AtmosISA(z[i], T, a_of_z[i], p_of_z[i], rho_of_z[i], g);
  }
}

void CBoom_AugBurgers::ScaleFactors(){

  p0 = p_inf;
  w0 = flt_M*a_inf/config->GetRefLength();
  beta = 1. + (atm_g - 1.)/2.;
  xbar = rho_inf*c_inf^3/(beta*w0*p0);

  scale_L = config->GetRefLength();
  scale_z = flt_h;
  scale_T = scale_L/(flt_M*a_inf);    // flow over aircraft [s]

}


void CBoom_AugBurgers:: Sph2Cart(su2double& nx, su2double& ny, su2double& nz, su2double az, su2double elev,
              su2double r){
  /*---Compute spherical coordinates from elevation, azimuth, and radius---*/
  nx = r*cos(elev)*cos(az);
  ny = r*cos(elev)*sin(az);
  nz = r*sin(elev);
}

void CBoom_AugBurgers::InitialWaveNormals(){

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
  a0 = a_of_z[n_prop-1];

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

su2double *derivs(su2double x, int m, su2double y[], CBoom_AugBurgers::RayData data){
  su2double *dydx;
  su2double a, rho, p, T;
  su2double g = 1.4;
  su2double xdim = x*data.L;

  AtmosISA(xdim, T, a, p, rho, g);

  su2double theta = acos(a/data.c0);
  su2double num = cos(theta);
  su2double denom = sin(theta);

  dydx = new su2double[3];
  dydx[0] = (num*sin(data.nu))/denom;
  dydx[1] = (num*cos(data.nu))/denom;
  dydx[2] = (data.L/(data.T*a))/denom;

  return dydx;
}

su2double *CBoom_AugBurgers::rk4(su2double x0, int m, su2double y0[], su2double dx, RayData data,
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

su2double *CBoom_AugBurgers::SplineGetDerivs(su2double x[], su2double y[], int N){
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

void CBoom_AugBurgers::RayTracer(unsigned short iPhi){

  /*---Scale factors---*/
  su2double L = flt_h;
  su2double T = scale_T;

  /*---Ambient conditions---*/
  su2double a0, rho0, p0;
  su2double a, rho, p;

  /*---Information for RK4 solver---*/
  su2double r0[3];
  su2double *f;
  su2double dz = (z[0] - z[1]);

  /*---Tolerances for initial ray tube---*/
  su2double tol_dr   = 1.0E-3;
  su2double tol_dphi = 1.0E-3;

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
    x_of_z[j] = new su2double[n_prop];
    y_of_z[j] = new su2double[n_prop];
    t_of_z[j] = new su2double[n_prop];
    ray_theta[j] = new su2double[n_prop];
  }

  /*---Primary ray---*/
  data.c0 = ray_c0[iPhi][0];
  data.nu = ray_nu[iPhi][0];
  r0[0] = r0[1] = r0[2] = 0.;
  x_of_z[0][n_prop-1] = r0[0];
  y_of_z[0][n_prop-1] = r0[1];
  t_of_z[0][n_prop-1] = r0[2];
  ray_theta[0][n_prop-1] = acos(a_of_z[n_prop-1]/data.c0);
  for(int j = n_prop-1; j > 0; j--){
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
  x_of_z[1][n_prop-1] = r0[0];
  y_of_z[1][n_prop-1] = r0[1];
  ray_theta[1][n_prop-1] = acos(a_of_z[n_prop-1]/data.c0);
  for(int j = n_prop-1; j > 0; j--){
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
  x_of_z[2][n_prop-1] = r0[0];
  y_of_z[2][n_prop-1] = r0[1];
  t_of_z[2][n_prop-1] = r0[2];
  ray_theta[2][n_prop-1] = acos(a_of_z[n_prop-1]/data.c0);
  for(int j = n_prop-1; j > 0; j--){
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
  x_of_z[3][n_prop-1] = r0[0];
  y_of_z[3][n_prop-1] = r0[1];
  ray_theta[3][n_prop-1] = acos(a_of_z[n_prop-1]/data.c0);
  for(int j = n_prop-1; j > 0; j--){
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

void CBoom_AugBurgers::RayTubeArea(unsigned short iPhi){

  su2double Ah, x_int, y_int, z_int;
  su2double corners[4][3];
  int M;

  ray_A = new su2double[n_prop];

  /*---Loop over rays---*/
  M = n_prop;
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

su2double CBoom_AugBurgers::matchr(int j, su2double h_L, su2double r0){
  su2double f;
  f = sqrt(pow(x_of_z[0][j],2) + pow(y_of_z[0][j],2) + pow(z[j]/scale_z-1.0,2))*h_L - r0/scale_L;
  return f;
}

void CBoom_AugBurgers::PropagateSignal(unsigned short iPhi){
  
  for(unsigned long iIter = 0; iIter < n_prop; iIter++){
  	Preprocessing(iPhi, iIter);
  	Nonlinearity(iPhi);
  	Attenuation(iPhi);
  	Relaxation(iPhi);
  	Spreading(iPhi);
  	Stratification(iPhi);
  	Iterate(iPhi);
  }
}

void CBoom_AugBurgers::Preprocessing(unsigned short iPhi, unsigned long iIter){
  
  /*---Preprocess signal for first iteration---*/
  if(iIter == 0){
  	signal.P       = new su2double[signal.len[iPhi]];
  	signal.t       = new su2double[signal.len[iPhi]];
  	signal.t_prime = new su2double[signal.len[iPhi]];
  	signal.tau.    = new su2double[signal.len[iPhi]];
  	for(unsigned long i = 0; i < signal.len[iPhi]; i++){
  		signal.P[i]       = signal.p_prime[iPhi][i]/p0;
  		signal.t[i]       = signal.x[iPhi][i]/(a_inf*flt_M);
  		signal.t_prime[i] = signal.t[i] - signal.x[iPhi][i]/a_inf;
  		signal.tau[i]     = w0*signal.t_prime[i];
  	}
  }

  /*---Compute other coefficients needed for solution of ABE---*/
  C_nu_O2 = 
  C_nu_N2 = 
}

void CBoom_AugBurgers::Nonlinearity(unsigned short iPhi){}

void CBoom_AugBurgers::Attenuation(unsigned short iPhi){}

void CBoom_AugBurgers::Relaxation(unsigned short iPhi){}

void CBoom_AugBurgers::Spreading(unsigned short iPhi){}

void CBoom_AugBurgers::Stratification(unsigned short iPhi){}

void CBoom_AugBurgers::Iterate(unsigned short iPhi){}

