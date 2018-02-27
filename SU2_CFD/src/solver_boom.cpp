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

void CBoom_AugBurgers::PropagateSignal(unsigned short iPhi){

  unsigned long iIter = 0;
  ground_flag = false;
  
  while(!ground_flag){
  	Preprocessing(iPhi, iIter);
  	Nonlinearity(iPhi);
  	Attenuation(iPhi);
  	Relaxation(iPhi);
  	Spreading(iPhi);
  	Stratification(iPhi);
  	Iterate(iPhi);
    iIter++;
  }
}

void CBoom_AugBurgers::Preprocessing(unsigned short iPhi, unsigned long iIter){
  
  unsigned long i_prop = n_prop-iIter-1; // Current index (since we iterate from N-1 to 0)
  su2double R = 287.058,         // Gas constant
            mu0 = 1.846E-5,      // Reference viscosity
            kappa0 = 2.624E-2,   // Reference thermal conduction coefficient
            T0  = 300.,          // Reference temperature (Sutherland)
            Ts  = 110.4,         // Reference temperature (Sutherland)
            Ta  = 245.4,         // Reference temperature (Sutherland)
            Tb  = 27.6,          // Reference temperature (Sutherland)
            T01 = 273.16,        // Reference temperature (Humidity)
            Tr  = 293.15,        // Reference temperature (Absorption)
            Pr  = 101325.;       // Reference pressure (Absorption)

  /*---Preprocess signal for first iteration---*/
  if(iIter == 0){

    /*---First create uniform grid since algorithm requires it---*/
    CreateUniformGridSignal(iPhi);

  	signal.P       = new su2double[signal.len[iPhi]];
    signal.dP_att  = new su2double[signal.len[iPhi]];
    signal.dP_rel  = new su2double[signal.len[iPhi]];
    signal.dP_spr  = new su2double[signal.len[iPhi]];
    signal.dP_str  = new su2double[signal.len[iPhi]];
  	signal.t       = new su2double[signal.len[iPhi]];
  	signal.t_prime = new su2double[signal.len[iPhi]];
  	signal.tau     = new su2double[signal.len[iPhi]];
    signal.taud    = new su2double[signal.len[iPhi]];
  	for(unsigned long i = 0; i < signal.len[iPhi]; i++){
  		signal.P[i]       = signal.p_prime[iPhi][i]/p0;
  		signal.t[i]       = (signal.x[iPhi][i]-signal.x[iPhi][0])/(a_inf);
  		signal.t_prime[i] = signal.t[i] - signal.x[iPhi][i]/a_inf;
  	}
    w0 = 2*M_PI/signal.t[iPhi][signal.len[iPhi]];
    beta = 1. + (atm_g - 1.)/2.;
    p0 = p_inf;
    for(unsigned long i = 0; i < signal.len[iPhi]; i++){
      signal.tau[i] = w0*signal.t[i];
    }
    dtau = signal.tau[1] - signal.tau[0];

    CreateInitialRayTube(iPhi);

    ground_flag = false;

  }

  /*---Compute other coefficients needed for solution of ABE---*/
  AtmosISA(ray_z, T_inf, c0, p_inf, rho0, atm_g);

  xbar = rho0*pow(c0,3)/(beta*w0*p0);

  mu        = mu0*pow(T_inf/T0,1.5)*(T0+Ts)/(T_inf+Ts);
  kappa     = kappa0*pow(T_inf/T0,1.5)*(T0+Ta*exp(-Tb/T0))/(T_inf+Ta*exp(-Tb/T_inf));
  delta     = mu/rho0*(4./3. + 0.6 + pow((atm_g-1.),2)*kappa/(atm_g*R*mu));
  alpha0_tv = delta*pow(w0,2)/(2.*pow(c0,3));
  Gamma     = 1./(alpha0_tv*xbar);

  su2double logPsat_Pr = -6.8346*pow(T01/T_inf,1.261) + 4.6151;  // ISO standard for saturation vapor pressure
  su2double hr = 20.;                                                     // Relative humidity (percent), hardcoding for now
  su2double h = hr * Pr/p_inf * pow(logPsat_Pr,10.);             // Absolute humidity (percent)

  su2double A_nu_O2 = 0.01275 * pow(T_inf/Tr,-2.5) * exp(-2339.1/T_inf);
  su2double A_nu_N2 = 0.1068  * pow(T_inf/Tr,-2.5) * exp(-3352./T_inf);

  m_nu_O2 = 2.*c0*A_nu_O2;
  m_nu_N2 = 2.*c0*A_nu_N2;

  su2double f_nu_O2 = p_inf/Pr * sqrt(T_inf/Tr) * (9. + 280.*h*exp(pow(Tr/T_inf, 1./3.) - 1.));
  su2double f_nu_N2 = p_inf/Pr * (24. + 4.04E4*h*(0.02+h)/(0.391+h));

  tau_nu_O2 = 1./(2.*M_PI*f_nu_O2);
  tau_nu_N2 = 1./(2.*M_PI*f_nu_N2);

  theta_nu_O2 = w0*tau_nu_O2;
  theta_nu_N2 = w0*tau_nu_N2;

  C_nu_O2 = m_nu_O2*tau_nu_O2*pow(w0,2)*xbar/(2*c0);
  C_nu_N2 = m_nu_N2*tau_nu_N2*pow(w0,2)*xbar/(2*c0);

}

void CBoom_AugBurgers::CreateUniformGridSignal(unsigned short iPhi){
  
  /*---Loop over signal and find smallest spacing---*/
  su2double dx, dx_min = 1.0E9;
  for(unsigned long i = 1; i < signal.len[iPhi]; i++){
    dx = signal.x[iPhi][i] - signal.x[iPhi][i-1];
    dx_min = min(dx, dx_min);
  }

  /*---Create new temp signal---*/
  unsigned long j = 0, len_new = ceil((signal.x[iPhi][signal.len[iPhi]-1]-signal.x[iPhi][0])/dx_min);
  su2double *xtmp = new su2double[len_new],
            *ptmp = new su2double[len_new];
  xtmp[0] = signal.x[iPhi][0];
  ptmp[0] = signal.p_prime[iPhi][0];
  for(unsigned long i = 1; i < len_new; i++){
    xtmp[i] = xtmp[i-1] + dx_min;
    if(xtmp[i] > signal.x[iPhi][j+1]) j++;

    if(j == (signal.len[iPhi]-1)){
      /*---Zero-pad end of signal, then decrement to avoid potential out-of-bounds access---*/
      ptmp[i] = 0.;
      j--;
    }
    else{
      /*---Interpolate signal---*/
      ptmp[i] = signal.p_prime[iPhi][j] + (xtmp[i] - signal.x[iPhi][j]) * (signal.p_prime[iPhi][j+1]-signal.p_prime[iPhi][j])/(signal.x[iPhi][j+1]-signal.x[iPhi][j]);
    }
  }

  /*---Store new signal---*/
  signal.len[iPhi] = len_new;
  signal.x[iPhi] = new su2double[len_new];
  signal.p_prime[iPhi] = new su2double[len_new];
  for(unsigned long i = 0; i < len_new; i++){
    signal.x[iPhi][i] = xtmp[i];
    signal.p_prime[iPhi][i] = ptmp[i];
  }
  delete [] xtmp;
  delete [] ptmp;

}

void CBoom_AugBurgers::CreateInitialRayTube(unsigned short iPhi){

  ray_x = new su2double[4];
  ray_y = new su2double[4];
  ray_gamma = new su2double[2];
  ray_theta = new su2double[2];

  /*---Ray tube origin---*/
  ray_x[0] = 0.;
  ray_y[0] = ray_r0*sin(ray_phi[iPhi]);
  ray_z    = flt_h - ray_r0*cos(ray_phi[iPhi]);

  /*---Ray tube corners---*/
  ray_y[1] = ray_y[0];
  ray_x[3] = ray_x[0];
  ray_x[1] = ray_x[2] = ray_x[0]+1.0E-3;
  ray_y[2] = ray_y[3] = ray_r0*sin(ray_phi[iPhi]+1.0E-3);

  ray_lambda = pow(flt_M*flt_M-1.,-0.5);
  ray_gamma[0]  = asin(ray_lambda*sin(ray_phi[iPhi])*pow(1. + pow(ray_lambda*sin(ray_phi[iPhi],2)), -0.5));
  ray_gamma[1]  = asin(ray_lambda*sin(ray_phi[iPhi]+1.0E-3)*pow(1. + pow(ray_lambda*sin(ray_phi[iPhi]+1.0E-3,2)), -0.5));
  ray_theta[0]  = acos(-1./(flt_M*cos(ray_gamma[0])));
  ray_theta[1]  = acos(-1./(flt_M*cos(ray_gamma[1])));

  su2double u[2] = {ray_x[2]-ray_x[0], ray_y[2]-ray_y[0]};
  su2double v[2] = {ray_x[3]-ray_x[1], ray_y[3]-ray_y[1]};
  su2double c    = u[0]*v[1] - u[1]*v[0];
  su2double Ah   = 0.5*sqrt(pow(c,2));

  ray_A = c0*A_h*sin(ray_theta[0])/(c0);  // TODO: Add wind contribution
}

void CBoom_AugBurgers::Nonlinearity(unsigned short iPhi){

  /*---Restrict dsigma to avoid multivalued waveforms---*/
  su2double dp_dtau, max_dp_dtau = -1.0E9;
  for(unsigned long i = 0; i < signal.len[iPhi]-1; i++){
    dp_dtau = (signal.P[i+1]-signal.P[i])/(signal.tau[i+1]-signal.tau[i]);
    max_dp_dtau = max(dp_dtau, max_dp_dtau);
  }

  dsigma = 0.5/max_dp_dtau; // dsigma < 1/max(dp/dtau)

  /*---Check for intersection with ground plane---*/
  su2double dxi_dz, deta_dz, ds_dz, ds;
  su2double uwind = 0., vwind = 0.;
  ds = xbar*dsigma;
  dxi_dz  = (c0*cos(ray_theta[0]) + uwind)/(c0*sin(ray_theta[0]));
  deta_dz = (vwind)/(c0*sin(ray_theta[0]));
  ds_dz   = sqrt(pow(dxi_dz,2)+pow(deta_dz,2)+1);
  dz = ds/ds_dz;

  if(ray_z-dz < 0.0){
    dz = ray_z;
    ds = ds_dz*dz;
    dsigma = ds/xbar;
    ground_flag = true;
  }

  /*---Now determine distorted time grid---*/
  for(unsigned long i = 0; i < signal.len[iPhi]; i++){
    signal.taud[i] = signal.tau[i] - signal.P[i]*dsigma;
  }
}

void CBoom_AugBurgers::Attenuation(unsigned short iPhi){

  su2double lambda = dsigma/(2.*Gamma*pow(dtau,2));
  su2double *BPk = new su2double[signal.len[iPhi]];

  /*---Tridiagonal matrix-vector multiplication B_tv * Pk (see Cleveland thesis)---*/
  BPk[0] = signal.P[0];
  BPk[signal.len[iPhi]-1] = signal.P[signal.len[iPhi]-1];
  for(unsigned long i = 1; i < signal.len[iPhi]-1; i++){
    BPk[i] = lambda*(signal.P[i+1]+signal.P[i-1]) + (1.-2.*lambda)*signal.P[i];
  }

  /*---Solve for Pk+1 via Thomas algorithm for tridiagonal matrix---*/
  su2double *ci = new su2double[signal.len[iPhi]-1],
            *di = new su2double[signal.len[iPhi]-1];

  ci[0] = 0.;
  di[0] = BPk[0];
  for(unsigned long i = 1; i < signal.len[iPhi]-1; i++){
    ci[i] = -lambda/((1.+2.*lambda) + lambda*ci[i-1]);
    di[i] = (BPk[i] + lambda*di[i-1])/((1.+2.*lambda) + lambda*ci[i-1]);
  }

  signal.dP_att[signal.len[iPhi]-1] = BPk[signal.len[iPhi]-1];
  for(unsigned long i = signal.len[iPhi]-2; i >= 0; i--){
    signal.dP_att[i] = di[i] - ci[i]*signal.dP_att[i+1];
  }

  /*---Get change in pressure---*/
  for(unsigned long i = 0; i < signal.len[iPhi]; i++){
    signal.dP_att[i] = signal.dP_att[i] - signal.P[i];
  }

  delete [] BPk;
  delete [] ci;
  delete [] di;
}

void CBoom_AugBurgers::Relaxation(unsigned short iPhi){

  su2double lambda[2] = {dsgima*C_nu_O2/pow(dtau,2), dsgima*C_nu_O2/pow(dtau,2)},
            mu[2] = {theta_nu_O2/(2.*dtau), theta_nu_N2/(2.*dtau)},
            alpha = 0.5,
            alpha_p = 1.-alpha;

  su2double *BPk = new su2double[signal.len[iPhi]],
            *dP  = new su2double[signal.len[iPhi]],
            *ci  = new su2double[signal.len[iPhi]-1],
            *di  = new su2double[signal.len[iPhi]-1];

  /*---Reset dP_rel---*/
  for(unsigned long i = 0; i < signal.len[iPhi]; i++){
    signal.dP_rel[i] = 0.;
  }

  /*---Compute effect of different relaxation modes---*/
  for(unsigned short j = 0; j < 2; j++){
    /*---Tridiagonal matrix-vector multiplication B_nu * Pk (see Cleveland thesis)---*/
    BPk[0] = signal.P[0];
    BPk[signal.len[iPhi]-1] = signal.P[signal.len[iPhi]-1];
    for(unsigned long i = 1; i < signal.len[iPhi]-1; i++){
      BPk[i] = (alpha_p*lambda[j]-mu[j])*signal.P[i-1] + (1.-2.*alpha_p*lambda[j])*signal.P[i] + (alpha_p*lambda[j]+mu[j])*signal.P[i+1];
    }

    /*---Solve for Pk+1 via Thomas algorithm for tridiagonal matrix---*/
    ci[0] = 0.;
    di[0] = BPk[0];
    for(unsigned long i = 1; i < signal.len[iPhi]-1; i++){
      ci[i] = -(alpha*lambda[j]-mu[j])/((1.+2.*alpha*lambda[j]) + (alpha*lambda[j]+mu[j])*ci[i-1]);
      di[i] = (BPk[i] + (alpha*lambda[j]+mu[j])*di[i-1])/((1.+2.*alpha*lambda[j]) + (alpha*lambda[j]+mu[j])*ci[i-1]);
    }

    dP[signal.len[iPhi]-1] = BPk[signal.len[iPhi]-1];
    for(unsigned long i = signal.len[iPhi]-2; i >= 0; i--){
      dP[i] = di[i] - ci[i]*dP[i+1];
    }

    /*---Get change in pressure---*/
    for(unsigned long i = 0; i < signal.len[iPhi]; i++){
      signal.dP_rel[i] += (dP[i]-signal.P[i]);
    }
  }

  delete [] BPk;
  delete [] dP;
  delete [] ci;
  delete [] di;
}

void CBoom_AugBurgers::Spreading(unsigned short iPhi){

  /*---Compute ray tube area at sigma+dsigma---*/
  su2double dxi_dz, deta_dz, ds_dz, ds;
  su2double xi_new, eta_new, x_new[4], y_new[4], z_new, A_new;
  su2double uwind = 0., vwind = 0.;

  ds = xbar*dsigma;
  dxi_dz  = (c0*cos(ray_theta[0]) + uwind)/(c0*sin(ray_theta[0]));
  deta_dz = (vwind)/(c0*sin(ray_theta[0]));
  ds_dz   = sqrt(pow(dxi_dz,2)+pow(deta_dz,2)+1);
  dz = ds/ds_dz;
  z_new = ray_z - dz;

  x_new[0] = ray_x[0] + dxi_dz*dz*cos(ray_gamma[0]) - deta_dz*dz*sin(ray_gamma[0]);
  y_new[0] = ray_y[0] + dxi_dz*dz*sin(ray_gamma[0]) + deta_dz*dz*cos(ray_gamma[0]);
  x_new[1] = ray_x[1] + dxi_dz*dz*cos(ray_gamma[0]) - deta_dz*dz*sin(ray_gamma[0]);
  y_new[1] = ray_y[1] + dxi_dz*dz*sin(ray_gamma[0]) + deta_dz*dz*cos(ray_gamma[0]);

  dxi_dz  = (c0*cos(ray_theta[0]) + uwind)/(c0*sin(ray_theta[0]));
  deta_dz = (vwind)/(c0*sin(ray_theta[0]));

  x_new[2] = ray_x[2] + dxi_dz*dz*cos(ray_gamma[1]) - deta_dz*dz*sin(ray_gamma[1]);
  y_new[2] = ray_y[2] + dxi_dz*dz*sin(ray_gamma[1]) + deta_dz*dz*cos(ray_gamma[1]);
  x_new[3] = ray_x[3] + dxi_dz*dz*cos(ray_gamma[1]) - deta_dz*dz*sin(ray_gamma[1]);
  y_new[3] = ray_y[3] + dxi_dz*dz*sin(ray_gamma[1]) + deta_dz*dz*cos(ray_gamma[1]);

  su2double u[2] = {x_new[2]-x_new[0], y_new[2]-y_new[0]};
  su2double v[2] = {x_new[3]-x_new[1], y_new[3]-y_new[1]};
  su2double c    = u[0]*v[1] - u[1]*v[0];
  su2double Ah   = 0.5*sqrt(pow(c,2));

  A_new = c0*A_h*sin(ray_theta[0])/(c0);  // TODO: Add wind contribution

  /*---Compute change in pressure from exact solution to general ray tube area eqn---*/
  for(unsigned long i = 0; i < signal.len[iPhi]; i++){
    signal.dP_spr[i] = (sqrt(ray_A/A_new)-1.0)*signal.P[i];
  }

  /*---Set new ray properties (except z, we'll do that later)---*/
  for(unsigned short i = 0; i < 4; i++){
    ray_x[i] = x_new[i];
    ray_y[i] = y_new[i];
  }
  ray_A = A_new;

}

void CBoom_AugBurgers::Stratification(unsigned short iPhi){

  /*---Compute atmospheric properties at sigma+dsigma---*/
  su2double Tp1, cp1, pp1, rhop1;
  AtmosISA(ray_z-dz, Tp1, cp1, pp1, rhop1, atm_g);

  /*---Compute change in pressure from exact solution to stratification eqn---*/
  for(unsigned long i = 0; i < signal.len[iPhi]; i++){
    signal.dP_str[i] = (sqrt((rhop1*cp1)/(rho0*c0))-1.0)*signal.P[i];
  }

}

void CBoom_AugBurgers::Iterate(unsigned short iPhi){

  /*---Propagate signal to sigma+dsigma---*/
  for(unsigned long i = 0; i < signal.len[iPhi]; i++){
    signal.P[i] += signal.dP_att[i] + signal.dP_rel[i] + signal.dP_spr[i] + signal.dP_str[i];
  }

  /*---Account for nonlinearity and interpolate new signal from distorted grid---*/
  su2double *Ptmp = new su2double[signal.len[iPhi]];
  unsigned long j, jp1;
  for(unsigned long i = 0; i < signal.len[iPhi]; i++){
    if(signal.tau[i] < signal.taud[0]){
      j   = 0;
      jp1 = 1;
    }
    else if(signal.tau[i] > signal.taud[signal.len[iPhi]-1]){
      j   = signal.len[iPhi]-2;
      jp1 = signal.len[iPhi]-1;
    }
    else{
      for(unsigned long k = 0; k < signal.len[iPhi]-1; k++){
        if((signal.tau[i] >= signal.taud[k]) && (signal.tau[i] <= signal.taud[k+1])){
          j   = k;
          jp1 = k+1;
          break;
        }
      }
    }

    Ptmp[i] = signal.P[j] + (signal.tau[i] - signal.taud[j]) * (signal.P[jp1]-signal.P[j])/(signal.taud[jp1]-signal.taud[j]);
  }

  for(unsigned long i = 0; i < signal.len[iPhi]; i++){
    signal.P[i] = Ptmp[i];
  }

  /*---Set new altitude---*/
  ray_z -= dz;

  delete [] Ptmp;

}

