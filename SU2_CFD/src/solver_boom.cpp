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

  Kind_Boom_Cost = config->GetKind_ObjFunc();
  AD_flag = false;
  if(config->GetAD_Mode()) AD_flag = true;
  CFL_reduce = config->GetBoom_cfl_reduce();
  Kind_Step = config->GetBoom_step_type();
  Step_size = config->GetBoom_step_size();
  Step_growth = config->GetBoom_step_growth();

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
  
  if(rank == MASTER_NODE)
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
  su2double *pp0, *pp1;

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
                    Coord_original = new su2double**[ray_N_phi];
                    Coord_original[iPhi] = new su2double*[1];
                    Coord_original[iPhi][0] = new su2double[nDim];

                    pointID_original[iPhi][0] = jElem;
                    Coord_original[iPhi][0][0] = (pp0[0] + pp1[0])/2.0;

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
                    Coord_original[iPhi][nPanel[iPhi]-1][0] = (pp0[0] + pp1[0])/2.0;

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
  delete [] pp0;
  delete [] pp1;

}

void CBoom_AugBurgers::ExtractLine(CGeometry *geometry, const su2double r0, unsigned short iPhi){
  bool inside, inside_iPanel, addPanel, end = false;
  unsigned short iElem, nElem;
  unsigned long jElem, jElem_m1, nElem_tot = geometry->GetnElem();
  su2double x_i, x_m1;

  unsigned long *pointID_tmp;
  su2double **Coord_tmp;
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
              Coord_original[iPhi][nPanel[iPhi]-1][0] = (pp0[0] + pp1[0])/2.0;
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

    for(iNode = 0; iNode < nNodeFace; iNode++){
      delete [] Coord_face[iNode];
    }
    delete [] Coord_face;
  }

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
    pp1[0] = 0.0;
    pp1[1] = y0;
    pp1[2] = z0;
    for(int iCoord = 0; iCoord < nCoord; iCoord++){
      pp1[0] += isoparams[iCoord]*Coord_i[iCoord][0];
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
  // su2double logPsat_Pr = 10.79586*(1. - T01/T_inf) - 5.02808 * log10(T_inf/T01) + 1.50474E-4*(1. - pow(10.,-8.29692*(T_inf/T01 - 1.))) - 4.2873E-4*(1. - pow(10.,-4.76955*(T_inf/T01 - 1.))) - 2.2195983;

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

void CBoom_AugBurgers::PropagateSignal(unsigned short iPhi){

int rank = 0;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  unsigned long iIter = 0;
  ground_flag = false;

  if(rank == MASTER_NODE){
  while(!ground_flag){

  	Preprocessing(iPhi, iIter);
    Scaling(iPhi);
    Relaxation(iPhi, iIter);
    Attenuation(iPhi);
  	Nonlinearity(iPhi);

    ray_z -= dz;

    iIter++;
  }
  }

  if(rank == MASTER_NODE){
    cout.width(5); cout << iIter;
    cout.width(12); cout.precision(6); cout << ray_z;
    cout.width(12); cout.precision(6); cout << p_peak << endl;
    cout << "Signal propagated in " << iIter << " iterations." << endl;
  }

  WriteGroundPressure(iPhi);
  if(Kind_Boom_Cost==BOOM_LOUD){
    if(rank == MASTER_NODE) PerceivedLoudness(iPhi);
  }
  else{
    if(rank == MASTER_NODE) AcousticEnergy(iPhi);
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
            T01 = 273.16,        // Reference temperature (Humidity)
            Tr  = 293.15,        // Reference temperature (Absorption)
            Pr  = 101325.;       // Reference pressure (Absorption)

  su2double c0_old, rho0_old, max_dp, h;

  /*---Preprocess signal for first iteration---*/
  if(iIter == 0){

    ray_z    = flt_h - ray_r0*cos(ray_phi[iPhi]);
    AtmosISA(ray_z, T_inf, c0, p_inf, rho0, atm_g);
    p0 = p_inf;
    flt_U = flt_M*c0;
    beta = (atm_g + 1.)/2.;
    HumidityISO(ray_z, p_inf, T_inf, h);
    c0 *= (1+0.0016*h); // Humidity correction to speed of sound

    /*---First create uniform grid since algorithm requires it---*/
    CreateUniformGridSignal(iPhi);
    dsigma_old = 0.001;
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
    cout << "Acoustic Mach number: " << M_a << endl;
  	for(unsigned long i = 0; i < signal.len[iPhi]; i++){
      signal.t[i]       = (signal.x[iPhi][i]-signal.x[iPhi][0])/flt_U ;
  	}
    f0 = flt_U/(signal.x[iPhi][signal.len[iPhi]-1]-signal.x[iPhi][0]); // Characteristic time governed by length of signal
    w0 = 2.*M_PI*f0; 
    for(unsigned long i = 0; i < signal.len[iPhi]; i++){
      signal.P[i] = signal.p_prime[iPhi][i]/p0;
      signal.tau[i] = w0*signal.t[i];
    }
    dtau = signal.tau[1] - signal.tau[0];

    CreateInitialRayTube(iPhi);

    /*---Initial signal---*/
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

  else{
    AtmosISA(ray_z, T_inf, c0, p_inf, rho0, atm_g);
    HumidityISO(ray_z, p_inf, T_inf, h);
    c0 *= (1+0.0016*h); // Humidity correction to speed of sound
  }

  /*---Compute other coefficients needed for solution of ABE---*/
  xbar = rho0*pow(c0,3)/(beta*w0*p0);
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

  DetermineStepSize(iPhi, iIter);

  if(iIter%10 == 0){
    cout.width(5); cout << iIter;
    cout.width(12); cout.precision(6); cout << ray_z;
    cout.width(12); cout.precision(6); cout << p_peak << endl;
  }

}

void CBoom_AugBurgers::CreateUniformGridSignal(unsigned short iPhi){

  /*---Clip signal---*/
  unsigned long istart = 0, iend = signal.len[iPhi]-1;
  for(unsigned long i = 0; i < signal.len[iPhi]; i++){
    if(abs(signal.p_prime[iPhi][i]/p_inf) > 1.0E-3){
      istart = i;
      break;
    }
  }

  for(unsigned long i = signal.len[iPhi]-1; i >= istart+1; i--){
    if(abs(signal.p_prime[iPhi][i]/p_inf) > 1.0E-3){
      iend = i;
      break;
    }
  }
  scale_L = signal.x[iPhi][iend] - signal.x[iPhi][istart];
  unsigned long len_new = signal.len[iPhi];
  su2double *xtmp, *ptmp;

  /*---Loop over signal and find smallest spacing---*/
  su2double dx, dx_min = 1.0E3; // TODO: add min spacing input scale_L/1001.; // sBOOM gets "sufficient" results in vicinity of 10000 pts
  for(unsigned long i = 1; i < signal.len[iPhi]; i++){
    dx = signal.x[iPhi][i] - signal.x[iPhi][i-1];
    dx_min = min(dx, dx_min);
  }

  /*---Create new temp signal---*/
  unsigned long j = 0;
  len_new = ceil((signal.x[iPhi][signal.len[iPhi]-1]-signal.x[iPhi][0])/dx_min);
  xtmp = new su2double[len_new];
  ptmp = new su2double[len_new];
  xtmp[0] = signal.x[iPhi][0];
  ptmp[0] = signal.p_prime[iPhi][0];
  for(unsigned long i = 1; i < len_new; i++){
    xtmp[i] = xtmp[i-1] + dx_min;
    if(xtmp[i] > signal.x[iPhi][j+1]) j++;

    if(j == (signal.len[iPhi]-1)){
      j--; // Decrement to avoid potential out-of-bounds access
    }
    /*---Interpolate signal---*/
    ptmp[i] = signal.p_prime[iPhi][j] + (xtmp[i] - signal.x[iPhi][j]) * (signal.p_prime[iPhi][j+1]-signal.p_prime[iPhi][j])/(signal.x[iPhi][j+1]-signal.x[iPhi][j]);
  }

  /*---Store new signal---*/
  su2double dp_dx_end = -ptmp[len_new-1]/(2.*scale_L);
  unsigned long len_recompress = ceil(2.*scale_L/dx_min);
  signal.len[iPhi] = ceil(len_new+len_recompress*8);
  signal.x[iPhi] = new su2double[signal.len[iPhi]];
  signal.p_prime[iPhi] = new su2double[signal.len[iPhi]];
  unsigned long i0 = floor(len_recompress*4), i1 = i0+len_new, i2 = signal.len[iPhi];
  /*---Zero-pad front of signal---*/
  for(unsigned long i = 0; i < i0; i++){
    signal.x[iPhi][i] = xtmp[0]-dx_min*su2double(i0-i);
    signal.p_prime[iPhi][i] = 0.;
  }
  /*---Interpolated signal---*/
  for(unsigned long i = i0; i < i1; i++){
    signal.x[iPhi][i] = xtmp[i-i0];
    signal.p_prime[iPhi][i] = ptmp[i-i0];
  }
  /*---Recompress aft of signal---*/
  for(unsigned long i = i1; i < i2; i++){
    signal.x[iPhi][i] = signal.x[iPhi][i1-1]+dx_min*su2double(i+1-i1);
    if(i-i1 < len_recompress){
      signal.p_prime[iPhi][i] = signal.p_prime[iPhi][i-1]+dp_dx_end*dx_min;
    }
    else{
      signal.p_prime[iPhi][i] = 0.;
    }
  }

  cout << "Signal refined and padded, now contains " << signal.len[iPhi] << " points." << endl;
  cout << "Length scale of waveform = " << scale_L << " m." << endl;
  cout << "Sample frequency = " << flt_U/dx_min << " Hz." << endl;

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
    su2double dp, max_dp = 1.0E-9, min_dp = -1.0E-9, max_p = -1E6, min_p = 1E6;
    p_peak = 0.;
    for(unsigned long i = 0; i < signal.len[iPhi]-1; i++){
      dp = signal.P[i+1]-signal.P[i];
      max_dp = max(dp, max_dp);
      p_peak = max(signal.P[i]*p0, p_peak);
    }

    dsigma_non = 0.9*dtau/max_dp; // dsigma < 1/max(dp/dtau)

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
    // dsigma = min(dsigma, dsigma_tv);
    // dsigma = min(dsigma, dsigma_relO);
    // dsigma = min(dsigma, dsigma_relN);
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

  su2double lambda = dsigma/(2.*Gamma*pow(dtau,2));
  su2double *y = new su2double[signal.len[iPhi]];

  /*---Tridiagonal matrix-vector multiplication B_tv * Pk (see Cleveland thesis)---*/
  y[0] = signal.P[0];
  y[signal.len[iPhi]-1] = signal.P[signal.len[iPhi]-1];
  for(unsigned long i = 1; i < signal.len[iPhi]-1; i++){
    y[i] = lambda*(signal.P[i+1]+signal.P[i-1]) + (1.-2.*lambda)*signal.P[i];
  }

  /*---Solve for Pk+1 via Thomas algorithm for tridiagonal matrix---*/
  su2double *ci = new su2double[signal.len[iPhi]-1],
            *di = new su2double[signal.len[iPhi]-1];

  ci[0] = 0.;
  di[0] = y[0];
  for(unsigned long i = 1; i < signal.len[iPhi]-1; i++){
    ci[i] = -lambda/((1.+2.*lambda) + lambda*ci[i-1]);
    di[i] = (y[i] + lambda*di[i-1])/((1.+2.*lambda) + lambda*ci[i-1]);
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
    y[0] = signal.P[0];
    y[signal.len[iPhi]-1] = signal.P[signal.len[iPhi]-1];
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
    di[0] = y[0];
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

  su2double *w, *p_of_w; // Frequency domain signal

  su2double w_min = 1., w_max = 20.E3; // [Hz]
  su2double p_ref = 20.E-6;             // [Pa]

  unsigned long n_sample;

  /*--- Compute frequency domain signal ---*/
  cout << "Performing Fourier Transform." << endl;
  FourierTransform(iPhi, w, p_of_w, n_sample);

  /*--- Compute 1/3-oct bands ---*/
  su2double fc[41], f_min[41], f_max[41], E_band[41], SPL_band[41];
  unsigned long band_edge_inds[41][2];

  cout << "Computing 1/3-oct band pressure levels." << endl;

  for(unsigned short i = 0; i < 41; i++){
    fc[i]     = pow(10.,3.)*pow(2.,-29./3.)*pow(2.,su2double(i)/3.);
    f_min[i]  = fc[i]/pow(2.,1./6.);
    f_max[i]  = fc[i]*pow(2.,1./6.);
    E_band[i] = 0.0;
  }

  /*--- Compute band energies and pressure levels ---*/
  for(unsigned short j = 0; j < 41; j++){
    for(unsigned long i = 0; i < n_sample; i++){
      if(f_min[40-j] > w[i]){
        band_edge_inds[j][0] = i;
        break;
      }
    }
    for(unsigned long i = 0; i < n_sample; i++){
      if(f_max[j] < w[i]){
        band_edge_inds[j][1] = i;
        break;
      }
    }
  }

  su2double ptmp;
  for(unsigned short j = 0; j < 41; j++){
    for(unsigned long i = band_edge_inds[j][0]; i < band_edge_inds[j][1]; i++){
      if(i == band_edge_inds[j][0]){ // Interpolate beginning of band
        ptmp = p_of_w[i] + (f_min[j]-w[i])*(p_of_w[i+1]-p_of_w[i])/(w[i+1]-w[i]);
        E_band[j] += 0.5*(ptmp*ptmp+p_of_w[i+1]*p_of_w[i+1])*(w[i+1]-f_min[j]);
      }
      else if(i == band_edge_inds[j][1]-1){ // Interpolate end of band
        ptmp = p_of_w[i] + (f_max[j]-w[i-1])*(p_of_w[i]-p_of_w[i-1])/(w[i]-w[i-1]);
        E_band[j] += 0.5*(ptmp*ptmp+p_of_w[i]*p_of_w[i])*(f_max[j]-w[i]);
      }
      else{
        E_band[j] += 0.5*(p_of_w[i]*p_of_w[i]+p_of_w[i+1]*p_of_w[i+1])*(w[i+1]-w[i]);
      }
    }
    E_band[j] /= 0.07;  // Critical time of human ear is 0.07 s (Shepherd, 1991)
    SPL_band[j] = 10.*log10(E_band[j]/pow(p_ref,2.));
  }

  /*--- Use band pressure levels to compute perceived loudness ---*/
  cout << "Computing perceived loudness (MarkVII)." << endl;
  MarkVII(iPhi, SPL_band);

  /*--- Clean up ---*/
  delete [] w;
  delete [] p_of_w;

}

void CBoom_AugBurgers::FourierTransform(unsigned short iPhi, su2double *w, su2double *p_of_w, unsigned long &n_sample){

  su2double t1, t2, y1, y2, m;
  su2double w_min = 1., w_max = 20.E3; // [Hz]
  su2double p_real, p_imag;
  unsigned long N = signal.len[iPhi];
  n_sample = 400001;  // w_max * 20 + 1

  w = new su2double[n_sample];
  p_of_w = new su2double[n_sample];

  /*---Initialize frequency domain signal ---*/

  for(unsigned long i = 0; i < n_sample; i++){
    w[i] = w_min + i/(n_sample-1)*(w_max-w_min);
    p_of_w[i] = 0.;
  }

  /*--- Transform ---*/

  for(unsigned long j = 0; j < n_sample; j++){

    p_real = 0.;
    p_imag = 0.;

    for(unsigned long i = 0; i < N-1; i++){
      t1 = signal.tau[i]/w0; t2 = signal.tau[i+1]/w0;
      y1 = signal.P[i]*p0;   y2 = signal.P[i+1]*p0;
      m = (y2-y1)/(t2-t1);

      p_real += (m*(cos(2*M_PI*w[j]*t2) - cos(2*M_PI*w[j]*t1)) + 2*M_PI*w[j]*((y1+m*(t2-t1))*sin(2*M_PI*w[j]*t2) - y1*sin(2*M_PI*w[j]*t1)))/pow(2*M_PI*w[j],2.);
      p_imag += (m*(-sin(2*M_PI*w[j]*t2) + sin(2*M_PI*w[j]*t1)) + 2*M_PI*w[j]*((y1+m*(t2-t1))*cos(2*M_PI*w[j]*t2) - y1*cos(2*M_PI*w[j]*t1)))/pow(2*M_PI*w[j],2.);
    }

    p_of_w[j] = sqrt(p_real*p_real + p_imag*p_imag);

  }

}

void CBoom_AugBurgers::MarkVII(unsigned short iPhi, su2double *SPL_band){

  su2double sonmax, sonsum, sonovr, sone[41];
  su2double uppt, dnpt, chkval, sonhi, sonlo;
  su2double f, f1, f2, s1, s2;

  /*--- Tables from Stevens, 1972 ---*/
  su2double ffactr[186] = {0.181,0.100,0.196,0.122,0.212,
                           10.140,0.230,0.158,0.248,0.174,0.269,0.187,0.290,0.200,0.314,0.212,
                           20.339,0.222,0.367,0.232,0.396,0.241,0.428,0.250,0.463,0.259,0.500,
                           30.267,0.540,0.274,0.583,0.281,0.630,0.287,0.680,0.293,0.735,0.298,
                           40.794,0.303,0.857,0.308,0.926,0.312,1.000,0.316,1.080,0.319,1.170,
                           50.320,1.260,0.322,1.360,0.322,1.470,0.320,1.590,0.319,1.720,0.317,
                           61.850,0.314,2.000,0.311,2.160,0.308,2.330,0.304,2.520,0.300,2.720,
                           70.296,2.940,0.292,3.180,0.288,3.430,0.284,3.700,0.279,4.000,0.275,
                           84.320,0.270,4.670,0.266,5.040,0.262,5.440,0.258,5.880,0.253,6.350,
                           90.248,6.860,0.244,7.410,0.240,8.000,0.235,8.640,0.230,9.330,0.226,
                           110.10,0.222,10.90,0.217,11.80,0.212,12.70,0.208,13.70,0.204,14.80,
                           20.200,16.00,0.197,17.30,0.195,18.70,0.194,20.20,0.193,21.80,0.192,
                           323.50,0.191,25.40,0.190,27.40,0.190,29.60,0.190,32.00,0.190,34.60,
                           40.190,37.30,0.190,40.30,0.191,43.50,0.191,47.00,0.192,50.80,0.193,
                           554.90,0.194,59.30,0.195,64.00,0.197,69.10,0.199,74.70,0.201,80.60,
                           60.203,87.10,0.205,94.10,0.208,102.0,0.210,110.0,0.212,119.0,0.215,
                           7128.0,0.217,138.0,0.219,149.0,0.221,161.0,0.223,174.0,0.224,188.0,
                           80.225,203.0,0.226,219.0,0.227},
            bpl[54]     = {48.78,24.42,51.50,28.13,55.33,
                           133.36,60.20,40.00,65.06,45.25,66.76,49.00,70.28,54.25,73.34,
                           258.00,77.18,63.25,79.93,67.00,82.04,69.87,85.87,75.17,89.38,
                           378.87,94.65,84.17,98.39,87.90,103.65,93.15,107.39,96.90,110.27,
                           499.75,115.54,105.08,119.27,108.80,124.57,114.04,128.29,117.80,
                           5133.00,123.06,135.72,126.79,137.84,129.68,141.96,134.95,144.65,
                           6138.68},
            sonval[27]  = {0.3,0.4,0.6,1.0,1.5,2.0,3.0,4.0,6.0,8.0,
                           110.0,15.0,20.0,30.0,40.0,60.0,80.0,100.0,150.0,200.0,300.0,400.0,
                           2600.0,800.0,1000.0,1500.0,2000.0},
            slope1[27],
            slope2[27],
            slope3      = -2.0,
            slope4      = 4.0,
            spl[135];

  for(unsigned short i = 0; i < 27; i++){
    spl[i*5]   = 160.;
    spl[i*5+1] = bpl[i*2];
    spl[i*5+2] = bpl[i*2+1];
    spl[i*5+3] = bpl[i*2+1]-8.0;
    spl[i*5+4] = bpl[i*2+1];

    /*--- Calculate the variable slopes for the segments of contours, ---
      --- defined as follows:                                         ---
      ---   Between bands 1-19  : Variable slope stored in slope1     ---
      ---   Between bands 19-26 : Variable slope stored in slope2     ---
      ---   Between bands 26-31 : Constant                            ---
      ---   Between bands 31-35 : -2 dB/band (slope3)                 ---
      ---   Between bands 35-39 : Constant                            ---
      ---   Between bands 39-41 : 4 dB/band (slope4)                  ---*/
    slope1[i] = (spl[i*5+1]-spl[i*5])/18.;
    slope2[i] = (spl[i*5+2]-spl[i*5+1])/7.;
  }

  /*--- Now find sone value for all 1/3 oct bands ---*/
  unsigned short l;
  for(unsigned short i = 0; i < 41; i++){
    sone[i] = 0.0;
    uppt = 0.0;
    dnpt = 0.0;
    chkval = SPL_band[i];
    l = 0;

    if(i+1 == 1){
      if(chkval >= 160.){}
      else{}
    }

    if(i+1 == 19){
      if(chkval > 144.5){}
      else{l = 1;}
    }

    if(i+1 == 41){
      if(chkval > 138.5){}
      else{l = 4;}
    }

    for(unsigned short j = 0; j < 26; j++){
      dnpt = spl[j*5+l];
      uppt = spl[(j+1)*5+l];
      if((chkval > dnpt) && (chkval < uppt)){
        sonlo = sonval[j]; sonhi = sonval[j+1];
        sone[i] = sonlo + (chkval - dnpt)*(sonhi-sonlo)/(uppt-dnpt);
      }
      else{}
    }

    if((i+1) > 1 && (i+1) < 19){
      for(unsigned short j = 0; j < 26; j++){
        dnpt = spl[j*5] + slope2[j]*(i-18);
        uppt = spl[(j+1)*5] + slope2[j+1]*(i-18);
        if((j+1) == 26 && chkval > uppt){}

        if((chkval >= dnpt) && (chkval < uppt)){
          sonlo = sonval[j]; sonhi = sonval[j+1];
          sone[i] = sonlo + (chkval - dnpt)*(sonhi-sonlo)/(uppt-dnpt);
        }
        else{}
      }
    }
    else{}

    if((i+1) > 19 && (i+1) < 26){
      for(unsigned short j = 0; j < 26; j++){
        dnpt = spl[j*5+1] + slope1[j]*i;
        uppt = spl[(j+1)*5+1] + slope1[j+1]*i;
        if((j+1) == 26 && chkval > uppt){}

        if((chkval >= dnpt) && (chkval < uppt)){
          sonlo = sonval[j]; sonhi = sonval[j+1];
          sone[i] = sonlo + (chkval - dnpt)*(sonhi-sonlo)/(uppt-dnpt);
        }
        else{}
      }
    }
    else{}

    /*--- BETWEEN BANDS 26-31, THE CONTOURS ARE CONSTANT AT A GIVEN BAND PRES-SURE LEVEL FOR A GIVEN SONE VALUE ---*/
    if((i+1) >= 26 && (i+1) <= 31){
      if(chkval > 138.5){}

      for(unsigned short j = 0; j < 26; j++){
        dnpt = spl[j*5+2];
        uppt = spl[(j+1)*5+2];
        if(chkval >= dnpt && chkval <= uppt){
          sonlo = sonval[j]; sonhi = sonval[j+1];
          sone[i] = sonlo + (chkval - dnpt)*(sonhi-sonlo)/(uppt-dnpt);
        }
        else{}
      }
    }
    else{}

    /*--- BETWEEN BANDS 31-35, THE CONTOURS SLOPE DOWN WITH A SLOPE VALUE OF 2 dB/BAND (SLOPE3) - THIS IS USED IN THE FOLLOWING LOGIC ---*/
    if((i+1) > 31 && (i+1) < 35){
      for(unsigned short j = 0; j < 26; j++){
        dnpt = spl[j*5+2] + slope3*(i-30);
        uppt = spl[(j+1)*5+2] + slope3*(i-30);
        if((j+1) == 26 && chkval > uppt){}

        if(chkval >= dnpt && chkval <= uppt){
          sonlo = sonval[j]; sonhi = sonval[j+1];
          sone[i] = sonlo + (chkval - dnpt)*(sonhi-sonlo)/(uppt-dnpt);
        }
        else{}
      }
    }
    else{}

    /*--- BETWEEN BANDS 35-39, THE CONTOURS ARE CONSTANT AT A GIVEN BAND PRESSURE LEVEL FOR A GIVEN SONE VALUE ---*/
    if((i+1) >= 35 && (i+1) <= 39){
      if(chkval > 130.5){}

      for(unsigned short j = 0; j < 26; j++){
        dnpt = spl[j*5+3];
        uppt = spl[(j+1)*5+3];
        if(chkval >= dnpt && chkval <= uppt){
          sonlo = sonval[j]; sonhi = sonval[j+1];
          sone[i] = sonlo + (chkval - dnpt)*(sonhi-sonlo)/(uppt-dnpt);
        }
        else{}
      }
    }
    else{}

    /*--- BETWEEN BANDS 39-41, THE CONTOURS RISE WITH A CONSTANT SLOPE OF 4 dB/BAND (SLOPE4) - THIS VALUE IS USED IN THE FOLLOWING LOGIC ---*/
    if((i+1) > 39 && (i+1) < 41){
      for(unsigned short j = 0; j < 26; j++){
        dnpt = spl[j*5+3] + slope4*(i-38);
        uppt = spl[(j+1)*5+3] + slope4*(i-38);
        if((j+1) == 26 && chkval > uppt){}

        if(chkval >= dnpt && chkval <= uppt){
          sonlo = sonval[j]; sonhi = sonval[j+1];
          sone[i] = sonlo + (chkval - dnpt)*(sonhi-sonlo)/(uppt-dnpt);
        }
        else{}
      }
    }
    else{}

  }

  /*--- Now compute overall loudness based on all 1/3-oct bands ---*/
  sonmax = 0.0;
  for(unsigned short i = 0; i < 41; i++){
    sonsum += sone[i];
    if(sone[i] > sonmax) sonmax = sone[i];
  }

  /*--- Find F factor based on maximally loud 1/3-oct band ---*/
  if(sonmax >= 219.) f = 0.227;

  for(unsigned short i = 0; i < 92; i++){
    s1 = ffactr[i*2]; s2 = ffactr[(i+1)*2];
    if(sonmax >= s1 && sonmax <= s2){
      f1 = ffactr[i*2+1]; f2 = ffactr[(i+1)*2+1];
      f = f1 + (sonmax - s1)*(f2-f1)/(s2-s1);
      break;
    }
  }

  sonovr = sonmax + f*(sonsum-sonmax);

  /*--- PLdB = 32 + 9 log_2(sonovr) = 32 + 9*3.322 log_10(sonovr) ---*/
  if(sonovr > 0){
    PLdB[iPhi] = 32. + 29.897*log10(sonovr);
  }
  else{
    PLdB[iPhi] = 0.;
  }

}

void CBoom_AugBurgers::AcousticEnergy(unsigned short iPhi){

  cout << "Comupting acoustic energy." << endl;

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
  /*---Final signal and boom strength---*/
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
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
#endif
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

