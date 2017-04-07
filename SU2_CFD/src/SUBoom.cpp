#include "../include/SUBoom.hpp"
#include <cmath>
#include <algorithm>
#include <complex>
#include <valarray>
#include <iostream>
#include <fstream>

typedef std::complex<su2double> Complex;
typedef std::valarray<Complex> CArray;

#define N_PROF 50001

SUBoom::SUBoom(){

}

SUBoom::SUBoom(CSolver *solver, CConfig *config, CGeometry *geometry){

#ifdef HAVE_MPI
  int rank, iProcessor, nProcessor;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
#endif

  /*---Make sure to read in hard-coded values in the future!---*/

  /*---Flight variables---*/
  flt_h = 15240; // altitude [m]
  flt_M = config->GetMach();
  flt_psi = 0.;  // heading angle [deg]
  flt_gamma = 0.; // flight path angle [deg]

  /*---Atmosphere variables---*/
  atm_g = config->GetGamma();
  atm_noise_flag = 0;

  /*---Ray variables---*/
  int N_phi = 1;
  ray_N_phi = N_phi;
  ray_phi = new su2double[ray_N_phi];
  for(int i = 0; i < N_phi; i++){
    //ray_phi[i] = phi[i];
    ray_phi[i] = 0.0;
  }
  ray_r0 = 2.0;

  /*---Tolerance variables---*/
//  char cstr [200];
  string str;
  ifstream tolfile;

  if(rank == MASTER_NODE){
  tolfile.open("tols.in", ios::in);
  if (tolfile.fail()) {
    cout << "There is no tol.in file. Using default tolerances for boom propagation. " << endl;
    tol_dphi = 1.0E-3;
    tol_dr = 1.0E-3;
    tol_m = 1.0E6;
    tol_dp = 1.0E-6;
    tol_l = 1.0E-6;
  }
  else{
    tolfile >> str >> tol_dphi;
    tolfile >> str >> tol_dr;
    tolfile >> str >> tol_m;
    tolfile >> str >> tol_dp;
    tolfile >> str >> tol_l;
  }
  }

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
  //p_inf = Pressure_FreeStream;
  if(rank == MASTER_NODE){
  cout << "Pressure_Ref = " << Pressure_Ref << ", Pressure_FreeStream = " << Pressure_FreeStream << endl;
  }

  /*---Loop over boundary markers to select those to extract pressure signature---*/
  unsigned long nMarker = config->GetnMarker_All();
  unsigned long nSig = 0, sigCount=0;
  unsigned long panelCount=0;
  unsigned long iMarker, iVertex, iPoint;
  su2double x, y, z;
  su2double rho, rho_ux, rho_uy, rho_uz, rho_E, TKE;
  su2double ux, uy, uz, StaticEnergy, p;
  su2double *Coord;

  nSig = 0;
  nDim = geometry->GetnDim();
  for(iMarker = 0; iMarker < nMarker; iMarker++){
    if(config->GetMarker_All_KindBC(iMarker) == INTERNAL_BOUNDARY){
      sigCount++;
      for(iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++){
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        if(geometry->node[iPoint]->GetDomain()){
//            Coord = geometry->node[iPoint]->GetCoord();
//            x = SU2_TYPE::GetValue(Coord[0]);
//            p = solver->node[iPoint]->GetSolution(2*nDim+2)*Pressure_Ref - Pressure_FreeStream;
            //if(p >= 0.1 || p <= -0.1){
                panelCount++;
            //}
        }
      }
      nSig = panelCount;
    }
  }
  cout<<"Rank= "<<rank<<", nSig= "<<nSig<<endl;
//  signal.original_len = sigCount;
//  nSig = sigCount;

  unsigned long Buffer_Send_sigCount[1], *Buffer_Recv_sigCount = NULL;
  unsigned long totSig, iPanel;
  if (rank == MASTER_NODE) Buffer_Recv_sigCount= new unsigned long [nProcessor];

  Buffer_Send_sigCount[0]=nSig;
#ifdef HAVE_MPI
   SU2_MPI::Gather(&Buffer_Send_sigCount, 1, MPI_UNSIGNED_LONG, Buffer_Recv_sigCount, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD); //send the number of vertices at each process to the master
   SU2_MPI::Allreduce(&nSig,&totSig,1,MPI_UNSIGNED_LONG,MPI_SUM,MPI_COMM_WORLD); //find the max num of vertices over all processes
//   SU2_MPI::Reduce(&nSig,&totSig,1,MPI_UNSIGNED_LONG,MPI_SUM,MASTER_NODE,MPI_COMM_WORLD); //find the total num of vertices (panels)
#endif

  su2double *Buffer_Send_Press = new su2double [totSig];
  su2double *Buffer_Send_x = new su2double [totSig];
  unsigned long *Buffer_Send_GlobalIndex = new unsigned long [totSig];
  //zero send buffers
  for (int i=0; i <totSig; i++){
   Buffer_Send_Press[i]=0.0;
   Buffer_Send_x[i] = 0.0;
  }
  for (int i=0; i <totSig; i++){
   Buffer_Send_GlobalIndex[i]=0;
  }
  su2double *Buffer_Recv_Press = NULL;
  su2double *Buffer_Recv_x = NULL;
  unsigned long *Buffer_Recv_GlobalIndex = NULL;

  if (rank == MASTER_NODE) {
   Buffer_Recv_Press = new su2double [nProcessor*totSig];
   Buffer_Recv_x = new su2double [nProcessor*totSig];
   Buffer_Recv_GlobalIndex = new unsigned long[nProcessor*totSig];
  }

  /*---Extract signature---*/
  for(iMarker = 0; iMarker < nMarker; iMarker++){
    if(config->GetMarker_All_KindBC(iMarker) == INTERNAL_BOUNDARY){
      for(iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++){
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        if(geometry->node[iPoint]->GetDomain()){
          Coord = geometry->node[iPoint]->GetCoord();
          x = SU2_TYPE::GetValue(Coord[0]);
          //p = solver->node[iPoint]->GetSolution(2*nDim+2)*Pressure_Ref - Pressure_FreeStream;
          //if(p >= 0.1 || p <= -0.1){
            panelCount++;
            y = SU2_TYPE::GetValue(Coord[1]);
            z = 0.0;
            if(nDim == 3) z = SU2_TYPE::GetValue(Coord[2]);

            /*---Extract conservative flow data---*/
            rho = solver->node[iPoint]->GetSolution(nDim);
            rho_ux = solver->node[iPoint]->GetSolution(nDim+1);
            rho_uy = solver->node[iPoint]->GetSolution(nDim+2);
            if(nDim == 3) rho_uz = solver->node[iPoint]->GetSolution(nDim+3);
            rho_E = solver->node[iPoint]->GetSolution(2*nDim+1);
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
            ux = rho_ux/rho;
            uy = rho_uy/rho;
            uz = 0.0;
            if(nDim == 3) uz= rho_uz/rho;
            StaticEnergy =  rho_E/rho-0.5*(ux*ux+uy*uy+uz*uz)-TKE;
            p = (config->GetGamma()-1)*rho*StaticEnergy - Pressure_FreeStream/Pressure_Ref;

            Buffer_Send_Press[panelCount-1] = p;
            Buffer_Send_x[panelCount-1] = x;
            Buffer_Send_GlobalIndex[panelCount-1] = geometry->node[iPoint]->GetGlobalIndex();
            //signal.original_p[panelCount-1] = p;
            //signal.x[panelCount-1] = x;
            //sigFile << x << "\t" << p << endl;

            //PointID[panelCount-1] = geometry->node[iPoint]->GetGlobalIndex();
          //}
        }

      }
    }
  }

  /*---Send signal to MASTER_NODE---*/
 /* for(iPanel=0; iPanel<nSig; iPanel++){
    Buffer_Send_Press[iPanel] = Press[iPanel];
    Buffer_Send_x[iPanel] =
  }

  for (iPanel=0; iPanel<nSig; iPanel++){
    Buffer_Send_GlobalIndex[iPanel] = PointID[iPanel];
  }*/

#ifdef HAVE_MPI
  SU2_MPI::Gather(Buffer_Send_Press, totSig, MPI_DOUBLE, Buffer_Recv_Press,  totSig , MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_x, totSig, MPI_DOUBLE, Buffer_Recv_x,  totSig , MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_GlobalIndex,totSig, MPI_UNSIGNED_LONG, Buffer_Recv_GlobalIndex, totSig , MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
#endif

  if (rank == MASTER_NODE){
  ofstream sigFile;
//  SPRINTF (cstr, "signal_original.dat");
  sigFile.precision(15);
//  sigFile.open(cstr, ios::out);
  sigFile.open("signal_original.dat");
  sigFile << "# x, p" << endl;
  panelCount = 0;
  nPanel = totSig;
  signal.x = new su2double[nSig];
  signal.original_p = new su2double[nSig];
  signal.original_T = new su2double[nSig];
  PointID = new unsigned long[nSig];

  int Total_Index;

  for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
    for (iPanel = 0; iPanel < Buffer_Recv_sigCount[iProcessor]; iPanel++) {
      PointID[panelCount] = Buffer_Recv_GlobalIndex[iProcessor*totSig+iPanel ];

      /*--- Current index position and global index ---*/
      Total_Index  = iProcessor*totSig  + iPanel;
      signal.x[panelCount] =
      signal.original_p[panelCount] = Buffer_Recv_Press[Total_Index];

      /*--- Write to file---*/
      sigFile << scientific << Buffer_Recv_x[Total_Index] << "\t";
      sigFile << scientific <<  Buffer_Recv_Press[Total_Index]   << "\t";
      sigFile << endl;

      panelCount++;
      }
      }

  sigFile.close();

  /*---Initialize sensitivities---*/
  if(config->GetAD_Mode()){
  dJdU = new su2double* [nDim+3];
  for(int iDim = 0; iDim < nDim+3 ; iDim++){
      dJdU[iDim] = new su2double[nSig];
      for(int iPanel = 0;  iPanel< nSig; iPanel++){
          dJdU[iDim][iPanel] = 0.0;
      }
  }
  }

  }
  else{
    ray_N_phi = NULL;
    nSig = NULL;
  }

  if(rank == MASTER_NODE) cout << "Tot_nSig = " << totSig << endl;

}

SUBoom::~SUBoom(void){

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
  z = new su2double[N_PROF];
  a_of_z = new su2double[N_PROF];
  rho_of_z = new su2double[N_PROF];
  for(int i = 0; i < N_PROF; i++){
    z[i] = h0*su2double(i)/(su2double(N_PROF)-1.);
    AtmosISA(z[i], T, a_of_z[i], p, rho_of_z[i], g);
  }
}

void SUBoom::ScaleFactors(){
  int len;
  su2double *x;
  su2double L;
  su2double h0 = flt_h;
  su2double M_inf = flt_M;

  len = signal.original_len;
  x = signal.x;

  su2double min_x = x[0], max_x = x[len-1];

  for(int i = 0; i < len; i++){
      max_x = x[i];
    if(x[i] < min_x)
      min_x = x[i];
  }
  // TODO: get a/c length
  L = max_x - min_x;
  //L = 1.0;

  scale_L = L;    // [m]
  scale_T = L/(M_inf*a_inf);    // flow over aircraft [s]
  scale_p = p_inf;    // ambient [Pa]
  scale_m = scale_p/scale_T;    // slope of boom signal [Pa/s]
  scale_z = h0;    // altitude [m]

  scale_C1 = scale_p;    // [Pa]
  scale_C2 = scale_T;    // [s]

}

void SUBoom::GetAtmosphericData(su2double& a, su2double& rho, su2double& p, su2double h){
  /*---Get speed of sound, density, and pressure at an altitude---*/
  int i_h = -1;
  for(int i = 0; i < N_PROF; i++){
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
  a0 = a_of_z[N_PROF-1];

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
  //  su2double theta[N_PROF];
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
  x = new su2double[N_PROF];
  y = new su2double[N_PROF];
  t = new su2double[N_PROF];
  kx = new su2double[N_PROF];
  ky = new su2double[N_PROF];
  kz = new su2double[N_PROF];

  x_of_z = new su2double**[ray_N_phi];
  y_of_z = new su2double**[ray_N_phi];
  t_of_z = new su2double**[ray_N_phi];
  dxdt = new su2double**[ray_N_phi];
  dydt = new su2double**[ray_N_phi];
  dzdt = new su2double**[ray_N_phi];
  theta = new su2double**[ray_N_phi];
  for(int i = 0; i < ray_N_phi; i++){
    x_of_z[i] = new su2double*[4];
    y_of_z[i] = new su2double*[4];
    t_of_z[i] = new su2double*[4];
    dxdt[i] = new su2double*[4];
    dydt[i] = new su2double*[4];
    dzdt[i] = new su2double*[4];
    theta[i] = new su2double*[4];
    for(int j = 0; j < 4; j++){
      x_of_z[i][j] = new su2double[N_PROF];
      y_of_z[i][j] = new su2double[N_PROF];
      t_of_z[i][j] = new su2double[N_PROF];
      dxdt[i][j] = new su2double[N_PROF];
      dydt[i][j] = new su2double[N_PROF];
      dzdt[i][j] = new su2double[N_PROF];
      theta[i][j] = new su2double[N_PROF];
    }
  }

   for(int i = 0; i < ray_N_phi; i++){

    /*---Primary ray---*/
    data.c0 = ray_c0[i][0];
    data.nu = ray_nu[i][0];
    r0[0] = r0[1] = r0[2] = 0.;

    for(int j = N_PROF-1; j > -1; j--){
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

    kx = SplineGetDerivs(t, x, N_PROF);
    ky = SplineGetDerivs(t, y, N_PROF);
    kz = SplineGetDerivs(t, z, N_PROF);
    for(int ik = 0; ik < N_PROF; ik++){
      dxdt[i][0][ik] = kx[ik];
      dydt[i][0][ik] = ky[ik];
      dzdt[i][0][ik] = kz[ik];
    }
    /*jj = 0;
    while(jj < N_PROF-1){
    for(int ik = 0; ik < 1001; ik++){
        xtmp[ik] = x[ik+jj];
        ytmp[ik] = y[ik+jj];
        ztmp[ik] = z[ik+jj];
        ttmp[ik] = t[ik+jj];
    }
    kx = SplineGetDerivs(ttmp, xtmp, 1001);
    ky = SplineGetDerivs(ttmp, ytmp, 1001);
    kz = SplineGetDerivs(ttmp, ztmp, 1001);
    for(int ik = 0; ik < 1001; ik++){
      dxdt[i][0][ik+jj] = kx[ik];
      dydt[i][0][ik+jj] = ky[ik];
      dzdt[i][0][ik+jj] = kz[ik];
    }
    jj = jj + 1001-1;
}*/

    /*---Ray tube corners: {0, +dheading}---*/
    r0[0] = flt_heading[0]*tol_dr;
    r0[1] = flt_heading[1]*tol_dr;
    r0[2] = flt_heading[2]*tol_dr;
    for(int j = N_PROF-1; j > -1; j--){
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

    kx = SplineGetDerivs(t, x, N_PROF);
    ky = SplineGetDerivs(t, y, N_PROF);
    for(int ik = 0; ik < N_PROF; ik++){
      dxdt[i][1][ik] = kx[ik];
      dydt[i][1][ik] = ky[ik];
    }
    /*jj = 0;
    while(jj < N_PROF-1){
    for(int ik = 0; ik < 1001; ik++){
        xtmp[ik] = x[ik+jj];
        ytmp[ik] = y[ik+jj];
        ttmp[ik] = t[ik+jj];
    }
    kx = SplineGetDerivs(ttmp, xtmp, 1001);
    ky = SplineGetDerivs(ttmp, ytmp, 1001);
    for(int ik = 0; ik < 1001; ik++){
      dxdt[i][1][ik+jj] = kx[ik];
      dydt[i][1][ik+jj] = ky[ik];
    }
    jj = jj + 1001-1;
}*/

    /*---Ray tube corners: {+dphi, 0}---*/
    data.c0 = ray_c0[i][1];
    data.nu = ray_nu[i][1];
    r0[0] = r0[1] = r0[2] = 0.;

    for(int j = N_PROF-1; j > -1; j--){
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

    kx = SplineGetDerivs(t, x, N_PROF);
    ky = SplineGetDerivs(t, y, N_PROF);
    kz = SplineGetDerivs(t, z, N_PROF);
    for(int ik = 0; ik < N_PROF; ik++){
      dxdt[i][2][ik] = kx[ik];
      dydt[i][2][ik] = ky[ik];
      dzdt[i][2][ik] = kz[ik];
    }
    /*jj = 0;
    while(jj < N_PROF-1){
    for(int ik = 0; ik < 1001; ik++){
      xtmp[ik] = x[ik+jj];
      ytmp[ik] = y[ik+jj];
      ztmp[ik] = z[ik+jj];
      ttmp[ik] = t[ik+jj];
    }
    kx = SplineGetDerivs(ttmp, xtmp, 1001);
    ky = SplineGetDerivs(ttmp, ytmp, 1001);
    kz = SplineGetDerivs(ttmp, ztmp, 1001);
    for(int ik = 0; ik < 1001; ik++){
      dxdt[i][2][ik+jj] = kx[ik];
      dydt[i][2][ik+jj] = ky[ik];
      dzdt[i][2][ik+jj] = kz[ik];
    }
    jj = jj + 1001-1;
}*/

    /*---Ray tube corners: {+dphi, +dheading}---*/
    r0[0] = flt_heading[0]*tol_dr;
    r0[1] = flt_heading[1]*tol_dr;
    r0[2] = flt_heading[2]*tol_dr;
    for(int j = N_PROF-1; j > -1; j--){
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

    kx = SplineGetDerivs(t, x, N_PROF);
    ky = SplineGetDerivs(t, y, N_PROF);
    for(int ik = 0; ik < N_PROF; ik++){
      dxdt[i][3][ik] = kx[ik];
      dydt[i][3][ik] = ky[ik];
    }
    /*jj = 0;
    while(jj < N_PROF-1){
    for(int ik = 0; ik < 1001; ik++){
      xtmp[ik] = x[ik+jj];
      ytmp[ik] = y[ik+jj];
      ttmp[ik] = t[ik+jj];
    }
    kx = SplineGetDerivs(ttmp, xtmp, 1001);
    ky = SplineGetDerivs(ttmp, ytmp, 1001);
    for(int ik = 0; ik < 1001; ik++){
      dxdt[i][3][ik+jj] = kx[ik];
      dydt[i][3][ik+jj] = ky[ik];
    }
    jj = jj + 1001-1;
}*/

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
    ray_A[i] = new su2double[N_PROF];
  }

  /*---Loop over rays---*/
  for(int i = 0; i < ray_N_phi; i++){
    M = N_PROF;
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

  su2double f[N_PROF], t[N_PROF];
  su2double f0, f1, t0, t1, tmp, trat;
  //int flag;

  su2double a, b;
  su2double *ks;

  ray_t0 = new su2double[ray_N_phi];

  ks = new su2double[N_PROF];

  for(int i = 0; i < ray_N_phi; i++){
    for(int j = 0; j < N_PROF; j++){
      t[j] = t_of_z[i][0][j];
      f[j] = matchr(i, j, h_L, ray_r0);
    }

    /*---Get spline---*/
    ks = SplineGetDerivs(t, f, N_PROF);

    /*---Secant method---*/
    f0 = f[(N_PROF-1)/3];
    f1 = f[2*(N_PROF-1)/3];
    t0 = t[(N_PROF-1)/3];
    t1 = t[2*(N_PROF-1)/3];
    eps = 1.E6;
    while(eps > 1.E-8){
      tmp = t1 - f1*(t1-t0)/(f1-f0);
      t0 = t1;
      t1 = tmp;
      f0 = f1;
      // Find interval which contains t1
      int j_sp = N_PROF-1;
      for(int j = N_PROF-1; j > -1; j--){
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

  su2double t[N_PROF], A[N_PROF];
  su2double *cn;
  su2double *dadt, *drhodt, *dAdt;
  su2double *dcndt;
  su2double g = atm_g;

  cn = new su2double[N_PROF];

  dadt = new su2double[N_PROF];
  drhodt = new su2double[N_PROF];
  dAdt = new su2double[N_PROF];
  dcndt = new su2double[N_PROF];

  ray_C1 = new su2double*[ray_N_phi];
  ray_C2 = new su2double*[ray_N_phi];
  ray_dC1 = new su2double*[ray_N_phi];
  ray_dC2 = new su2double*[ray_N_phi];
  for(int i = 0; i < ray_N_phi; i++){
    for (int j = 0; j < N_PROF; j++){
      t[j] = t_of_z[i][0][j];
      A[j] = ray_A[i][j];
      cn[j] = ray_c0[i][0]*cos(theta[i][0][j]);
    }

    /*--- Spline interpolation of a, rho, A---*/
    dadt = SplineGetDerivs(t, a_of_z, N_PROF);
    drhodt = SplineGetDerivs(t, rho_of_z, N_PROF);
    dAdt = SplineGetDerivs(t, A, N_PROF);
    dcndt = SplineGetDerivs(t, cn, N_PROF);

    ray_C1[i] = new su2double[N_PROF];
    ray_C2[i] = new su2double[N_PROF];
    for(int j = 0; j < N_PROF; j++){
      ray_C1[i][j] = ((g+1.)/(2.*g))*a_of_z[j]/cn[j];
      //if(A[j] > 1E-16){
          ray_C2[i][j] = 0.5*((3./a_of_z[j])*dadt[j] + drhodt[j]/rho_of_z[j] - (2./cn[j])*dcndt[j] - (1./A[j])*dAdt[j]);
      //}
      //else{
        //ray_C2[i][j] = 0.5*((3./a_of_z[j])*dadt[j] + drhodt[j]/rho_of_z[j] - (2./cn[j])*dcndt[j]);
      //}
      //cout << "A = " << A[j] << ", dAdt = " << dAdt[j] << ", C2 = " << ray_C2[i][j] << endl;
    }
    ray_dC1[i] = SplineGetDerivs(t, ray_C1[i], N_PROF);
    ray_dC2[i] = SplineGetDerivs(t, ray_C2[i], N_PROF);

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
    if(mm[i] > tol_m/scale_m || mm[i] < -tol_m/scale_m){  // shock present
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
    //else if(mm[i] < -tol_m){  // "expansion shock" present
    //  /*---Remove segment i---*/
    //  ll[i] = tol_l;
    //  mm[i] = (pp[1][i] - pp[0][i])/ll[i];
    //}
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

su2double **WaveformToPressureSignal(su2double fvec[], int M, int &Msig){
  su2double TT[2][M], pp[2][M];
  su2double m[M], dp[M], l[M];
  su2double **sig;

  /*---Get m, dp, l from fvec---*/
  for(int i = 0; i < M; i++){
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
  for(int seg = 1; seg < M; seg++){
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
    for(int j = 0; j < 2; j++){sig[j] = new su2double[M];}

    for(int j = 0; j < M; j++){
        sig[0][j] = pp[0][j];
        sig[1][j] = pp[1][j];

    }
  }
  else{
    /*---Build 2-D vector---*/
    sig = new su2double*[2];
    for(int j = 0; j < 2; j++){sig[j] = new su2double[2*M];}

    /*---First segment---*/
    sig[0][0] = TT[0][0];
    sig[0][1] = TT[1][0];
    sig[1][0] = pp[0][0];
    sig[1][1] = pp[1][0];

    int i = 2;
    for(int seg = 1; seg < M; seg++){
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
    if(l[i] <= tol_l/scale_T || m[i] >= tol_m/scale_m || m[i] <= -tol_m/scale_m){
      /*---Record pressure gap---*/
      current_signal = WaveformToPressureSignal(fvec, N, Msig);
      dp_seg = dp[i] + (current_signal[1][i] - current_signal[0][i]);
      /*---Add to next segment if needed---*/
      if(dp_seg > tol_dp || dp_seg < -tol_dp){
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
//  int len = signal.original_len;
  su2double t0, tf, dt;
  int j0;
  RayData data;
  su2double *fvec;
  su2double *f;
  su2double **ground_signal;
  int M = signal.M;

  for(int i = 0; i < ray_N_phi; i++){
    t0 = ray_t0[i];

    /*---Assemble f vector and ray data for integration---*/
    signal.fvec = new su2double[3*signal.M];
    for(int j = 0; j < 3*signal.M; j++){
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
    for(int j = N_PROF-1; j > 0; j--){
        if(t_of_z[i][0][j] >= ray_t0[i]){
            j0 = j+1;
            break;
        }
    }
    for(int j = j0; j > 0; j--){
      /*---Integrate---*/
      int j0_dat;
      if(j >= N_PROF-4) j0_dat = N_PROF-4;
      else if(j < 4) j0_dat = 4;
      else j0_dat = j;
      for(int jj = 0; jj < 9; jj++){
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
    for(int ii = 0; ii < 2; ii++){
      ground_signal[ii] = new su2double[2*M];
    }
    int Msig = M;
    ground_signal = WaveformToPressureSignal(fvec, M, Msig);

    signal.final_T = new su2double[Msig];
    signal.final_p = new su2double[Msig];
    signal.final_M = Msig;

    /*---Final signal and boom strength---*/
//    char cstr [200];
    ofstream sigFile;
//    SPRINTF (cstr, "signal_final.dat");
    sigFile.precision(15);
//    sigFile.open(cstr, ios::out);
    sigFile.open("signal_final.dat", ios::out);
    sigFile << "# T, p" << endl;
    p_max = -1E10;
    p_int2 = 0;
    for(int j = 0; j < Msig; j++){
      signal.final_T[j] = ground_signal[0][j]*scale_T;
      signal.final_p[j] = ground_signal[1][j]*scale_p;
      sigFile << signal.final_T[j] << "\t" << signal.final_p[j] << endl;
      if(signal.final_p[j] > p_max) p_max = signal.final_p[j];
      if(j > 0) p_int2 = p_int2 + 0.5*(signal.final_p[j]+signal.final_p[j-1])/(signal.final_T[j]-signal.final_T[j]);
    }
    sigFile.close();
    p_rise = signal.final_p[0];
    if(signal.final_p[0] > -signal.final_p[Msig-1]) p_rise2 = signal.final_p[0];
    else p_rise2 = signal.final_p[Msig-1];
    cout << "Scale_T = " << scale_T << ", Scale_p = " << scale_p << endl;
    cout << "p_rise = " << p_rise << ", p_max = " << p_max << endl;
  }

  /*---Clean up---*/
  for(int i = 0; i < ray_N_phi; i++){
  for(int j = 0; j < 4; j++){
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

void SUBoom::WriteSensitivities(CSolver *solver, CConfig *config, CGeometry *geometry){
  unsigned long iVar, iPanel, Max_nPanel, Tot_nPanel, nVar, Global_Index ;
  ofstream Boom_AdjointFile;

  int rank, iProcessor, nProcessor;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
  unsigned long Buffer_Send_nPanel[1], *Buffer_Recv_nPanel = NULL;

  if (rank == MASTER_NODE) Buffer_Recv_nPanel= new unsigned long [nProcessor];

  Buffer_Send_nPanel[0]=nPanel;
#ifdef HAVE_MPI
  SU2_MPI::Gather(&Buffer_Send_nPanel, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nPanel, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD); //send the number of vertices at each process to the master
  SU2_MPI::Allreduce(&nPanel,&Max_nPanel,1,MPI_UNSIGNED_LONG,MPI_MAX,MPI_COMM_WORLD); //find the max num of vertices over all processes
  SU2_MPI::Reduce(&nPanel,&Tot_nPanel,1,MPI_UNSIGNED_LONG,MPI_SUM,MASTER_NODE,MPI_COMM_WORLD); //find the total num of vertices (panels)
#endif

  nVar = nDim+3;

  /* pack sensitivity values in each processor and send to root */
  su2double *Buffer_Send_dJdU = new su2double [Max_nPanel*nVar];
  unsigned long *Buffer_Send_GlobalIndex = new unsigned long [Max_nPanel];
  //zero send buffers
  for (int i=0; i <Max_nPanel*nVar; i++){
    Buffer_Send_dJdU[i]=0.0;
  }
  for (int i=0; i <Max_nPanel; i++){
    Buffer_Send_GlobalIndex[i]=0;
  }
  su2double *Buffer_Recv_dJdU = NULL;
  unsigned long *Buffer_Recv_GlobalIndex = NULL;

  if (rank == MASTER_NODE) {
    Buffer_Recv_dJdU = new su2double [nProcessor*Max_nPanel*nVar];
    Buffer_Recv_GlobalIndex = new unsigned long[nProcessor*Max_nPanel];
  }

  for (iVar=0; iVar<nVar; iVar++){
      for(iPanel=0; iPanel<nPanel; iPanel++){
          Buffer_Send_dJdU[iVar*nPanel+iPanel] = dJdU[iVar][iPanel];
        }
    }

  for (iPanel=0; iPanel<nPanel; iPanel++){
     Buffer_Send_GlobalIndex[iPanel] = PointID[iPanel];
  }


#ifdef HAVE_MPI
  SU2_MPI::Gather(Buffer_Send_dJdU, Max_nPanel*nVar, MPI_DOUBLE, Buffer_Recv_dJdU,  Max_nPanel*nVar , MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_GlobalIndex,Max_nPanel, MPI_UNSIGNED_LONG, Buffer_Recv_GlobalIndex, Max_nPanel , MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
#endif

  /* root opens a file at each time step and write out the merged dJdU values at that time step into the file */
  if (rank == MASTER_NODE){
//  char cstr [200];

//  SPRINTF (cstr, "Adj_Boom.dat");
  Boom_AdjointFile.precision(15);
//  Boom_AdjointFile.open(cstr, ios::out);
  Boom_AdjointFile.open("Adj_Boom.dat", ios::out);

  /*--- Loop through all of the collected data and write each node's values ---*/
  unsigned long Total_Index;
  for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
    for (iPanel = 0; iPanel < Buffer_Recv_nPanel[iProcessor]; iPanel++) {
        Global_Index = Buffer_Recv_GlobalIndex[iProcessor*Max_nPanel+iPanel ];
        Boom_AdjointFile  << scientific << Global_Index << "\t";

       for (iVar = 0; iVar < nVar; iVar++){
         /*--- Current index position and global index ---*/
         Total_Index  = iProcessor*Max_nPanel*nVar + iVar*Buffer_Recv_nPanel[iProcessor]  + iPanel;

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
/*  for (int i=0; i<nDim+3; i++){
    delete [] dJdU[i];
  }
  delete [] dJdU;
  delete [] PointID;*/

  cout << "\nFinished writing boom adjoint file." << endl;
}
