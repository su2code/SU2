/*!
 * \file postprocessing_structure.cpp
 * \brief Source file of the post-processing structure.
 * \author T. Albring, Beckett Y. Zhou
 * \version 4.3.0 "Cardinal"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2016 SU2, the open-source CFD code.
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

#include "../include/postprocessing_structure.hpp"
#include "../../Common/include/math_op_structure.hpp"
#include <fstream>
#include <iostream>
//#include "../include/gsl_sf_bessel.h"
//#include <cmath>
#include <complex>
//#include "../include/gsl_errno.h"
//#include "../include/gsl_fft_real.h"
//#include "../include/gsl_fft_halfcomplex.h"
#include <valarray>
//#include "fft_r2.cpp"
//#include "bessel_new.c"
#include <string.h>
typedef std::complex<su2double> Complex;
typedef std::valarray<Complex> CArray;

FWHSolver::FWHSolver(CConfig *config,CGeometry *geometry) {

    unsigned long  i, nMarker,iMarker,panelCount, end_iter, start_iter, iVertex,iPoint,  iSample,iPanel,iObserver, iDim;
  // su2double FreeStreamPressure;
   su2double FreeStreamTemperature;
    su2double pi=3.141592653589793;
    su2double *Coord,   *Normal;
    su2double  x, y, z, nx, ny, nz,  Area;
    su2double R = 287.058;
    su2double CheckNormal=0.0;
    unsigned long nFWH, maxFWHcount;
    nDim = geometry->GetnDim();
#ifdef HAVE_MPI
  int rank, nProcessor;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
     MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
#endif

    TimeDomain3D = config->GetTimeDomain3D();
    UseAnalytic = false;
    if (rank==MASTER_NODE){
    if (TimeDomain3D){
        if(nDim==2){
        cout<<endl<<"***************  WARNING!!! Time domain FWH implementation NOT available in 2D!!! Use frequency domain!!!  ***************"<<endl;
          }else{
         cout<<endl<<"-------------- Initiating 3D FWH Solver in Time Domain --------------"<<endl;
          }

      }else{
        if(nDim==2){
        cout<<endl<<"-------------- Initiating 2D FWH Solver in Frequency Domain --------------"<<endl;
          }else{
         cout<<endl<<"-------------- Initiating 3D FWH Solver in Frequency Domain --------------"<<endl;
          }
      }
          cout<<endl<<"-----------------------------------------VARIABLE SAMPLING FREQUENCY TEST --------------------------------------"<<endl;
      }


    totFWH = 0;

    SPL = 0.0;
    end_iter  = config->GetnExtIter();
    start_iter  =  config->GetUnst_RestartIter();
    SamplingFreq = config->GetWrt_Sol_Freq_DualTime();    //Sampling Freq: defined as the number of dt's between each snapsho (restart file)
    nSample  =  config->GetIter_Avg_Objective()/SamplingFreq;


    /* Setting Observer locations -- change to config option later! */
    string  text_line;
    ifstream Observer_LocationFile;
    string filename = "Observer_Locations.dat";
    Observer_LocationFile.open(filename.data() , ios::in);
    if (Observer_LocationFile.fail()) {
        cout << "There is no file!!! " <<  filename.data()  << "."<< endl;
      exit(EXIT_FAILURE);
    }
    getline (Observer_LocationFile, text_line);
    istringstream point_line(text_line);
     point_line >> nObserver ;
    //nObserver = 128;
    Observer_Locations = new su2double* [nObserver];
    closet_coord_AllObs = new su2double* [nObserver];
    for(iObserver = 0;  iObserver <nObserver  ;  iObserver++)
    {
       Observer_Locations[iObserver] = new su2double[nDim];
       closet_coord_AllObs[iObserver] = new su2double[nDim];
       for (iDim=0; iDim < nDim; iDim++){
         Observer_Locations[ iObserver][iDim]= 0.0;
         closet_coord_AllObs[ iObserver][iDim]= 0.0;
       }
    }

    iObserver=0;
  while (getline (Observer_LocationFile, text_line)) {
        istringstream point_line(text_line);
        if (nDim==2){
        point_line >> Observer_Locations[iObserver][0]>> Observer_Locations[iObserver][1];
          }
        if (nDim==3){
        point_line >> Observer_Locations[iObserver][0]>> Observer_Locations[iObserver][1]>> Observer_Locations[iObserver][2];
          }
        iObserver++;
    }





    M = config->GetMach();
    beta_sq = 1-M*M;

    FreeStreamPressure=config->GetPressure_FreeStream();
    FreeStreamTemperature = config->GetTemperature_FreeStream();
    FreeStreamDensity = FreeStreamPressure/R/FreeStreamTemperature;

    a_inf = sqrt(config->GetGamma()*FreeStreamPressure / FreeStreamDensity);
    AOA = config->GetAoA();

      U1 = M*a_inf*cos(AOA*pi/180) ;
      U2 = M*a_inf*sin(AOA*pi/180) ;
      U3 = 0.0;    //FIX THIS LATER!!!


    if (UseAnalytic){
        FreeStreamPressure=0.0;
        FreeStreamDensity=1.225;
        a_inf=340.0;
        Amp_a = 1.0;
//        T_a = (nSample-1)*config->GetDelta_UnstTime();
        T_a = 0.1;
        Freq_a = 2*M_PI/T_a;
        U1=0.0; U2=0.0; U3=0.0; M=0.0;
      }


  if (rank==MASTER_NODE && UseAnalytic) cout<<endl<<std::setprecision(8)<<"dt= "<< config->GetDelta_UnstTime()<<", T_a= "<< T_a<<", Freq_a= "<< Freq_a<<endl;
  if (rank==MASTER_NODE) cout<<U1<<",  "<<U2<<",  "<<M<<",  "<<a_inf<<",  "<<FreeStreamPressure<<", "<<FreeStreamDensity<<endl;

        nPanel =  0;
        panelCount = 0;
        maxFWHcount = 0;

    nMarker      = config->GetnMarker_All();
   for (iMarker = 0; iMarker < nMarker; iMarker++){

     /* --- Loop over boundary markers to select those on the FWH surface --- */
       if (config->GetMarker_All_KindBC(iMarker) == INTERNAL_BOUNDARY) {
      //  cout<<"Internal Boundary Detected"<<endl;
    //            cout<<"Rank: "<<rank<<", "<< config->GetMarker_All_TagBound(iMarker)<<endl;
           size_t last_index = config->GetMarker_All_TagBound(iMarker).find_last_not_of("0123456789");
         string result = config->GetMarker_All_TagBound(iMarker).substr(last_index + 1);
                 int f = std::strtol(result.c_str(),NULL, 10);
   //     cout<<"Rank: "<<rank<<", "<< config->GetMarker_All_TagBound(iMarker)<<", trailing digit="<< result <<endl;

        cout<<"Rank: "<<rank<<", "<< config->GetMarker_All_TagBound(iMarker)<<", trailing digit="<< f <<", coeff= "<<1.0/f<<endl;

              if (f>maxFWHcount){
                maxFWHcount = f;
              }
                //FWHcount++;

//       if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY) {

         for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++){
             iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
               if( geometry->node[iPoint]->GetDomain()){

                   Coord = geometry->node[iPoint]->GetCoord();
                   Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
                   Area  = 0.0; for ( iDim = 0; iDim < nDim; iDim++)   Area += Normal[iDim]*Normal[iDim];  Area  = sqrt( Area );

                   x  = SU2_TYPE::GetValue(Coord[0]);                                                                     // x
                   y  = SU2_TYPE::GetValue(Coord[1]);                                                                     // y
                   z  = 0.0;
                   if (nDim==3) z = SU2_TYPE::GetValue(Coord[2]);

                   nx = Normal[0]/Area  ;                                                                                 // n_x
                   ny = Normal[1]/Area  ;                                                                                 // n_y
                   nz = 0.0;
                   if (nDim==3)  nz = Normal[2]/Area  ;

        //           CheckNormal =  x*nx+y*ny+z*nz;
               //    if(CheckNormal<0) nx=-nx; ny=-ny;

                //  if (CheckNormal>0){
                     panelCount++;
                //   }
  //   cout<<"Normal="<<CheckNormal<<", x="<<x<<", y="<<y<<", nx="<<nx<<", ny="<<ny<<endl;
                 }
           }
          nPanel = panelCount;
     }

    }
 nFWH = maxFWHcount;
 cout<<"Process "<<rank<<" contains "<<nPanel <<" Panels."<<endl;
 //cout<<"Rank= "<<rank<<", nFWH= "<< nFWH<<endl;
 //communicate the global total count of the FWH surfaces to all processors
 //unsigned long Buffer_Send_nFWH[1], *Buffer_Recv_nFWH = NULL;
 //if (rank == MASTER_NODE) Buffer_Recv_nFWH= new unsigned long [nProcessor];

//  Buffer_Send_nFWH[0]=nFWH;
#ifdef HAVE_MPI
 //  SU2_MPI::Gather(&Buffer_Send_nFWH, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nFWH, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD); //send the number of vertices at each process to the master
   SU2_MPI::Allreduce(&nFWH,&totFWH,1,MPI_UNSIGNED_LONG,MPI_MAX,MPI_COMM_WORLD); //find the max num of vertices over all processes
 //  SU2_MPI::Reduce(&nPanel,&Tot_nPanel,1,MPI_UNSIGNED_LONG,MPI_SUM,MASTER_NODE,MPI_COMM_WORLD); //find the total num of vertices (panels)
#endif

  cout<<"Highest FWH Surface Index Found in Process "<<rank<<": "<< nFWH<<", Total Number of FWH Surfaces for All Processes= "<< totFWH <<endl;



// //  Interpolation* Interp_container = new Interpolation() ;


//  //Interp_container->spline_pchip_set ( int n, double x[], double f[], double d[] );
//  unsigned long n_x = 9;
//  unsigned long n_xe = 33;

////  double del_nx = M_PI/4.0;
////  double del_nx_e = M_PI/16.0;
////  double sinx[n_x], vec_x[n_x], derivative[n_x], sinx_e[n_xe], vec_x_e[n_xe],sinx_exact[n_xe];
////  ofstream data_file ;
////  data_file.open("data");
////  for (int i = 0; i < n_x; i++)
////    {
////      vec_x[i] = i*del_nx;
////      sinx[i]=sin(i*del_nx);
////      derivative[i] = 0.0;
////      cout<<i<<", "<<vec_x[i]<<", "<<sinx[i]<<endl;
////      data_file<<std::setprecision(15) << vec_x[i]<<", "<<sinx[i]<<endl;
////    }
////      cout<<"*****************************************"<<endl;
////  for (int i = 0; i < n_xe; i++)
////    {
////      vec_x_e[i] = i*del_nx_e;
////      sinx_e[i]=0.0;
////      sinx_exact[i]=sin(vec_x_e[i]);
////      cout<<i<<", "<<vec_x_e[i]<<", "<<sinx_e[i]<<endl;
////    }

////  Interp_container->spline_pchip_set ( n_x , vec_x , sinx ,  derivative );
////  Interp_container->spline_pchip_val ( n_x , vec_x , sinx ,  derivative  ,n_xe,vec_x_e,sinx_e);

////  cout<<"*****************************************"<<endl;
////  ofstream interp_file ;
////  interp_file.open("interp_PCHIP");
////for (int i = 0; i < n_xe; i++)
////{
////     interp_file << std::setprecision(15) << vec_x_e[i]<<", "<<sinx_e[i]<<", "<<sinx_exact[i]<<endl;

////}


//su2double del_nx = M_PI/4.0;
//su2double del_nx_e = M_PI/16.0;
//vector <su2double> sinx(n_x), vec_x(n_x), derivative(n_x),sinx_exact(n_xe);
//su2double sinx_e=0.0, vec_x_e=0.0;

//ofstream data_file ;
//su2double yp1=10.0e31;su2double ypn=10.0e31;
//data_file.open("data");
//for (int i = 0; i < n_x; i++)
//  {
//    vec_x[i] = i*del_nx;
//    sinx[i]=sin(i*del_nx);
//    derivative[i] = 0.0;
//    cout<<i<<", "<<vec_x[i]<<", "<<sinx[i]<<endl;
//    data_file<<std::setprecision(15) << vec_x[i]<<", "<<sinx[i]<<endl;
//  }
//    cout<<"*****************************************"<<endl;
//for (int i = 0; i < n_xe; i++)
//  {
//   // vec_x_e[i] = i*del_nx_e;
//  //  sinx_e[i]=0.0;
//    sinx_exact[i]=sin(i*del_nx_e);
////    cout<<i<<", "<<vec_x_e[i]<<", "<<sinx_e[i]<<endl;
//  }


////void CGeometry::SetSpline(vector<su2double> &x, vector<su2double> &y, unsigned long n, su2double yp1, su2double ypn, vector<su2double> &y2) {
//geometry->SetSpline(vec_x, sinx, n_x, yp1, ypn, derivative);
////GetSpline(vector<su2double>&xa, vector<su2double>&ya, vector<su2double>&y2a, unsigned long n, su2double x) {


////Interp_container->spline_pchip_set ( n_x , vec_x , sinx ,  derivative );
////Interp_container->spline_pchip_val ( n_x , vec_x , sinx ,  derivative  ,n_xe,vec_x_e,sinx_e);

//cout<<"*****************************************"<<endl;
//ofstream interp_file ;
//interp_file.open("interp_Cubic");
//for (int i = 0; i < n_xe; i++)
//{
//    vec_x_e = i*del_nx_e;
//    sinx_e=geometry->GetSpline(vec_x, sinx, derivative, n_x, vec_x_e);
//   interp_file << std::setprecision(15) << vec_x_e<<", "<<sinx_e<<", "<<sinx_exact[i]<<endl;

//}







   dJdU = new su2double** [nDim+3];
   for(int iDim = 0; iDim < nDim+3 ; iDim++)
   {
       dJdU[iDim] = new su2double*[nPanel];
       for(int iPanel = 0;  iPanel< nPanel; iPanel++)
       {
           dJdU[iDim][iPanel] = new su2double [nSample];
           for (int iSample=0; iSample <nSample; iSample++){
              dJdU[iDim][iPanel][iSample]= 0.0;
            }
       }
   }

   if(!TimeDomain3D){
   G = new complex <su2double>** [nObserver ];
   dGdy1 = new complex <su2double>** [nObserver ];
   dGdy2 = new complex <su2double>** [nObserver ];
   dGdy3 = new complex <su2double>** [nObserver ];
   for(iObserver = 0; iObserver < nObserver ; iObserver++)
   {
       G[iObserver] = new complex <su2double>*[nPanel];
       dGdy1[iObserver] = new complex <su2double>*[nPanel];
       dGdy2[iObserver] = new complex <su2double>*[nPanel];
       dGdy3[iObserver] = new complex <su2double>*[nPanel];
       for(iPanel = 0;  iPanel< nPanel; iPanel++)
       {
           G[iObserver][iPanel] = new complex <su2double>[nSample];
           dGdy1[iObserver][iPanel] = new complex <su2double>[nSample];
           dGdy2[iObserver][iPanel] = new complex <su2double>[nSample];
           dGdy3[iObserver][iPanel] = new complex <su2double>[nSample];
           for (iSample=0; iSample < nSample; iSample++){
              G[iObserver][iPanel][iSample]= 0.0;
              dGdy1[iObserver][iPanel][iSample]= 0.0;
              dGdy2[iObserver][iPanel][iSample]= 0.0;
              dGdy3[iObserver][iPanel][iSample]= 0.0;
            }
       }
   }

   fpp = new complex <su2double>* [nObserver];
   fpp_r = new su2double* [nObserver];
   fpp_i = new su2double* [nObserver];
   fpp_r_root = new su2double* [nObserver];
   fpp_i_root = new su2double* [nObserver];
   pp = new complex <su2double>* [nObserver];
   pp_CFD = new su2double* [nObserver];
   pp_CFD_mean = new su2double [nObserver];
   //Fr = new su2double* [nObserver];
   for(iObserver = 0; iObserver < nObserver ; iObserver++)
   {
           fpp[iObserver] = new complex <su2double>[nSample];
           fpp_r_root[iObserver] = new su2double [nSample];
           fpp_i_root[iObserver] = new su2double [nSample];
           fpp_r[iObserver] = new su2double [nSample];
           fpp_i[iObserver] = new su2double [nSample];
           pp[iObserver] = new complex <su2double> [nSample];
           pp_CFD[iObserver] = new su2double [nSample];        
           pp_CFD_mean [iObserver]=0.0;
          // Fr[iObserver] = new su2double [nSample];
           for (iSample=0; iSample < nSample; iSample++){
              fpp[iObserver][iSample]= 0.0;
              fpp_r_root[iObserver][iSample]= 0.0;
              fpp_i_root[iObserver][iSample]= 0.0;
              fpp_r[iObserver][iSample]= 0.0;
              fpp_i[iObserver][iSample]= 0.0;
//              fpp_r[iObserver] = new su2double [nSample];
//              fpp_i[iObserver] = new su2double [nSample];
              pp[iObserver][iSample]= 0.0;
              pp_CFD[iObserver][iSample]= 0.0;
            //  Fr[iObserver][iSample]= 0.0;
            }

   }
   }

   pp_CFD = new su2double* [nObserver];
   pp_CFD_mean = new su2double [nObserver];
   for(iObserver = 0; iObserver < nObserver ; iObserver++)
   {

           pp_CFD[iObserver] = new su2double [nSample];
           pp_CFD_mean [iObserver]=0.0;
           for (iSample=0; iSample < nSample; iSample++){
              pp_CFD[iObserver][iSample]= 0.0;
            }

   }





       surface_geo = new su2double* [nPanel];
       for(iPanel = 0;  iPanel< nPanel; iPanel++)
       {
           if (TimeDomain3D){
               surface_geo[iPanel] = new su2double[2*nDim+1+4*nObserver];
               for (i=0; i < 2*nDim+1+4*nObserver; i++){
                  surface_geo[iPanel][i]= 0.0;
                 }
             }else{
           surface_geo[iPanel] = new su2double[2*nDim+1];
           for (i=0; i < 2*nDim+1; i++){
              surface_geo[iPanel][i]= 0.0;
             }
             }
       }


       if(TimeDomain3D){

          t_Obs = new su2double** [nObserver ];
          pp_interp = new su2double** [nObserver ];
          for(iObserver = 0; iObserver<nObserver ; iObserver++)
               {
                   t_Obs[iObserver] = new su2double*[nPanel];
                   pp_interp[iObserver] = new su2double*[nPanel];
                   for(iPanel = 0;  iPanel< nPanel; iPanel++)
                   {
                       t_Obs[iObserver][iPanel] = new su2double[nSample];
                       pp_interp[iObserver][iPanel] = new su2double[nSample];

                       for (iSample=0; iSample < nSample; iSample++){

                          t_Obs[iObserver][iPanel][iSample]= 0.0;
                          pp_interp[iObserver][iPanel][iSample]= 0.0;
                        }
                   }
               }
          pp_TimeDomain = new su2double* [nObserver];
          pp_TimeDomain_root = new su2double* [nObserver];
          t_interp = new su2double* [nObserver];
          r_minmax = new su2double* [nObserver];
          for(iObserver = 0; iObserver<nObserver ; iObserver++)
               {
                       pp_TimeDomain[iObserver]= new su2double[nSample];
                       pp_TimeDomain_root[iObserver]= new su2double[nSample];
                       t_interp[iObserver]= new su2double[nSample];
                       r_minmax[iObserver]= new su2double[2];
                       r_minmax[iObserver][0]= 0.0; r_minmax[iObserver][1]= 0.0;
                       for (iSample=0; iSample < nSample; iSample++){
                          pp_TimeDomain[iObserver][iSample]= 0.0;
                          pp_TimeDomain_root[iObserver][iSample]= 0.0;
                          t_interp[iObserver][iSample]= 0.0;
                        }

               }

//          t_interp = new su2double [nSample];
//          for (iSample=0; iSample < nSample; iSample++){
//             t_interp[iSample] = 0.0;
//           }
          Fr = new su2double** [nObserver ];
          pp_ret = new su2double** [nObserver ];
          Fr_mean = new su2double* [nObserver ];
          for(iObserver = 0; iObserver<nObserver ; iObserver++)
          {
              Fr[iObserver] = new su2double*[nPanel];
              pp_ret[iObserver] = new su2double*[nPanel];
              Fr_mean[iObserver] = new su2double[nPanel];
              for(iPanel = 0;  iPanel< nPanel; iPanel++)
              {
                  Fr[iObserver][iPanel] = new su2double[nSample];
                  pp_ret[iObserver][iPanel] = new su2double[nSample];
                  Fr_mean[iObserver][iPanel] = 0.0;
                  for (iSample=0; iSample < nSample; iSample++){
                     Fr[iObserver][iPanel][iSample]= 0.0;
                     pp_ret[iObserver][iPanel][iSample]= 0.0;
                   }
              }
          }

          }

      if(!TimeDomain3D){
       F = new su2double** [nDim ];
       F_mean = new su2double* [nDim];
       fF = new complex <su2double>**  [nDim ];
       for(iDim = 0; iDim < nDim ; iDim++)
       {
           F[iDim] = new su2double*[nPanel];
           F_mean[iDim] = new su2double[nPanel];
           fF[iDim] = new complex <su2double>*[nPanel];
           for(iPanel = 0;  iPanel< nPanel; iPanel++)
           {
               F[iDim][iPanel] = new su2double[nSample];
               F_mean[iDim][iPanel]=0.0;
               fF[iDim][iPanel] = new complex <su2double>[nSample];
               for (iSample=0; iSample < nSample; iSample++){
                  F[iDim][iPanel][iSample]= 0.0;
                  fF[iDim][iPanel][iSample]= 0.0;
                }
           }
       }
        }



       if (!TimeDomain3D){
       fQ = new complex <su2double>* [nPanel];
       fF1 = new complex <su2double>* [nPanel];
       fF2 = new complex <su2double>* [nPanel];
       for(iPanel = 0;  iPanel< nPanel; iPanel++)
       {
          fQ [iPanel] = new complex <su2double>[nSample];
          fF1 [iPanel] = new complex <su2double>[nSample];
          fF2 [iPanel] = new complex <su2double>[nSample];
          for (iSample=0; iSample < nSample; iSample++){
            fQ[iPanel][iSample]= 0.0;
            fF1[iPanel][iSample]= 0.0;
            fF2[iPanel][iSample]= 0.0;
          }
       }
       }


       Q = new su2double* [nPanel];
       for(iPanel = 0;  iPanel< nPanel; iPanel++)
       {
          Q [iPanel] = new su2double[nSample];
          for (iSample=0; iSample < nSample; iSample++){
            Q[iPanel][iSample]= 0.0;
          }
       }


       Q_mean = new su2double [nPanel];
       F1_mean = new su2double [nPanel];
       F2_mean = new su2double [nPanel];
       PointID = new unsigned long [nPanel];

       for(iPanel = 0;  iPanel< nPanel; iPanel++)
       {
          Q_mean [iPanel] = 0.0;
          F1_mean [iPanel] = 0.0;
          F2_mean [iPanel] = 0.0;
          PointID [iPanel] = 0;
       }

      /* Hanning windowing scales and weight */
       Hanning_W = new su2double [nSample];
       idx_window_l = floor(nSample/8);
       idx_window_r =  nSample - floor(nSample/8);
     //  idx_window_l=0; idx_window_r=nSample+1;
       Hanning_Scale = 0.0;
       for(iSample = 0;  iSample< nSample; iSample ++)
       {
          if (iSample<idx_window_l || iSample>=idx_window_r){
              Hanning_W [iSample]= 0.5*(1-cos(8*pi*(iSample)/(nSample-1)));
            }
          else{
              Hanning_W [iSample] = 1.0;
            }
          Hanning_Scale = Hanning_Scale+Hanning_W [iSample]*Hanning_W [iSample];
        }
       Hanning_Scale = sqrt(nSample/Hanning_Scale);

}






FWHSolver::~FWHSolver(void) {
  unsigned long iSample, iObserver, iPanel;

  if (Q_mean != NULL)         delete [] Q_mean;
  if (F1_mean != NULL)         delete [] F1_mean;
  if (F2_mean != NULL)         delete [] F2_mean;
  if (Hanning_W != NULL)         delete [] Hanning_W;

  if (G != NULL){
  for(iObserver = 0; iObserver < nObserver; iObserver++)
  {
      for(iPanel  = 0;  iPanel< nPanel; iPanel++)
      {
        if (G[iObserver][iPanel] != NULL)  delete[] G[iObserver][iPanel];
      }
     if ( G[iObserver] != NULL ) delete[] G[iObserver];

  }
  delete[] G;
  }

  if (dGdy1 != NULL){
  for(iObserver = 0; iObserver < nObserver; iObserver++)
  {
      for(iPanel  = 0;  iPanel< nPanel; iPanel++)
      {
        if (dGdy1[iObserver][iPanel] != NULL)  delete[] dGdy1[iObserver][iPanel];
      }
     if ( dGdy1[iObserver] != NULL ) delete[] dGdy1[iObserver];

  }
  delete[] dGdy1;
  }

  if (dGdy2 != NULL){
  for(iObserver = 0; iObserver < nObserver; iObserver++)
  {
      for(iPanel  = 0;  iPanel< nPanel; iPanel++)
      {
        if (dGdy2[iObserver][iPanel] != NULL)  delete[] dGdy2[iObserver][iPanel];
      }
     if ( dGdy2[iObserver] != NULL ) delete[] dGdy2[iObserver];

  }
  delete[] dGdy2;
  }

  if (Q != NULL) {
    for (iPanel = 0; iPanel < nPanel; iPanel ++)
    if (Q[iPanel] != NULL)  delete [] Q[iPanel];
    delete [] Q;
  }
  if (F1 != NULL) {
    for (iPanel = 0; iPanel < nPanel; iPanel ++)
     if (F1[iPanel] != NULL) delete [] F1[iPanel];
    delete [] F1;
  }
  if (F2 != NULL) {
    for (iPanel = 0; iPanel < nPanel; iPanel ++)
     if (F2[iPanel] != NULL) delete [] F2[iPanel];
    delete [] F2;
  }
  if (fQ != NULL) {
    for (iPanel = 0; iPanel < nPanel; iPanel ++)
    if (fQ[iPanel] != NULL)  delete [] fQ[iPanel];
    delete [] fQ;
  }
  if (fF1 != NULL) {
    for (iPanel = 0; iPanel < nPanel; iPanel ++)
     if (fF1[iPanel] != NULL) delete [] fF1[iPanel];
    delete [] fF1;
  }
  if (fF2 != NULL) {
    for (iPanel = 0; iPanel < nPanel; iPanel ++)
     if (fF2[iPanel] != NULL) delete [] fF2[iPanel];
    delete [] fF2;
  }
  if (surface_geo != NULL) {
    for (iPanel = 0; iPanel < nPanel; iPanel ++)
     if (surface_geo[iPanel] != NULL) delete []surface_geo[iPanel];
    delete [] surface_geo;
  }
  if (pp != NULL) {
    for (iObserver = 0; iObserver < nPanel; iObserver++)
     if (pp[iObserver] != NULL) delete []pp[iObserver];
    delete [] pp;
  }
  if (fpp != NULL) {
    for (iObserver = 0; iObserver < nPanel; iObserver++)
     if (fpp[iObserver] != NULL) delete []fpp[iObserver];
    delete [] fpp;
  }
  if (fpp_r != NULL) {
    for (iObserver = 0; iObserver < nPanel; iObserver++)
     if (fpp_r[iObserver] != NULL) delete []fpp_r[iObserver];
    delete [] fpp_r;
  }
  if (fpp_i != NULL) {
    for (iObserver = 0; iObserver < nPanel; iObserver++)
     if (fpp_i[iObserver] != NULL) delete []fpp_i[iObserver];
    delete [] fpp_i;
  }
  if (fpp_r_root != NULL) {
    for (iObserver = 0; iObserver < nPanel; iObserver++)
     if (fpp_r_root[iObserver] != NULL) delete []fpp_r_root[iObserver];
    delete [] fpp_r_root;
  }
  if (fpp_i_root != NULL) {
    for (iObserver = 0; iObserver < nPanel; iObserver++)
     if (fpp_i_root[iObserver] != NULL) delete []fpp_i_root[iObserver];
    delete [] fpp_i_root;
  }
  if (Observer_Locations != NULL) {
    for (iObserver = 0; iObserver < nPanel; iObserver++)
     if (Observer_Locations[iObserver] != NULL) delete []Observer_Locations[iObserver];
    delete []Observer_Locations ;
  }
}


void FWHSolver::SetAeroacoustic_Analysis(CSolver *solver, CConfig *config,CGeometry *geometry,
                              ofstream &CFD_pressure_file){

    unsigned long  nMarker, iObserver;
    unsigned short iDim;
    su2double Physical_dt,Physical_t;
    su2double  CFD_PressureFluctuation=0.0;
     unsigned long ExtIter = config->GetExtIter();
     unsigned long start_iter  =  config->GetUnst_RestartIter();
     unsigned long iSample = (ExtIter-start_iter)/SamplingFreq;
     Physical_dt = config->GetDelta_UnstTime();
     Physical_t  = ExtIter*Physical_dt;
     nMarker      = config->GetnMarker_All();


//      Extract_NoiseSources(solver, config, geometry );
      Extract_NoiseSources(solver, config, geometry );
  //   SetCAA_PressureFluctuation(solver, config, geometry,SRC_p_file, SRC_ux_file, SRC_uy_file, SRC_rho_file,time_file);

      //Extract static pressure data from observer locations (for validation only)//
#ifdef HAVE_MPI
  int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
      if (rank == MASTER_NODE){
          CFD_pressure_file  <<  ExtIter  << " ";
          CFD_pressure_file << std::setprecision(15) << Physical_t  << " ";
     //     t_data  <<  ExtIter  << " ";
     //     t_data << std::setprecision(15) << Physical_t  << endl;
     //     t_data << endl;

      }


     for (iObserver=0; iObserver<nObserver; iObserver++){


      SetCFD_PressureFluctuation(solver, config, geometry, iObserver);
      CFD_PressureFluctuation=GetCFD_PressureFluctuation();


      if (rank == MASTER_NODE){

             CFD_pressure_file << std::setprecision(15) << CFD_PressureFluctuation << "\t";
            // cout<<iSample<<",  "<<iObserver<<endl;
             pp_CFD[iObserver][iSample] = CFD_PressureFluctuation;
             pp_CFD_mean[iObserver] = (pp_CFD_mean[iObserver]*iSample+ pp_CFD[iObserver][iSample])/(iSample+1);

      }
   }
      if (rank == MASTER_NODE) CFD_pressure_file << "\n";

      if (rank == MASTER_NODE){

             if (iSample==0){
                 ofstream closest_obs_coords ;
                 closest_obs_coords.open("Closest_Coords.dat");
                 closest_obs_coords <<nObserver<<endl;
                 for(iObserver = 0;  iObserver <nObserver  ;  iObserver++)
                 {
                   for (iDim=0; iDim < nDim; iDim++){
                        closest_obs_coords << std::setprecision(15) << closet_coord_AllObs[iObserver][iDim] <<"\t";
                     }
                    closest_obs_coords<< endl;
                 }
               }



      }






}




void FWHSolver::SetCFD_PressureFluctuation(CSolver *solver, CConfig *config, CGeometry *geometry, unsigned long iObserver){
  unsigned long iPoint, ClosestIndex,iSample;
  unsigned short iDim;
  su2double  Local_AvgPressureFluctuation;
  double this_distance, shortest_distance;
  su2double *Coord, *Observer_Location, *closest_coord;
  unsigned long ExtIter = config->GetExtIter();
  Observer_Location  = new su2double [3];
  closest_coord  = new su2double [nDim];
  closest_coord[0]=0.0; closest_coord[1]=0.0;
  if (nDim==3) closest_coord[2]=0.0;
  double ain, aout;
    int  ind;
   ind = 0;
    struct {
        double val;
        int rank;
    } in, out;
   int myrank;

  Observer_Location[0] = Observer_Locations[iObserver][0];
  Observer_Location[1] = Observer_Locations[iObserver][1];
  Observer_Location[2] = 0.0;
  if (nDim==3)  Observer_Location[2] = Observer_Locations[iObserver][2];
  shortest_distance = 10000000;
  CFD_PressureFluctuation= 0.0;
  Local_AvgPressureFluctuation= 0.0;

  iSample = (ExtIter - config->GetUnst_RestartIter())/SamplingFreq;

  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
          if( geometry->node[iPoint]->GetDomain()){
      Coord = geometry->node[iPoint]->GetCoord();
//     this_distance = (SU2_TYPE::GetValue(Coord[0])-SU2_TYPE::GetValue(Observer_Location[0]))*(SU2_TYPE::GetValue(Coord[0])-SU2_TYPE::GetValue(Observer_Location[0]));
//     this_distance = this_distance+(SU2_TYPE::GetValue(Coord[1])-SU2_TYPE::GetValue(Observer_Location[1]))*(SU2_TYPE::GetValue(Coord[1])-SU2_TYPE::GetValue(Observer_Location[1]));
       this_distance =0.0;
     for ( iDim = 0; iDim < nDim; iDim++){
         this_distance += (SU2_TYPE::GetValue(Coord[iDim])-SU2_TYPE::GetValue(Observer_Location[iDim]))*(SU2_TYPE::GetValue(Coord[iDim])-SU2_TYPE::GetValue(Observer_Location[iDim]));
       }

     this_distance = sqrt(this_distance);
     if (this_distance<shortest_distance){
         shortest_distance = this_distance;
        ClosestIndex = iPoint;
     }
    }
   }
       ain=shortest_distance;
     MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
        in.val = ain;
        in.rank = myrank;
  MPI_Reduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MINLOC, MASTER_NODE, MPI_COMM_WORLD );
   if (myrank == MASTER_NODE) {
            aout = out.val;
            ind = out.rank;
    }
  MPI_Bcast( &ind, 1,MPI_INT , MASTER_NODE, MPI_COMM_WORLD );
     if (myrank==ind && !UseAnalytic){
         if (config->GetKind_Solver()==NAVIER_STOKES){
//         Local_AvgPressureFluctuation = solver->node[ClosestIndex]->GetSolution(6);
           Local_AvgPressureFluctuation = solver->node[ClosestIndex]->GetSolution(2*nDim+2);
         }
        if (config->GetKind_Solver()==RANS){
//         Local_AvgPressureFluctuation = solver->node[ClosestIndex]->GetSolution(7);
           Local_AvgPressureFluctuation = solver->node[ClosestIndex]->GetSolution(2*nDim+3);
        }
     if (nDim==2){
        cout<<std::setprecision(9)<<"Obs_x= "<<  SU2_TYPE::GetValue(Observer_Location[0]) << ", Obs_y= "<< SU2_TYPE::GetValue(Observer_Location[1]) <<", x= "<<solver->node[ClosestIndex]->GetSolution(0)<<", y="<<solver->node[ClosestIndex]->GetSolution(1)<<", p="<<Local_AvgPressureFluctuation<<endl;
     }else{
        cout<<std::setprecision(9)<<"Obs_x= "<<  SU2_TYPE::GetValue(Observer_Location[0]) << ", Obs_y= "<< SU2_TYPE::GetValue(Observer_Location[1]) << ", Obs_z= "<< SU2_TYPE::GetValue(Observer_Location[2])<<", x= "<<solver->node[ClosestIndex]->GetSolution(0)<<", y="<<solver->node[ClosestIndex]->GetSolution(1)<<", z="<<solver->node[ClosestIndex]->GetSolution(2)<<", p="<<Local_AvgPressureFluctuation<<endl;
       }


    if (iSample == 0){
     closest_coord[0]= solver->node[ClosestIndex]->GetSolution(0);
     closest_coord[1]= solver->node[ClosestIndex]->GetSolution(1);
     if (nDim==3) closest_coord[2]= solver->node[ClosestIndex]->GetSolution(2);

       }
        //        cout<<"Obs_x= "<<  SU2_TYPE::GetValue(Observer_Location[0]) << ", Obs_y= "<< SU2_TYPE::GetValue(Observer_Location[1]) <<", x= "<<solver->node[ClosestIndex]->GetSolution(0)<<", y="<<solver->node[ClosestIndex]->GetSolution(1)<<", p="<<Local_AvgPressureFluctuation<<endl;

     }
#ifdef HAVE_MPI
  SU2_MPI::Reduce(&Local_AvgPressureFluctuation,&CFD_PressureFluctuation,1,MPI_DOUBLE,MPI_SUM,MASTER_NODE,MPI_COMM_WORLD);
  SU2_MPI::Reduce(closest_coord,closet_coord_AllObs[iObserver],nDim,MPI_DOUBLE,MPI_SUM,MASTER_NODE,MPI_COMM_WORLD);
#endif
//cout<<closest_coord[0]<<" "<<closest_coord[1]<<" "<<closest_coord[2]<<" "<<myrank<<endl;
  // if (myrank==MASTER_NODE) cout<<endl<<closet_coord_AllObs[iObserver][0]<<" "<<closet_coord_AllObs[iObserver][1]<<" "<<closet_coord_AllObs[iObserver][2]<<endl;

}


void FWHSolver::Extract_NoiseSources(CSolver *solver, CConfig *config, CGeometry *geometry){
  unsigned short iDim ;
  unsigned long iObserver, iPoint ,iMarker,iVertex,  nMarker,iSample, iPanel;
  su2double *Coord,   *Normal;
  su2double  x, y, z, nx, ny, nz, dS, Area;
  su2double rho, rho_ux, rho_uy, rho_uz, rho_E, TKE, ux, uy, uz, StaticEnergy, p;
  su2double CheckNormal=0.0;
  su2double Area_Factor=1.0;
  su2double F1, F2, F3;
  su2double x1,x2,x3,y1,y2,y3,r1,r2,r3,r_mag;
  su2double theta, r_panel, u_mag;


  nMarker      = config->GetnMarker_All();
  unsigned long ExtIter = config->GetExtIter();
  iSample = (ExtIter - config->GetUnst_RestartIter())/SamplingFreq;

#ifdef HAVE_MPI
  int rank, nProcessor;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
     MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
#endif
  iPanel = 0;
  for (iMarker = 0; iMarker < nMarker; iMarker++){

    /* --- Loop over boundary markers to select those on the FWH surface --- */
      if (config->GetMarker_All_KindBC(iMarker) == INTERNAL_BOUNDARY) {
//      if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY) {
  //      cout<<"INTERNAL BOUNDARY DETECTED!!! RANK="<<rank <<endl;
         /* Determine the area factor for outflow boundary averaging */
        size_t last_index = config->GetMarker_All_TagBound(iMarker).find_last_not_of("0123456789");
        string result = config->GetMarker_All_TagBound(iMarker).substr(last_index + 1);
        int iFWH = std::strtol(result.c_str(),NULL, 10);
        //For the case of an open FWH surface
        if (totFWH==1){
            Area_Factor = 1.0;
          }
        //Closed FWH surface with at least one outflow boundary
        else if (totFWH>1){
            if (iFWH>0 && iFWH<totFWH){
                Area_Factor = 1.0*(totFWH-iFWH)/(totFWH-1.0);
              }
            else if (iFWH == totFWH){
                Area_Factor = 1.0/(totFWH-1.0);
              }
            else{
                cout<<"WARNING!!! iFWH > totFWH or iFWH < 0 !!!! Check FWH Markers!!!"<<endl;
              }


          }
        else{
            cout<<"WARNING!!! totFWH < 1 !!!! Check FWH Markers!!!"<<endl;
          }


        for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++){
            iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
              if( geometry->node[iPoint]->GetDomain()){

                      Coord = geometry->node[iPoint]->GetCoord();
                      Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
                      Area  = 0.0; for ( iDim = 0; iDim < nDim; iDim++)   Area += Normal[iDim]*Normal[iDim];  Area  = sqrt( Area );

                      x  = SU2_TYPE::GetValue(Coord[0]);                                                                            // x
                      y  = SU2_TYPE::GetValue(Coord[1]);                                                                            // y
                      z  = 0.0;
                      if (nDim==3) z = SU2_TYPE::GetValue(Coord[2]);
                      nx = Normal[0]/Area  ;                                                                                        // n_x
                      ny = Normal[1]/Area  ;                                                                                        // n_y
                      nz = 0.0;
                      if (nDim==3)  nz = Normal[2]/Area  ;
                      dS = Area ;                                                                                                   // dS

                    //only extract flow data for points on the boundary that has surface normal pointing AWAY from the body.
                     CheckNormal =  x*nx+y*ny+z*nz;
                     if(CheckNormal<0) {
                  //    cout<<"Inward Pointing Normal Detected!!! Flipping Normals"<<", x="<<x<<", y="<<y<<", nx="<<nx<<", ny="<<ny<<endl;
                      nx=-nx; ny=-ny; nz=-nz;
                      }

                //  if (CheckNormal>0){

                      /*write out geometric info of the permeable surface only once*/
                      if (iSample==0){
                            surface_geo[iPanel][0] = x;
                            surface_geo[iPanel][1] = y;
                            surface_geo[iPanel][2] = nx;
                            surface_geo[iPanel][3] = ny;
                            surface_geo[iPanel][4] = Area_Factor*dS;
                            if (nDim==3){
                                surface_geo[iPanel][5] = z;
                                surface_geo[iPanel][6] = nz;
                                if (TimeDomain3D){
                                    for(iObserver = 0;  iObserver <nObserver  ;  iObserver++)
                                    {
                                        r1 = Observer_Locations[iObserver][0]-x;
                                        r2 = Observer_Locations[iObserver][1]-y;
                                        r3 = Observer_Locations[iObserver][2]-z;
                                        r_mag = sqrt(r1*r1+r2*r2+r3*r3);
                                        r1 = r1/r_mag; r2 = r2/r_mag;r3 = r3/r_mag;
                                        surface_geo[iPanel][7+iObserver*4] = r1;
                                        surface_geo[iPanel][8+iObserver*4] = r2;
                                        surface_geo[iPanel][9+iObserver*4] = r3;
                                        surface_geo[iPanel][10+iObserver*4] = r_mag;
                                    }
                                  }


                              }
                            PointID[iPanel] = geometry->node[iPoint]->GetGlobalIndex();
                         //   if (rank == MASTER_NODE) cout<<PointID[iPanel]<<",  "<<x<<", "<<y <<",  "<<sqrt(x*x+y*y)<<endl;

             //          cout << std::setprecision(15) << x<<' '<<y<<' '<<z<<' '<<nx<<' '<<ny<<' '<<nz<<' '<< Observer_Locations[0][0]-x<<' '<<Observer_Locations[0][1]-y<<' '<< Observer_Locations[0][2]-z<<' '<<dS<<endl;




                      }

                      if(!UseAnalytic){
                      //extract CONSERVATIVE flow data from a particular panel on the FWH surface
                      rho = solver->node[iPoint]->GetSolution(nDim);
                      rho_ux = solver->node[iPoint]->GetSolution(nDim+1);
                      rho_uy = solver->node[iPoint]->GetSolution(nDim+2);


                      if (nDim==3)  rho_uz = solver->node[iPoint]->GetSolution(nDim+3);
                      rho_E = solver->node[iPoint]->GetSolution(2*nDim+1);
                      TKE = 0.0;
                        }

                      //Register CONSERVATIVE variables as input for adjoint computation
                      if (config->GetAD_Mode()){
                      AD::RegisterInput(rho );
                      AD::RegisterInput(rho_ux );
                      AD::RegisterInput(rho_uy );
                      if (nDim==3) AD::RegisterInput(rho_uz );
                      AD::RegisterInput(rho_E );
                      AD::RegisterInput(TKE);
                        }



                      if (!UseAnalytic){
                      //compute primitive variables from conservative variables
                      ux = rho_ux/rho;
                      uy = rho_uy/rho;
                      uz = 0.0;
                      if (nDim==3) uz= rho_uz/rho;
                      StaticEnergy =  rho_E/rho-0.5*(ux*ux+uy*uy+uz*uz)-TKE;
                      p = (config->GetGamma()-1)*rho*StaticEnergy;
                        }else{
                          //use analytic solution on the porous surface
                          r_panel = sqrt(x*x+y*y+z*z);
                          theta = (iSample*SamplingFreq*config->GetDelta_UnstTime()+r_panel/a_inf)*Freq_a;
                          p=Amp_a/4/M_PI/r_panel*cos(theta);
                          u_mag = Amp_a/4/M_PI/r_panel/FreeStreamDensity*(sin(theta)/r_panel/Freq_a+cos(theta)/a_inf);
                          ux = u_mag*nx;
                          uy = u_mag*ny;
                          uz = u_mag*nz;
                          rho = p/a_inf/a_inf+FreeStreamDensity;
                       //   if (iPanel==0) cout<<p<<' '<<u_mag<<' '<<r_panel<<endl;

                        }
//                      rho_data << std::setprecision(15) << rho << ' ';
//                      p_data << std::setprecision(15) << p << ' ';
//                      ux_data << std::setprecision(15) << ux << ' ';
//                      uy_data << std::setprecision(15) << uy << ' ';
//                      uz_data << std::setprecision(15) << uz << ' ';


                      //compute monopole and dipole source terms on the FWH surface
                      if (TimeDomain3D){
                      Q[iPanel][iSample] = rho*(ux*nx+uy*ny+uz*nz)/FreeStreamDensity;
                      F1= rho*(ux*nx+uy*ny+uz*nz)*(ux)+(p-FreeStreamPressure)*nx;
                      F2= rho*(ux*nx+uy*ny+uz*nz)*(uy)+(p-FreeStreamPressure)*ny;
                      F3= rho*(ux*nx+uy*ny+uz*nz)*(uz)+(p-FreeStreamPressure)*nz;
                      for (iObserver=0; iObserver<nObserver; iObserver++){
                          Fr [iObserver][iPanel][iSample]= F1*surface_geo[iPanel][7+iObserver*4]+F2*surface_geo[iPanel][8+iObserver*4] +F3*surface_geo[iPanel][9+iObserver*4] ;
                        }
                        }
                      else{
                      Q[iPanel][iSample] = rho*(ux*nx+uy*ny+uz*nz);
                      F[0][iPanel][iSample]  = rho*(ux*nx+uy*ny+uz*nz)*(ux-2*U1)+p*nx;
                      F[1][iPanel][iSample]  = rho*(ux*nx+uy*ny+uz*nz)*(uy-2*U2)+p*ny;
                      if (nDim==3) F[2][iPanel][iSample]  = rho*(ux*nx+uy*ny+uz*nz)*(uz-2*U3)+p*nz;
                        }

                      //Compute mean source terms
                      Q_mean[iPanel] = (Q_mean[iPanel]*iSample+ Q[iPanel][iSample])/(iSample+1);
                      if (TimeDomain3D){
                       for (iObserver=0; iObserver<nObserver; iObserver++){
                      Fr_mean[iObserver][iPanel] = (Fr_mean[iObserver][iPanel]*iSample+ Fr[iObserver][iPanel][iSample])/(iSample+1);
                         }
                        }else{
                      for (iDim=0; iDim<nDim; iDim++){
                      F_mean[iDim][iPanel] = (F_mean[iDim][iPanel]*iSample+ F[iDim][iPanel][iSample])/(iSample+1);
                        }
                        }
                      iPanel++;

                }
          }

    }
   }

//  cout<<"loc-1, after extract, Rank= "<<rank<<endl;
//  p_data <<endl;
//  ux_data <<endl;
//  uy_data <<endl;
//  uz_data <<endl;
//  rho_data<<endl;
}




void FWHSolver::Window_SourceTerms(){
  unsigned long iPanel, iSample,iDim, iObserver;

   if (TimeDomain3D){
       for (iPanel=0; iPanel<nPanel; iPanel++){
            for (iSample=0; iSample<nSample; iSample++){
               Q[iPanel][iSample]=  Q[iPanel][iSample] - Q_mean[iPanel];
          //     cout<< Q_mean[iPanel]<< endl;
                for (iObserver=0; iObserver<nObserver; iObserver++){
               Fr[iObserver][iPanel][iSample] =    Fr[iObserver][iPanel][iSample] - Fr_mean[iObserver][iPanel] ;
             //  cout<<Fr_mean[iObserver][iPanel] <<endl;
                  }
              }
         }
     }else{
     for (iPanel=0; iPanel<nPanel; iPanel++){
          for (iSample=0; iSample<nSample; iSample++){
              fQ[iPanel][iSample]=  (Q[iPanel][iSample] - Q_mean[iPanel])*Hanning_W[iSample]*Hanning_Scale ;
              for (iDim=0; iDim<nDim; iDim++){
              fF[iDim][iPanel][iSample] =    (F[iDim][iPanel][iSample] - F_mean[iDim][iPanel])*Hanning_W[iSample]*Hanning_Scale ;
                }
            }
       }
     }


//   ofstream Fr_file ;
//   ofstream Un_file ;
//   Fr_file.open("Fr2"); Un_file.open("Un2");
//   for (iPanel=0; iPanel<nPanel; iPanel++){
//       for (iSample=0; iSample<nSample; iSample++){
//           Fr_file  << std::setprecision(15) << Fr[0][iPanel][iSample]<< ' ';
//           Un_file  << std::setprecision(15) << Q[iPanel][iSample]<< ' ';
//         }
//       Fr_file  << endl;Un_file  << endl;
//     }


}


void FWHSolver::Compute_TimeDomainPanelSignal(CConfig *config){
su2double Un_dot, Fr_dot, dt, dS;
unsigned long iPanel, iSample, iObserver;

  dt = config->GetDelta_UnstTime()*SamplingFreq;
  for (iObserver=0; iObserver<nObserver; iObserver++){
      for (iPanel=0; iPanel<nPanel; iPanel++){
          dS = surface_geo[iPanel][4];
          for (iSample=0; iSample<nSample; iSample++){
          if (iSample==0){
              Un_dot = (-Q[iPanel][iSample+2]+4.0*Q[iPanel][iSample+1]-3.0*Q[iPanel][iSample])/2.0/dt;
              Fr_dot = (-Fr[iObserver][iPanel][iSample+2]+4.0*Fr[iObserver][iPanel][iSample+1]-3.0*Fr[iObserver][iPanel][iSample])/2.0/dt;

            }else if(iSample==nSample-1){
              Un_dot = (3.0*Q[iPanel][iSample]-4.0*Q[iPanel][iSample-1]+Q[iPanel][iSample-2])/2.0/dt;
              Fr_dot = (3.0*Fr[iObserver][iPanel][iSample]-4.0*Fr[iObserver][iPanel][iSample-1]+Fr[iObserver][iPanel][iSample-2])/2.0/dt;
            }else{
              Un_dot = (Q[iPanel][iSample+1]-Q[iPanel][iSample-1])/2.0/dt;
              Fr_dot = (Fr[iObserver][iPanel][iSample+1]-Fr[iObserver][iPanel][iSample-1])/2.0/dt;
            }
            pp_ret[iObserver][iPanel][iSample] = FreeStreamDensity*Un_dot/surface_geo[iPanel][10+iObserver*4]+Fr_dot/surface_geo[iPanel][10+iObserver*4]/a_inf+ Fr[iObserver][iPanel][iSample]/surface_geo[iPanel][10+iObserver*4]/surface_geo[iPanel][10+iObserver*4];
            pp_ret[iObserver][iPanel][iSample] = pp_ret[iObserver][iPanel][iSample]*dS/4/M_PI;
            }

        }
    }
//  ofstream pp_ret_file ;
//  pp_ret_file.open("pp_ret");
//  for (iPanel=0; iPanel<nPanel; iPanel++){
//      for (iSample=0; iSample<nSample; iSample++){
//          pp_ret_file  << std::setprecision(15) << pp_ret[0][iPanel][iSample]<< ' ';
//        }
//      pp_ret_file  << endl;
//    }
}


void FWHSolver::Compute_ObserverTime(CConfig *config){
  su2double r, t_src,t_interp_start,t_interp_end, dt_interp, dt;
unsigned long iPanel, iSample, iObserver;
unsigned long start_iter  =  config->GetUnst_RestartIter();
su2double r_min=10.0e31;su2double r_max=0.0;


  dt = config->GetDelta_UnstTime()*SamplingFreq;
  for (iObserver=0; iObserver<nObserver; iObserver++){
      r_min=10.0e31; r_max=0.0;
      for (iPanel=0; iPanel<nPanel; iPanel++){
          r = surface_geo[iPanel][10+iObserver*4];
          if (r>r_max){
              r_max = r;
            }
          if (r<r_min){
              r_min = r;
            }

          for (iSample=0; iSample<nSample; iSample++){
              t_src = config->GetDelta_UnstTime()*(start_iter+iSample*SamplingFreq);
              if(UseAnalytic) t_src = t_src + sqrt(surface_geo[iPanel][0]*surface_geo[iPanel][0]+surface_geo[iPanel][1]*surface_geo[iPanel][1]+surface_geo[iPanel][5]*surface_geo[iPanel][5])/a_inf;
              t_Obs[iObserver][iPanel][iSample]=t_src + r/a_inf;
            }

        }

      //send rmax and rmin to root
    SU2_MPI::Allreduce(&r_min, &r_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
     SU2_MPI::Allreduce(&r_max, &r_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
     r_minmax[iObserver][0]= r_min; r_minmax[iObserver][1]= r_max;


    }


  //cout<<"Time Shift INFO: "<< r_max<<", "<< r_min<<", "<< ;
  for (iObserver=0; iObserver<nObserver; iObserver++){
  t_interp_start = config->GetDelta_UnstTime()*(start_iter)+r_minmax[iObserver][1]/a_inf;
  t_interp_end = config->GetDelta_UnstTime()*(start_iter+nSample*SamplingFreq-1)+r_minmax[iObserver][0]/a_inf;
  dt_interp = (t_interp_end - t_interp_start)/(nSample-1);
  //cout<<t_interp_start<<", "<<t_interp_end<<", "<< dt_interp;
  for (iSample=0; iSample<nSample; iSample++){
      t_interp[iObserver][iSample] = t_interp_start + dt_interp*iSample;
    }
    }
//  ofstream t_obs_file ;
//  t_obs_file.open("t_obs2");
//  for (iPanel=0; iPanel<nPanel; iPanel++){
//      for (iSample=0; iSample<nSample; iSample++){
//           t_obs_file<< std::setprecision(15) << t_Obs[0][iPanel][iSample]<< ' ';
//        }
//       t_obs_file << endl;
//    }
//  ofstream t_interp_file ;
//  t_interp_file.open("t_interp2");
//  for (iSample=0; iSample<nSample; iSample++){
//       t_interp_file<< std::setprecision(15) << t_interp[iSample]<< ' ';
//    }
}

void FWHSolver::Interpolate_PressureSignal(CGeometry *geometry){
  vector <su2double> x(nSample), t(nSample), derivative(nSample);
  unsigned long iPanel, iSample, iObserver;
//  su2double sinx_e=0.0, vec_x_e=0.0;

  su2double yp1=10.0e31;su2double ypn=10.0e31;
  for (iObserver=0; iObserver<nObserver; iObserver++){
      for (iPanel=0; iPanel<nPanel; iPanel++){
          for (iSample=0; iSample<nSample; iSample++){
               x[iSample]= pp_ret[iObserver][iPanel][iSample];
               t[iSample]= t_Obs[iObserver][iPanel][iSample];
            }
          geometry->SetSpline(t, x, nSample, yp1, ypn, derivative);
          for (iSample=0; iSample<nSample; iSample++){
            pp_interp[iObserver][iPanel][iSample] =geometry->GetSpline(t, x, derivative, nSample, t_interp[iObserver][iSample]);
            }
        }
    }

//  ofstream pp_interp_file ;
//  pp_interp_file.open("pp_interp2");
//  for (iPanel=0; iPanel<nPanel; iPanel++){
//      for (iSample=0; iSample<nSample; iSample++){
//           pp_interp_file << std::setprecision(15) << pp_interp[0][iPanel][iSample]<< ' ';
//        }
//       pp_interp_file << endl;
//    }
}


void FWHSolver::FFT_SourceTermsR2(){
  unsigned long  iPanel, iSample;

 FFT* FFT_container = new FFT() ;



       //perform FFT on each row (panel) of Q, F1 and F2
       for (iPanel=0; iPanel<nPanel; iPanel++){

           // Writing the complex array data
            CArray dataQ(fQ[iPanel] ,nSample);
//            CArray dataF1(fF1[iPanel] ,nSample);
//            CArray dataF2(fF2[iPanel] ,nSample);
            CArray dataF1(fF[0][iPanel] ,nSample);
            CArray dataF2(fF[1][iPanel] ,nSample);
            CArray dataF3(fF[1][iPanel] ,nSample);
            if (nDim==3) CArray dataF3(fF[2][iPanel] ,nSample);

            // Call the FFT function for performing the forward fft
            FFT_container->fft_r2(dataQ);
            FFT_container->fft_r2(dataF1);
            FFT_container->fft_r2(dataF2);
            if (nDim==3)  FFT_container->fft_r2(dataF3);

            for (iSample=0; iSample<nSample; iSample++){
                fQ[iPanel][iSample] = dataQ[iSample];
//                fF1[iPanel][iSample] = dataF1[iSample];
//                fF2[iPanel][iSample] = dataF2[iSample];
                fF[0][iPanel][iSample] = dataF1[iSample];
                fF[1][iPanel][iSample] = dataF2[iSample];
                if (nDim==3)  fF[2][iPanel][iSample] = dataF3[iSample];
//              cout<<fQ[iPanel][iSample]<<endl;
              }
         }

}

//void FWHSolver::RegisterVariables(){
//unsigned long iPanel, iSample;
//  for (int i =0; i< nDim+3; i++){
//     for (iPanel=0; iPanel<nPanel; iPanel++){
//       for (iSample=0; iSample<nSample; iSample++){
//               AD::RegisterInput(U[i][iPanel][iSample]);
//         }w
//    }
//  }
//}



void FWHSolver::Compute_FarfieldNoise(CSolver *solver, CConfig *config, CGeometry *geometry){

#ifdef HAVE_MPI
 int rank, nProcessor;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
  // cout<<"START OF Compute_FarfieldNoise, Rank= "<<rank<<endl;
#endif

   //  cout<<"loc0a, Rank= "<<rank;
  if (!TimeDomain3D){
  if (nDim==2)
  Compute_GreensFunction2D(config);
  if (nDim==3)
  Compute_GreensFunction3D(config);
    }


 // cout<<"loc1, Rank= "<<rank;
  Window_SourceTerms();

 // cout<<"loc2, Rank= "<<rank;

  if (!TimeDomain3D){
  FFT_SourceTermsR2();
    }

  if (TimeDomain3D){
  Compute_TimeDomainPanelSignal(config);

//  cout<<"loc3, Rank= "<<rank;
  Compute_ObserverTime(config);

 // cout<<"loc4, Rank= "<<rank;
  Interpolate_PressureSignal(geometry);
    }


 // cout<<"loc5, Rank= "<<rank;

  Integrate_Sources(config);


//  ofstream pp_FWH_file ;
//  pp_FWH_file.open("pp_FWH2");
//  for (int iObserver=0; iObserver<nObserver; iObserver++){
//      for (int iSample=0; iSample<nSample; iSample++){
//           pp_FWH_file << std::setprecision(15) << pp_TimeDomain[iObserver][iSample]<< ' ';
//        }
//      pp_FWH_file  << endl;
//    }


//cout<<"before MPI Reduce, Rank= "<<rank;

for (int iObserver = 0; iObserver<nObserver; iObserver++){
#ifdef HAVE_MPI
    if(TimeDomain3D){
       SU2_MPI::Reduce( pp_TimeDomain[iObserver],pp_TimeDomain_root[iObserver], nSample,MPI_DOUBLE,MPI_SUM,MASTER_NODE,MPI_COMM_WORLD);
   //     SU2_MPI::Reduce( pp_TimeDomain[iObserver],pp[iObserver], nSample,MPI_DOUBLE,MPI_SUM,MASTER_NODE,MPI_COMM_WORLD);
      }else{
    SU2_MPI::Reduce( fpp_r[iObserver],fpp_r_root[iObserver], nSample,MPI_DOUBLE,MPI_SUM,MASTER_NODE,MPI_COMM_WORLD);
    SU2_MPI::Reduce( fpp_i[iObserver], fpp_i_root[iObserver], nSample,MPI_DOUBLE,MPI_SUM,MASTER_NODE,MPI_COMM_WORLD);
      }
#endif
}

if (rank==MASTER_NODE){



if (TimeDomain3D){

 //   for (int iObserver=0; iObserver<nObserver; iObserver++){
 //       cout << std::setprecision(15) <<"Obs No.: "<<iObserver<<", r_min= "<<r_minmax[iObserver][0]<<", r_max= "<< r_minmax[iObserver][1] <<endl;
 //    }

Compute_SPL();
    ofstream pp_FWH_file ;
    pp_FWH_file.open("pp_FWH");

for (int iSample=0; iSample<nSample; iSample++){
   //pp_FWH_file << std::setprecision(15)<< t_interp[iSample]<<' ';

     for (int iObserver=0; iObserver<nObserver; iObserver++){
        pp_FWH_file << std::setprecision(15) <<t_interp[iObserver][iSample] <<' '<<pp_TimeDomain_root[iObserver][iSample]<< ' ';
    //     pp_FWH_file << std::setprecision(15) <<t_interp[iObserver][iSample] <<' '<<pp[iObserver][iSample]<< ' ';
      }
    pp_FWH_file  << endl;
  }

  }


if (!TimeDomain3D){

iFFT_SignalR2();

Compute_SPL();

ofstream ppa_file ;
ppa_file.open("pp_FWH");
ofstream ppb_file ;
ppb_file.open("ppbSU2");
//ofstream pp_CFD_file ;
//pp_CFD_file.open("pp_CFD");

unsigned long start_iter  =  config->GetUnst_RestartIter();
for (int i=0;i<nSample; i++){
     ppa_file <<start_iter+i<<", ";
     ppa_file<< std::setprecision(15) <<  SPL<< ", ";
 //   pp_CFD_file<<start_iter+i<<", "<< config->GetDelta_UnstTime()*(start_iter+i)<<", ";
    for (int j=0; j<nObserver; j++){
      ppa_file << std::setprecision(15) << real(pp[j][i]);
      if (j != nObserver - 1) ppa_file << ", ";
      ppb_file << std::setprecision(15) << imag(pp[j][i]);
      if (j != nObserver - 1) ppb_file << ", ";
//      pp_CFD_file<< std::setprecision(15) << pp_CFD[j][i]-pp_CFD_mean[j];
 //     if (j != nObserver - 1) pp_CFD_file << ", ";

      }
    ppa_file << endl;
 //   pp_CFD_file <<endl;
  }
 }




cout<<endl<<std::setprecision(15)<<"****** RMS(p') averaged over "<<nObserver<<" observer locations = "<<SPL <<"  ******"<<endl;


//write out static pressure fluctuation
ofstream pp_CFD_file ;
pp_CFD_file.open("pp_CFD");
unsigned long start_iter  =  config->GetUnst_RestartIter();
for (int i=0;i<nSample; i++){
    pp_CFD_file<<start_iter+i<<", "<< config->GetDelta_UnstTime()*(start_iter+i)<<", ";
    for (int j=0; j<nObserver; j++){
      pp_CFD_file<< std::setprecision(15) << pp_CFD[j][i]-pp_CFD_mean[j];
      if (j != nObserver - 1) pp_CFD_file << ", ";
      }
    pp_CFD_file <<endl;
  }

  }
}




void FWHSolver::Compute_GreensFunction2D(CConfig *config){
  unsigned long iObserver, iPanel, iSample;
  su2double r1, r2, r_beta, x1, x2, y1, y2, w, dt, arg;
  su2double J0, J1, Y0, Y1;
  complex <su2double> H_2_0, H_2_1, H_2_0xi, Amp, Amp1, exp_arg;
  su2double pi =  acos(-1);
  complex <su2double> I = complex<su2double> (0, 1);


       dt = config->GetDelta_UnstTime()*SamplingFreq;


       for (iObserver=0; iObserver<nObserver; iObserver++){
           for (iPanel=0; iPanel<nPanel; iPanel++){
               x1 = Observer_Locations[iObserver][0];
               x2 = Observer_Locations[iObserver][1];
               y1 = surface_geo[iPanel][0];
               y2 = surface_geo[iPanel][1];
               r1 = (x1-y1)*cos(AOA*pi/180)+(x2-y2)*sin(AOA*pi/180);
               r2 = -(x1-y1)*sin(AOA*pi/180)+(x2-y2)*cos(AOA*pi/180);
               r_beta = sqrt(r1*r1+beta_sq*r2*r2);

               for (iSample=0; iSample<nSample; iSample++){
                  w = 2*pi*iSample/(nSample-1)/dt;
                  arg = w/a_inf/beta_sq*r_beta;
                  J0 = bessj0(SU2_TYPE::GetValue(arg));
                  Y0 = bessy0(SU2_TYPE::GetValue(arg));
                  J1 = bessj1(SU2_TYPE::GetValue(arg));
                  Y1 = bessy1(SU2_TYPE::GetValue(arg));
                  H_2_0 = J0 - Y0*I;
                  H_2_1 = J1 - Y1*I;
                  H_2_0xi = Y0 + J0*I;
                  exp_arg = I*M*w/a_inf*r1/beta_sq;

                  Amp1=  su2double(1.0/4.0/sqrt(beta_sq))*I;
                  Amp = Amp1*(su2double(cos(exp_arg.imag()))+su2double(sin(exp_arg.imag()))*I);

                   if(iSample != 0){
                     G[iObserver][iPanel][iSample] = Amp*H_2_0;
                     dGdy1[iObserver][iPanel][iSample]= -Amp*M*w/a_inf/beta_sq*su2double(cos(AOA*pi/180))*H_2_0xi +Amp*w/a_inf/beta_sq/r_beta*su2double((r1*cos(AOA*pi/180))-beta_sq*r2*sin(AOA*pi/180))*H_2_1;
                     dGdy2[iObserver][iPanel][iSample]= -Amp*M*w/a_inf/beta_sq*su2double(sin(AOA*pi/180))*H_2_0xi +Amp*w/a_inf/beta_sq/r_beta*su2double((r1*sin(AOA*pi/180))+beta_sq*r2*cos(AOA*pi/180))*H_2_1;
                     }
                  else{
                      G[iObserver][iPanel][iSample] = 0;
                      dGdy1[iObserver][iPanel][iSample]= 0;
                      dGdy2[iObserver][iPanel][iSample]= 0;
                    }

                 }
             }
         }

}

void FWHSolver::Compute_GreensFunction3D(CConfig *config){
  unsigned long iObserver, iPanel, iSample;
  su2double  x1, x2, x3, y1, y2, y3, w, dt, arg;
  su2double r, r_plus, r_star, M_r;
  complex <su2double>   Amp, Amp1, exp_arg;
  su2double pi =  acos(-1);
  complex <su2double> I = complex<su2double> (0, 1);
  su2double fac1_y1, fac2_y1, fac3_y1, fac1_y2, fac2_y2, fac3_y2, fac1_y3, fac2_y3, fac3_y3;


       dt = config->GetDelta_UnstTime()*SamplingFreq;


       for (iObserver=0; iObserver<nObserver; iObserver++){
           for (iPanel=0; iPanel<nPanel; iPanel++){
               x1 = Observer_Locations[iObserver][0];
               x2 = Observer_Locations[iObserver][1];
               x3 = Observer_Locations[iObserver][2];
               y1 = surface_geo[iPanel][0];
               y2 = surface_geo[iPanel][1];
               y3 = surface_geo[iPanel][5];
               r = sqrt((x1-y1)*(x1-y1)+(x2-y2)*(x2-y2)+(x3-y3)*(x3-y3));
               M_r = (U1*(x1-y1)+U2*(x2-y2)+U3*(x3-y3))/a_inf/r;
               r_star = r*sqrt(M_r*M_r+1-M*M);
               r_plus = (-r*M_r+r_star)/beta_sq;
               fac1_y1=(x1-y1)/r/r;
               fac2_y1=((x1-y1)*M_r*M_r -(U1*M_r*r)/a_inf)/r_star/r_star;
               fac3_y1=(U1/a_inf - (x1-y1)*r_star/r/r + ((x1-y1)*M_r*M_r-(U1*M_r*r)/a_inf)/r_star);
               fac1_y2=(x2-y2)/r/r;
               fac2_y2=((x2-y2)*M_r*M_r -(U2*M_r*r)/a_inf)/r_star/r_star;
               fac3_y2=(U2/a_inf - (x2-y2)*r_star/r/r + ((x2-y2)*M_r*M_r-(U2*M_r*r)/a_inf)/r_star);
               fac1_y3=(x3-y3)/r/r;
               fac2_y3=((x3-y3)*M_r*M_r -(U3*M_r*r)/a_inf)/r_star/r_star;
               fac3_y3=(U3/a_inf - (x3-y3)*r_star/r/r + ((x3-y3)*M_r*M_r-(U3*M_r*r)/a_inf)/r_star);

               for (iSample=0; iSample<nSample; iSample++){
                  w = 2*pi*iSample/(nSample-1)/dt;
                  arg = -w/a_inf*r_plus;
                  exp_arg = I*arg;
                  Amp1=  su2double(-1.0/4.0/pi/r_star);
                  Amp = Amp1*(su2double(cos(exp_arg.imag()))+su2double(sin(exp_arg.imag()))*I);
                  G[iObserver][iPanel][iSample]= Amp;
//                  dGdy1[iObserver][iPanel][iSample]= Amp*fac1_y1- Amp*fac2_y1 -w*Amp*I/a_inf/beta_sq*fac3_y1  ;
                  dGdy1[iObserver][iPanel][iSample]= Amp*fac1_y1- Amp*fac2_y1 -w*Amp*I/a_inf/beta_sq*su2double(U3/a_inf - (x3-y3)*r_star/r/r + ((x3-y3)*M_r*M_r-(U3*M_r*r)/a_inf)/r_star);
                  dGdy2[iObserver][iPanel][iSample]= Amp*fac1_y2- Amp*fac2_y2 -w*Amp*I/a_inf/beta_sq*fac3_y2  ;
                  dGdy3[iObserver][iPanel][iSample]= Amp*fac1_y3- Amp*fac2_y3 -w*Amp*I/a_inf/beta_sq*fac3_y3  ;
                 }
             }
         }



}
//dGdy1 =G*(x1-y1)/r^2 ...
//   - G*((x1-y1)*M_r^2 -(U1*M_r*r)/a_inf)/r_star^2 ...
//   - w*G*(U1/a_inf - (x1-y1)*r_star/r^2 + ((x1-y1)*M_r^2-(U1*M_r*r)/a_inf)/r_star)*1i/(a_inf*beta_sq)


su2double FWHSolver::bessj0 (su2double x){
  su2double ax,z;
  su2double xx,y,ans,ans1,ans2;

  if ((ax=fabs(x)) < 8.0) {
     y=x*x;
     ans1=57568490574.0+y*(-13362590354.0+y*(651619640.7
        +y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))));
     ans2=57568490411.0+y*(1029532985.0+y*(9494680.718
        +y*(59272.64853+y*(267.8532712+y*1.0))));
     ans=ans1/ans2;
  } else {
     z=8.0/ax;
     y=z*z;
     xx=ax-0.785398164;
     ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4
        +y*(-0.2073370639e-5+y*0.2093887211e-6)));
     ans2 = -0.1562499995e-1+y*(0.1430488765e-3
        +y*(-0.6911147651e-5+y*(0.7621095161e-6
        -y*0.934935152e-7)));
     ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
  }
  return ans;

}

su2double FWHSolver::bessj1 (su2double x){
  su2double ax,z;
  su2double xx,y,ans,ans1,ans2;

  if ((ax=fabs(x)) < 8.0) {
     y=x*x;
     ans1=x*(72362614232.0+y*(-7895059235.0+y*(242396853.1
        +y*(-2972611.439+y*(15704.48260+y*(-30.16036606))))));
     ans2=144725228442.0+y*(2300535178.0+y*(18583304.74
        +y*(99447.43394+y*(376.9991397+y*1.0))));
     ans=ans1/ans2;
  } else {
     z=8.0/ax;
     y=z*z;
     xx=ax-2.356194491;
     ans1=1.0+y*(0.183105e-2+y*(-0.3516396496e-4
        +y*(0.2457520174e-5+y*(-0.240337019e-6))));
     ans2=0.04687499995+y*(-0.2002690873e-3
        +y*(0.8449199096e-5+y*(-0.88228987e-6
        +y*0.105787412e-6)));
     ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
     if (x < 0.0) ans = -ans;
  }
  return ans;

}

su2double FWHSolver::bessy0 (su2double x){
  su2double z;
  su2double xx,y,ans,ans1,ans2;

  if (x < 8.0) {
     y=x*x;
     ans1 = -2957821389.0+y*(7062834065.0+y*(-512359803.6
        +y*(10879881.29+y*(-86327.92757+y*228.4622733))));
     ans2=40076544269.0+y*(745249964.8+y*(7189466.438
        +y*(47447.26470+y*(226.1030244+y*1.0))));
     ans=(ans1/ans2)+0.636619772*bessj0(x)*log(x);
  } else {
     z=8.0/x;
     y=z*z;
     xx=x-0.785398164;
     ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4
        +y*(-0.2073370639e-5+y*0.2093887211e-6)));
     ans2 = -0.1562499995e-1+y*(0.1430488765e-3
        +y*(-0.6911147651e-5+y*(0.7621095161e-6
        +y*(-0.934945152e-7))));
     ans=sqrt(0.636619772/x)*(sin(xx)*ans1+z*cos(xx)*ans2);
  }
  return ans;

}

su2double FWHSolver::bessy1 (su2double x){
  su2double z;
  su2double xx,y,ans,ans1,ans2;

  if (x < 8.0) {
     y=x*x;
     ans1=x*(-0.4900604943e13+y*(0.1275274390e13
        +y*(-0.5153438139e11+y*(0.7349264551e9
        +y*(-0.4237922726e7+y*0.8511937935e4)))));
     ans2=0.2499580570e14+y*(0.4244419664e12
        +y*(0.3733650367e10+y*(0.2245904002e8
        +y*(0.1020426050e6+y*(0.3549632885e3+y)))));
     ans=(ans1/ans2)+0.636619772*(bessj1(x)*log(x)-1.0/x);
  } else {
     z=8.0/x;
     y=z*z;
     xx=x-2.356194491;
     ans1=1.0+y*(0.183105e-2+y*(-0.3516396496e-4
        +y*(0.2457520174e-5+y*(-0.240337019e-6))));
     ans2=0.04687499995+y*(-0.2002690873e-3
        +y*(0.8449199096e-5+y*(-0.88228987e-6
        +y*0.105787412e-6)));
     ans=sqrt(0.636619772/x)*(sin(xx)*ans1+z*cos(xx)*ans2);
  }
  return ans;

}




void FWHSolver::Integrate_Sources(CConfig *config){
  unsigned long iObserver, iPanel, iSample;
  su2double dS, dt, w;
  complex <su2double> wxi;
  su2double pi =  acos(-1);
  complex <su2double> I = complex<su2double> (0, 1);

     if (TimeDomain3D){
         for (iObserver=0; iObserver<nObserver; iObserver++){
             for (iSample=0; iSample<nSample; iSample++){
                 for (iPanel=0; iPanel<nPanel; iPanel++){
                   pp_TimeDomain[iObserver][iSample] = pp_TimeDomain[iObserver][iSample] + pp_interp[iObserver][iPanel][iSample];
                   }
               }
           }

       }else{

       dt = config->GetDelta_UnstTime()*SamplingFreq;
       for (iObserver=0; iObserver<nObserver; iObserver++){
           for (iSample=0; iSample<nSample; iSample++){

               w = 2*pi*iSample/(nSample-1)/dt;
               wxi = w*I;
               for (iPanel=0; iPanel<nPanel; iPanel++){
                    dS = surface_geo[iPanel][4];
//                    fpp[iObserver][iSample] = fpp[iObserver][iSample] - dS*(wxi*fQ[iPanel][iSample]*G[iObserver][iPanel][iSample] +  fF1[iPanel][iSample]*dGdy1[iObserver][iPanel][iSample] + fF2[iPanel][iSample]*dGdy2[iObserver][iPanel][iSample]);
                    fpp[iObserver][iSample] = fpp[iObserver][iSample] - dS*(wxi*fQ[iPanel][iSample]*G[iObserver][iPanel][iSample] +  fF[0][iPanel][iSample]*dGdy1[iObserver][iPanel][iSample] + fF[1][iPanel][iSample]*dGdy2[iObserver][iPanel][iSample]);
                    if(nDim==3)  fpp[iObserver][iSample] -= dS*fF[2][iPanel][iSample]*dGdy3[iObserver][iPanel][iSample];
                 }
               fpp_r[iObserver][iSample] = real(fpp[iObserver][iSample]);
               fpp_i[iObserver][iSample] = imag(fpp[iObserver][iSample]);

 //              cout<<fpp[iObserver][iSample] <<",  "<<fpp_r[iObserver][iSample] <<",  "<<  fpp_i[iObserver][iSample]<<endl;

             }
         }
       }

}



void FWHSolver::iFFT_SignalR2(){
  unsigned long iObserver,   iSample;
  complex <su2double> I = complex<su2double> (0, 1);

     FFT* FFT_container = new FFT() ;

       //perform inverse FFT on each row (observer) of frequency-domain pp data
       for (iObserver=0; iObserver < nObserver;  iObserver++){
           cout<<"iFFT at Observer Location: "<<iObserver<<endl;

           //repack data in a conjugate symmetric fashion for inverse
           for (iSample= 0; iSample < nSample; iSample++)
             {
               pp[iObserver][iSample].real(fpp_r_root[iObserver][iSample]);
               pp[iObserver][iSample].imag(fpp_i_root[iObserver][iSample]);
//               pp[iObserver][iSample]= fpp_r_root[iObserver][iSample]+fpp_i_root[iObserver][iSample]*I;
             }

           for (iSample = 1; iSample < nSample/2; iSample++) {
            //   cout<<pp[iObserver][iSample]<<endl;
               pp[iObserver][nSample-iSample] = conj(pp[iObserver][iSample]);
             }
           CArray dataQ(pp[iObserver] ,nSample);


           FFT_container->ifft_r2(dataQ);

           //extract pressure fluctuation data from buffer
           for (iSample= 0; iSample < nSample; iSample++)
             {
               pp[iObserver][iSample].real(dataQ[iSample].real());
               pp[iObserver][iSample].imag(dataQ[iSample].imag());
 //               pp[iObserver][iSample]= real(dataQ[iSample])+imag(dataQ[iSample])*I;
             }


         }

}


void FWHSolver::Compute_SPL(){
  unsigned long iObserver, iSample;
  su2double SPL_iObserver;

       for (iObserver=0; iObserver<nObserver; iObserver++){
           SPL_iObserver = 0.0;
//           for (iSample=idx_window_l; iSample<idx_window_r; iSample++){
           for (iSample=0; iSample<nSample; iSample++){
//               cout<<pp[iObserver][iSample].real()<<endl;
               if (TimeDomain3D){
               SPL_iObserver  =  SPL_iObserver + pp_TimeDomain_root[iObserver][iSample]*pp_TimeDomain_root[iObserver][iSample];
                 }else{
               SPL_iObserver  =  SPL_iObserver + real(pp[iObserver][iSample])*real(pp[iObserver][iSample]);
                 }
             }
            SPL_iObserver = sqrt(SPL_iObserver/nSample);
            SPL = SPL + SPL_iObserver;
         }
      SPL = SPL/nObserver;

}





void FWHSolver::Write_Sensitivities(CSolver *solver, CConfig *config, CGeometry *geometry){
    unsigned long iVar, iPanel, iSample, start_iter,iExtIter,Max_nPanel, Tot_nPanel,nVar,Global_Index ;
      ofstream CAA_AdjointFile;
    char buffer [50];
    start_iter  =  config->GetUnst_RestartIter();

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
//cout<<Max_nPanel<<",  "<<nPanel<<",   " <<Tot_nPanel<<",   "<<rank<<endl;
      nVar = nDim+3;



      /* Loop through all time steps */
      for (iSample=0; iSample<nSample; iSample++){
      iExtIter = start_iter + iSample*SamplingFreq;

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
              Buffer_Send_dJdU[iVar*nPanel+iPanel] = dJdU[iVar][iPanel][iSample];
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
      char cstr [200];

      SPRINTF (cstr, "Adj_CAA");
      if ((SU2_TYPE::Int(iExtIter) >= 0)    && (SU2_TYPE::Int(iExtIter) < 10))    SPRINTF (buffer, "_0000%d.dat", SU2_TYPE::Int(iExtIter));
      if ((SU2_TYPE::Int(iExtIter) >= 10)   && (SU2_TYPE::Int(iExtIter) < 100))   SPRINTF (buffer, "_000%d.dat",  SU2_TYPE::Int(iExtIter));
      if ((SU2_TYPE::Int(iExtIter) >= 100)  && (SU2_TYPE::Int(iExtIter) < 1000))  SPRINTF (buffer, "_00%d.dat",   SU2_TYPE::Int(iExtIter));
      if ((SU2_TYPE::Int(iExtIter) >= 1000) && (SU2_TYPE::Int(iExtIter) < 10000)) SPRINTF (buffer, "_0%d.dat",    SU2_TYPE::Int(iExtIter));
      if (SU2_TYPE::Int(iExtIter) >= 10000) SPRINTF (buffer, "_%d.dat", SU2_TYPE::Int(iExtIter));
      strcat (cstr, buffer);
      cout<<cstr<<endl;
      CAA_AdjointFile.precision(15);
      CAA_AdjointFile.open(cstr, ios::out);

      /*--- Loop through all of the collected data and write each node's values ---*/
      unsigned long Total_Index;
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        for (iPanel = 0; iPanel < Buffer_Recv_nPanel[iProcessor]; iPanel++) {
            Global_Index = Buffer_Recv_GlobalIndex[iProcessor*Max_nPanel+iPanel ];
            CAA_AdjointFile  << scientific << Global_Index << "\t";

           for (iVar = 0; iVar < nVar; iVar++){

          /*--- Current index position and global index ---*/
          Total_Index  = iProcessor*Max_nPanel*nVar + iVar*Buffer_Recv_nPanel[iProcessor]  + iPanel;

          /*--- Write to file---*/
          CAA_AdjointFile << scientific <<  Buffer_Recv_dJdU[Total_Index]   << "\t";
           }
          CAA_AdjointFile  << endl;

         }
      }

      CAA_AdjointFile.close();
      delete [] Buffer_Recv_dJdU;
      delete [] Buffer_Recv_GlobalIndex;
       }

      delete [] Buffer_Send_dJdU;
      delete [] Buffer_Send_GlobalIndex;

    }



}






















SNG::SNG(CConfig *config,CGeometry *geometry) {

    unsigned long  i, iDim, iPoint, iSNGPt, SNGPtCount, nF, iT;
    su2double *Coord;
    su2double FreeStreamTemperature,FreeStreamPressure,FreeStreamDensity,M, AOA;
    su2double pi=3.141592653589793;
    su2double R = 287.058;

    nDim = geometry->GetnDim();

#ifdef HAVE_MPI
  int rank, nProcessor;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
     MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
#endif

//     cout<<"Entered SNG Constructor"<<endl;

    /* Reading in SNG boundary points -- change to config option later! */
    string  text_line;
    ifstream SNG_BoundaryFile;
    string filename = "NoiseSourceRegionDef.dat";
    SNG_BoundaryFile.open(filename.data() , ios::in);
    if (SNG_BoundaryFile.fail()) {
        cout << "There is no file!!! " <<  filename.data()  << "."<< endl;
      exit(EXIT_FAILURE);
    }

    NoiseSourceZone = new su2double* [2];

    for(i = 0;  i <2  ;  i++)
    {
       NoiseSourceZone[i] = new su2double[nDim];
       for (iDim=0; iDim < nDim; iDim++){
         NoiseSourceZone[i][iDim]= 0.0;
       }
    }

    i=0;
    getline (SNG_BoundaryFile, text_line);
    istringstream point_line(text_line);
    point_line >> f_min >> f_max >> NF >> GenNewRand;
    getline (SNG_BoundaryFile, text_line);
    istringstream point_line2(text_line);
    point_line2 >> dt >> NT>> N_Tij_Out;
    getline (SNG_BoundaryFile, text_line);
    istringstream point_line3(text_line);
    point_line3 >> Type_JBBN;
    while (getline (SNG_BoundaryFile, text_line)) {
        istringstream point_line(text_line);
        if (nDim==2){
        point_line >> NoiseSourceZone[i][0]>> NoiseSourceZone[i][1];
          }
        if (nDim==3){
        point_line >> NoiseSourceZone[i][0]>> NoiseSourceZone[i][1]>> NoiseSourceZone[i][2];
          }
        i++;
    }



  //   cout<<"Finished Reading Bounds"<<endl;

   if (rank==MASTER_NODE)  cout<<"f_min= "<<f_min<<", f_max= "<<f_max<<", N_F= "<<NF<<endl;
   if (rank==MASTER_NODE)  cout<<"dt= "<<dt<<", NT= "<<NT<<", Tij Ouput= "<<N_Tij_Out<<endl;

   T_ij_OutputIdx = new unsigned long [N_Tij_Out];
   for (int i=0; i<N_Tij_Out; i++){
    T_ij_OutputIdx[i] = NT/N_Tij_Out*(i);
//    if (rank==MASTER_NODE)  cout<<"i= "<<i<<", T_ij_OutputIdx[i]= "<<T_ij_OutputIdx[i]<<endl;
   }

    nSNGPts = 0;
    SNGPtCount = 0;
    //Go over the whole mesh and identify the points within the SNG zone
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
            if( geometry->node[iPoint]->GetDomain()){
        Coord = geometry->node[iPoint]->GetCoord();

        if (SU2_TYPE::GetValue(Coord[0])>NoiseSourceZone[0][0] &&  SU2_TYPE::GetValue(Coord[0])<NoiseSourceZone[1][0] && SU2_TYPE::GetValue(Coord[1])>NoiseSourceZone[0][1] && SU2_TYPE::GetValue(Coord[1])<NoiseSourceZone[1][1]){
//            cout<<"x= "<<SU2_TYPE::GetValue(Coord[0])<<"; y= "<<SU2_TYPE::GetValue(Coord[1]) <<endl;
//          SNG_Coords[SNGPtCount][0] = SU2_TYPE::GetValue(Coord[0]);
//          SNG_Coords[SNGPtCount][1] = SU2_TYPE::GetValue(Coord[1]);

          SNGPtCount++;
        }

      }
     }
    nSNGPts = SNGPtCount;
    cout<<"Process "<<rank<<" contains "<<nSNGPts <<" SNG Points."<<endl;


    SNG_Coords = new su2double* [nSNGPts];
    for(iSNGPt = 0; iSNGPt<nSNGPts; iSNGPt++)
    {
       SNG_Coords[iSNGPt] = new su2double[3];
       SNG_Coords[iSNGPt][0]= 0.0;
       SNG_Coords[iSNGPt][1]= 0.0;
       SNG_Coords[iSNGPt][2]= 0.0;
    }

    TKE = new su2double [nSNGPts];
    omega = new su2double [nSNGPts];
    SNG_CellVol = new su2double [nSNGPts];
    PointID = new unsigned long [nSNGPts];
    for(iSNGPt = 0; iSNGPt<nSNGPts; iSNGPt++)
    {
       TKE [iSNGPt] = 0.0;
       omega [iSNGPt] = 0.0;
       SNG_CellVol [iSNGPt] = 0.0;
       PointID[iSNGPt] = 0;
    }


    k_n = new su2double* [NF];
    sigma_n = new su2double* [NF];
    Psi_n = new su2double [NF];
    for(nF = 0; nF<NF; nF++)
    {
       k_n[nF] = new su2double[3];
       sigma_n[nF] = new su2double[3];
       Psi_n[nF] = 0.0;
       for (iDim=0; iDim < 3; iDim++){
         k_n[nF][iDim]= 0.0;
         sigma_n[nF][iDim]= 0.0;
       }
    }



    M = config->GetMach();
    FreeStreamPressure=config->GetPressure_FreeStream();
    FreeStreamTemperature = config->GetTemperature_FreeStream();
    FreeStreamDensity = FreeStreamPressure/R/FreeStreamTemperature;
    a_inf = sqrt(config->GetGamma()*FreeStreamPressure / FreeStreamDensity);
    AOA = config->GetAoA();
    U1 = M*a_inf*cos(AOA*pi/180) ;
    U2 = M*a_inf*sin(AOA*pi/180) ;
    U3 = 0.0;    //FIX THIS LATER!!!

    //Multiplicative factors to re-dimensionalize TKE and omega
    //Missing Gamma=1.4 factor but this is how TKE and omega are nondimensionalized in SU2 RANS
    TKE_ReDimFac = R*FreeStreamTemperature;
    omega_ReDimFac = sqrt(TKE_ReDimFac);



    if (rank==MASTER_NODE) cout<<U1<<",  "<<U2<<",  "<<M<<",  "<<a_inf<<",  "<<FreeStreamPressure<<", "<<FreeStreamDensity<<", "<<FreeStreamTemperature<<endl;

//    if (rank==MASTER_NODE)  cout<<"FS Turb Intensity= "<<config->GetTurbulenceIntensity_FreeStream()<<endl;
//    if (rank==MASTER_NODE)  cout<<"FS Ux= "<<config->GetVelocity_FreeStream()[0]<<", FS Uy= "<<config->GetVelocity_FreeStream()[1]<<endl;
//    if (rank==MASTER_NODE)  cout<<"FS Ux_ND= "<<config->GetVelocity_FreeStreamND()[0]<<", FS Uy_ND= "<<config->GetVelocity_FreeStreamND()[1]<<endl;
//    if (rank==MASTER_NODE)  cout<<"FS p_ND= "<<config->GetPressure_FreeStreamND()<<", FS p_ref= "<<config->GetPressure_Ref()<<endl;
//    if (rank==MASTER_NODE)  cout<<"FS mu= "<<config->GetViscosity_FreeStream()<<", FS mu_ND= "<<config->GetViscosity_FreeStreamND()<<", FS mu_ref= "<<config->GetViscosity_Ref()<<endl;


    u_turb = new su2double**[nSNGPts];
    for (iSNGPt = 0; iSNGPt<nSNGPts; iSNGPt++)
    {
        u_turb[iSNGPt] = new su2double*[nDim];
        for(iDim = 0; iDim<nDim; iDim++)
        {
            u_turb[iSNGPt][iDim] = new su2double [NT];
            for(iT=0; iT<NT; iT++){
                u_turb[iSNGPt][iDim][iT] = 0.0;
            }
        }
    }


    T_tilda = new su2double**[nSNGPts];
    T_tilda_mean = new su2double*[nSNGPts];
    for (iSNGPt = 0; iSNGPt<nSNGPts; iSNGPt++)
    {
        T_tilda[iSNGPt] = new su2double*[3*(nDim-1)];
        T_tilda_mean[iSNGPt] = new su2double[3*(nDim-1)];
        for(iDim = 0; iDim<3*(nDim-1); iDim++)
        {
            T_tilda_mean[iSNGPt][iDim] = 0.0;
            T_tilda[iSNGPt][iDim] = new su2double [NT];
            for(iT=0; iT<NT; iT++){
                T_tilda[iSNGPt][iDim][iT] = 0.0;
            }
        }
    }

    J_BBN = 0.0;

    dJBBN_dU = new su2double* [nDim+2];
    for(int iDim = 0; iDim < nDim+2 ; iDim++)
    {
        dJBBN_dU[iDim] = new su2double[nSNGPts];
        for (iSNGPt = 0; iSNGPt<nSNGPts; iSNGPt++)
        {
               dJBBN_dU[iDim][iSNGPt]= 0.0;
        }
    }



//        dJBBN_dU = new su2double[nSNGPts];
//        for (iSNGPt = 0; iSNGPt<nSNGPts; iSNGPt++)
//        {
//               dJBBN_dU[iSNGPt]= 0.0;
//        }





}






SNG::~SNG(void) {

}



void SNG::SetSNG_Analysis(CSolver *solver, CConfig *config,CGeometry *geometry ){
#ifdef HAVE_MPI
  int rank, nProcessor;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
     MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
#endif
//    cout<<"In SetSNG_Analysis"<<endl;

    //Extract TKE, omega, cell volume and cell coordinates within the SNG zone from the RANS solution
    Extract_RANS(solver, config, geometry );
    if (rank == MASTER_NODE) cout<<"Finished RANS Extraction"<<endl;

    //Set all NF random Fourier modes
    SetRandomFourierModes();



}



void SNG::Perform_SNG_Analysis(){

#ifdef HAVE_MPI
  int rank, nProcessor;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
     MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
#endif
   if (rank == MASTER_NODE) cout<<"Performing SNG Analysis"<<endl;


    //Synthesize turbulent velocity fields
    Compute_TurbVelocity();
    if (rank == MASTER_NODE) cout<<"u_turb synthesized"<<endl;

    //Compute broadband noise source
    Compute_BroadBandNoiseSource();
    if (rank == MASTER_NODE) cout<<"Tij Computed"<<endl;

    //Compute broadband noise objective function
     Compute_BBN_ObjFunc();
   if (rank == MASTER_NODE) cout<<"J_BBN Computed"<<endl;

    //Merge broadband noise source to root and write out required snapshots
   if (rank == MASTER_NODE) cout<<"Writing BBN Source to file"<<endl;
    Write_BroadBandNoiseSource();

}





void SNG::SetRandomFourierModes(){
  unsigned long  nF;
  su2double pi =  acos(-1.0);
  su2double f_n, k_mod_n, Phi_n, Theta_n,Phi_n1, Theta_n1, kt0, kt1, kt2, kcoef, sig1, sig2, sig3, sig_mag;
  ofstream RandNum_File_OUT; ifstream RandNum_File_IN;
  string  text_line;

#ifdef HAVE_MPI
  int rank, nProcessor;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
     MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
#endif


    srand(time(NULL));
    if (rank == MASTER_NODE){
    if (GenNewRand){
      cout<<"Generating a new set of random angles for SNG" << endl;
      RandNum_File_OUT.open("RandNumFile", ios::out);
    }
    else{
       cout<<"Reading pre-generated random angles for SNG" << endl;
       RandNum_File_IN.open("RandNumFile", ios::in);
       if (RandNum_File_IN.fail()) {
           cout << "There is no RandNumFile!!! "<< endl;
         exit(EXIT_FAILURE);
      }
    }
     for(nF = 0; nF<NF; nF++)
     {
       f_n = f_min + (f_max-f_min)/(NF-1)*nF;
       k_mod_n = 2*pi*f_n/a_inf;

       //Read random angles from a pre-generated file.
       if (!GenNewRand){
          getline (RandNum_File_IN, text_line);
          istringstream point_line(text_line);
          point_line >> Psi_n[nF] >> Phi_n >>Theta_n >> Phi_n1 >>Theta_n1;
       }

       //Random generation of Psi_n
       if (GenNewRand)Psi_n[nF] = ((double)rand()/(double)RAND_MAX)*2.0*pi - pi;


       //Random generation of Phi_n and Theta_n
       if (GenNewRand){
       Phi_n = ((double)rand()/(double)RAND_MAX)*pi;
       Theta_n = ((double)rand()/(double)RAND_MAX)*2*pi;
       }

       //Convert k_mod_n, Phi_n and Theta_n in spherical coord to cartesian coord for k_n
       k_n[nF][0]= k_mod_n*sin(Theta_n)*cos(Phi_n);
       k_n[nF][1]= k_mod_n*sin(Theta_n)*sin(Phi_n);
       k_n[nF][2]= k_mod_n*cos(Theta_n);

       //Random generation of Phi_n1 and Theta_n1
       if (GenNewRand){
       Phi_n1 = ((double)rand()/(double)RAND_MAX)*pi;
       Theta_n1 = ((double)rand()/(double)RAND_MAX)*2*pi;
       }

       //Compute sigma_n by projecting the unit vector defined by Phi_n1 and Theta_n1 onto a plane orthogonal to k_n
       kt0 = sin(Theta_n1)*cos(Phi_n1);
       kt1 = sin(Theta_n1)*sin(Phi_n1);
       kt2 = cos(Theta_n1);
       kcoef = (k_n[nF][0]*kt0 + k_n[nF][1]*kt1 + k_n[nF][2]*kt2)/k_mod_n/k_mod_n;
       sig1 = kt0 - k_n[nF][0]*kcoef;
       sig2 = kt1 - k_n[nF][1]*kcoef;
       sig3 = kt2 - k_n[nF][2]*kcoef;
       sig_mag = sqrt(sig1*sig1+sig2*sig2+sig3*sig3);
       //scale to unit vector
       sigma_n[nF][0]= sig1/sig_mag;
       sigma_n[nF][1]= sig2/sig_mag;
       sigma_n[nF][2]= sig3/sig_mag;

       if (GenNewRand) RandNum_File_OUT<< std::setprecision(15) << Psi_n[nF]<<"  "<<Phi_n<<"  "<<Theta_n<<"  "<<Phi_n1<<"  "<<Theta_n1<<endl;
     }
      if (GenNewRand) RandNum_File_OUT.close();
      else RandNum_File_IN.close();

    }


    //Broadcast to all other processes
#ifdef HAVE_MPI
    SU2_MPI::Bcast(Psi_n, NF, MPI_DOUBLE,  MASTER_NODE, MPI_COMM_WORLD);
    for (int i = 0; i<NF; i++) SU2_MPI::Bcast(k_n[i], 3, MPI_DOUBLE,  MASTER_NODE, MPI_COMM_WORLD);
    for (int i = 0; i<NF; i++) SU2_MPI::Bcast(sigma_n[i], 3, MPI_DOUBLE,  MASTER_NODE, MPI_COMM_WORLD);
#endif


//    for (int i=0; i<nProcessor; i++){
//            if (rank==i){
//    //            cout<<"Printing Psi_n of Processor "<<rank<<endl;
//                cout<<"Printing k_n of Processor "<<rank<<endl;
//    for(nF = 0; nF<NF; nF++){
//        for(int iDim=0; iDim<3; iDim++){
//            cout<< std::setprecision(15)<<k_n[nF][iDim]<<" ";
//        }
//        cout<<endl;
//    }
//     cout<<"Printing sigma_n of Processor "<<rank<<endl;
//    for(nF = 0; nF<NF; nF++){
//        for(int iDim=0; iDim<3; iDim++){
//            cout<< std::setprecision(15)<<sigma_n[nF][iDim]<<" ";
//        }
//        cout<<endl;
//    }
//            }
//    }



}







void SNG::Extract_RANS(CSolver *solver, CConfig *config, CGeometry *geometry){
    unsigned long   i, iPoint,  SNGPtCount;
    su2double *Coord;
    su2double Turb_Kinetic_Energy, x_coord, y_coord, z_coord, Specific_Dissip_Rate;

#ifdef HAVE_MPI
  int rank, nProcessor;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
     MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
#endif
    Turb_Kinetic_Energy = 0.0;
    Specific_Dissip_Rate=0.0;
    SNGPtCount = 0;
    //Go over the whole mesh and identify the points within the SNG zone
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
            if( geometry->node[iPoint]->GetDomain()){
        Coord = geometry->node[iPoint]->GetCoord();

        if (SU2_TYPE::GetValue(Coord[0])>NoiseSourceZone[0][0] &&  SU2_TYPE::GetValue(Coord[0])<NoiseSourceZone[1][0] && SU2_TYPE::GetValue(Coord[1])>NoiseSourceZone[0][1] && SU2_TYPE::GetValue(Coord[1])<NoiseSourceZone[1][1]){
          x_coord = SU2_TYPE::GetValue(Coord[0]);
          y_coord = SU2_TYPE::GetValue(Coord[1]);
          if(nDim==3) z_coord = SU2_TYPE::GetValue(Coord[2]);
      //    SNG_CellVol[SNGPtCount] = geometry->node[iPoint]->GetVolume();
          SNG_CellVol[SNGPtCount] = 1.0;
     //     cout<<"i= "<< SNGPtCount<<", Cell Vol= "<<geometry->node[iPoint]->GetVolume()<<",  "<<geometry->node[iPoint]->GetVolume_n() <<endl;
          Turb_Kinetic_Energy =  solver->node[iPoint]->GetSolution(2*nDim+2);
          Specific_Dissip_Rate =  solver->node[iPoint]->GetSolution(2*nDim+3);


          if (config->GetAD_Mode()){
          AD::RegisterInput(x_coord);
          AD::RegisterInput(y_coord);
          if(nDim==3)AD::RegisterInput(z_coord);
          AD::RegisterInput(Turb_Kinetic_Energy);
          AD::RegisterInput(Specific_Dissip_Rate);
          }

          SNG_Coords[SNGPtCount][0] = x_coord;
          SNG_Coords[SNGPtCount][1] = y_coord;
          if(nDim==3)SNG_Coords[SNGPtCount][2] = z_coord;

          TKE[SNGPtCount] = Turb_Kinetic_Energy;
          omega[SNGPtCount] = Specific_Dissip_Rate;

          PointID[SNGPtCount] = geometry->node[iPoint]->GetGlobalIndex();
//          cout<<"Registering TKE and omega, pass: "<< SNGPtCount<<", TKE= "<<Turb_Kinetic_Energy<<endl;
          SNGPtCount++;
        }

      }
     }


    if (nProcessor==1){
    ofstream SNG_file;
    SNG_file.open("SNG");

    for (i=0;i<nSNGPts; i++){

        SNG_file<< std::setprecision(15) << SNG_Coords[i][0]  <<"  "<<SNG_Coords[i][1]<<"  "<<TKE[i]<<"  "<<omega[i]<<endl;

     }
    SNG_file.close();
    }

    Write_ExtractedRANS();



}




void SNG::Compute_TurbVelocity(){
    su2double A = 1.453;
    su2double c1 = 1.0;
    su2double beta = 0.09;
    su2double nu = 1.57E-5; //kinematic viscosity of air at 25 deg Celsius
    su2double pi =  acos(-1.0);
    unsigned long nF, iSNGPt, iT;
    su2double u_n, E_kn, del_kn, f_n, k_mod_n, k_e, k_eta, epsilon, L_T, u_prime;
    su2double *convection_term, dot_prod, t;

    convection_term = new su2double[nDim];
    for (int iDim=0; iDim<nDim; iDim++) convection_term[iDim] = 0.0;

    del_kn = (f_max-f_min)/(NF-1);

    for(nF=0; nF<NF; nF++){
        f_n = f_min + (f_max-f_min)/(NF-1)*nF;
        k_mod_n = 2*pi*f_n/a_inf;
    //    cout<<k_mod_n<<", "<<pi<<", "<<f_n<<", "<<a_inf<<endl;
        for(iSNGPt=0; iSNGPt<nSNGPts; iSNGPt++){

            //Compute epsilon from TKE and omega
            epsilon = beta*TKE[iSNGPt]*TKE_ReDimFac*omega[iSNGPt]*omega_ReDimFac;

            //Compute k_e and k_eta
            u_prime = sqrt(2.0*TKE[iSNGPt]*TKE_ReDimFac/3.0);
            L_T = c1*pow(u_prime,3.0)/epsilon;
            k_e = 0.747/L_T;
            k_eta = pow(epsilon,0.25)*pow(nu,-0.75);

            //Compute E at k_n
            E_kn = 2.0*A/3.0*TKE[iSNGPt]*TKE_ReDimFac/k_e*pow(k_mod_n/k_e,4.0)*exp(-2*pow(k_mod_n/k_eta,2.0))*pow(1.0+pow(k_mod_n/k_e,2.0),-17.0/6.0);

            //Compute velocity magnitude u_n
            u_n = sqrt(E_kn*del_kn);

            for(iT=0; iT<NT; iT++){
                t = dt*iT;
                convection_term[0] = SNG_Coords[iSNGPt][0] - U1*t;
                convection_term[1] = SNG_Coords[iSNGPt][1] - U2*t;
                if (nDim==3) convection_term[2] = SNG_Coords[iSNGPt][2] - U3*t;
                dot_prod = 0.0;
                for (int iDim=0; iDim<nDim; iDim++) dot_prod += k_n[nF][iDim]*convection_term[iDim];
                u_turb[iSNGPt][0][iT] += 2*u_n*cos(dot_prod+Psi_n[nF])*sigma_n[nF][0];
                u_turb[iSNGPt][1][iT] += 2*u_n*cos(dot_prod+Psi_n[nF])*sigma_n[nF][1];
                if (nDim==3) u_turb[iSNGPt][2][iT] += 2*u_n*cos(dot_prod+Psi_n[nF])*sigma_n[nF][2];
            }
        //    cout<<dot_prod<<", "<<u_n<<", "<<epsilon<<", "<<k_mod_n<<", "<< convection_term[0]<<", "<<convection_term[1]<<", "<< t <<endl;

        }

    }


}



void SNG::Compute_BroadBandNoiseSource(){
    unsigned long iSNGPt, iT,idx_out,iVar, idx;
    ofstream T_tilde_File;
        char buffer [50];
#ifdef HAVE_MPI
    int rank, nProcessor;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
#endif

    for(iT=0; iT<NT; iT++){
        for (iSNGPt=0; iSNGPt<nSNGPts; iSNGPt++){
            T_tilda[iSNGPt][0][iT] = u_turb[iSNGPt][0][iT]*u_turb[iSNGPt][0][iT];  //T_11 = u1*u1
            T_tilda[iSNGPt][1][iT] = u_turb[iSNGPt][0][iT]*u_turb[iSNGPt][1][iT];  //T_12 = u1*u2
            T_tilda[iSNGPt][2][iT] = u_turb[iSNGPt][1][iT]*u_turb[iSNGPt][1][iT];  //T_22 = u2*u2
            T_tilda_mean[iSNGPt][0] += T_tilda[iSNGPt][0][iT]/NT;
            T_tilda_mean[iSNGPt][1] += T_tilda[iSNGPt][1][iT]/NT;
            T_tilda_mean[iSNGPt][2] += T_tilda[iSNGPt][2][iT]/NT;
            if (nDim == 3){
                T_tilda[iSNGPt][3][iT] = u_turb[iSNGPt][0][iT]*u_turb[iSNGPt][2][iT];  //T_13 = u1*u3
                T_tilda[iSNGPt][4][iT] = u_turb[iSNGPt][1][iT]*u_turb[iSNGPt][2][iT];  //T_23 = u2*u3
                T_tilda[iSNGPt][5][iT] = u_turb[iSNGPt][2][iT]*u_turb[iSNGPt][2][iT];  //T_33 = u3*u3
                T_tilda_mean[iSNGPt][3] += T_tilda[iSNGPt][3][iT]/NT;
                T_tilda_mean[iSNGPt][4] += T_tilda[iSNGPt][4][iT]/NT;
                T_tilda_mean[iSNGPt][5] += T_tilda[iSNGPt][5][iT]/NT;
            }
        }
    }


    //Write out T_ij files only if it is in serial
    if (nProcessor==1){
        for (idx=0; idx<N_Tij_Out; idx++){
            idx_out = T_ij_OutputIdx[idx];
                  char cstr [200];
            SPRINTF (cstr, "Tij_Master");
            if ((SU2_TYPE::Int(idx_out) >= 0)    && (SU2_TYPE::Int(idx_out) < 10))    SPRINTF (buffer, "_0000%d.dat", SU2_TYPE::Int(idx_out));
            if ((SU2_TYPE::Int(idx_out) >= 10)   && (SU2_TYPE::Int(idx_out) < 100))   SPRINTF (buffer, "_000%d.dat",  SU2_TYPE::Int(idx_out));
            if ((SU2_TYPE::Int(idx_out) >= 100)  && (SU2_TYPE::Int(idx_out) < 1000))  SPRINTF (buffer, "_00%d.dat",   SU2_TYPE::Int(idx_out));
            if ((SU2_TYPE::Int(idx_out) >= 1000) && (SU2_TYPE::Int(idx_out) < 10000)) SPRINTF (buffer, "_0%d.dat",    SU2_TYPE::Int(idx_out));
            if (SU2_TYPE::Int(idx_out) >= 10000) SPRINTF (buffer, "_%d.dat", SU2_TYPE::Int(idx_out));
            strcat (cstr, buffer);
            cout<<cstr<<endl;
            T_tilde_File.precision(15);
            T_tilde_File.open(cstr, ios::out);
            for (iSNGPt=0; iSNGPt<nSNGPts; iSNGPt++){
                for(iVar=0; iVar<3*(nDim-1); iVar++){
                 T_tilde_File <<scientific<< T_tilda[iSNGPt][iVar][idx_out] << "  ";
                }
                 T_tilde_File <<endl;
            }
                 T_tilde_File.close();
        }


        ofstream Tij_Mean_file;
        Tij_Mean_file.open("Tij_Master_Mean");
        Tij_Mean_file.precision(15);
        for (iSNGPt=0; iSNGPt<nSNGPts; iSNGPt++){
            for(iVar=0; iVar<3*(nDim-1); iVar++){
             Tij_Mean_file <<scientific<< T_tilda_mean[iSNGPt][iVar] << "  ";
            }
             Tij_Mean_file <<endl;
        }
             Tij_Mean_file.close();


    }
}



void SNG::Compute_BBN_ObjFunc(){
    unsigned long iSNGPt, nVar,iVar;
    su2double TotCellVol_Local;
    su2double TotCellVol_Global;
    TotCellVol_Local=0.0;
    TotCellVol_Global = 0.0;
    su2double Max_T22_Local=-1e16;
    su2double Max_T22_Global=0.0;

#ifdef HAVE_MPI
   int rank, nProcessor;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
#endif




    nVar = 3*(nDim-1);
    su2double *T_tilda_mean_VolAvg_Local = new su2double [nVar];
    su2double *T_tilda_mean_VolAvg_Global = new su2double [nVar];
    for (iVar=0; iVar <nVar; iVar++){
      T_tilda_mean_VolAvg_Local[iVar] = 0.0;
      T_tilda_mean_VolAvg_Global[iVar] = 0.0;
    }

    for (iSNGPt=0; iSNGPt<nSNGPts; iSNGPt++){
        for (iVar=0; iVar <nVar; iVar++){
            T_tilda_mean_VolAvg_Local[iVar] +=  T_tilda_mean[iSNGPt][iVar]*SNG_CellVol [iSNGPt];
        }
        TotCellVol_Local += SNG_CellVol [iSNGPt];
        if(T_tilda_mean[iSNGPt][2]>Max_T22_Local) Max_T22_Local = T_tilda_mean[iSNGPt][2];
    }


#ifdef HAVE_MPI
    SU2_MPI::Reduce(T_tilda_mean_VolAvg_Local ,T_tilda_mean_VolAvg_Global, nVar,MPI_DOUBLE,MPI_SUM,MASTER_NODE,MPI_COMM_WORLD);
    SU2_MPI::Reduce( &TotCellVol_Local, &TotCellVol_Global, 1, MPI_DOUBLE,MPI_SUM,MASTER_NODE,MPI_COMM_WORLD);
    SU2_MPI::Reduce( &Max_T22_Local, &Max_T22_Global, 1, MPI_DOUBLE,MPI_MAX,MASTER_NODE,MPI_COMM_WORLD);
#endif

    if (rank == MASTER_NODE){
    for (iVar=0; iVar <nVar; iVar++){
        T_tilda_mean_VolAvg_Global[iVar] = T_tilda_mean_VolAvg_Global[iVar]/TotCellVol_Global;
    }

       cout<<"Time-averaged over "<<NT<<" samples & spatially-averaged: T11= "<<  T_tilda_mean_VolAvg_Global[0] <<", T12= "<<  T_tilda_mean_VolAvg_Global[1] <<", T22= "<<  T_tilda_mean_VolAvg_Global[2]<<endl;
       cout<<"T22_Max= "<< Max_T22_Global<<endl;

       //BBN Objective is the MAX (in space) of T22 component of the time-averaged Lighthill's stress tensor
       if (Type_JBBN==-1){
           J_BBN = Max_T22_Global;
           cout<<"J_BBN set to spatial max of T22. J_BBN= "<<J_BBN<<endl;
       }

       //BBN Objective is the Frobenius norm of the time and spatially-averaged Lighthill's stress tensor
       if (Type_JBBN== 0) {
           for (iVar=0; iVar <nVar; iVar++) J_BBN += T_tilda_mean_VolAvg_Global[iVar];
           J_BBN = sqrt(J_BBN);
           cout<<"J_BBN set to Frobenius norm of the time and spatially-averaged Lighthill's stress tensor. J_BBN= "<<J_BBN<<endl;
       }

       //BBN Objective is the magnitude of one component of the time and spatially-averaged Lighthill's stress tensor
       //Type_JBBN=1 -> J=|T11|, Type_JBBN=2 -> J=|T12|, Type_JBBN=3 -> J=|T22|
       else{
           J_BBN = abs(T_tilda_mean_VolAvg_Global[Type_JBBN-1]);
           cout<<"J_BBN set to magnitude of component "<<  Type_JBBN-1<<  " of the time and spatially-averaged Lighthill's stress tensor. J_BBN= " <<J_BBN<<endl;
       }


       ofstream J_BBN_file;
       J_BBN_file.open("J_BBN");
       J_BBN_file.precision(15);
            J_BBN_file <<scientific<< J_BBN<<endl;
            J_BBN_file.close();

    }
}




void SNG::Write_ExtractedRANS(){
    unsigned long  iSNGPts, Max_nSNGPts, Tot_nSNGPts;
      ofstream ExtractedRANS_Merged;


    int rank, iProcessor, nProcessor;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
    unsigned long Buffer_Send_nSNGPts[1], *Buffer_Recv_nSNGPts = NULL;



    if (rank == MASTER_NODE) Buffer_Recv_nSNGPts= new unsigned long [nProcessor];

     Buffer_Send_nSNGPts[0]= nSNGPts;
#ifdef HAVE_MPI
      SU2_MPI::Gather(&Buffer_Send_nSNGPts, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nSNGPts, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD); //send the number of vertices at each process to the master
      SU2_MPI::Allreduce(&nSNGPts,&Max_nSNGPts,1,MPI_UNSIGNED_LONG,MPI_MAX,MPI_COMM_WORLD); //find the max num of vertices over all processes
      SU2_MPI::Reduce(&nSNGPts,&Tot_nSNGPts,1,MPI_UNSIGNED_LONG,MPI_SUM,MASTER_NODE,MPI_COMM_WORLD); //find the total num of vertices (SNGPtss)
#endif

cout<<"Max nSNGPts for all processes: "<<Max_nSNGPts<<",  nSNGPTs of Process "<<rank<<": "<<nSNGPts<<", Total Num of SNG Pts: " <<Tot_nSNGPts<<endl;


      /* pack sensitivity values in each processor and send to root */
      su2double *Buffer_Send_xCoord = new su2double [Max_nSNGPts];
      su2double *Buffer_Send_yCoord = new su2double [Max_nSNGPts];
      su2double *Buffer_Send_TKE = new su2double [Max_nSNGPts];
      su2double *Buffer_Send_omega = new su2double [Max_nSNGPts];

      //zero send buffers
      for (int i=0; i <Max_nSNGPts; i++){
       Buffer_Send_xCoord[i]=0.0;
       Buffer_Send_yCoord[i]=0.0;
       Buffer_Send_TKE[i]=0.0;
       Buffer_Send_omega[i]=0.0;
      }

      su2double *Buffer_Recv_xCoord = NULL;
      su2double *Buffer_Recv_yCoord = NULL;
      su2double *Buffer_Recv_TKE = NULL;
      su2double *Buffer_Recv_omega = NULL;

      if (rank == MASTER_NODE) {
       Buffer_Recv_xCoord = new su2double [nProcessor*Max_nSNGPts];
       Buffer_Recv_yCoord = new su2double [nProcessor*Max_nSNGPts];
       Buffer_Recv_TKE = new su2double [nProcessor*Max_nSNGPts];
       Buffer_Recv_omega = new su2double [nProcessor*Max_nSNGPts];
      }

      for(iSNGPts=0; iSNGPts<nSNGPts; iSNGPts++){
          Buffer_Send_xCoord[iSNGPts] = SNG_Coords[iSNGPts][0];
          Buffer_Send_yCoord[iSNGPts] = SNG_Coords[iSNGPts][1];
          Buffer_Send_TKE[iSNGPts] = TKE[iSNGPts];
          Buffer_Send_omega[iSNGPts] = omega[iSNGPts];
      }


#ifdef HAVE_MPI
     SU2_MPI::Gather(Buffer_Send_xCoord, Max_nSNGPts, MPI_DOUBLE, Buffer_Recv_xCoord,  Max_nSNGPts , MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
     SU2_MPI::Gather(Buffer_Send_yCoord, Max_nSNGPts, MPI_DOUBLE, Buffer_Recv_yCoord,  Max_nSNGPts , MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
     SU2_MPI::Gather(Buffer_Send_TKE, Max_nSNGPts, MPI_DOUBLE, Buffer_Recv_TKE,  Max_nSNGPts , MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
     SU2_MPI::Gather(Buffer_Send_omega, Max_nSNGPts, MPI_DOUBLE, Buffer_Recv_omega,  Max_nSNGPts , MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#endif

     /* root opens a file at each time step and write out the merged dJdU values at that time step into the file */
      if (rank == MASTER_NODE){
      char cstr [200];

      SPRINTF (cstr, "ExtractedRANS_Merged");

      ExtractedRANS_Merged.precision(15);
      ExtractedRANS_Merged.open(cstr, ios::out);

      /*--- Loop through all of the collected data and write each node's values ---*/
      unsigned long Total_Index;
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        for (iSNGPts = 0; iSNGPts < Buffer_Recv_nSNGPts[iProcessor]; iSNGPts++) {



          /*--- Current index position and global index ---*/
          Total_Index  = iProcessor*Max_nSNGPts + iSNGPts;

          /*--- Write to file---*/
          ExtractedRANS_Merged << scientific <<  Buffer_Recv_xCoord[Total_Index]   << "  "<<  Buffer_Recv_yCoord[Total_Index]   << "  "<<  Buffer_Recv_TKE[Total_Index]   << "  "<<  Buffer_Recv_omega[Total_Index]   << "  "  << endl;

         }
      }

      ExtractedRANS_Merged.close();
      delete [] Buffer_Recv_xCoord;
      delete [] Buffer_Recv_yCoord;
      delete [] Buffer_Recv_TKE;
      delete [] Buffer_Recv_omega;
       }
      delete [] Buffer_Send_xCoord;
      delete [] Buffer_Send_yCoord;
      delete [] Buffer_Send_TKE;
      delete [] Buffer_Send_omega;

}





void SNG::Write_BroadBandNoiseSource(){

    unsigned long iVar, idx, idx_out, iSNGPts, Max_nSNGPts, Tot_nSNGPts,nVar;
      ofstream T_tilde_File;  ofstream T_tilde_Mean_File;
    char buffer [50];
    int rank, iProcessor, nProcessor;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
    unsigned long Buffer_Send_nSNGPts[1], *Buffer_Recv_nSNGPts = NULL;



    if (rank == MASTER_NODE) Buffer_Recv_nSNGPts= new unsigned long [nProcessor];

      Buffer_Send_nSNGPts[0]=nSNGPts;
#ifdef HAVE_MPI
      SU2_MPI::Gather(&Buffer_Send_nSNGPts, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nSNGPts, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD); //send the number of vertices at each process to the master
      SU2_MPI::Allreduce(&nSNGPts,&Max_nSNGPts,1,MPI_UNSIGNED_LONG,MPI_MAX,MPI_COMM_WORLD); //find the max num of vertices over all processes
      SU2_MPI::Reduce(&nSNGPts,&Tot_nSNGPts,1,MPI_UNSIGNED_LONG,MPI_SUM,MASTER_NODE,MPI_COMM_WORLD); //find the total num of vertices (panels)
#endif

      cout<<"Max nSNGPts for all processes: "<<Max_nSNGPts<<",  nSNGPTs of Process "<<rank<<": "<<nSNGPts<<", Total Num of SNG Pts: " <<Tot_nSNGPts<<endl;


      nVar = 3*(nDim-1);



      /* Loop the list of Tij Output Indices */
      for (idx=0; idx<N_Tij_Out; idx++){

      idx_out = T_ij_OutputIdx[idx];

      /* pack sensitivity values in each processor and send to root */
      su2double *Buffer_Send_Tij = new su2double [Max_nSNGPts*nVar];
      su2double *Buffer_Send_Tij_Mean = new su2double [Max_nSNGPts*nVar];

      //zero send buffers
      for (int i=0; i <Max_nSNGPts*nVar; i++){
       Buffer_Send_Tij[i]=0.0;
       if (idx==0)  Buffer_Send_Tij_Mean[i]=0.0;
      }

      su2double *Buffer_Recv_Tij = NULL;
      su2double *Buffer_Recv_Tij_Mean = NULL;

      if (rank == MASTER_NODE) {
       Buffer_Recv_Tij = new su2double [nProcessor*Max_nSNGPts*nVar];
       if (idx==0) Buffer_Recv_Tij_Mean = new su2double [nProcessor*Max_nSNGPts*nVar];
      }

      for (iVar=0; iVar<nVar; iVar++){
          for(iSNGPts=0; iSNGPts<nSNGPts; iSNGPts++){
              Buffer_Send_Tij[iVar*nSNGPts+iSNGPts] = T_tilda[iSNGPts][iVar][idx_out];
              if (idx==0)Buffer_Send_Tij_Mean[iVar*nSNGPts+iSNGPts]= T_tilda_mean[iSNGPts][iVar];
            }
      }



#ifdef HAVE_MPI
     SU2_MPI::Gather(Buffer_Send_Tij, Max_nSNGPts*nVar, MPI_DOUBLE, Buffer_Recv_Tij,  Max_nSNGPts*nVar , MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
     if (idx==0) SU2_MPI::Gather(Buffer_Send_Tij_Mean, Max_nSNGPts*nVar, MPI_DOUBLE, Buffer_Recv_Tij_Mean,  Max_nSNGPts*nVar , MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);

#endif

     /* root opens a file at each time step and write out the merged dJdU values at that time step into the file */
      if (rank == MASTER_NODE){
      char cstr [200];
      if (idx==0){
          char cstr2 [200];
          SPRINTF (cstr2, "Tij_Merged_Mean");
          cout<<"Writing Merged File: "<<cstr2<<endl;
          T_tilde_Mean_File.precision(15);
          T_tilde_Mean_File.open(cstr2, ios::out);
      }
      SPRINTF (cstr, "Tij_Merged");
      if ((SU2_TYPE::Int(idx_out) >= 0)    && (SU2_TYPE::Int(idx_out) < 10))    SPRINTF (buffer, "_0000%d.dat", SU2_TYPE::Int(idx_out));
      if ((SU2_TYPE::Int(idx_out) >= 10)   && (SU2_TYPE::Int(idx_out) < 100))   SPRINTF (buffer, "_000%d.dat",  SU2_TYPE::Int(idx_out));
      if ((SU2_TYPE::Int(idx_out) >= 100)  && (SU2_TYPE::Int(idx_out) < 1000))  SPRINTF (buffer, "_00%d.dat",   SU2_TYPE::Int(idx_out));
      if ((SU2_TYPE::Int(idx_out) >= 1000) && (SU2_TYPE::Int(idx_out) < 10000)) SPRINTF (buffer, "_0%d.dat",    SU2_TYPE::Int(idx_out));
      if (SU2_TYPE::Int(idx_out) >= 10000) SPRINTF (buffer, "_%d.dat", SU2_TYPE::Int(idx_out));
      strcat (cstr, buffer);
      cout<<"Writing Merged File: "<<cstr<<endl;
      T_tilde_File.precision(15);
      T_tilde_File.open(cstr, ios::out);


      /*--- Loop through all of the collected data and write each node's values ---*/
      unsigned long Total_Index;
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        for (iSNGPts = 0; iSNGPts < Buffer_Recv_nSNGPts[iProcessor]; iSNGPts++) {

           for (iVar = 0; iVar < nVar; iVar++){

          /*--- Current index position and global index ---*/
          Total_Index  = iProcessor*Max_nSNGPts*nVar + iVar*Buffer_Recv_nSNGPts[iProcessor]  + iSNGPts;

          /*--- Write to file---*/
          T_tilde_File << scientific <<  Buffer_Recv_Tij[Total_Index]   << "\t";
         if (idx==0) T_tilde_Mean_File << scientific <<  Buffer_Recv_Tij_Mean[Total_Index]   << "\t";
           }
           T_tilde_File  << endl;
           if (idx==0)T_tilde_Mean_File  << endl;
         }
      }

      T_tilde_File.close();
      if (idx==0) T_tilde_Mean_File.close();
      delete [] Buffer_Recv_Tij;
      if (idx==0) delete [] Buffer_Recv_Tij_Mean;
       }
      delete [] Buffer_Send_Tij;
      if (idx==0) delete [] Buffer_Send_Tij_Mean;
    }

}






void SNG::Write_SNGSensitivities(){


    unsigned long iVar, iSNGPts, Max_nSNGPts, Tot_nSNGPts,nVar,Global_Index;
       ofstream dJBBN_dU_File;

    int rank, iProcessor, nProcessor;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
    unsigned long Buffer_Send_nSNGPts[1], *Buffer_Recv_nSNGPts = NULL;



    if (rank == MASTER_NODE) Buffer_Recv_nSNGPts= new unsigned long [nProcessor];

      Buffer_Send_nSNGPts[0]=nSNGPts;
#ifdef HAVE_MPI
      SU2_MPI::Gather(&Buffer_Send_nSNGPts, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nSNGPts, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD); //send the number of vertices at each process to the master
      SU2_MPI::Allreduce(&nSNGPts,&Max_nSNGPts,1,MPI_UNSIGNED_LONG,MPI_MAX,MPI_COMM_WORLD); //find the max num of vertices over all processes
      SU2_MPI::Reduce(&nSNGPts,&Tot_nSNGPts,1,MPI_UNSIGNED_LONG,MPI_SUM,MASTER_NODE,MPI_COMM_WORLD); //find the total num of vertices (panels)
#endif

      cout<<"Max nSNGPts for all processes: "<<Max_nSNGPts<<",  nSNGPTs of Process "<<rank<<": "<<nSNGPts<<", Total Num of SNG Pts: " <<Tot_nSNGPts<<endl;


      nVar = nDim+2;

      /* pack sensitivity values in each processor and send to root */
      su2double *Buffer_Send_dJBBN_dU = new su2double [Max_nSNGPts*nVar];
      unsigned long *Buffer_Send_GlobalIndex = new unsigned long [Max_nSNGPts];

      //zero send buffers
      for (int i=0; i <Max_nSNGPts*nVar; i++){
             Buffer_Send_dJBBN_dU[i]=0.0;
      }
      for (int i=0; i <Max_nSNGPts; i++){
       Buffer_Send_GlobalIndex[i]=0;
      }


      su2double *Buffer_Recv_dJBBN_dU = NULL;
      unsigned long *Buffer_Recv_GlobalIndex = NULL;


      if (rank == MASTER_NODE) {
          Buffer_Recv_dJBBN_dU = new su2double [nProcessor*Max_nSNGPts*nVar];
          Buffer_Recv_GlobalIndex = new unsigned long[nProcessor*Max_nSNGPts];
      }

      for (iVar=0; iVar<nVar; iVar++){
          for(iSNGPts=0; iSNGPts<nSNGPts; iSNGPts++){
               Buffer_Send_dJBBN_dU[iVar*nSNGPts+iSNGPts]= dJBBN_dU[iVar][iSNGPts];
            }
      }


      for (iSNGPts=0; iSNGPts<nSNGPts; iSNGPts++){
         Buffer_Send_GlobalIndex[iSNGPts] = PointID[iSNGPts];
        }

#ifdef HAVE_MPI
     SU2_MPI::Gather(Buffer_Send_dJBBN_dU, Max_nSNGPts*nVar, MPI_DOUBLE, Buffer_Recv_dJBBN_dU,  Max_nSNGPts*nVar , MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
     SU2_MPI::Gather(Buffer_Send_GlobalIndex,Max_nSNGPts, MPI_UNSIGNED_LONG, Buffer_Recv_GlobalIndex, Max_nSNGPts , MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
#endif

     /* root opens a file at each time step and write out the merged dJdU values at that time step into the file */
      if (rank == MASTER_NODE){

          char cstr [200];
          SPRINTF (cstr, "dJBBN_dU");
          cout<<"Writing Merged File: "<<cstr<<endl;
          dJBBN_dU_File.precision(15);
          dJBBN_dU_File.open(cstr, ios::out);


      /*--- Loop through all of the collected data and write each node's values ---*/
      unsigned long Total_Index;
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        for (iSNGPts = 0; iSNGPts < Buffer_Recv_nSNGPts[iProcessor]; iSNGPts++) {
            Global_Index = Buffer_Recv_GlobalIndex[iProcessor*Max_nSNGPts+iSNGPts ];
            dJBBN_dU_File  << scientific << Global_Index << "\t";
           for (iVar = 0; iVar < nVar; iVar++){

          /*--- Current index position and global index ---*/
          Total_Index  = iProcessor*Max_nSNGPts*nVar + iVar*Buffer_Recv_nSNGPts[iProcessor]  + iSNGPts;

          /*--- Write to file---*/
            dJBBN_dU_File << scientific <<  Buffer_Recv_dJBBN_dU[Total_Index]   << "\t";
           }
            dJBBN_dU_File  << endl;
         }
      }

         dJBBN_dU_File.close();
         delete [] Buffer_Recv_dJBBN_dU;
         delete [] Buffer_Recv_GlobalIndex;
       }
        delete [] Buffer_Send_dJBBN_dU;
        delete [] Buffer_Send_GlobalIndex;

}






