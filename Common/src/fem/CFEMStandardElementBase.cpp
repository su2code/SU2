/*!
 * \file CFEMStandardElementBase.cpp
 * \brief Functions for the base class CFEMStandardElementBase.
 * \author E. van der Weide
 * \version 7.0.8 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/fem/CFEMStandardElementBase.hpp"
#include "../../include/toolboxes/CSquareMatrixCM.hpp"
#include "../../include/omp_structure.hpp"

/*----------------------------------------------------------------------------------*/
/*          Constructor and destructor of CFEMStandardElementBase.                  */
/*----------------------------------------------------------------------------------*/

CFEMStandardElementBase::~CFEMStandardElementBase() {

#if defined(PRIMAL_SOLVER) && defined(HAVE_MKL)
  if( jitter ) {
    mkl_jit_destroy(jitter);
    jitter = nullptr;
  }
#endif
}

/*----------------------------------------------------------------------------------*/
/*          Public member functions of CFEMStandardElementBase.                     */
/*----------------------------------------------------------------------------------*/

unsigned short CFEMStandardElementBase::GetNIntStatic(unsigned short VTK_Type,
                                                      unsigned short orderExact) {

  /*--- Determine the number of integration points in 1D for a line. ---*/
  const unsigned short nInt1D = orderExact/2 + 1;

  /*--- Determine the element type and set the number of integration points. ---*/
  unsigned short nInt;
  switch(VTK_Type) {

    case LINE:                 // Tensor products elements for
    case QUADRILATERAL:        // which the number of integration
    case HEXAHEDRON:           // points in 1D is returned.
      nInt = nInt1D;
      break;

    case TRIANGLE:
      nInt = GetNIntTriangleStatic(orderExact);
      break;

    case TETRAHEDRON:
      nInt = GetNIntTetrahedronStatic(orderExact);
      break;

    case PRISM:
      nInt = nInt1D*GetNIntTetrahedronStatic(orderExact);
      break;

    case PYRAMID:
      nInt = nInt1D*nInt1D*nInt1D;
      break;

    default:
      nInt = 0;  // Indicates an invalid element.
      break;
  }

  return nInt;
}

unsigned short CFEMStandardElementBase::GetNDOFsStatic(unsigned short VTK_Type,
                                                       unsigned short nPoly) {
  unsigned short nDOFsEdge = nPoly + 1;
  unsigned short nDOFs;

  switch(VTK_Type) {

    case LINE:
      nDOFs = nDOFsEdge;
      break;

    case TRIANGLE:
      nDOFs = nDOFsEdge*(nDOFsEdge+1)/2;
      break;

    case QUADRILATERAL:
      nDOFs = nDOFsEdge*nDOFsEdge;
      break;

    case TETRAHEDRON:
      nDOFs = nDOFsEdge*(nDOFsEdge+1)*(nDOFsEdge+2)/6;
      break;

    case HEXAHEDRON:
      nDOFs = nDOFsEdge*nDOFsEdge*nDOFsEdge;
      break;

    case PRISM:
      nDOFs = nDOFsEdge*nDOFsEdge*(nDOFsEdge+1)/2;
      break;

    case PYRAMID:
      nDOFs = nDOFsEdge*(nDOFsEdge+1)*(2*nDOFsEdge+1)/6;
      break;

    default:
      nDOFs = 0;  // Indicates an invalid element.
      break;
  }

  return nDOFs;
}

unsigned short CFEMStandardElementBase::GetNDOFsPerSubElem(unsigned short val_VTK_Type) const {

  /*--- Distinguish between the possible element types for a linear
        sub-element and set the nDOFs accordingly. ---*/
  unsigned short nDOFsSubElem = 0;   // To avoid a compiler warning.
  switch( val_VTK_Type ) {
    case NONE:          nDOFsSubElem = 0; break;
    case LINE:          nDOFsSubElem = 2; break;
    case TRIANGLE:      nDOFsSubElem = 3; break;
    case QUADRILATERAL: nDOFsSubElem = 4; break;
    case TETRAHEDRON:   nDOFsSubElem = 4; break;
    case PYRAMID:       nDOFsSubElem = 5; break;
    case PRISM:         nDOFsSubElem = 6; break;
    case HEXAHEDRON:    nDOFsSubElem = 8; break;
    default:
      ostringstream message;
      message << "Impossible FEM sub element type, " << val_VTK_Type
              << ", encountered.";
      SU2_MPI::Error(message.str(), CURRENT_FUNCTION);
  }

  /*--- Return the number of DOFs for a subface. ---*/
  return nDOFsSubElem;
}

void CFEMStandardElementBase::MetricTermsVolume(vector<ColMajorMatrix<su2double> > &matMetricTerms,
                                                su2activevector                    &Jacobians) {

  /*--- Determine the number of items for which the metric terms
        must be computed. ---*/
  const unsigned short nItems = Jacobians.rows();

  /*--- Check the dimension of the problem. ---*/
  if(matMetricTerms.size() == 2) {

    /*--- 2D problem. Set the references for the derivatives of the Cartesian
          coordinates w.r.t. the parametric coordinates, which are currently
          stored in matMetricTerms. ---*/
    ColMajorMatrix<su2double> &dCoorDr = matMetricTerms[0];
    ColMajorMatrix<su2double> &dCoorDs = matMetricTerms[1];

    /*--- Set the references to store the derivatives of the parametric coordinates
          w.r.t. the Cartesian coordinates. ---*/
    ColMajorMatrix<su2double> &dParDx = matMetricTerms[0];
    ColMajorMatrix<su2double> &dParDy = matMetricTerms[1];

    /*--- Loop over the padded number of integration points
          to compute the 2D metric terms. ---*/
    SU2_OMP_SIMD_IF_NOT_AD
    for(unsigned short i=0; i<nItems; ++i) {

      /*--- Easier storage of the derivatives of the Cartesian coordinates
            w.r.t. the parametric coordinates. ---*/
      const su2double dxdr = dCoorDr(i,0);
      const su2double dydr = dCoorDr(i,1);

      const su2double dxds = dCoorDs(i,0);
      const su2double dyds = dCoorDs(i,1);

      /*--- Compute the Jacobian. ---*/
      Jacobians(i) = dxdr*dyds - dxds*dydr;

      /*--- Compute the derivatives of the parametric coordinates w.r.t.
            to the Cartesian coordinates. They are multiplied by the Jacobian
            to avoid a division by this Jacobian at this stage. ---*/
      dParDx(i,0) =  dyds;   // J drdx
      dParDx(i,1) = -dydr;   // J dsdx

      dParDy(i,0) = -dxds;   // J drdy
      dParDy(i,1) =  dxdr;   // J dsdy
    }
  }
  else {

    /*--- 3D problem. Set the references for the derivatives of the Cartesian
          coordinates w.r.t. the parametric coordinates, which are currently
          stored in matMetricTerms. ---*/
    ColMajorMatrix<su2double> &dCoorDr = matMetricTerms[0];
    ColMajorMatrix<su2double> &dCoorDs = matMetricTerms[1];
    ColMajorMatrix<su2double> &dCoorDt = matMetricTerms[2];

    /*--- Set the references to store the derivatives of the parametric coordinates
          w.r.t. the Cartesian coordinates. ---*/
    ColMajorMatrix<su2double> &dParDx = matMetricTerms[0];
    ColMajorMatrix<su2double> &dParDy = matMetricTerms[1];
    ColMajorMatrix<su2double> &dParDz = matMetricTerms[2];

    /*--- Loop over the padded number of integration points
          to compute the 3D metric terms. ---*/
    SU2_OMP_SIMD_IF_NOT_AD
    for(unsigned short i=0; i<nItems; ++i) {

      /*--- Easier storage of the derivatives of the Cartesian coordinates
            w.r.t. the parametric coordinates. ---*/
      const su2double dxdr = dCoorDr(i,0);
      const su2double dydr = dCoorDr(i,1);
      const su2double dzdr = dCoorDr(i,2);

      const su2double dxds = dCoorDs(i,0);
      const su2double dyds = dCoorDs(i,1);
      const su2double dzds = dCoorDs(i,2);

      const su2double dxdt = dCoorDt(i,0);
      const su2double dydt = dCoorDt(i,1);
      const su2double dzdt = dCoorDt(i,2);

      /*--- Compute the Jacobian. ---*/
      Jacobians(i) = dxdr*(dyds*dzdt - dzds*dydt)
                   - dxds*(dydr*dzdt - dzdr*dydt)
                   + dxdt*(dydr*dzds - dzdr*dyds);

      /*--- Compute the derivatives of the parametric coordinates w.r.t.
            to the Cartesian coordinates. They are multiplied by the Jacobian
            to avoid a division by this Jacobian at this stage. ---*/
      dParDx(i,0) = dyds*dzdt - dzds*dydt;  // J drdx
      dParDx(i,1) = dzdr*dydt - dydr*dzdt;  // J dsdx
      dParDx(i,2) = dydr*dzds - dzdr*dyds;  // J dtdx

      dParDy(i,0) = dzds*dxdt - dxds*dzdt;  // J drdy
      dParDy(i,1) = dxdr*dzdt - dzdr*dxdt;  // J dsdy
      dParDy(i,2) = dzdr*dxds - dxdr*dzds;  // J dtdy

      dParDz(i,0) = dxds*dydt - dyds*dxdt;  // J drdz
      dParDz(i,1) = dydr*dxdt - dxdr*dydt;  // J dsdz
      dParDz(i,2) = dxdr*dyds - dydr*dxds;  // J dtdz
    }
  }
}

void CFEMStandardElementBase::MetricTermsVolumeIntPoints(const bool                         LGLDistribution,
                                                         ColMajorMatrix<su2double>          &matCoor,
                                                         vector<ColMajorMatrix<su2double> > &matMetricTerms,
                                                         su2activevector                    &Jacobians) {

  /*--- Compute the derivatives of the coordinates in the volume integration points.
        The matrix matMetricTerms is used as storage for this data. ---*/
  DerivativesCoorIntPoints(LGLDistribution, matCoor, matMetricTerms);

  /*--- Compute the actual metric terms. ---*/
  MetricTermsVolume(matMetricTerms, Jacobians);
}

void CFEMStandardElementBase::MetricTerms2ndDerVolumeIntPoints(const bool                         LGLDistribution,
                                                               ColMajorMatrix<su2double>          &matCoor,
                                                               vector<ColMajorMatrix<su2double> > &matMetricTerms,
                                                               su2activevector                    &Jacobians,
                                                               vector<ColMajorMatrix<su2double> > &matMetricTerms2ndDer) {

  /*--- Compute the 2nd derivatives of the coordinates in the volume integration points.
        The matrix matMetricTerms2ndDer is used as storage for this data. ---*/
  Derivatives2ndCoorIntPoints(LGLDistribution, matCoor, matMetricTerms2ndDer);

  /*--- Determine the number of items for which the metric terms
        must be computed. ---*/
  const unsigned short nItems = Jacobians.rows();

  /*--- Check the dimension of the problem. ---*/
  if(matMetricTerms.size() == 2) {

    /*--- 2D problem. Set the references for the derivatives of the
          parametric coordinates w.r.t. the Cartesian coordinates. ---*/
    ColMajorMatrix<su2double> &dParDx = matMetricTerms[0];
    ColMajorMatrix<su2double> &dParDy = matMetricTerms[1];

    /*--- Set the references for the second derivatives. ---*/
    ColMajorMatrix<su2double> &d2dr2 = matMetricTerms2ndDer[0];
    ColMajorMatrix<su2double> &d2ds2 = matMetricTerms2ndDer[1];
    ColMajorMatrix<su2double> &d2drs = matMetricTerms2ndDer[2];

    /*--- Loop over the padded number of integration points
          to compute the 2D metric terms. ---*/
    SU2_OMP_SIMD_IF_NOT_AD
    for(unsigned short i=0; i<nItems; ++i) {

      /*--- Compute the inverse of the Jacobian and the square of it. ---*/
      const su2double JInv  = 1.0/Jacobians(i);
      const su2double JInv2 = JInv*JInv;

      /*--- Compute the derivatives of the parametric coordinates w.r.t.
            the Cartesian coordinates. These are the standard metric terms,
            but they must be corrected with the inverse of the Jacobian. ---*/
      const su2double drdx = JInv*dParDx(i,0);
      const su2double dsdx = JInv*dParDx(i,1);
      const su2double drdy = JInv*dParDy(i,0);
      const su2double dsdy = JInv*dParDy(i,1);

      /*--- Retrieve the first derivatives of the Cartesian coordinates
            from the metric terms. ---*/
      const su2double dxdr =  dParDy(i,1);
      const su2double dxds = -dParDy(i,0);
      const su2double dydr = -dParDx(i,1);
      const su2double dyds =  dParDx(i,0);

      /*--- Store the second derivatives of the Cartesian coordinates. ---*/
      const su2double d2xdr2 = d2dr2(i,0), d2ydr2 = d2dr2(i,1);
      const su2double d2xds2 = d2ds2(i,0), d2yds2 = d2ds2(i,1);
      const su2double d2xdrs = d2drs(i,0), d2ydrs = d2drs(i,1);

      /*--- Compute the derivatives of the Jacobian. ---*/
      const su2double dJdr = d2xdr2*dyds + d2ydrs*dxdr - d2xdrs*dydr - d2ydr2*dxds;
      const su2double dJds = d2xdrs*dyds + d2yds2*dxdr - d2xds2*dydr - d2ydrs*dxds;

      /*--- Compute the derivatives of the standard metric terms w.r.t. the
            parametric coordinates. ---*/
      const su2double ddrdx_dr =  JInv*d2ydrs - dyds*JInv2*dJdr;
      const su2double ddrdx_ds =  JInv*d2yds2 - dyds*JInv2*dJds;
      const su2double ddrdy_dr = -JInv*d2xdrs + dxds*JInv2*dJdr;
      const su2double ddrdy_ds = -JInv*d2xds2 + dxds*JInv2*dJds;

      const su2double ddsdx_dr = -JInv*d2ydr2 + dydr*JInv2*dJdr;
      const su2double ddsdx_ds = -JInv*d2ydrs + dydr*JInv2*dJds;
      const su2double ddsdy_dr =  JInv*d2xdr2 - dxdr*JInv2*dJdr;
      const su2double ddsdy_ds =  JInv*d2xdrs - dxdr*JInv2*dJds;

      /*--- Compute the metric terms needed to compute the Cartesian 2nd
            derivatives. Note that for the cross derivatives an average of
            the two possibilities is taken. ---*/
      d2dr2(i,0) = drdx*ddrdx_dr + dsdx*ddrdx_ds;
      d2dr2(i,1) = drdx*ddsdx_dr + dsdx*ddsdx_ds;

      d2ds2(i,0) = drdy*ddrdy_dr + dsdy*ddrdy_ds;
      d2ds2(i,1) = drdy*ddsdy_dr + dsdy*ddsdy_ds;

      d2drs(i,0) = 0.5*(drdx*ddrdy_dr + dsdx*ddrdy_ds + drdy*ddrdx_dr + dsdy*ddrdx_ds);
      d2drs(i,1) = 0.5*(drdx*ddsdy_dr + dsdx*ddsdy_ds + drdy*ddsdx_dr + dsdy*ddsdx_ds);
    }
  }
  else {

    /*--- 3D problem. Set the references for the derivatives of the
          parametric coordinates w.r.t. the Cartesian coordinates. ---*/
    ColMajorMatrix<su2double> &dParDx = matMetricTerms[0];
    ColMajorMatrix<su2double> &dParDy = matMetricTerms[1];
    ColMajorMatrix<su2double> &dParDz = matMetricTerms[2];

    /*--- Set the references for the second derivatives. ---*/
    ColMajorMatrix<su2double> &d2dr2 = matMetricTerms2ndDer[0];
    ColMajorMatrix<su2double> &d2ds2 = matMetricTerms2ndDer[1];
    ColMajorMatrix<su2double> &d2dt2 = matMetricTerms2ndDer[2];
    ColMajorMatrix<su2double> &d2drs = matMetricTerms2ndDer[3];
    ColMajorMatrix<su2double> &d2drt = matMetricTerms2ndDer[4];
    ColMajorMatrix<su2double> &d2dst = matMetricTerms2ndDer[5];

    /*--- Loop over the padded number of integration points
          to compute the 3D metric terms. ---*/
    SU2_OMP_SIMD_IF_NOT_AD
    for(unsigned short i=0; i<nItems; ++i) {

      /*--- Compute the inverse of the Jacobian and the square of it. ---*/
      const su2double JInv  = 1.0/Jacobians(i);
      const su2double JInv2 = JInv*JInv;

      /*--- Compute the derivatives of the parametric coordinates w.r.t.
            the Cartesian coordinates. These are the standard metric terms,
            but they must be corrected with the inverse of the Jacobian. ---*/
      const su2double drdx = JInv*dParDx(i,0);
      const su2double dsdx = JInv*dParDx(i,1);
      const su2double dtdx = JInv*dParDx(i,2);

      const su2double drdy = JInv*dParDy(i,0);
      const su2double dsdy = JInv*dParDy(i,1);
      const su2double dtdy = JInv*dParDy(i,2);

      const su2double drdz = JInv*dParDz(i,0);
      const su2double dsdz = JInv*dParDz(i,1);
      const su2double dtdz = JInv*dParDz(i,2);

      /*--- Retrieve the first derivatives of the Cartesian coordinates
            from the metric terms. ---*/
      const su2double dxdr = (dsdy*dtdz - dsdz*dtdy)*Jacobians(i);
      const su2double dxds = (drdz*dtdy - drdy*dtdz)*Jacobians(i);
      const su2double dxdt = (drdy*dsdz - drdz*dsdy)*Jacobians(i);

      const su2double dydr = (dsdz*dtdx - dsdx*dtdz)*Jacobians(i);
      const su2double dyds = (dtdz*drdx - dtdz*dtdx)*Jacobians(i);
      const su2double dydt = (drdz*dsdx - drdx*dsdz)*Jacobians(i);

      const su2double dzdr = (dsdx*dtdy - dsdy*dtdx)*Jacobians(i);
      const su2double dzds = (drdy*dtdx - drdx*dtdy)*Jacobians(i);
      const su2double dzdt = (drdz*dsdy - drdy*dsdx)*Jacobians(i);

       /*--- Store the second derivatives of the Cartesian coordinates. ---*/
      const su2double d2xdr2 = d2dr2(i,0), d2ydr2 = d2dr2(i,1), d2zdr2 = d2dr2(i,2);
      const su2double d2xds2 = d2ds2(i,0), d2yds2 = d2ds2(i,1), d2zds2 = d2ds2(i,2);
      const su2double d2xdt2 = d2dt2(i,0), d2ydt2 = d2dt2(i,1), d2zdt2 = d2dt2(i,2);

      const su2double d2xdrs = d2drs(i,0), d2ydrs = d2drs(i,1), d2zdrs = d2drs(i,2);
      const su2double d2xdrt = d2drt(i,0), d2ydrt = d2drt(i,1), d2zdrt = d2drt(i,2);
      const su2double d2xdst = d2dst(i,0), d2ydst = d2dst(i,1), d2zdst = d2dst(i,2);

      /*--- Compute the derivatives of the Jacobian. ---*/
      const su2double dJdr = d2xdr2*(dyds*dzdt - dzds*dydt)
                           + dxdr  *(d2ydrs*dzdt + d2zdrt*dyds - d2zdrs*dydt - d2ydrt*dzds)
                           - d2xdrs*(dydr*dzdt - dzdr*dydt)
                           - dxds  *(d2ydr2*dzdt + d2zdrt*dydr - d2zdr2*dydt - d2ydrt*dzdr)
                           + d2xdrt*(dydr*dzds - dzdr*dyds)
                           + dxdt  *(d2ydr2*dzds + d2zdrs*dydr - d2zdr2*dyds - d2ydrs*dzdr);

      const su2double dJds = d2xdrs*(dyds*dzdt - dzds*dydt)
                           + dxdr  *(d2yds2*dzdt + d2zdst*dyds - d2zds2*dydt - d2ydst*dzds)
                           - d2xds2*(dydr*dzdt - dzdr*dydt)
                           - dxds  *(d2ydrs*dzdt + d2zdst*dydr - d2zdrs*dydt - d2ydst*dzdr)
                           + d2xdst*(dydr*dzds - dzdr*dyds)
                           + dxdt  *(d2ydrs*dzds + d2zds2*dydr - d2zdrs*dyds - d2yds2*dzdr);

      const su2double dJdt = d2xdrt*(dyds*dzdt - dzds*dydt)
                           + dxdr  *(d2ydst*dzdt + d2zdt2*dyds - d2zdst*dydt - d2ydt2*dzds)
                           - d2xdst*(dydr*dzdt - dzdr*dydt)
                           - dxds  *(d2ydrt*dzdt + d2zdt2*dydr - d2zdrt*dydt - d2ydt2*dzdr)
                           + d2xdt2*(dydr*dzds - dzdr*dyds)
                           + dxdt  *(d2ydrt*dzds + d2zdst*dydr - d2zdrt*dyds - d2ydst*dzdr);

      /*--- Compute the derivatives of the standard metric terms w.r.t. the
            parametric coordinates. ---*/
      const su2double ddrdx_dr =  JInv*(d2ydrs*dzdt + d2zdrt*dyds - d2zdrs*dydt - d2ydrt*dzds)
                               - JInv2*(dyds*dzdt - dzds*dydt)*dJdr;
      const su2double ddrdx_ds =  JInv*(d2yds2*dzdt + d2zdst*dyds - d2zds2*dydt - d2ydst*dzds)
                               - JInv2*(dyds*dzdt - dzds*dydt)*dJds;
      const su2double ddrdx_dt =  JInv*(d2ydst*dzdt + d2zdt2*dyds - d2zdst*dydt - d2ydt2*dzds)
                               - JInv2*(dyds*dzdt - dzds*dydt)*dJdt;

      const su2double ddrdy_dr =  JInv*(d2zdrs*dxdt + d2xdrt*dzds - d2xdrs*dzdt - d2zdrt*dxds)
                               - JInv2*(dzds*dxdt - dxds*dzdt)*dJdr;
      const su2double ddrdy_ds =  JInv*(d2zds2*dxdt + d2xdst*dzds - d2xds2*dzdt - d2zdst*dxds)
                               - JInv2*(dzds*dxdt - dxds*dzdt)*dJds;
      const su2double ddrdy_dt =  JInv*(d2zdst*dxdt + d2xdt2*dzds - d2xdst*dzdt - d2zdt2*dxds)
                               - JInv2*(dzds*dxdt - dxds*dzdt)*dJdt;

      const su2double ddrdz_dr =  JInv*(d2xdrs*dydt + d2ydrt*dxds - d2ydrs*dxdt - d2xdrt*dyds)
                               - JInv2*(dxds*dydt - dyds*dxdt)*dJdr;
      const su2double ddrdz_ds =  JInv*(d2xds2*dydt + d2ydst*dxds - d2yds2*dxdt - d2xdst*dyds)
                               - JInv2*(dxds*dydt - dyds*dxdt)*dJds;
      const su2double ddrdz_dt =  JInv*(d2xdst*dydt + d2ydt2*dxds - d2ydst*dxdt - d2xdt2*dyds)
                               - JInv2*(dxds*dydt - dyds*dxdt)*dJdt;

      const su2double ddsdx_dr =  JInv*(d2zdr2*dydt + d2ydrt*dzdr - d2ydr2*dzdt - d2zdrt*dydr)
                               - JInv2*(dzdr*dydt - dydr*dzdt)*dJdr;
      const su2double ddsdx_ds =  JInv*(d2zdrs*dydt + d2ydst*dzdr - d2ydrs*dzdt - d2zdst*dydr)
                               - JInv2*(dzdr*dydt - dydr*dzdt)*dJds;
      const su2double ddsdx_dt =  JInv*(d2zdrt*dydt + d2ydt2*dzdr - d2ydrt*dzdt - d2zdt2*dydr)
                               - JInv2*(dzdr*dydt - dydr*dzdt)*dJdt;

      const su2double ddsdy_dr =  JInv*(d2xdr2*dzdt + d2zdrt*dxdr - d2zdr2*dxdt - d2xdrt*dzdr)
                               - JInv2*(dxdr*dzdt - dzdr*dxdt)*dJdr;
      const su2double ddsdy_ds =  JInv*(d2xdrs*dzdt + d2zdst*dxdr - d2zdrs*dxdt - d2xdst*dzdr)
                               - JInv2*(dxdr*dzdt - dzdr*dxdt)*dJds;
      const su2double ddsdy_dt =  JInv*(d2xdrt*dzdt + d2zdt2*dxdr - d2zdrt*dxdt - d2xdt2*dzdr)
                               - JInv2*(dxdr*dzdt - dzdr*dxdt)*dJdt;

      const su2double ddsdz_dr =  JInv*(d2ydr2*dxdt + d2xdrt*dydr - d2xdr2*dydt - d2ydrt*dxdr)
                               - JInv2*(dydr*dxdt - dxdr*dydt)*dJdr;
      const su2double ddsdz_ds =  JInv*(d2ydrs*dxdt + d2xdst*dydr - d2xdrs*dydt - d2ydst*dxdr)
                               - JInv2*(dydr*dxdt - dxdr*dydt)*dJds;
      const su2double ddsdz_dt =  JInv*(d2ydrt*dxdt + d2xdt2*dydr - d2xdrt*dydt - d2ydt2*dxdr)
                               - JInv2*(dydr*dxdt - dxdr*dydt)*dJdt;

      const su2double ddtdx_dr =  JInv*(d2ydr2*dzds + d2zdrs*dydr - d2zdr2*dyds - d2ydrs*dzdr)
                               - JInv2*(dydr*dzds - dzdr*dyds)*dJdr;
      const su2double ddtdx_ds =  JInv*(d2ydrs*dzds + d2zds2*dydr - d2zdrs*dyds - d2yds2*dzdr)
                               - JInv2*(dydr*dzds - dzdr*dyds)*dJds;
      const su2double ddtdx_dt =  JInv*(d2ydrt*dzds + d2zdst*dydr - d2zdrt*dyds - d2ydst*dzdr)
                               - JInv2*(dydr*dzds - dzdr*dyds)*dJdt;

      const su2double ddtdy_dr =  JInv*(d2zdr2*dxds + d2xdrs*dzdr - d2xdr2*dzds - d2zdrs*dxdr)
                               - JInv2*(dzdr*dxds - dxdr*dzds)*dJdr;
      const su2double ddtdy_ds =  JInv*(d2zdrs*dxds + d2xds2*dzdr - d2xdrs*dzds - d2zds2*dxdr)
                               - JInv2*(dzdr*dxds - dxdr*dzds)*dJds;
      const su2double ddtdy_dt =  JInv*(d2zdrt*dxds + d2xdst*dzdr - d2xdrt*dzds - d2zdst*dxdr)
                               - JInv2*(dzdr*dxds - dxdr*dzds)*dJdt;

      const su2double ddtdz_dr =  JInv*(d2xdr2*dyds + d2ydrs*dxdr - d2ydr2*dxds - d2xdrs*dydr)
                               - JInv2*(dxdr*dyds - dydr*dxds)*dJdr;
      const su2double ddtdz_ds =  JInv*(d2xdrs*dyds + d2yds2*dxdr - d2ydrs*dxds - d2xds2*dydr)
                               - JInv2*(dxdr*dyds - dydr*dxds)*dJds;
      const su2double ddtdz_dt =  JInv*(d2xdrt*dyds + d2ydst*dxdr - d2ydrt*dxds - d2xdst*dydr)
                               - JInv2*(dxdr*dyds - dydr*dxds)*dJdt;

      /*--- Compute the metric terms needed to compute the Cartesian 2nd
            derivatives. Note that for the cross derivatives an average of
            the two possibilities is taken. ---*/
      d2dr2(i,0) = drdx*ddrdx_dr + dsdx*ddrdx_ds + dtdx*ddrdx_dt;
      d2dr2(i,1) = drdx*ddsdx_dr + dsdx*ddsdx_ds + dtdx*ddsdx_dt;
      d2dr2(i,2) = drdx*ddtdx_dr + dsdx*ddtdx_ds + dtdx*ddtdx_dt;

      d2ds2(i,0) = drdy*ddrdy_dr + dsdy*ddrdy_ds + dtdy*ddrdy_dt;
      d2ds2(i,1) = drdy*ddsdy_dr + dsdy*ddsdy_ds + dtdy*ddsdy_dt;
      d2ds2(i,2) = drdy*ddtdy_dr + dsdy*ddtdy_ds + dtdy*ddtdy_dt;

      d2dt2(i,0) = drdz*ddrdz_dr + dsdz*ddrdz_ds + dtdz*ddrdz_dt;
      d2dt2(i,1) = drdz*ddsdz_dr + dsdz*ddsdz_ds + dtdz*ddsdz_dt;
      d2dt2(i,2) = drdz*ddtdz_dr + dsdz*ddtdz_ds + dtdz*ddtdz_dt;

      d2drs(i,0) = 0.5*(drdx*ddrdy_dr + dsdx*ddrdy_ds + dtdx*ddrdy_dt
                 +      drdy*ddrdx_dr + dsdy*ddrdx_ds + dtdy*ddrdx_dt);
      d2drs(i,1) = 0.5*(drdx*ddsdy_dr + dsdx*ddsdy_ds + dtdx*ddsdy_dt
                 +      drdy*ddsdx_dr + dsdy*ddsdx_ds + dtdy*ddsdx_dt);
      d2drs(i,2) = 0.5*(drdx*ddtdy_dr + dsdx*ddtdy_ds + dtdx*ddtdy_dt
                 +      drdy*ddtdx_dr + dsdy*ddtdx_ds + dtdy*ddtdx_dt);

      d2drt(i,0) = 0.5*(drdx*ddrdz_dr + dsdx*ddrdz_ds + dtdx*ddrdz_dt
                 +      drdz*ddrdx_dr + dsdz*ddrdx_ds + dtdz*ddrdx_dt);
      d2drt(i,1) = 0.5*(drdx*ddsdz_dr + dsdx*ddsdz_ds + dtdx*ddsdz_dt
                 +      drdz*ddsdx_dr + dsdz*ddsdx_ds + dtdz*ddsdx_dt);
      d2drt(i,2) = 0.5*(drdx*ddtdz_dr + dsdx*ddtdz_ds + dtdx*ddtdz_dt
                 +      drdz*ddtdx_dr + dsdz*ddtdx_ds + dtdz*ddtdx_dt);

      d2dst(i,0) = 0.5*(drdy*ddrdz_dr + dsdy*ddrdz_ds + dtdy*ddrdz_dt
                 +      drdz*ddrdy_dr + dsdz*ddrdy_ds + dtdz*ddrdy_dt);
      d2dst(i,1) = 0.5*(drdy*ddsdz_dr + dsdy*ddsdz_ds + dtdy*ddsdz_dt
                 +      drdz*ddsdy_dr + dsdz*ddsdy_ds + dtdz*ddsdy_dt);
      d2dst(i,2) = 0.5*(drdy*ddtdz_dr + dsdy*ddtdz_ds + dtdy*ddtdz_dt
                 +      drdz*ddtdy_dr + dsdz*ddtdy_ds + dtdz*ddtdy_dt);
    }
  }
}

void CFEMStandardElementBase::MetricTermsSolDOFs(ColMajorMatrix<su2double>          &matCoor,
                                                 vector<ColMajorMatrix<su2double> > &matMetricTerms,
                                                 su2activevector                    &Jacobians) {

  /*--- Compute the derivatives of the coordinates in the solution DOFs.
        The matrix matMetricTerms is used as storage for this data. ---*/
  DerivativesCoorSolDOFs(matCoor, matMetricTerms);

  /*--- Compute the actual metric terms. ---*/
  MetricTermsVolume(matMetricTerms, Jacobians);
}

void CFEMStandardElementBase::MinMaxFaceJacobians(const bool                         LGLDistribution,
                                                  ColMajorMatrix<su2double>          &matCoor,
                                                  vector<ColMajorMatrix<su2double> > &matDerCoor,
                                                  ColMajorMatrix<su2double>          &unitNormals,
                                                  su2activevector                    &Jacobians,
                                                  su2double                          &jacMin,
                                                  su2double                          &jacMax,
                                                  su2double                          &cosAngleMin) {

  /*--- Compute the unit face normals, which will be stored in matCoor on return. ---*/
  UnitFaceNormals(LGLDistribution, matCoor, matDerCoor, unitNormals, Jacobians);

  /*--- Loop over the integration points and determine the minimum
        and maximum value of the Jacobian. ---*/
  jacMin = jacMax = Jacobians(0);
  for(unsigned short i=0; i<nIntegration; ++i) {
    jacMin = min(jacMin, Jacobians(i));
    jacMax = max(jacMax, Jacobians(i));
  }

  /*--- Initialize the cosine of the minimum angle to 1.0, i.e. 90 degrees. ---*/
  cosAngleMin = 1.0;

  /*--- Check the dimension of the problem. ---*/
  if(matDerCoor.size() == 1) {

    /*--- Face of a 2D element. Double loop over the integration points to determine
          the minimum value of the cosine of the angle between the normals. ---*/
    for(unsigned short i=0; i<nIntegration; ++i) {
      for( unsigned short j=i+1; j<nIntegration; ++j) {
        const su2double dot = unitNormals(i,0)*unitNormals(j,0)
                            + unitNormals(i,1)*unitNormals(j,1);
        cosAngleMin = min(cosAngleMin, dot);
      }
    }
  }
  else {

    /*--- Face of a 3D element. Double loop over the integration points to determine
          the minimum value of the cosine of the angle between the normals. ---*/
    for(unsigned short i=0; i<nIntegration; ++i) {
      for( unsigned short j=i+1; j<nIntegration; ++j) {
        const su2double dot = unitNormals(i,0)*unitNormals(j,0)
                            + unitNormals(i,1)*unitNormals(j,1)
                            + unitNormals(i,2)*unitNormals(j,2);
        cosAngleMin = min(cosAngleMin, dot);
      }
    }
  }
}

void CFEMStandardElementBase::MinMaxJacobians(const bool                         LGLDistribution,
                                              ColMajorMatrix<su2double>          &matCoor,
                                              vector<ColMajorMatrix<su2double> > &matMetricTerms,
                                              su2activevector                    &Jacobians,
                                              su2double                          &jacMin,
                                              su2double                          &jacMax) {

  /*--- Compute the metric terms in the volume integration points. ---*/
  MetricTermsVolumeIntPoints(LGLDistribution, matCoor, matMetricTerms, Jacobians);

  /*--- Loop over the integration points and determine the minimum
        and maximum value of the Jacobian. ---*/
  jacMin = jacMax = Jacobians(0);
  for(unsigned short i=0; i<nIntegration; ++i) {
    jacMin = min(jacMin, Jacobians(i));
    jacMax = max(jacMax, Jacobians(i));
  }
}

void CFEMStandardElementBase::UnitFaceNormals(const bool                         LGLDistribution,
                                              ColMajorMatrix<su2double>          &matCoor,
                                              vector<ColMajorMatrix<su2double> > &matDerCoor,
                                              ColMajorMatrix<su2double>          &unitNormals,
                                              su2activevector                    &Jacobians) {

  /*--- Compute the derivatives of the coordinates in the surface integration points. ---*/
  DerivativesCoorIntPoints(LGLDistribution, matCoor, matDerCoor);

  /*--- Check the dimension of the problem. ---*/
  if(matDerCoor.size() == 1) {

    /*--- Face of a 2D element. Set the reference for the derivative of the
          Cartesian coordinates w.r.t. the parametric coordinate. ---*/
    const ColMajorMatrix<su2double> &dCoorDr = matDerCoor[0];

    /*--- Loop over the padded number of integration points
          to compute the normals. ---*/
    SU2_OMP_SIMD_IF_NOT_AD
    for(unsigned short i=0; i<nIntegrationPad; ++i) {

      /*--- Abbreviate dxdr and dydr for readability. ---*/
      const su2double dxdr = dCoorDr(i,0);
      const su2double dydr = dCoorDr(i,1);

      /*--- Determine the Jacobian and its inverse. There is a clipping to
            avoid a division by zero, but this should not be active. ---*/
      Jacobians[i]           = sqrt(dxdr*dxdr + dydr*dydr);
      const su2double invJac = Jacobians[i] < su2double(1.e-35) ? su2double(1.e+35) : 1.0/Jacobians[i];

      /*--- Compute and store the outward pointing unit normal vector. ---*/
      unitNormals(i,0) =  dydr*invJac;
      unitNormals(i,1) = -dxdr*invJac;
    }

  }
  else {

    /*--- Face of a 3D element. Set the references for the derivatives of the Cartesian
          coordinates w.r.t. the parametric coordinates. ---*/
    const ColMajorMatrix<su2double> &dCoorDr = matDerCoor[0];
    const ColMajorMatrix<su2double> &dCoorDs = matDerCoor[1];

    /*--- Loop over the padded number of integration points
          to compute the normals. ---*/
    SU2_OMP_SIMD_IF_NOT_AD
    for(unsigned short i=0; i<nIntegrationPad; ++i) {

      /*--- Easier storage of the derivatives of the Cartesian coordinates. ---*/
      const su2double dxdr = dCoorDr(i,0);
      const su2double dydr = dCoorDr(i,1);
      const su2double dzdr = dCoorDr(i,2);

      const su2double dxds = dCoorDs(i,0);
      const su2double dyds = dCoorDs(i,1);
      const su2double dzds = dCoorDs(i,2);

      /*--- Compute the vector product dxdr X dxds, where x is the coordinate
            vector (x,y,z). Compute the length of this vector, which is the
            Jacobian as well as the inverse.  Make sure that a division by zero
            is avoided, although this is most likely never active. ---*/
      const su2double nx = dydr*dzds - dyds*dzdr;
      const su2double ny = dxds*dzdr - dxdr*dzds;
      const su2double nz = dxdr*dyds - dxds*dydr;

      Jacobians[i]           = sqrt(nx*nx + ny*ny + nz*nz);
      const su2double invJac = Jacobians[i] < su2double(1.e-35) ? su2double(1.e+35) : 1.0/Jacobians[i];

      /*--- Store the unit normal. Note that the current direction is inward
            pointing. However, the outward pointing normal is needed, so the
            direction  must be reversed. ---*/
      unitNormals(i,0) = -nx*invJac;
      unitNormals(i,1) = -ny*invJac;
      unitNormals(i,2) = -nz*invJac;
    }
  }
}

/*----------------------------------------------------------------------------------*/
/*          Protected member functions of CFEMStandardElementBase.                  */
/*----------------------------------------------------------------------------------*/

void CFEMStandardElementBase::ChangeDirectionQuadConn(vector<unsigned short> &connQuad,
                                                      const unsigned short   vert0,
                                                      const unsigned short   vert1,
                                                      const unsigned short   vert2,
                                                      const unsigned short   vert3) {

  /*--- Determine the indices of the 4 corner vertices of the quad. ---*/
  const unsigned short ind0 = 0;
  const unsigned short ind1 = nPoly;
  const unsigned short ind2 = (nPoly+1)*(nPoly+1) -1;
  const unsigned short ind3 = ind2 - nPoly;

  /*--- There exists a linear mapping from the indices of the numbering used in the
        connectivity of this face to the indices of the target numbering. This
        mapping is of the form ii = a + b*i + c*j and jj = d + e*i + f*j, where
        ii,jj are the indices of the target numbering and i,j the indices of the
        numbering used for this face. The values of the coefficients a,b,c,d,e,f
        depend on how the corner points coincide with each other. This is
        determined below. The bool verticesDontMatch is there to check if vertices
        do not match. This should not happen, but it is checked for security. ---*/

  signed short a=0, b=0, c=0, d=0, e=0, f=0;  // Initialization to avoid a compiler warning.
  bool verticesDontMatch = false;

  if(vert0 == connQuad[ind0]) {

    /*--- Vert0 coincides with the first vertex of the face connectivity.
          Set the coefficients a and d accordingly.   ---*/
    a = d = 0;
    if(vert2 != connQuad[ind2]) verticesDontMatch = true;

    /*--- Check the situation for the neighboring vertices. ---*/
    if(vert1 == connQuad[ind1]) {

      /*--- The vertex numbering is the same for both faces. ---*/
      if(vert3 != connQuad[ind3]) verticesDontMatch = true;

      b = f = 1; c = e = 0;
    }
    else if(vert1 == connQuad[ind3]) {

      /*--- The i and j numbering are swapped. ---*/
      if(vert3 != connQuad[ind1]) verticesDontMatch = true;

      b = f = 0; c = e = 1;
    }
    else {
      verticesDontMatch = true;
    }
  }
  else if(vert1 == connQuad[ind0]) {

    /*--- Vert1 coincides with the first vertex of the face connectivity.
          Set the coefficients a and d accordingly.  ---*/
    a = nPoly; d = 0;
    if(vert3 != connQuad[ind2]) verticesDontMatch = true;

    /*--- Check the situation for the neighboring vertices. ---*/
    if(vert0 == connQuad[ind1]) {

      /*--- The i-direction is negated while the j-direction coincides. ---*/
      if(vert2 != connQuad[ind3]) verticesDontMatch = true;

      b = -1; f = 1; c = e = 0;
    }
    else if(vert0 == connQuad[ind3]) {

      /*--- The j-direction of the current face corresponds with the negative
            i-direction of the target, while the i-direction coincides with
            the j-direction of the target.     ---*/
      if(vert2 != connQuad[ind1]) verticesDontMatch = true;

      b = f = 0; c = -1; e = 1;
    }
    else {
      verticesDontMatch = true;
    }
  }
  else if(vert2 == connQuad[ind0]) {

    /*--- Vert2 coincides with the first vertex of the face connectivity.
          Set the coefficients a and d accordingly.  ---*/
    a = d = nPoly;
    if(vert0 != connQuad[ind2]) verticesDontMatch = true;

    /*--- Check the situation for the neighboring vertices. ---*/
    if(vert1 == connQuad[ind3]) {

      /*--- Both the i and j-direction are negated. ---*/
      if(vert3 != connQuad[ind1]) verticesDontMatch = true;

      b = f = -1; c = e = 0;
    }
    else if(vert1 == connQuad[ind1]) {

      /*--- The i and j-direction are negated and swapped. ---*/
      if(vert3 != connQuad[ind3]) verticesDontMatch = true;

      b = f = 0; c = e = -1;
    }
    else {
      verticesDontMatch = true;
    }
  }
  else if(vert3 == connQuad[ind0]) {

    /*--- Vert3 coincides with the first vertex of the face connectivity.
          Set the coefficients a and d accordingly.  ---*/
    a = 0; d = nPoly;
    if(vert1 != connQuad[ind2]) verticesDontMatch = true;

    /*--- Check the situation for the neighboring vertices. ---*/
    if(vert0 == connQuad[ind3]) {

      /*--- The i-direction coincides while the j-direction is negated. ---*/
      if(vert2 != connQuad[ind1]) verticesDontMatch = true;

      b = 1; f = -1; c = e = 0;
    }
    else if(vert0 == connQuad[ind1]) {

      /*--- The j-direction of the current face corresponds with the i-direction
            of the target, while the i-direction coincides with the negative
            j-direction of the target.    ---*/
      if(vert2 != connQuad[ind3]) verticesDontMatch = true;

      b = f = 0; c = 1; e = -1;
    }
    else {
      verticesDontMatch = true;
    }
  }
  else {
    verticesDontMatch = true;
  }

  /*--- If non-matching vertices have been found, terminate with an error message. ---*/
  if( verticesDontMatch )
    SU2_MPI::Error("Corner vertices do not match. This should not happen.",
                   CURRENT_FUNCTION);

  /*--- Copy the connectivity, such that things works out correctly when carrying
        out the renumbering.      ---*/
  vector<unsigned short> connQuadOr = connQuad;

  /*--- Loop over the vertices of the original face to copy the connectivity data. ---*/
  unsigned short ind = 0;
  for(unsigned short j=0; j<=nPoly; ++j) {
    for(unsigned short i=0; i<=nPoly; ++i, ++ind) {

      /*--- Determine the ii and jj indices of the target, convert it to
            a 1D index and shore the modified index in connQuad. ---*/
      const unsigned short ii = a + i*b + j*c;
      const unsigned short jj = d + i*e + j*f;

      const unsigned short iind = jj*(nPoly+1) + ii;

      connQuad[iind] = connQuadOr[ind];
    }
  }
}

void CFEMStandardElementBase::ChangeDirectionTriangleConn(vector<unsigned short> &connTriangle,
                                                          const unsigned short   vert0,
                                                          const unsigned short   vert1,
                                                          const unsigned short   vert2) {

  /*--- Determine the indices of the 3 corner vertices of the triangle. ---*/
  const unsigned short ind0 = 0;
  const unsigned short ind1 = nPoly;
  const unsigned short ind2 = (nPoly+1)*(nPoly+2)/2 -1;

  /*--- There exists a linear mapping from the indices of the numbering used in the
        connectivity of this face to the indices of the target numbering. This
        mapping is of the form ii = a + b*i + c*j and jj = d + e*i + f*j, where
        ii,jj are the indices of the target numbering and i,j the indices of the
        numbering used for this face. The values of the coefficients a,b,c,d,e,f
        depend on how the corner points coincide with each other. This is
        determined below. The bool verticesDontMatch is there to check if vertices
        do not match. This should not happen, but it is checked for security. ---*/
  signed short a=0, b=0, c=0, d=0, e=0, f=0;
  bool verticesDontMatch = false;

  if(vert0 == connTriangle[ind0]) {

    /*--- Vert0 coincides with the first vertex of the face connectivity.
          Check the situation for the neighboring vertices.  ---*/
    if(vert1 == connTriangle[ind1]) {

      /*--- The vertex numbering is the same for both faces. ---*/
      if(vert2 != connTriangle[ind2]) verticesDontMatch = true;

      a = 0; b = 1; c = 0; d = 0; e = 0; f = 1;
    }
    else if(vert1 == connTriangle[ind2]) {

      /*--- The ii-index corresponds to the j-index and the
            jj-index corresponds to i-index.    ---*/
      if(vert2 != connTriangle[ind1]) verticesDontMatch = true;

      a = 0; b = 0; c = 1; d = 0; e = 1; f = 0;
    }
    else {
      verticesDontMatch = true;
    }
  }
  else if(vert0 == connTriangle[ind1]) {

    /*--- Vert0 coincides with the second vertex of the face connectivity.
          Check the situation for the neighboring vertices.    ---*/
    if(vert1 == connTriangle[ind0]) {

      /*--- The ii-index corresponds to a combination of the i and j index
            and the jj-index corresponds to the j-index.   ---*/
      if(vert2 != connTriangle[ind2]) verticesDontMatch = true;

      a = nPoly; b = -1; c = -1; d = 0; e = 0; f = 1;
    }
    else if(vert1 == connTriangle[ind2]) {

      /*--- The jj-index corresponds to a combination of the i and j index
            and the ii-index corresponds to the j-index.  ---*/
      if(vert2 != connTriangle[ind0]) verticesDontMatch = true;

      a = 0; b = 0; c = 1; d = nPoly; e = -1; f = -1;
    }
    else {
      verticesDontMatch = true;
    }
  }
  else if(vert0 == connTriangle[ind2]) {

    /*--- Vert0 coincides with the third vertex of the face connectivity.
          Check the situation for the neighboring vertices.    ---*/
    if(vert1 == connTriangle[ind0]) {

      /*--- The ii-index corresponds to a combination of the i and j index
            and the jj-index corresponds with the i-index.  ---*/
      if(vert2 != connTriangle[ind1]) verticesDontMatch = true;

      a = nPoly; b = -1; c = -1; d = 0; e = 1; f = 0;
    }
    else if(vert1 == connTriangle[ind1]) {

      /*--- The jj-index corresponds to a combination of the i and j index
            and the ii-index corresponds with the i-index.   ---*/
      if(vert2 != connTriangle[ind0]) verticesDontMatch = true;

      a = 0; b = 1; c = 0; d = nPoly; e = -1; f = -1;
    }
    else {
      verticesDontMatch = true;
    }
  }
  else {
    verticesDontMatch = true;
  }

  /*--- If non-matching vertices have been found, terminate with an error message. ---*/
  if( verticesDontMatch )
    SU2_MPI::Error("Corner vertices do not match. This should not happen.",
                   CURRENT_FUNCTION);

  /*--- Copy the connectivity, such that things works out correctly when carrying
        out the renumbering.  ---*/
  vector<unsigned short> connTriangleOr = connTriangle;

  /*--- Loop over the vertices of the original face to copy the connectivity data. ---*/
  unsigned short ind = 0;
  for(unsigned short j=0; j<=nPoly; ++j) {
    for(unsigned short i=0; i<=(nPoly-j); ++i, ++ind) {

      /*--- Determine the ii and jj indices of the target, convert it to
            a 1D index and shore the modified index in connTriangle. ---*/
      const unsigned short ii = a + i*b + j*c;
      const unsigned short jj = d + i*e + j*f;

      const unsigned short iind = jj*(nPoly+1) + ii - jj*(jj-1)/2;

      connTriangle[iind] = connTriangleOr[ind];
    }
  }
}

void CFEMStandardElementBase::Location1DGridDOFsEquidistant(vector<passivedouble> &r) {

  /*--- Allocate the memory and set the location of the DOFs using
        equidistant spacing. ---*/
  r.resize(nPoly+1);
  const passivedouble dh = 2.0/nPoly;

  for(unsigned short i=0; i<=nPoly; ++i)
    r[i] = -1.0 + i*dh;
}

void CFEMStandardElementBase::Location1DGridDOFsLGL(vector<passivedouble> &r) {

  /*--- Allocate the memory. ---*/
  const unsigned short nPoints = nPoly+1;
  r.resize(nPoints);

  /*--- The distribution of points is symmetric. Hence only half the number
        of points must be computed. The first and last are at the end of
        the interval and for and add number the mid point is zero. ---*/
  const unsigned short nn = nPoints/2;

  r[0]     = -1.0;
  r[nPoly] =  1.0;

  if(2*nn < nPoints) r[nn] = 0.0;

  /*--- Constants used in the initial guess of the roots in the loop below. ---*/
  const passivedouble t1 = 1.0 - 3.0*(nPoly-1)/(8.0*nPoly*nPoly*nPoly);
  const passivedouble t2 = SU2_TYPE::GetValue(PI_NUMBER)/(4.0*nPoly + 1.0);

  /*--- The remaing points must be computed. These are the roots of P'_{n-1}(x),
        P_n is the classic Legendre polynomial of order n. Loop over roots to
        be computed. ---*/
  unsigned short ii = nPoly-1;
  for(unsigned short i=1; i<nn; ++i, --ii) {

    /*--- Initial guess of this root. ---*/
    passivedouble x = t1*cos(t2*(4*i+1));

    /*--- Determine the Legendre Polynomials P_{n-2} and P_{n-1} and
          the value f = P'_{n-1}(x). ---*/
    passivedouble Pnm1 = Legendre(nPoly,   x);
    passivedouble Pnm2 = Legendre(nPoly-1, x);
    passivedouble f    = nPoly*(Pnm2 - x*Pnm1)/(1.0-x*x);

    /*--- Solve the root using Halley's method.
          Loop until machine precision has been reached. ---*/
    for(;;) {

      /*--- Determine the value of the first and second derivative of f. ---*/
      const passivedouble df  = (2.0*x*f - nPoints*nPoly*Pnm1)/(1.0-x*x);
      const passivedouble d2f = (2.0*x*df - (nPoints*nPoly-2)*f)/(1.0-x*x);

      /*--- Compute the new value of the root. ---*/
      x = x - 2.0*f*df/(2.0*df*df - f*d2f);

      /*--- Determine the new value of the Legendre polynomials and
            compute the new value of f. Store the old value. ---*/
      const passivedouble fOld = f;

      Pnm1 = Legendre(nPoly,   x);
      Pnm2 = Legendre(nPoly-1, x);
      f    = nPoly*(Pnm2 - x*Pnm1)/(1.0-x*x);

      /*--- Convergence criterion. ---*/
      if(fabs(fOld) <= fabs(f)) break; 
    }

    /*--- Store the value as well as the symmetric equivalent. ---*/
    r[ii] =  x;
    r[i]  = -x;
  }
}

void CFEMStandardElementBase::LocationTriangleGridDOFsEquidistant(vector<passivedouble> &r,
                                                                  vector<passivedouble> &s) {

  /*--- Determine the number of DOFs of the triangle. As this function is also
        called for prisms, the member variable nDOFs cannot be used for this.
        Allocate the memory for r and s afterwards. ---*/
  const unsigned short nD = (nPoly+1)*(nPoly+2)/2;
  r.resize(nD);
  s.resize(nD);

  /*--- Determine the equidistant spacing in parametric space. ---*/
  const passivedouble dh = 2.0/nPoly;

  /*--- Double loop to determine the parametric coordinates. ---*/
  unsigned short ii = 0;
  for(unsigned short j=0; j<=nPoly; ++j) {
    const unsigned short uppBoundI = nPoly - j;
    for(unsigned short i=0; i<=uppBoundI; ++i, ++ii) {
      r[ii] = -1.0 + i*dh;
      s[ii] = -1.0 + j*dh;
    }
  }
}

void CFEMStandardElementBase::LocationTriangleGridDOFsLGL(vector<passivedouble> &r,
                                                          vector<passivedouble> &s) {

  /*--- The code to determine the parametric coordinates of the DOFs of the
        triangle is a translation of the Matlab code belonging to the book
        Nodal Discontinuous Galerkin Methods, Algorithms, Analysis and Applications,
        written by Jan S. Hesthaven and Tim Warburton. ---*/

  /*--- Local parameters. ---*/
  const passivedouble sqrt3 = sqrt(3.0);
  const passivedouble alphaOpt[] = {0.0000, 0.0000, 0.0000, 1.4152, 0.1001, 0.2751,
                                    0.9800, 1.0999, 1.2832, 1.3648, 1.4773, 1.4959,
                                    1.5743, 1.5770, 1.6223, 1.6258};

  /*--- Determine the number of DOFs of the triangle. As this function is also
        called for prisms, the member variable nDOFs cannot be used for this. ---*/
  const unsigned short nD = (nPoly+1)*(nPoly+2)/2;

  /*--- Allocate the memory for the help variables, as well as r and s. ---*/
  vector<passivedouble> L1(nD), L2(nD), L3(nD), rout(nD), warpf1(nD), warpf2(nD), warpf3(nD);

  r.resize(nD);
  s.resize(nD);

  /*--- Determine the uniform spacing. ---*/
  /*--- Determine the equidistant spacing in parametric space. ---*/
  const passivedouble dh = 1.0/nPoly;

  /*-----------------------------------------------------------------------------*/
  /*--- Step 1: Create the equidistributed nodes on the equilateral triangle. ---*/
  /*-----------------------------------------------------------------------------*/

  unsigned short ii = 0;
  for(unsigned short j=0; j<=nPoly; ++j) {
    const unsigned short uppBoundI = nPoly - j;
    for(unsigned short i=0; i<=uppBoundI; ++i, ++ii) {
      L1[ii] = j*dh;
      L3[ii] = i*dh;
      L2[ii] = 1.0 - L1[ii] - L3[ii];

      r[ii] = L3[ii] - L2[ii];
      s[ii] = (2.0*L1[ii] - L2[ii] - L3[ii])/sqrt3;
    }
  }

  /*---------------------------------------------*/
  /*--- Step 2: Modify the node distribution. ---*/
  /*---------------------------------------------*/

  /*--- Set the optimal value of alp for this case. ---*/
  const passivedouble alp = (nPoly < 16) ? alphaOpt[nPoly] : 5.0/3.0;

  /*--- Compute the warp factors for the DOFs for each edge. ---*/
  for(unsigned short i=0; i<nD; ++i) rout[i] = L3[i] - L2[i];
  WarpFactor(rout, warpf1);

  for(unsigned short i=0; i<nD; ++i) rout[i] = L1[i] - L3[i];
  WarpFactor(rout, warpf2);

  for(unsigned short i=0; i<nD; ++i) rout[i] = L2[i] - L1[i];
  WarpFactor(rout, warpf3);

  /*--- Loop over the DOFs to correct the coordinates. ---*/
  for(unsigned short i=0; i<nD; ++i) {

    /*--- Compute the blending functions for each edge. ---*/
    const passivedouble blend1 = 4.0*L2[i]*L3[i];
    const passivedouble blend2 = 4.0*L1[i]*L3[i];
    const passivedouble blend3 = 4.0*L1[i]*L2[i];

    /*--- Compute the combined blending and warp factor. ---*/
    const passivedouble warp1 = blend1*warpf1[i]*(1.0 + alp*alp*L1[i]*L1[i]);
    const passivedouble warp2 = blend2*warpf2[i]*(1.0 + alp*alp*L2[i]*L2[i]);
    const passivedouble warp3 = blend3*warpf3[i]*(1.0 + alp*alp*L3[i]*L3[i]);

    /*--- Compute the new coordinates. The multiplication factors in front
          of the warp factors corresponds to the cosine and sine of 0,
          2*pi/3 and 4*pi/3 respectively. ---*/
    r[i] = r[i] + warp1 - 0.5*(warp2 + warp3);
    s[i] = s[i] +   0.5*sqrt3*(warp2 - warp3);
  }
}

void CFEMStandardElementBase::DerLagBasisIntPointsLine(const vector<passivedouble>   &rDOFs,
                                                       const vector<passivedouble>   &rInt,
                                                       ColMajorMatrix<passivedouble> &derLag) {

  /*--- Determine the number of integration points along the line
        and its padded value. ---*/
  const unsigned short nIntLine    = rInt.size();
  const unsigned short nIntPadLine = ((nIntLine+baseVectorLen-1)/baseVectorLen)*baseVectorLen;

  /*--- Determine the inverse of the Vandermonde matrix of rDOFs. ---*/
  CSquareMatrixCM VInv(nPoly+1);
  Vandermonde1D(rDOFs, VInv.GetMat());
  VInv.Invert();

  /*--- Determine the gradient of the Vandermonde matrix of rInt. Make sure to
        allocate the number of rows to nIntPadLine and initialize them to zero. ---*/ 
  ColMajorMatrix<passivedouble> gradV(nIntPadLine,nPoly+1);
  gradV.setConstant(0.0);
  GradVandermonde1D(rInt, gradV);

  /*--- The derivatives of the Lagrangian basis functions can be obtained
        by multiplying gradV and VInv. ---*/
  VInv.MatMatMult('R', gradV, derLag);

  /*--- Check if the sum of the elements of the relevant rows of derLag is 0. ---*/
  for(unsigned short i=0; i<nIntLine; ++i) {
    passivedouble rowsum = 0.0;
    for(unsigned short j=0; j<=nPoly; ++j) rowsum += derLag(i,j);
    assert(fabs(rowsum) < 1.e-6);
  }
}

void CFEMStandardElementBase::HesLagBasisIntPointsLine(const vector<passivedouble>   &rDOFs,
                                                       const vector<passivedouble>   &rInt,
                                                       ColMajorMatrix<passivedouble> &hesLag) {

  /*--- Determine the number of integration points along the line
        and its padded value. ---*/
  const unsigned short nIntLine    = rInt.size();
  const unsigned short nIntPadLine = ((nIntLine+baseVectorLen-1)/baseVectorLen)*baseVectorLen;

  /*--- Determine the inverse of the Vandermonde matrix of rDOFs. ---*/
  CSquareMatrixCM VInv(nPoly+1);
  Vandermonde1D(rDOFs, VInv.GetMat());
  VInv.Invert();

  /*--- Determine the Hessian of the Vandermonde matrix of rInt. Make sure to
        allocate the number of rows to nIntPadLine and initialize them to zero. ---*/ 
  ColMajorMatrix<passivedouble> hesV(nIntPadLine,nPoly+1);
  hesV.setConstant(0.0);
  HesVandermonde1D(rInt, hesV);

  /*--- The Hessian of the Lagrangian basis functions can be obtained
        by multiplying hesV and VInv. ---*/
  VInv.MatMatMult('R', hesV, hesLag);

  /*--- Check if the sum of the elements of the relevant rows of hesLag is 0. ---*/
  for(unsigned short i=0; i<nIntLine; ++i) {
    passivedouble rowsum = 0.0;
    for(unsigned short j=0; j<=nPoly; ++j) rowsum += hesLag(i,j);
    assert(fabs(rowsum) < 1.e-6);
  }
}

void CFEMStandardElementBase::LagBasisIntPointsLine(const vector<passivedouble>   &rDOFs,
                                                    const vector<passivedouble>   &rInt,
                                                    ColMajorMatrix<passivedouble> &lag) {

  /*--- Determine the number of integration points along the line
        and its padded value. ---*/
  const unsigned short nIntLine    = rInt.size();
  const unsigned short nIntPadLine = ((nIntLine+baseVectorLen-1)/baseVectorLen)*baseVectorLen;

  /*--- Determine the inverse of the Vandermonde matrix of rDOFs. ---*/
  CSquareMatrixCM VInv(nPoly+1);
  Vandermonde1D(rDOFs, VInv.GetMat());
  VInv.Invert();

  /*--- Determine the Vandermonde matrix of rInt. Make sure to allocate
        the number of rows to nIntPadLine and initialize them to zero. ---*/ 
  ColMajorMatrix<passivedouble> V(nIntPadLine,nPoly+1);
  V.setConstant(0.0);
  Vandermonde1D(rInt, V);

  /*--- The Lagrangian basis functions can be obtained by multiplying
        V and VInv. ---*/
  VInv.MatMatMult('R', V, lag);

  /*--- Check if the sum of the elements of the relevant rows of lag is 1. ---*/
  for(unsigned short i=0; i<nIntLine; ++i) {
    passivedouble rowsum = -1.0;
    for(unsigned short j=0; j<=nPoly; ++j) rowsum += lag(i,j);
    assert(fabs(rowsum) < 1.e-6);
  }
}

passivedouble CFEMStandardElementBase::GradNormJacobi(unsigned short n,
                                                      unsigned short alpha,
                                                      unsigned short beta,
                                                      passivedouble  x) {

  /*--- Make a distinction for n == 0 and n > 0. For n == 0 the derivative is
        zero, because the polynomial itself is constant. ---*/
  passivedouble grad;
  if(n == 0) grad = 0.0;
  else
  {
    passivedouble tmp = n*(n+alpha+beta+1.0);
    grad              = sqrt(tmp)*NormJacobi(n-1, alpha+1, beta+1, x);
  }

  /*--- Return the gradient. ---*/
  return grad;
}

passivedouble CFEMStandardElementBase::HesNormJacobi(unsigned short n,
                                                     unsigned short alpha,
                                                     unsigned short beta,
                                                     passivedouble  x) {

  /*--- Make a distinction for n < 2 and n >= 2. For n < 2 the 2nd
        derivative is zero, because the polynomial itself is either
        constant or linear. ---*/
  passivedouble hes;
  if(n < 2) hes = 0.0;
  else
  {
    passivedouble tmp = n*(n+alpha+beta+1.0);
    hes               = sqrt(tmp)*GradNormJacobi(n-1, alpha+1, beta+1, x);
  }

  /*--- Return the gradient. ---*/
  return hes;
}


passivedouble CFEMStandardElementBase::NormJacobi(unsigned short n,
                                                  unsigned short alpha,
                                                  unsigned short beta,
                                                  passivedouble  x) {
  /*--- Some abbreviations. ---*/
  const passivedouble ap1   = alpha + 1;
  const passivedouble bp1   = beta  + 1;
  const passivedouble apb   = alpha + beta;
  const passivedouble apbp1 = apb + 1;
  const passivedouble apbp2 = apb + 2;
  const passivedouble apbp3 = apb + 3;
  const passivedouble b2ma2 = beta*beta - alpha*alpha;

  /*--- Determine the term, which involves the gamma function for P0. As the
        arguments are integers, this term can be computed easily, because
        Gamma(n+1) = n!. ---*/
  passivedouble Gamap1 = 1.0, Gambp1 = 1.0, Gamapbp2 = 1.0;
  for(unsigned short i=2; i<=alpha; ++i)          Gamap1   *= i;
  for(unsigned short i=2; i<=beta; ++i)           Gambp1   *= i;
  for(unsigned short i=2; i<=(alpha+beta+1); ++i) Gamapbp2 *= i;

  /*--- Initialize the normalized polynomials. ---*/
  passivedouble Pnm1 = sqrt(pow(0.5,apbp1)*Gamapbp2/(Gamap1*Gambp1));
  passivedouble Pn   = 0.5*Pnm1*(apbp2*x + alpha - beta)*sqrt(apbp3/(ap1*bp1));

  /*--- Take care of the special situation of n == 0. ---*/
  if(n == 0) Pn = Pnm1;
  else
  {
    /*--- The value of the normalized Jacobi polynomial must be obtained via recursion. ---*/
    for(unsigned short i=2; i<=n; ++i)
    {
      /*--- Compute the coefficients a for i and i-1 and the coefficient bi. ---*/
      unsigned short j = i-1;
      passivedouble tmp  = 2*j + apb;
      passivedouble aim1 = 2.0*sqrt(j*(j+apb)*(j+alpha)*(j+beta)/((tmp-1.0)*(tmp+1.0)))
                         / tmp;

      passivedouble bi = b2ma2/(tmp*(tmp+2.0));

      tmp              = 2*i + apb;
      passivedouble ai = 2.0*sqrt(i*(i+apb)*(i+alpha)*(i+beta)/((tmp-1.0)*(tmp+1.0)))
                       / tmp;

      /*--- Compute the new value of Pn and make sure to store Pnm1 correctly. ---*/
      tmp  = Pnm1;
      Pnm1 = Pn;

      Pn = ((x-bi)*Pn - aim1*tmp)/ai;
    }
  }

  /*--- Return Pn. ---*/
  return Pn;
}

void CFEMStandardElementBase::OwnGemm(dgemm_jit_kernel_t            &gemm,
                                      const int                     M,
                                      const int                     N,
                                      const int                     K,
                                      ColMajorMatrix<passivedouble> &A,
                                      ColMajorMatrix<su2double>     &B,
                                      ColMajorMatrix<su2double>     &C,
                                      const CConfig                 *config) {

  /*--- Check if the jitted gemm call can be used. ---*/
#if defined(PRIMAL_SOLVER) && defined(HAVE_MKL)

  /*--- Carry out the timing, if desired and call the jitted gemm function. ---*/
#ifdef PROFILE
  double timeGemm;
  if( config ) config->GEMM_Tick(&timeGemm);
#endif

  gemm(jitter, A.data(), B.data(), C.data());

#ifdef PROFILE
  if( config ) config->GEMM_Tock(timeGemm, M, N, K);
#endif

#else

  /*--- Use the interface to the more standard BLAS functionality. ---*/
  blasFunctions.gemm(M, N, K, A.data(), B.data(), C.data(), config);

#endif
}

void CFEMStandardElementBase::SetUpJittedGEMM(const int          M,
                                              const int          N,
                                              const int          K,
                                              void               *&val_jitter,
                                              dgemm_jit_kernel_t &val_gemm) {

  /*--- Check if the jitted GEMM call can be created. ---*/
#if defined(PRIMAL_SOLVER) && defined(HAVE_MKL)

  /*--- Create the GEMM Kernel and check if it went okay. ---*/
  passivedouble alpha = 1.0;
  passivedouble beta  = 0.0;
  mkl_jit_status_t status = mkl_jit_create_dgemm(&val_jitter, MKL_COL_MAJOR, MKL_NOTRANS,
                                                 MKL_NOTRANS, M, N, K, alpha, M, K, beta, M);

  if(status == MKL_JIT_ERROR)
    SU2_MPI::Error(string("Jitted gemm kernel could not be created"), CURRENT_FUNCTION);

  /*--- Retrieve the function pointer to the DGEMM kernel
        void val_gemm(void*, double*, double*, double*). ---*/
  val_gemm = mkl_jit_get_dgemm_ptr(val_jitter);

#endif
}

/*----------------------------------------------------------------------------------*/
/*          Private member functions of CFEMStandardElementBase.                    */
/*----------------------------------------------------------------------------------*/

passivedouble CFEMStandardElementBase::Legendre(unsigned short n,
                                                passivedouble  x) {

  /*--- Initialization of the polynomials Pnm1 and Pn. ---*/
  passivedouble Pnm1 = 1.0;
  passivedouble Pn   = x;

  /*--- Take care of the special situation of n == 0. ---*/
  if(n == 0) Pn = Pnm1;
  else {

    /*--- Recursive definition of Pn. ---*/
    for(unsigned short i=2; i<=n; ++i)
    {
      const passivedouble tmp = Pnm1;
      Pnm1 = Pn;

      Pn = ((2*i-1)*x*Pn - (i-1)*tmp)/i;
    }
  }

  /*--- Return Pn. ---*/
  return Pn;
}

void CFEMStandardElementBase::GradVandermonde1D(const vector<passivedouble>   &r,
                                                ColMajorMatrix<passivedouble> &VDr) {

  /*--- Compute the gradient of the 1D Vandermonde matrix. ---*/
  for(unsigned short j=0; j<=nPoly; ++j)
    for(unsigned short i=0; i<r.size(); ++i)
      VDr(i,j) = GradNormJacobi(j, 0, 0, r[i]);
}

void CFEMStandardElementBase::HesVandermonde1D(const vector<passivedouble>   &r,
                                               ColMajorMatrix<passivedouble> &VD2r) {

  /*--- Compute the Hessian of the 1D Vandermonde matrix. ---*/
  for(unsigned short j=0; j<=nPoly; ++j)
    for(unsigned short i=0; i<r.size(); ++i)
      VD2r(i,j) = HesNormJacobi(j, 0, 0, r[i]);
}

void CFEMStandardElementBase::Vandermonde1D(const vector<passivedouble>   &r,
                                            ColMajorMatrix<passivedouble> &V) {

  /*--- Compute the 1D Vandermonde matrix. ---*/
  for(unsigned short j=0; j<=nPoly; ++j)
    for(unsigned short i=0; i<r.size(); ++i)
      V(i,j) = NormJacobi(j, 0, 0, r[i]);
}

void CFEMStandardElementBase::WarpFactor(const vector<passivedouble> &rout,
                                         vector<passivedouble>       &warp) {

  /*--- The code to determine the warp factor is a translation of the Matlab code
        belonging to the book
        Nodal Discontinuous Galerkin Methods, Algorithms, Analysis and Applications,
        written by Jan S. Hesthaven and Tim Warburton. ---*/

  /*--- Determine the number of DOFs for which the warp factor
        must be computed. ---*/
  const unsigned short nD = rout.size();

  /*--- Compute the 1D Gauss-Lobato and an equidistant node distribution. ---*/
  vector<passivedouble> r1DEqui, r1DLGL;
  Location1DGridDOFsEquidistant(r1DEqui);
  Location1DGridDOFsLGL(r1DLGL);

  /*--- Determine the 1D Vandermonde matrix based on the equidistant node
        distribution. Invert the matrix afterwards. ---*/
  CSquareMatrixCM V1D(nPoly+1);
  Vandermonde1D(r1DEqui, V1D.GetMat());
  V1D.Invert();

  /*--- Determine the Lagrange polynomials at rout. This is accomplished by
        first computing the Legendre polynomials at rout and multiply
        this value by V1D. The result is stored in Lmat. ---*/
  ColMajorMatrix<passivedouble> Pmat(nD,nPoly+1), Lmat;
  Vandermonde1D(rout, Pmat);
  V1D.MatMatMult('R', Pmat, Lmat);

  /*--- Determine the difference between the 1D LGL points and
        the equidistant points. The result is stored in r1DLGL. ---*/
  for(unsigned short i=0; i<=nPoly; ++i) r1DLGL[i] -= r1DEqui[i];

  /*--- Compute the unscaled warp factors, which is Lmat*r1DLGL. ---*/
  for(unsigned short j=0; j<nD; ++j) {
    warp[j] = 0.0;
    for(unsigned short i=0; i<=nPoly; ++i)
      warp[j] += Lmat(j,i)*r1DLGL[i];
  }

  /*--- Scale the warp factor by 1-r*r. ---*/
  for(unsigned short i=0; i<nD; ++i) {

    /*--- Take care of the exceptional case that r = 1. ---*/
    passivedouble zerof = 0.0;
    if(fabs(rout[i]) < (1.0-1.e-10)) zerof = 1.0;

    /*--- Compute 1-r*r. ---*/
    const passivedouble sf = 1.0 - zerof*rout[i]*rout[i];

    /*--- Compute the scaled warp factor. ---*/
    warp[i] = warp[i]/sf + warp[i]*(zerof-1.0);
    if(fabs(warp[i]) < 1.e-10) warp[i] = 0.0;
  }
}
