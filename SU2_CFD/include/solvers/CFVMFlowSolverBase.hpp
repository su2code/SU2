/*!
 * \file CFVMFlowSolverBase.hpp
 * \brief Base class template for all FVM flow solvers.
 * \version 7.0.5 "Blackbird"
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

#pragma once

#include "CSolver.hpp"
#include "../../../Common/include/omp_structure.hpp"


template<class VariableType, ENUM_REGIME FlowRegime>
class CFVMFlowSolverBase : public CSolver {
protected:
  enum : size_t {MAXNDIM = 3};    /*!< \brief Max number of space dimensions, used in some static arrays. */

  enum : size_t {OMP_MAX_SIZE = 512};  /*!< \brief Max chunk size for light point loops. */
  enum : size_t {OMP_MIN_SIZE = 32};   /*!< \brief Min chunk size for edge loops (max is color group size). */

  unsigned long omp_chunk_size;  /*!< \brief Chunk size used in light point loops. */

  su2double
  Mach_Inf = 0.0,            /*!< \brief Mach number at the infinity. */
  Density_Inf = 0.0,         /*!< \brief Density at the infinity. */
  Energy_Inf = 0.0,          /*!< \brief Energy at the infinity. */
  Temperature_Inf = 0.0,     /*!< \brief Energy at the infinity. */
  Pressure_Inf = 0.0,        /*!< \brief Pressure at the infinity. */
  *Velocity_Inf = nullptr;   /*!< \brief Flow Velocity vector at the infinity. */

  su2double Viscosity_Inf;   /*!< \brief Viscosity at the infinity. */
  su2double Tke_Inf;         /*!< \brief Turbulent kinetic energy at the infinity. */

  su2double StrainMag_Max;   /*!< \brief Maximum Strain Rate magnitude. */
  su2double Omega_Max;       /*!< \brief Maximum Omega. */

  /*!
   * \brief Auxilary types to store common aero coefficients (avoids repeating oneself so much).
   */
  struct AeroCoeffsArray {
    su2double* CD = nullptr;      /*!< \brief Drag coefficient. */
    su2double* CL = nullptr;      /*!< \brief Lift coefficient. */
    su2double* CSF = nullptr;     /*!< \brief Sideforce coefficient. */
    su2double* CEff = nullptr;    /*!< \brief Efficiency (Cl/Cd). */
    su2double* CFx = nullptr;     /*!< \brief x Force coefficient. */
    su2double* CFy = nullptr;     /*!< \brief y Force coefficient. */
    su2double* CFz = nullptr;     /*!< \brief z Force coefficient. */
    su2double* CMx = nullptr;     /*!< \brief x Moment coefficient. */
    su2double* CMy = nullptr;     /*!< \brief y Moment coefficient. */
    su2double* CMz = nullptr;     /*!< \brief z Moment coefficient. */
    su2double* CoPx = nullptr;    /*!< \brief x Moment coefficient. */
    su2double* CoPy = nullptr;    /*!< \brief y Moment coefficient. */
    su2double* CoPz = nullptr;    /*!< \brief z Moment coefficient. */
    su2double* CT = nullptr;      /*!< \brief Thrust coefficient. */
    su2double* CQ = nullptr;      /*!< \brief Torque coefficient. */
    su2double* CMerit = nullptr;  /*!< \brief Rotor Figure of Merit. */
    int _size = 0;                /*!< \brief Array size. */

    void allocate(int size);      /*!< \brief Allocates arrays. */

    void setZero(int i);          /*!< \brief Sets all values to zero at a particular index. */
    void setZero() {              /*!< \brief Sets all values to zero for all indices. */
      for(int i=0; i<_size; ++i) setZero(i);
    }

    AeroCoeffsArray(int size = 0) : _size(size) { if(size) allocate(size); }

    ~AeroCoeffsArray();
  };

  /*!
   * \brief Scalar version of the coefficients type.
   */
  struct AeroCoeffs {
    su2double CD,CL,CSF,CEff,CFx,CFy,CFz,CMx,CMy,CMz,CoPx,CoPy,CoPz,CT,CQ,CMerit;

    void setZero() {
      CD=CL=CSF=CEff=CFx=CFy=CFz=CMx=CMy=CMz=CoPx=CoPy=CoPz=CT=CQ=CMerit=0.0;
    }

    AeroCoeffs() { setZero(); }
  };

  AeroCoeffsArray InvCoeff;          /*!< \brief Inviscid pressure contributions for each boundary. */
  AeroCoeffsArray SurfaceInvCoeff;   /*!< \brief Inviscid pressure contributions for each monitoring boundary. */
  AeroCoeffs AllBoundInvCoeff;       /*!< \brief Total pressure contribution for all the boundaries. */

  AeroCoeffsArray MntCoeff;          /*!< \brief Inviscid momentum contributions for each boundary. */
  AeroCoeffsArray SurfaceMntCoeff;   /*!< \brief Inviscid momentum contributions for each monitoring boundary. */
  AeroCoeffs AllBoundMntCoeff;       /*!< \brief Total momentum contribution for all the boundaries. */

  AeroCoeffsArray ViscCoeff;         /*!< \brief Viscous contributions for each boundary. */
  AeroCoeffsArray SurfaceViscCoeff;  /*!< \brief Viscous contributions for each monitoring boundary. */
  AeroCoeffs AllBoundViscCoeff;      /*!< \brief Total pressure contribution for all the boundaries. */

  AeroCoeffsArray SurfaceCoeff;      /*!< \brief Totals for each monitoring surface. */
  AeroCoeffs TotalCoeff;             /*!< \brief Totals for all boundaries. */

  su2double
  *CNearFieldOF_Inv = nullptr,      /*!< \brief Near field pressure (inviscid contribution) for each boundary. */
  Total_CNearFieldOF = 0.0,         /*!< \brief Total Near-Field Pressure coefficient for all the boundaries. */
  AllBound_CNearFieldOF_Inv = 0.0,  /*!< \brief Near-Field press coefficient (inviscid contribution) for all the boundaries. */
  Total_Heat = 0.0,                 /*!< \brief Total heat load for all the boundaries. */
  Total_MaxHeat = 0.0,              /*!< \brief Maximum heat flux on all boundaries. */
  *Surface_HF_Visc = nullptr,       /*!< \brief Total (integrated) heat flux for each monitored surface. */
  *Surface_MaxHF_Visc = nullptr,    /*!< \brief Maximum heat flux for each monitored surface. */
  *HF_Visc = nullptr,               /*!< \brief Heat load (viscous contribution) for each boundary. */
  *MaxHF_Visc = nullptr,            /*!< \brief Maximum heat flux (viscous contribution) for each boundary. */
  AllBound_HF_Visc,                 /*!< \brief Heat load (viscous contribution) for all the boundaries. */
  AllBound_MaxHF_Visc;              /*!< \brief Maximum heat flux (viscous contribution) for all boundaries. */

  su2double
  **HeatFlux = nullptr,          /*!< \brief Heat transfer coefficient for each boundary and vertex. */
  ***CSkinFriction = nullptr,    /*!< \brief Skin friction coefficient for each boundary and vertex. */
  **CPressure = nullptr,         /*!< \brief Pressure coefficient for each boundary and vertex. */
  **YPlus = nullptr;             /*!< \brief Yplus for each boundary and vertex. */

  bool space_centered,       /*!< \brief True if space centered scheeme used. */
  euler_implicit,            /*!< \brief True if euler implicit scheme used. */
  least_squares;             /*!< \brief True if computing gradients by least squares. */
  su2double Gamma;           /*!< \brief Fluid's Gamma constant (ratio of specific heats). */
  su2double Gamma_Minus_One; /*!< \brief Fluids's Gamma - 1.0  . */

  /*--- Sliding meshes variables ---*/

  su2double ****SlidingState = nullptr;
  int **SlidingStateNodes = nullptr;

  /*--- Shallow copy of grid coloring for OpenMP parallelization. ---*/

#ifdef HAVE_OMP
  vector<GridColor<> > EdgeColoring;   /*!< \brief Edge colors. */
  bool ReducerStrategy = false;        /*!< \brief If the reducer strategy is in use. */
#else
  array<DummyGridColor<>,1> EdgeColoring;
  /*--- Never use the reducer strategy if compiling for MPI-only. ---*/
  static constexpr bool ReducerStrategy = false;
#endif

  /*--- Edge fluxes, for OpenMP parallelization off difficult-to-color grids.
   * We first store the fluxes and then compute the sum for each cell.
   * This strategy is thread-safe but lower performance than writting to both
   * end points of each edge, so we only use it when necessary, i.e. when the
   * coloring does not allow "enough" parallelism. ---*/

  CSysVector<su2double> EdgeFluxes; /*!< \brief Flux across each edge. */


  VariableType* nodes = nullptr;  /*!< \brief The highest level in the variable hierarchy this solver can safely use. */

  /*!
   * \brief Return nodes to allow CSolver::base_nodes to be set.
   */
  inline CVariable* GetBaseClassPointerToNodes() final { return nodes; }

  /*!
   * \brief Default constructor.
   */
  CFVMFlowSolverBase() = default;

public:

  /*!
   * \brief Compute the density at the infinity.
   * \return Value of the density at the infinity.
   */
  inline su2double GetDensity_Inf(void) const final { return Density_Inf; }

  /*!
   * \brief Compute 2-norm of the velocity at the infinity.
   * \return Value of the 2-norm of the velocity at the infinity.
   */
  inline su2double GetModVelocity_Inf(void) const final {
    su2double Vel2 = 0;
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      Vel2 += Velocity_Inf[iDim]*Velocity_Inf[iDim];
    return sqrt(Vel2);
  }


  /*!
   * \brief Compute the pressure at the infinity.
   * \return Value of the pressure at the infinity.
   */
  inline su2double GetPressure_Inf(void) const final { return Pressure_Inf; }

    /*!
   * \brief Get the temperature value at infinity.
   * \return Value of the temperature at infinity.
   */
  inline su2double GetTemperature_Inf(void) const { return Temperature_Inf; }

  /*!
   * \brief Compute the density multiply by velocity at the infinity.
   * \param[in] val_dim - Index of the velocity vector.
   * \return Value of the density multiply by the velocity at the infinity.
   */
  inline su2double GetDensity_Velocity_Inf(unsigned short val_dim) const final { return Density_Inf*Velocity_Inf[val_dim]; }

  /*!
   * \brief Get the velocity at the infinity.
   * \param[in] val_dim - Index of the velocity vector.
   * \return Value of the velocity at the infinity.
   */
  inline su2double GetVelocity_Inf(unsigned short val_dim) const final { return Velocity_Inf[val_dim]; }

  /*!
   * \brief Get the velocity at the infinity.
   * \return Value of the velocity at the infinity.
   */
  inline su2double *GetVelocity_Inf(void) final { return Velocity_Inf; }

  /*!
   * \brief Set the velocity at infinity.
   * \param[in] val_dim - Index of the velocity vector.
   * \param[in] val_velocity - Value of the velocity.
   */
  inline void SetVelocity_Inf(unsigned short val_dim, su2double val_velocity) final {
    Velocity_Inf[val_dim] = val_velocity;
  }

  /*!
   * \brief Compute the density multiply by energy at the infinity.
   * \return Value of the density multiply by  energy at the infinity.
   */
  inline su2double GetDensity_Energy_Inf(void) const final { return Density_Inf*Energy_Inf; }

  /*!
   * \brief Set the freestream pressure.
   * \param[in] Value of freestream pressure.
   */
  inline void SetPressure_Inf(su2double p_inf) final { Pressure_Inf = p_inf; }

  /*!
   * \brief Set the freestream temperature.
   * \param[in] Value of freestream temperature.
   */
  inline void SetTemperature_Inf(su2double t_inf) final { Temperature_Inf = t_inf; }

  /*!
   * \brief Set the freestream temperature.
   * \param[in] Value of freestream temperature.
   */
  inline void SetDensity_Inf(su2double rho_inf) final { Density_Inf = rho_inf; }

  /*!
   * \brief Compute the viscosity at the infinity.
   * \return Value of the viscosity at the infinity.
   */
  inline su2double GetViscosity_Inf(void) const final { return Viscosity_Inf; }

  /*!
   * \brief Get the turbulent kinetic energy at the infinity.
   * \return Value of the turbulent kinetic energy at the infinity.
   */
  inline su2double GetTke_Inf(void) const final { return Tke_Inf; }

  /*!
   * \brief Get the max Omega.
   * \return Value of the max Omega.
   */
  inline su2double GetOmega_Max(void) const final { return Omega_Max; }

  /*!
   * \brief Get the max Strain rate magnitude.
   * \return Value of the max Strain rate magnitude.
   */
  inline su2double GetStrainMag_Max(void) const final { return StrainMag_Max; }

  /*!
   * \brief Provide the non dimensional lift coefficient (inviscid contribution).
   * \param val_marker Surface where the coefficient is going to be computed.
   * \return Value of the lift coefficient (inviscid contribution) on the surface <i>val_marker</i>.
   */
  inline su2double GetCL_Inv(unsigned short val_marker) const final { return InvCoeff.CL[val_marker]; }

  /*!
   * \brief Provide the non dimensional drag coefficient (inviscid contribution).
   * \param val_marker Surface where the coeficient is going to be computed.
   * \return Value of the drag coefficient (inviscid contribution) on the surface <i>val_marker</i>.
   */
  inline su2double GetCD_Inv(unsigned short val_marker) const final { return InvCoeff.CD[val_marker]; }

  /*!
   * \brief Provide the non dimensional lift coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the lift coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CL(unsigned short val_marker) const final { return SurfaceCoeff.CL[val_marker]; }

  /*!
   * \brief Provide the non dimensional drag coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the drag coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CD(unsigned short val_marker) const final { return SurfaceCoeff.CD[val_marker]; }

  /*!
   * \brief Provide the non dimensional side-force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side-force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CSF(unsigned short val_marker) const final { return SurfaceCoeff.CSF[val_marker]; }

  /*!
   * \brief Provide the non dimensional side-force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side-force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CEff(unsigned short val_marker) const final { return SurfaceCoeff.CEff[val_marker]; }

  /*!
   * \brief Provide the non dimensional x force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CFx(unsigned short val_marker) const final { return SurfaceCoeff.CFx[val_marker]; }

  /*!
   * \brief Provide the non dimensional y force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CFy(unsigned short val_marker) const final { return SurfaceCoeff.CFy[val_marker]; }

  /*!
   * \brief Provide the non dimensional z force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CFz(unsigned short val_marker) const final { return SurfaceCoeff.CFz[val_marker]; }

  /*!
   * \brief Provide the non dimensional x moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x moment coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CMx(unsigned short val_marker) const final { return SurfaceCoeff.CMx[val_marker]; }

  /*!
   * \brief Provide the non dimensional y moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y moment coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CMy(unsigned short val_marker) const final { return SurfaceCoeff.CMy[val_marker]; }

  /*!
   * \brief Provide the non dimensional z moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z moment coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CMz(unsigned short val_marker) const final { return SurfaceCoeff.CMz[val_marker]; }

  /*!
   * \brief Provide the non dimensional lift coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the lift coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CL_Inv(unsigned short val_marker) const final { return SurfaceInvCoeff.CL[val_marker]; }

  /*!
   * \brief Provide the non dimensional drag coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the drag coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CD_Inv(unsigned short val_marker) const final { return SurfaceInvCoeff.CD[val_marker]; }

  /*!
   * \brief Provide the non dimensional side-force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side-force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CSF_Inv(unsigned short val_marker) const final { return SurfaceInvCoeff.CSF[val_marker]; }

  /*!
   * \brief Provide the non dimensional side-force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side-force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CEff_Inv(unsigned short val_marker) const final { return SurfaceInvCoeff.CEff[val_marker]; }

  /*!
   * \brief Provide the non dimensional x force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CFx_Inv(unsigned short val_marker) const final { return SurfaceInvCoeff.CFx[val_marker]; }

  /*!
   * \brief Provide the non dimensional y force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CFy_Inv(unsigned short val_marker) const final { return SurfaceInvCoeff.CFy[val_marker]; }

  /*!
   * \brief Provide the non dimensional z force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CFz_Inv(unsigned short val_marker) const final { return SurfaceInvCoeff.CFz[val_marker]; }

  /*!
   * \brief Provide the non dimensional x moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x moment coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CMx_Inv(unsigned short val_marker) const final { return SurfaceInvCoeff.CMx[val_marker]; }

  /*!
   * \brief Provide the non dimensional y moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y moment coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CMy_Inv(unsigned short val_marker) const final { return SurfaceInvCoeff.CMy[val_marker]; }

  /*!
   * \brief Provide the non dimensional z moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z moment coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CMz_Inv(unsigned short val_marker) const final { return SurfaceInvCoeff.CMz[val_marker]; }

  /*!
   * \brief Provide the non dimensional lift coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the lift coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CL_Mnt(unsigned short val_marker) const final { return SurfaceMntCoeff.CL[val_marker]; }

  /*!
   * \brief Provide the non dimensional drag coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the drag coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CD_Mnt(unsigned short val_marker) const final { return SurfaceMntCoeff.CD[val_marker]; }

  /*!
   * \brief Provide the non dimensional side-force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side-force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CSF_Mnt(unsigned short val_marker) const final { return SurfaceMntCoeff.CSF[val_marker]; }

  /*!
   * \brief Provide the non dimensional side-force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side-force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CEff_Mnt(unsigned short val_marker) const final { return SurfaceMntCoeff.CEff[val_marker]; }

  /*!
   * \brief Provide the non dimensional x force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CFx_Mnt(unsigned short val_marker) const final { return SurfaceMntCoeff.CFx[val_marker]; }

  /*!
   * \brief Provide the non dimensional y force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CFy_Mnt(unsigned short val_marker) const final { return SurfaceMntCoeff.CFy[val_marker]; }

  /*!
   * \brief Provide the non dimensional z force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CFz_Mnt(unsigned short val_marker) const final { return SurfaceMntCoeff.CFz[val_marker]; }

  /*!
   * \brief Provide the non dimensional x moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x moment coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CMx_Mnt(unsigned short val_marker) const final { return SurfaceMntCoeff.CMx[val_marker]; }

  /*!
   * \brief Provide the non dimensional y moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y moment coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CMy_Mnt(unsigned short val_marker) const final { return SurfaceMntCoeff.CMy[val_marker]; }

  /*!
   * \brief Provide the non dimensional z moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z moment coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CMz_Mnt(unsigned short val_marker) const final { return SurfaceMntCoeff.CMz[val_marker]; }

  /*!
   * \brief Provide the non dimensional sideforce coefficient (inviscid contribution).
   * \param val_marker Surface where the coeficient is going to be computed.
   * \return Value of the sideforce coefficient (inviscid contribution) on the surface <i>val_marker</i>.
   */
  inline su2double GetCSF_Inv(unsigned short val_marker) const final { return InvCoeff.CSF[val_marker]; }

  /*!
   * \brief Provide the non dimensional efficiency coefficient (inviscid contribution).
   * \param val_marker Surface where the coeficient is going to be computed.
   * \return Value of the efficiency coefficient (inviscid contribution) on the surface <i>val_marker</i>.
   */
  inline su2double GetCEff_Inv(unsigned short val_marker) const final { return InvCoeff.CEff[val_marker]; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional sideforce coefficient.
   * \return Value of the sideforce coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CSF() const final { return TotalCoeff.CSF; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CEff() const final { return TotalCoeff.CEff; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional lift coefficient.
   * \return Value of the lift coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CL() const final { return TotalCoeff.CL; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CD() const final { return TotalCoeff.CD; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional x moment coefficient.
   * \return Value of the moment x coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CMx() const final { return TotalCoeff.CMx; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional y moment coefficient.
   * \return Value of the moment y coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CMy() const final { return TotalCoeff.CMy; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional z moment coefficient.
   * \return Value of the moment z coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CMz() const final { return TotalCoeff.CMz; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional x moment coefficient.
   * \return Value of the moment x coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CoPx() const final { return TotalCoeff.CoPx; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional y moment coefficient.
   * \return Value of the moment y coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CoPy() const final { return TotalCoeff.CoPy; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional z moment coefficient.
   * \return Value of the moment z coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CoPz() const final { return TotalCoeff.CoPz; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional x force coefficient.
   * \return Value of the force x coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CFx() const final { return TotalCoeff.CFx; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional y force coefficient.
   * \return Value of the force y coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CFy() const final { return TotalCoeff.CFy; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional z force coefficient.
   * \return Value of the force z coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CFz() const final { return TotalCoeff.CFz; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional thrust coefficient.
   * \return Value of the rotor efficiency coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CT() const final { return TotalCoeff.CT; }

  /*!
   * \brief Store the total (inviscid + viscous) non dimensional thrust coefficient.
   * \param[in] val_Total_CT - Value of the total thrust coefficient.
   */
  inline void SetTotal_CT(su2double val_Total_CT) final { TotalCoeff.CT = val_Total_CT; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional torque coefficient.
   * \return Value of the rotor efficiency coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CQ() const final { return TotalCoeff.CQ; }

  /*!
   * \brief Store the total (inviscid + viscous) non dimensional torque coefficient.
   * \param[in] val_Total_CQ - Value of the total torque coefficient.
   */
  inline void SetTotal_CQ(su2double val_Total_CQ) final { TotalCoeff.CQ = val_Total_CQ; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional rotor Figure of Merit.
   * \return Value of the rotor efficiency coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CMerit() const final { return TotalCoeff.CMerit; }

  /*!
   * \brief Store the total (inviscid + viscous) non dimensional drag coefficient.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  inline void SetTotal_CD(su2double val_Total_CD) final { TotalCoeff.CD = val_Total_CD; }

  /*!
   * \brief Store the total (inviscid + viscous) non dimensional lift coefficient.
   * \param[in] val_Total_CL - Value of the total lift coefficient.
   */
  inline void SetTotal_CL(su2double val_Total_CL) final { TotalCoeff.CL = val_Total_CL; }

  /*!
   * \brief Get the inviscid contribution to the lift coefficient.
   * \return Value of the lift coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CL_Inv() const final { return AllBoundInvCoeff.CL; }

  /*!
   * \brief Get the inviscid contribution to the drag coefficient.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CD_Inv() const final { return AllBoundInvCoeff.CD; }

  /*!
   * \brief Get the inviscid contribution to the sideforce coefficient.
   * \return Value of the sideforce coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CSF_Inv() const final { return AllBoundInvCoeff.CSF; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CEff_Inv() const final { return AllBoundInvCoeff.CEff; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CMx_Inv() const final { return AllBoundInvCoeff.CMx; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CMy_Inv() const final { return AllBoundInvCoeff.CMy; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CMz_Inv() const final { return AllBoundInvCoeff.CMz; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CoPx_Inv() const final { return AllBoundInvCoeff.CoPx; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CoPy_Inv() const final { return AllBoundInvCoeff.CoPy; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CoPz_Inv() const final { return AllBoundInvCoeff.CoPz; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CFx_Inv() const final { return AllBoundInvCoeff.CFx; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CFy_Inv() const final { return AllBoundInvCoeff.CFy; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CFz_Inv() const final { return AllBoundInvCoeff.CFz; }

  /*!
   * \brief Get the inviscid contribution to the lift coefficient.
   * \return Value of the lift coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CL_Mnt() const final { return AllBoundMntCoeff.CL; }

  /*!
   * \brief Get the inviscid contribution to the drag coefficient.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CD_Mnt() const final { return AllBoundMntCoeff.CD; }

  /*!
   * \brief Get the inviscid contribution to the sideforce coefficient.
   * \return Value of the sideforce coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CSF_Mnt() const final { return AllBoundMntCoeff.CSF; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CEff_Mnt() const final { return AllBoundMntCoeff.CEff; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CMx_Mnt() const final { return AllBoundMntCoeff.CMx; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CMy_Mnt() const final { return AllBoundMntCoeff.CMy; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CMz_Mnt() const final { return AllBoundMntCoeff.CMz; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CoPx_Mnt() const final { return AllBoundMntCoeff.CoPx; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CoPy_Mnt() const final { return AllBoundMntCoeff.CoPy; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CoPz_Mnt() const final { return AllBoundMntCoeff.CoPz; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CFx_Mnt() const final { return AllBoundMntCoeff.CFx; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CFy_Mnt() const final { return AllBoundMntCoeff.CFy; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CFz_Mnt() const final { return AllBoundMntCoeff.CFz; }

  /*!
   * \brief Provide the non dimensional lift coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the lift coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CL_Visc(unsigned short val_marker) const final { return SurfaceViscCoeff.CL[val_marker]; }

  /*!
   * \brief Provide the non dimensional drag coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the drag coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CD_Visc(unsigned short val_marker) const final { return SurfaceViscCoeff.CD[val_marker]; }

  /*!
   * \brief Provide the non dimensional side-force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side-force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CSF_Visc(unsigned short val_marker) const final { return SurfaceViscCoeff.CSF[val_marker]; }

  /*!
   * \brief Provide the non dimensional side-force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side-force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CEff_Visc(unsigned short val_marker) const final { return SurfaceViscCoeff.CEff[val_marker]; }

  /*!
   * \brief Provide the non dimensional x force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CFx_Visc(unsigned short val_marker) const final { return SurfaceViscCoeff.CFx[val_marker]; }

  /*!
   * \brief Provide the non dimensional y force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CFy_Visc(unsigned short val_marker) const final { return SurfaceViscCoeff.CFy[val_marker]; }

  /*!
   * \brief Provide the non dimensional z force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CFz_Visc(unsigned short val_marker) const final { return SurfaceViscCoeff.CFz[val_marker]; }

  /*!
   * \brief Provide the non dimensional x moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x moment coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CMx_Visc(unsigned short val_marker) const final { return SurfaceViscCoeff.CMx[val_marker]; }

  /*!
   * \brief Provide the non dimensional y moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y moment coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CMy_Visc(unsigned short val_marker) const final { return SurfaceViscCoeff.CMy[val_marker]; }

  /*!
   * \brief Provide the non dimensional z moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z moment coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CMz_Visc(unsigned short val_marker) const final { return SurfaceViscCoeff.CMz[val_marker]; }

  /*!
   * \brief Get the inviscid contribution to the lift coefficient.
   * \return Value of the lift coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CL_Visc() const final { return AllBoundViscCoeff.CL; }

  /*!
   * \brief Get the inviscid contribution to the drag coefficient.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CD_Visc() const final { return AllBoundViscCoeff.CD; }

  /*!
   * \brief Get the inviscid contribution to the sideforce coefficient.
   * \return Value of the sideforce coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CSF_Visc() const final { return AllBoundViscCoeff.CSF; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CEff_Visc() const final { return AllBoundViscCoeff.CEff; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CMx_Visc() const final { return AllBoundViscCoeff.CMx; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CMy_Visc() const final { return AllBoundViscCoeff.CMy; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CMz_Visc() const final { return AllBoundViscCoeff.CMz; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CoPx_Visc() const final { return AllBoundViscCoeff.CoPx; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CoPy_Visc() const final { return AllBoundViscCoeff.CoPy; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CoPz_Visc() const final { return AllBoundViscCoeff.CoPz; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CFx_Visc() const final { return AllBoundViscCoeff.CFx; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CFy_Visc() const final { return AllBoundViscCoeff.CFy; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CFz_Visc() const final { return AllBoundViscCoeff.CFz; }

  /*!
   * \brief Get the non dimensional lift coefficient (viscous contribution).
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the lift coefficient (viscous contribution) on the surface <i>val_marker</i>.
   */
  inline su2double GetCL_Visc(unsigned short val_marker) const final { return ViscCoeff.CL[val_marker]; }

  /*!
   * \brief Get the non dimensional sideforce coefficient (viscous contribution).
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the sideforce coefficient (viscous contribution) on the surface <i>val_marker</i>.
   */
  inline su2double GetCSF_Visc(unsigned short val_marker) const final { return ViscCoeff.CSF[val_marker]; }

  /*!
   * \brief Get the non dimensional drag coefficient (viscous contribution).
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the drag coefficient (viscous contribution) on the surface <i>val_marker</i>.
   */
  inline su2double GetCD_Visc(unsigned short val_marker) const final { return ViscCoeff.CD[val_marker]; }

  /*!
   * \brief Get the total heat flux.
   * \param[in] val_marker - Surface marker where the heat flux is computed.
   * \return Value of the integrated heat flux (viscous contribution) on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_HF_Visc(unsigned short val_marker) const final { return Surface_HF_Visc[val_marker]; }

  /*!
   * \brief Get the maximum (per surface) heat flux.
   * \param[in] val_marker - Surface marker where the heat flux is computed.
   * \return Value of the maximum heat flux (viscous contribution) on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_MaxHF_Visc(unsigned short val_marker) const final { return Surface_MaxHF_Visc[val_marker]; }

  /*!
   * \brief Compute the pressure forces and all the adimensional coefficients.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Pressure_Forces(const CGeometry* geometry, const CConfig* config) final;

  /*!
   * \brief Compute the pressure forces and all the adimensional coefficients.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Momentum_Forces(const CGeometry* geometry, const CConfig* config) final;

  /*!
   * \brief Compute the viscous forces and all the addimensional coefficients.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Friction_Forces(const CGeometry* geometry, const CConfig* config) final;

  /*!
   * \brief Allocates the final pointer of SlidingState depending on how many donor vertex donate to it. That number is stored in SlidingStateNodes[val_marker][val_vertex].
   * \param[in] val_marker   - marker index
   * \param[in] val_vertex   - vertex index
   */
  inline void SetSlidingStateStructure(unsigned short val_marker, unsigned long val_vertex) final {
    for (int iVar = 0; iVar < nPrimVar+1; iVar++){
      if( SlidingState[val_marker][val_vertex][iVar] != nullptr )
        delete [] SlidingState[val_marker][val_vertex][iVar];
    }

    for (int iVar = 0; iVar < nPrimVar+1; iVar++)
      SlidingState[val_marker][val_vertex][iVar] = new su2double[ GetnSlidingStates(val_marker, val_vertex) ];
  }

  /*!
   * \brief Set the outer state for fluid interface nodes.
   * \param[in] val_marker   - marker index
   * \param[in] val_vertex   - vertex index
   * \param[in] val_state    - requested state component
   * \param[in] donor_index  - index of the donor node to set
   * \param[in] component    - set value
   */
  inline void SetSlidingState(unsigned short val_marker,
                              unsigned long val_vertex,
                              unsigned short val_state,
                              unsigned long donor_index,
                              su2double component) final {
    SlidingState[val_marker][val_vertex][val_state][donor_index] = component;
  }

  /*!
   * \brief Set the number of outer state for fluid interface nodes.
   * \param[in] val_marker - marker index
   * \param[in] val_vertex - vertex index
   * \param[in] value - number of outer states
   */
  inline void SetnSlidingStates(unsigned short val_marker,
                                unsigned long val_vertex,
                                int value) final {
    SlidingStateNodes[val_marker][val_vertex] = value;
  }

  /*!
   * \brief Get the number of outer state for fluid interface nodes.
   * \param[in] val_marker - marker index
   * \param[in] val_vertex - vertex index
   */
  inline int GetnSlidingStates(unsigned short val_marker, unsigned long val_vertex) const final{
    return SlidingStateNodes[val_marker][val_vertex];
  }

  /*!
   * \brief Get the outer state for fluid interface nodes.
   * \param[in] val_marker - marker index
   * \param[in] val_vertex - vertex index
   * \param[in] val_state  - requested state component
   * \param[in] donor_index- index of the donor node to get
   */
  inline su2double GetSlidingState(unsigned short val_marker,
                                   unsigned long val_vertex,
                                   unsigned short val_state,
                                   unsigned long donor_index) const final {
    return SlidingState[val_marker][val_vertex][val_state][donor_index];
  }

};

template<class V, ENUM_REGIME R>
void CFVMFlowSolverBase<V,R>::AeroCoeffsArray::allocate(int size) {
  _size = size;
  CD = new su2double[size]; CL = new su2double[size]; CSF = new su2double[size]; CEff = new su2double[size];
  CFx = new su2double[size]; CFy = new su2double[size]; CFz = new su2double[size]; CMx = new su2double[size];
  CMy = new su2double[size]; CMz = new su2double[size]; CoPx = new su2double[size]; CoPy = new su2double[size];
  CoPz = new su2double[size]; CT = new su2double[size]; CQ = new su2double[size]; CMerit = new su2double[size];
  setZero();
}

template<class V, ENUM_REGIME R>
CFVMFlowSolverBase<V,R>::AeroCoeffsArray::~AeroCoeffsArray() {
  delete [] CD; delete [] CL; delete [] CSF; delete [] CEff;
  delete [] CFx; delete [] CFy; delete [] CFz; delete [] CMx;
  delete [] CMy; delete [] CMz; delete [] CoPx; delete [] CoPy;
  delete [] CoPz; delete [] CT; delete [] CQ; delete [] CMerit;
}

template<class V, ENUM_REGIME R>
void CFVMFlowSolverBase<V,R>::AeroCoeffsArray::setZero(int i) {
  CD[i] = CL[i] = CSF[i] = CEff[i] = 0.0;
  CFx[i] = CFy[i] = CFz[i] = CMx[i] = 0.0;
  CMy[i] = CMz[i] = CoPx[i] = CoPy[i] = 0.0;
  CoPz[i] = CT[i] = CQ[i] = CMerit[i] = 0.0;
}

template<class V, ENUM_REGIME FlowRegime>
void CFVMFlowSolverBase<V,FlowRegime>::Pressure_Forces(const CGeometry* geometry, const CConfig* config) {

  unsigned long iVertex, iPoint;
  unsigned short iDim, iMarker, Boundary, Monitoring, iMarker_Monitoring;
  su2double Pressure = 0.0, factor, NFPressOF, RefVel2 = 0.0,
  RefTemp, RefDensity = 0.0, RefPressure, Mach2Vel, Mach_Motion;
  const su2double *Normal = nullptr, *Coord = nullptr;
  string Marker_Tag, Monitoring_Tag;
  su2double AxiFactor;

  su2double Alpha = config->GetAoA()*PI_NUMBER/180.0;
  su2double Beta = config->GetAoS()*PI_NUMBER/180.0;
  su2double RefArea = config->GetRefArea();
  su2double RefLength = config->GetRefLength();
  su2double Gas_Constant = config->GetGas_ConstantND();
  const su2double *Origin = nullptr;
  if (config->GetnMarker_Monitoring() != 0){
    Origin = config->GetRefOriginMoment(0);
  }
  bool axisymmetric = config->GetAxisymmetric();

  /// TODO: Move these ifs to specialized functions.

  if (FlowRegime == COMPRESSIBLE) {
    /*--- Evaluate reference values for non-dimensionalization.
     For dynamic meshes, use the motion Mach number as a reference value
     for computing the force coefficients. Otherwise, use the freestream values,
     which is the standard convention. ---*/

    RefTemp = Temperature_Inf;
    RefDensity = Density_Inf;
    if (dynamic_grid) {
      Mach2Vel = sqrt(Gamma*Gas_Constant*RefTemp);
      Mach_Motion = config->GetMach_Motion();
      RefVel2 = (Mach_Motion*Mach2Vel)*(Mach_Motion*Mach2Vel);
    } else {
      RefVel2 = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        RefVel2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
    }
  }

  if (FlowRegime == INCOMPRESSIBLE) {
    /*--- Evaluate reference values for non-dimensionalization.
     For dimensional or non-dim based on initial values, use
     the far-field state (inf). For a custom non-dim based
     on user-provided reference values, use the ref values
     to compute the forces. ---*/

    if ((config->GetRef_Inc_NonDim() == DIMENSIONAL) ||
        (config->GetRef_Inc_NonDim() == INITIAL_VALUES)) {
      RefDensity  = Density_Inf;
      RefVel2 = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        RefVel2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
    }
    else if (config->GetRef_Inc_NonDim() == REFERENCE_VALUES) {
      RefDensity = config->GetInc_Density_Ref();
      RefVel2    = config->GetInc_Velocity_Ref()*config->GetInc_Velocity_Ref();
    }
  }

  factor = 1.0 / (0.5*RefDensity*RefArea*RefVel2);

  /*--- Reference pressure is always the far-field value. ---*/

  RefPressure = Pressure_Inf;

  /*-- Variables initialization ---*/

  TotalCoeff.setZero();

  Total_CNearFieldOF = 0.0; Total_Heat = 0.0;  Total_MaxHeat = 0.0;

  AllBoundInvCoeff.setZero();

  AllBound_CNearFieldOF_Inv = 0.0;

  SurfaceInvCoeff.setZero();
  SurfaceCoeff.setZero();

  /*--- Loop over the Euler and Navier-Stokes markers ---*/

  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    Boundary   = config->GetMarker_All_KindBC(iMarker);
    Monitoring = config->GetMarker_All_Monitoring(iMarker);

    /*--- Obtain the origin for the moment computation for a particular marker ---*/

    if (Monitoring == YES) {
      for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
        Monitoring_Tag = config->GetMarker_Monitoring_TagBound(iMarker_Monitoring);
        Marker_Tag = config->GetMarker_All_TagBound(iMarker);
        if (Marker_Tag == Monitoring_Tag)
          Origin = config->GetRefOriginMoment(iMarker_Monitoring);
      }
    }

    if (config->GetSolid_Wall(iMarker) || (Boundary == NEARFIELD_BOUNDARY) ||
        (Boundary == INLET_FLOW) || (Boundary == OUTLET_FLOW) ||
        (Boundary == ACTDISK_INLET) || (Boundary == ACTDISK_OUTLET)||
        (Boundary == ENGINE_INFLOW) || (Boundary == ENGINE_EXHAUST)) {

      /*--- Forces initialization at each Marker ---*/

      InvCoeff.setZero(iMarker);

      CNearFieldOF_Inv[iMarker] = 0.0;

      su2double ForceInviscid[MAXNDIM] = {0.0}, MomentInviscid[MAXNDIM] = {0.0};
      su2double MomentX_Force[MAXNDIM] = {0.0}, MomentY_Force[MAXNDIM] = {0.0}, MomentZ_Force[MAXNDIM] = {0.0};

      NFPressOF = 0.0;

      /*--- Loop over the vertices to compute the forces ---*/

      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {

        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

        Pressure = nodes->GetPressure(iPoint);

        CPressure[iMarker][iVertex] = (Pressure - RefPressure)*factor*RefArea;

        /*--- Note that the pressure coefficient is computed at the
         halo cells (for visualization purposes), but not the forces ---*/

        if ( (geometry->nodes->GetDomain(iPoint)) && (Monitoring == YES) ) {

          Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
          Coord = geometry->nodes->GetCoord(iPoint);

          /*--- Quadratic objective function for the near-field.
           This uses the infinity pressure regardless of Mach number. ---*/

          NFPressOF += 0.5*(Pressure - Pressure_Inf)*(Pressure - Pressure_Inf)*Normal[nDim-1];

          su2double MomentDist[MAXNDIM] = {0.0};
          for (iDim = 0; iDim < nDim; iDim++) {
            MomentDist[iDim] = Coord[iDim] - Origin[iDim];
          }

          /*--- Axisymmetric simulations ---*/

          if (axisymmetric) AxiFactor = 2.0*PI_NUMBER*geometry->nodes->GetCoord(iPoint, 1);
          else AxiFactor = 1.0;

          /*--- Force computation, note the minus sign due to the
           orientation of the normal (outward) ---*/

          su2double Force[MAXNDIM] = {0.0};
          for (iDim = 0; iDim < nDim; iDim++) {
            Force[iDim] = -(Pressure - Pressure_Inf) * Normal[iDim] * factor * AxiFactor;
            ForceInviscid[iDim] += Force[iDim];
          }

          /*--- Moment with respect to the reference axis ---*/

          if (nDim == 3) {
            MomentInviscid[0] += (Force[2]*MomentDist[1]-Force[1]*MomentDist[2])/RefLength;
            MomentX_Force[1]  += (-Force[1]*Coord[2]);
            MomentX_Force[2]  += (Force[2]*Coord[1]);

            MomentInviscid[1] += (Force[0]*MomentDist[2]-Force[2]*MomentDist[0])/RefLength;
            MomentY_Force[2]  += (-Force[2]*Coord[0]);
            MomentY_Force[0]  += (Force[0]*Coord[2]);
          }
          MomentInviscid[2] += (Force[1]*MomentDist[0]-Force[0]*MomentDist[1])/RefLength;
          MomentZ_Force[0]  += (-Force[0]*Coord[1]);
          MomentZ_Force[1]  += (Force[1]*Coord[0]);
        }

      }

      /*--- Project forces and store the non-dimensional coefficients ---*/

      if (Monitoring == YES) {

        if (Boundary != NEARFIELD_BOUNDARY) {
          if (nDim == 2) {
            InvCoeff.CD[iMarker]     =  ForceInviscid[0]*cos(Alpha) + ForceInviscid[1]*sin(Alpha);
            InvCoeff.CL[iMarker]     = -ForceInviscid[0]*sin(Alpha) + ForceInviscid[1]*cos(Alpha);
            InvCoeff.CEff[iMarker]   = InvCoeff.CL[iMarker] / (InvCoeff.CD[iMarker]+EPS);
            InvCoeff.CMz[iMarker]    = MomentInviscid[2];
            InvCoeff.CoPx[iMarker]   = MomentZ_Force[1];
            InvCoeff.CoPy[iMarker]   = -MomentZ_Force[0];
            InvCoeff.CFx[iMarker]    = ForceInviscid[0];
            InvCoeff.CFy[iMarker]    = ForceInviscid[1];
            InvCoeff.CT[iMarker]     = -InvCoeff.CFx[iMarker];
            InvCoeff.CQ[iMarker]     = -InvCoeff.CMz[iMarker];
            InvCoeff.CMerit[iMarker] = InvCoeff.CT[iMarker] / (InvCoeff.CQ[iMarker] + EPS);
          }
          if (nDim == 3) {
            InvCoeff.CD[iMarker]      =  ForceInviscid[0]*cos(Alpha)*cos(Beta) + ForceInviscid[1]*sin(Beta) + ForceInviscid[2]*sin(Alpha)*cos(Beta);
            InvCoeff.CL[iMarker]      = -ForceInviscid[0]*sin(Alpha) + ForceInviscid[2]*cos(Alpha);
            InvCoeff.CSF[iMarker]     = -ForceInviscid[0]*sin(Beta)*cos(Alpha) + ForceInviscid[1]*cos(Beta) - ForceInviscid[2]*sin(Beta)*sin(Alpha);
            InvCoeff.CEff[iMarker]    = InvCoeff.CL[iMarker] / (InvCoeff.CD[iMarker] + EPS);
            InvCoeff.CMx[iMarker]     = MomentInviscid[0];
            InvCoeff.CMy[iMarker]     = MomentInviscid[1];
            InvCoeff.CMz[iMarker]     = MomentInviscid[2];
            InvCoeff.CoPx[iMarker]    = -MomentY_Force[0];
            InvCoeff.CoPz[iMarker]    = MomentY_Force[2];
            InvCoeff.CFx[iMarker]     = ForceInviscid[0];
            InvCoeff.CFy[iMarker]     = ForceInviscid[1];
            InvCoeff.CFz[iMarker]     = ForceInviscid[2];
            InvCoeff.CT[iMarker]      = -InvCoeff.CFz[iMarker];
            InvCoeff.CQ[iMarker]      = -InvCoeff.CMz[iMarker];
            InvCoeff.CMerit[iMarker]  = InvCoeff.CT[iMarker] / (InvCoeff.CQ[iMarker] + EPS);
          }

          AllBoundInvCoeff.CD           += InvCoeff.CD[iMarker];
          AllBoundInvCoeff.CL           += InvCoeff.CL[iMarker];
          AllBoundInvCoeff.CSF          += InvCoeff.CSF[iMarker];
          AllBoundInvCoeff.CEff          = AllBoundInvCoeff.CL / (AllBoundInvCoeff.CD + EPS);
          AllBoundInvCoeff.CMx          += InvCoeff.CMx[iMarker];
          AllBoundInvCoeff.CMy          += InvCoeff.CMy[iMarker];
          AllBoundInvCoeff.CMz          += InvCoeff.CMz[iMarker];
          AllBoundInvCoeff.CoPx         += InvCoeff.CoPx[iMarker];
          AllBoundInvCoeff.CoPy         += InvCoeff.CoPy[iMarker];
          AllBoundInvCoeff.CoPz         += InvCoeff.CoPz[iMarker];
          AllBoundInvCoeff.CFx          += InvCoeff.CFx[iMarker];
          AllBoundInvCoeff.CFy          += InvCoeff.CFy[iMarker];
          AllBoundInvCoeff.CFz          += InvCoeff.CFz[iMarker];
          AllBoundInvCoeff.CT           += InvCoeff.CT[iMarker];
          AllBoundInvCoeff.CQ           += InvCoeff.CQ[iMarker];
          AllBoundInvCoeff.CMerit        = AllBoundInvCoeff.CT / (AllBoundInvCoeff.CQ + EPS);

          /*--- Compute the coefficients per surface ---*/

          for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
            Monitoring_Tag = config->GetMarker_Monitoring_TagBound(iMarker_Monitoring);
            Marker_Tag = config->GetMarker_All_TagBound(iMarker);
            if (Marker_Tag == Monitoring_Tag) {
              SurfaceInvCoeff.CL[iMarker_Monitoring]      += InvCoeff.CL[iMarker];
              SurfaceInvCoeff.CD[iMarker_Monitoring]      += InvCoeff.CD[iMarker];
              SurfaceInvCoeff.CSF[iMarker_Monitoring]     += InvCoeff.CSF[iMarker];
              SurfaceInvCoeff.CEff[iMarker_Monitoring]     = InvCoeff.CL[iMarker] / (InvCoeff.CD[iMarker] + EPS);
              SurfaceInvCoeff.CFx[iMarker_Monitoring]     += InvCoeff.CFx[iMarker];
              SurfaceInvCoeff.CFy[iMarker_Monitoring]     += InvCoeff.CFy[iMarker];
              SurfaceInvCoeff.CFz[iMarker_Monitoring]     += InvCoeff.CFz[iMarker];
              SurfaceInvCoeff.CMx[iMarker_Monitoring]     += InvCoeff.CMx[iMarker];
              SurfaceInvCoeff.CMy[iMarker_Monitoring]     += InvCoeff.CMy[iMarker];
              SurfaceInvCoeff.CMz[iMarker_Monitoring]     += InvCoeff.CMz[iMarker];
            }
          }

        }

        /*--- At the Nearfield SU2 only cares about the pressure coeffient ---*/

        else {
          CNearFieldOF_Inv[iMarker] = NFPressOF;
          AllBound_CNearFieldOF_Inv += CNearFieldOF_Inv[iMarker];
        }

      }

    }
  }

#ifdef HAVE_MPI

  /*--- Add AllBound information using all the nodes ---*/

  if (config->GetComm_Level() == COMM_FULL) {

    auto Allreduce = [](su2double x) {
      su2double tmp = x; x = 0.0;
      SU2_MPI::Allreduce(&tmp, &x, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      return x;
    };
    AllBoundInvCoeff.CD = Allreduce(AllBoundInvCoeff.CD);
    AllBoundInvCoeff.CL = Allreduce(AllBoundInvCoeff.CL);
    AllBoundInvCoeff.CSF = Allreduce(AllBoundInvCoeff.CSF);
    AllBoundInvCoeff.CEff = AllBoundInvCoeff.CL / (AllBoundInvCoeff.CD + EPS);

    AllBoundInvCoeff.CMx = Allreduce(AllBoundInvCoeff.CMx);
    AllBoundInvCoeff.CMy = Allreduce(AllBoundInvCoeff.CMy);
    AllBoundInvCoeff.CMz = Allreduce(AllBoundInvCoeff.CMz);

    AllBoundInvCoeff.CoPx = Allreduce(AllBoundInvCoeff.CoPx);
    AllBoundInvCoeff.CoPy = Allreduce(AllBoundInvCoeff.CoPy);
    AllBoundInvCoeff.CoPz = Allreduce(AllBoundInvCoeff.CoPz);

    AllBoundInvCoeff.CFx = Allreduce(AllBoundInvCoeff.CFx);
    AllBoundInvCoeff.CFy = Allreduce(AllBoundInvCoeff.CFy);
    AllBoundInvCoeff.CFz = Allreduce(AllBoundInvCoeff.CFz);

    AllBoundInvCoeff.CT = Allreduce(AllBoundInvCoeff.CT);
    AllBoundInvCoeff.CQ = Allreduce(AllBoundInvCoeff.CQ);
    AllBoundInvCoeff.CMerit = AllBoundInvCoeff.CT / (AllBoundInvCoeff.CQ + EPS);
    AllBound_CNearFieldOF_Inv = Allreduce(AllBound_CNearFieldOF_Inv);

  }

  /*--- Add the forces on the surfaces using all the nodes ---*/

  if (config->GetComm_Level() == COMM_FULL) {

    int nMarkerMon = config->GetnMarker_Monitoring();

    /*--- Use the same buffer for all reductions. We could avoid the copy back into
     *    the original variable by swaping pointers, but it is safer this way... ---*/

    su2double* buffer = new su2double [nMarkerMon];

    auto Allreduce_inplace = [buffer](int size, su2double* x) {
      SU2_MPI::Allreduce(x, buffer, size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      for(int i=0; i<size; ++i) x[i] = buffer[i];
    };

    Allreduce_inplace(nMarkerMon, SurfaceInvCoeff.CL);
    Allreduce_inplace(nMarkerMon, SurfaceInvCoeff.CD);
    Allreduce_inplace(nMarkerMon, SurfaceInvCoeff.CSF);

    for (iMarker_Monitoring = 0; iMarker_Monitoring < nMarkerMon; iMarker_Monitoring++)
      SurfaceInvCoeff.CEff[iMarker_Monitoring] = SurfaceInvCoeff.CL[iMarker_Monitoring] / (SurfaceInvCoeff.CD[iMarker_Monitoring] + EPS);

    Allreduce_inplace(nMarkerMon, SurfaceInvCoeff.CFx);
    Allreduce_inplace(nMarkerMon, SurfaceInvCoeff.CFy);
    Allreduce_inplace(nMarkerMon, SurfaceInvCoeff.CFz);

    Allreduce_inplace(nMarkerMon, SurfaceInvCoeff.CMx);
    Allreduce_inplace(nMarkerMon, SurfaceInvCoeff.CMy);
    Allreduce_inplace(nMarkerMon, SurfaceInvCoeff.CMz);

    delete [] buffer;

  }

#endif

  /*--- Update the total coefficients (note that all the nodes have the same value) ---*/

  TotalCoeff.CD            = AllBoundInvCoeff.CD;
  TotalCoeff.CL            = AllBoundInvCoeff.CL;
  TotalCoeff.CSF           = AllBoundInvCoeff.CSF;
  TotalCoeff.CEff          = TotalCoeff.CL / (TotalCoeff.CD + EPS);
  TotalCoeff.CFx           = AllBoundInvCoeff.CFx;
  TotalCoeff.CFy           = AllBoundInvCoeff.CFy;
  TotalCoeff.CFz           = AllBoundInvCoeff.CFz;
  TotalCoeff.CMx           = AllBoundInvCoeff.CMx;
  TotalCoeff.CMy           = AllBoundInvCoeff.CMy;
  TotalCoeff.CMz           = AllBoundInvCoeff.CMz;
  TotalCoeff.CoPx          = AllBoundInvCoeff.CoPx;
  TotalCoeff.CoPy          = AllBoundInvCoeff.CoPy;
  TotalCoeff.CoPz          = AllBoundInvCoeff.CoPz;
  TotalCoeff.CT            = AllBoundInvCoeff.CT;
  TotalCoeff.CQ            = AllBoundInvCoeff.CQ;
  TotalCoeff.CMerit        = TotalCoeff.CT / (TotalCoeff.CQ + EPS);
  Total_CNearFieldOF       = AllBound_CNearFieldOF_Inv;

  /*--- Update the total coefficients per surface (note that all the nodes have the same value)---*/

  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    SurfaceCoeff.CL[iMarker_Monitoring]      = SurfaceInvCoeff.CL[iMarker_Monitoring];
    SurfaceCoeff.CD[iMarker_Monitoring]      = SurfaceInvCoeff.CD[iMarker_Monitoring];
    SurfaceCoeff.CSF[iMarker_Monitoring]     = SurfaceInvCoeff.CSF[iMarker_Monitoring];
    SurfaceCoeff.CEff[iMarker_Monitoring]    = SurfaceInvCoeff.CL[iMarker_Monitoring] / (SurfaceInvCoeff.CD[iMarker_Monitoring] + EPS);
    SurfaceCoeff.CFx[iMarker_Monitoring]     = SurfaceInvCoeff.CFx[iMarker_Monitoring];
    SurfaceCoeff.CFy[iMarker_Monitoring]     = SurfaceInvCoeff.CFy[iMarker_Monitoring];
    SurfaceCoeff.CFz[iMarker_Monitoring]     = SurfaceInvCoeff.CFz[iMarker_Monitoring];
    SurfaceCoeff.CMx[iMarker_Monitoring]     = SurfaceInvCoeff.CMx[iMarker_Monitoring];
    SurfaceCoeff.CMy[iMarker_Monitoring]     = SurfaceInvCoeff.CMy[iMarker_Monitoring];
    SurfaceCoeff.CMz[iMarker_Monitoring]     = SurfaceInvCoeff.CMz[iMarker_Monitoring];
  }

}

template<class V, ENUM_REGIME FlowRegime>
void CFVMFlowSolverBase<V,FlowRegime>::Momentum_Forces(const CGeometry* geometry, const CConfig* config) {

  unsigned long iVertex, iPoint;
  unsigned short iDim, iMarker, Boundary, Monitoring, iMarker_Monitoring;
  su2double Area, factor, RefVel2 = 0.0, RefTemp, RefDensity = 0.0,  Mach2Vel, Mach_Motion, MassFlow, Density;
  const su2double *Normal = nullptr, *Coord = nullptr;
  string Marker_Tag, Monitoring_Tag;
  su2double AxiFactor;

  su2double Alpha = config->GetAoA()*PI_NUMBER/180.0;
  su2double Beta = config->GetAoS()*PI_NUMBER/180.0;
  su2double RefArea = config->GetRefArea();
  su2double RefLength = config->GetRefLength();
  su2double Gas_Constant = config->GetGas_ConstantND();
  const su2double *Origin = nullptr;
  if (config->GetnMarker_Monitoring() != 0){
    Origin = config->GetRefOriginMoment(0);
  }
  bool axisymmetric          = config->GetAxisymmetric();

  /// TODO: Move these ifs to specialized functions.

  if (FlowRegime == COMPRESSIBLE) {
    /*--- Evaluate reference values for non-dimensionalization.
     For dynamic meshes, use the motion Mach number as a reference value
     for computing the force coefficients. Otherwise, use the freestream values,
     which is the standard convention. ---*/

    RefTemp = Temperature_Inf;
    RefDensity = Density_Inf;
    if (dynamic_grid) {
      Mach2Vel = sqrt(Gamma*Gas_Constant*RefTemp);
      Mach_Motion = config->GetMach_Motion();
      RefVel2 = (Mach_Motion*Mach2Vel)*(Mach_Motion*Mach2Vel);
    } else {
      RefVel2 = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        RefVel2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
    }
  }

  if (FlowRegime == INCOMPRESSIBLE) {
    /*--- Evaluate reference values for non-dimensionalization.
     For dimensional or non-dim based on initial values, use
     the far-field state (inf). For a custom non-dim based
     on user-provided reference values, use the ref values
     to compute the forces. ---*/

    if ((config->GetRef_Inc_NonDim() == DIMENSIONAL) ||
        (config->GetRef_Inc_NonDim() == INITIAL_VALUES)) {
      RefDensity  = Density_Inf;
      RefVel2 = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        RefVel2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
    }
    else if (config->GetRef_Inc_NonDim() == REFERENCE_VALUES) {
      RefDensity = config->GetInc_Density_Ref();
      RefVel2    = config->GetInc_Velocity_Ref()*config->GetInc_Velocity_Ref();
    }
  }

  factor = 1.0 / (0.5*RefDensity*RefArea*RefVel2);

  /*-- Variables initialization ---*/

  AllBoundMntCoeff.setZero();
  SurfaceMntCoeff.setZero();

  /*--- Loop over the Inlet -Outlet Markers  ---*/

  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    Boundary   = config->GetMarker_All_KindBC(iMarker);
    Monitoring = config->GetMarker_All_Monitoring(iMarker);

    /*--- Obtain the origin for the moment computation for a particular marker ---*/

    if (Monitoring == YES) {
      for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
        Monitoring_Tag = config->GetMarker_Monitoring_TagBound(iMarker_Monitoring);
        Marker_Tag = config->GetMarker_All_TagBound(iMarker);
        if (Marker_Tag == Monitoring_Tag)
          Origin = config->GetRefOriginMoment(iMarker_Monitoring);
      }
    }

    if ((Boundary == INLET_FLOW) || (Boundary == OUTLET_FLOW) ||
        (Boundary == ACTDISK_INLET) || (Boundary == ACTDISK_OUTLET)||
        (Boundary == ENGINE_INFLOW) || (Boundary == ENGINE_EXHAUST)) {

      /*--- Forces initialization at each Marker ---*/

      MntCoeff.setZero(iMarker);

      su2double ForceMomentum[MAXNDIM] = {0.0}, MomentMomentum[MAXNDIM] = {0.0};
      su2double MomentX_Force[3] = {0.0}, MomentY_Force[3] = {0.0}, MomentZ_Force[3] = {0.0};

      /*--- Loop over the vertices to compute the forces ---*/

      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {

        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

        /*--- Note that the pressure coefficient is computed at the
         halo cells (for visualization purposes), but not the forces ---*/

        if ( (geometry->nodes->GetDomain(iPoint)) && (Monitoring == YES) ) {

          Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
          Coord = geometry->nodes->GetCoord(iPoint);
          Density   = nodes->GetDensity(iPoint);

          Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);

          MassFlow = 0.0;
          su2double Velocity[MAXNDIM] = {0.0}, MomentDist[MAXNDIM] = {0.0};
          for (iDim = 0; iDim < nDim; iDim++) {
            Velocity[iDim]  = nodes->GetVelocity(iPoint,iDim);
            MomentDist[iDim] = Coord[iDim] - Origin[iDim];
            MassFlow -= Normal[iDim]*Velocity[iDim]*Density;
          }

          /*--- Axisymmetric simulations ---*/

          if (axisymmetric) AxiFactor = 2.0*PI_NUMBER*geometry->nodes->GetCoord(iPoint, 1);
          else AxiFactor = 1.0;

          /*--- Force computation, note the minus sign due to the
           orientation of the normal (outward) ---*/

          su2double Force[MAXNDIM] = {0.0};
          for (iDim = 0; iDim < nDim; iDim++) {
            Force[iDim] = MassFlow * Velocity[iDim] * factor * AxiFactor;
            ForceMomentum[iDim] += Force[iDim];
          }

          /*--- Moment with respect to the reference axis ---*/

          if (iDim == 3) {
            MomentMomentum[0] += (Force[2]*MomentDist[1]-Force[1]*MomentDist[2])/RefLength;
            MomentX_Force[1]  += (-Force[1]*Coord[2]);
            MomentX_Force[2]  += (Force[2]*Coord[1]);

            MomentMomentum[1] += (Force[0]*MomentDist[2]-Force[2]*MomentDist[0])/RefLength;
            MomentY_Force[2]  += (-Force[2]*Coord[0]);
            MomentY_Force[0]  += (Force[0]*Coord[2]);
          }
          MomentMomentum[2] += (Force[1]*MomentDist[0]-Force[0]*MomentDist[1])/RefLength;
          MomentZ_Force[0]  += (-Force[0]*Coord[1]);
          MomentZ_Force[1]  += (Force[1]*Coord[0]);

        }

      }

      /*--- Project forces and store the non-dimensional coefficients ---*/

      if (Monitoring == YES) {

        if (nDim == 2) {
          MntCoeff.CD[iMarker]     =  ForceMomentum[0]*cos(Alpha) + ForceMomentum[1]*sin(Alpha);
          MntCoeff.CL[iMarker]     = -ForceMomentum[0]*sin(Alpha) + ForceMomentum[1]*cos(Alpha);
          MntCoeff.CEff[iMarker]   = MntCoeff.CL[iMarker] / (MntCoeff.CD[iMarker]+EPS);
          MntCoeff.CFx[iMarker]    = ForceMomentum[0];
          MntCoeff.CFy[iMarker]    = ForceMomentum[1];
          MntCoeff.CMz[iMarker]    = MomentMomentum[2];
          MntCoeff.CoPx[iMarker]   = MomentZ_Force[1];
          MntCoeff.CoPy[iMarker]   = -MomentZ_Force[0];
          MntCoeff.CT[iMarker]     = -MntCoeff.CFx[iMarker];
          MntCoeff.CQ[iMarker]     = -MntCoeff.CMz[iMarker];
          MntCoeff.CMerit[iMarker] = MntCoeff.CT[iMarker] / (MntCoeff.CQ[iMarker] + EPS);
        }
        if (nDim == 3) {
          MntCoeff.CD[iMarker]         =  ForceMomentum[0]*cos(Alpha)*cos(Beta) + ForceMomentum[1]*sin(Beta) + ForceMomentum[2]*sin(Alpha)*cos(Beta);
          MntCoeff.CL[iMarker]         = -ForceMomentum[0]*sin(Alpha) + ForceMomentum[2]*cos(Alpha);
          MntCoeff.CSF[iMarker]        = -ForceMomentum[0]*sin(Beta)*cos(Alpha) + ForceMomentum[1]*cos(Beta) - ForceMomentum[2]*sin(Beta)*sin(Alpha);
          MntCoeff.CEff[iMarker]       = MntCoeff.CL[iMarker] / (MntCoeff.CD[iMarker] + EPS);
          MntCoeff.CFx[iMarker]        = ForceMomentum[0];
          MntCoeff.CFy[iMarker]        = ForceMomentum[1];
          MntCoeff.CFz[iMarker]        = ForceMomentum[2];
          MntCoeff.CMx[iMarker]        = MomentMomentum[0];
          MntCoeff.CMy[iMarker]        = MomentMomentum[1];
          MntCoeff.CMz[iMarker]        = MomentMomentum[2];
          MntCoeff.CoPx[iMarker]       = -MomentY_Force[0];
          MntCoeff.CoPz[iMarker]       =  MomentY_Force[2];
          MntCoeff.CT[iMarker]         = -MntCoeff.CFz[iMarker];
          MntCoeff.CQ[iMarker]         = -MntCoeff.CMz[iMarker];
          MntCoeff.CMerit[iMarker]     = MntCoeff.CT[iMarker] / (MntCoeff.CQ[iMarker] + EPS);
        }

        AllBoundMntCoeff.CD           += MntCoeff.CD[iMarker];
        AllBoundMntCoeff.CL           += MntCoeff.CL[iMarker];
        AllBoundMntCoeff.CSF          += MntCoeff.CSF[iMarker];
        AllBoundMntCoeff.CEff          = AllBoundMntCoeff.CL / (AllBoundMntCoeff.CD + EPS);
        AllBoundMntCoeff.CFx          += MntCoeff.CFx[iMarker];
        AllBoundMntCoeff.CFy          += MntCoeff.CFy[iMarker];
        AllBoundMntCoeff.CFz          += MntCoeff.CFz[iMarker];
        AllBoundMntCoeff.CMx          += MntCoeff.CMx[iMarker];
        AllBoundMntCoeff.CMy          += MntCoeff.CMy[iMarker];
        AllBoundMntCoeff.CMx          += MntCoeff.CMz[iMarker];
        AllBoundMntCoeff.CoPx         += MntCoeff.CoPx[iMarker];
        AllBoundMntCoeff.CoPy         += MntCoeff.CoPy[iMarker];
        AllBoundMntCoeff.CoPz         += MntCoeff.CoPz[iMarker];
        AllBoundMntCoeff.CT           += MntCoeff.CT[iMarker];
        AllBoundMntCoeff.CQ           += MntCoeff.CQ[iMarker];
        AllBoundMntCoeff.CMerit       += AllBoundMntCoeff.CT / (AllBoundMntCoeff.CQ + EPS);

        /*--- Compute the coefficients per surface ---*/

        for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
          Monitoring_Tag = config->GetMarker_Monitoring_TagBound(iMarker_Monitoring);
          Marker_Tag = config->GetMarker_All_TagBound(iMarker);
          if (Marker_Tag == Monitoring_Tag) {
            SurfaceMntCoeff.CL[iMarker_Monitoring]      += MntCoeff.CL[iMarker];
            SurfaceMntCoeff.CD[iMarker_Monitoring]      += MntCoeff.CD[iMarker];
            SurfaceMntCoeff.CSF[iMarker_Monitoring]     += MntCoeff.CSF[iMarker];
            SurfaceMntCoeff.CEff[iMarker_Monitoring]     = MntCoeff.CL[iMarker] / (MntCoeff.CD[iMarker] + EPS);
            SurfaceMntCoeff.CFx[iMarker_Monitoring]     += MntCoeff.CFx[iMarker];
            SurfaceMntCoeff.CFy[iMarker_Monitoring]     += MntCoeff.CFy[iMarker];
            SurfaceMntCoeff.CFz[iMarker_Monitoring]     += MntCoeff.CFz[iMarker];
            SurfaceMntCoeff.CMx[iMarker_Monitoring]     += MntCoeff.CMx[iMarker];
            SurfaceMntCoeff.CMy[iMarker_Monitoring]     += MntCoeff.CMy[iMarker];
            SurfaceMntCoeff.CMz[iMarker_Monitoring]     += MntCoeff.CMz[iMarker];
          }
        }

      }

    }
  }

#ifdef HAVE_MPI

  /*--- Add AllBound information using all the nodes ---*/

  if (config->GetComm_Level() == COMM_FULL) {

    auto Allreduce = [](su2double x) {
      su2double tmp = x; x = 0.0;
      SU2_MPI::Allreduce(&tmp, &x, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      return x;
    };

    AllBoundMntCoeff.CD = Allreduce(AllBoundMntCoeff.CD);
    AllBoundMntCoeff.CL = Allreduce(AllBoundMntCoeff.CL);
    AllBoundMntCoeff.CSF = Allreduce(AllBoundMntCoeff.CSF);
    AllBoundMntCoeff.CEff = AllBoundMntCoeff.CL / (AllBoundMntCoeff.CD + EPS);

    AllBoundMntCoeff.CFx = Allreduce(AllBoundMntCoeff.CFx);
    AllBoundMntCoeff.CFy = Allreduce(AllBoundMntCoeff.CFy);
    AllBoundMntCoeff.CFz = Allreduce(AllBoundMntCoeff.CFz);

    AllBoundMntCoeff.CMx = Allreduce(AllBoundMntCoeff.CMx);
    AllBoundMntCoeff.CMy = Allreduce(AllBoundMntCoeff.CMy);
    AllBoundMntCoeff.CMz = Allreduce(AllBoundMntCoeff.CMz);

    AllBoundMntCoeff.CoPx = Allreduce(AllBoundMntCoeff.CoPx);
    AllBoundMntCoeff.CoPy = Allreduce(AllBoundMntCoeff.CoPy);
    AllBoundMntCoeff.CoPz = Allreduce(AllBoundMntCoeff.CoPz);

    AllBoundMntCoeff.CT = Allreduce(AllBoundMntCoeff.CT);
    AllBoundMntCoeff.CQ = Allreduce(AllBoundMntCoeff.CQ);
    AllBoundMntCoeff.CMerit = AllBoundMntCoeff.CT / (AllBoundMntCoeff.CQ + EPS);

  }

  /*--- Add the forces on the surfaces using all the nodes ---*/

  if (config->GetComm_Level() == COMM_FULL) {

    int nMarkerMon = config->GetnMarker_Monitoring();

    /*--- Use the same buffer for all reductions. We could avoid the copy back into
     *    the original variable by swaping pointers, but it is safer this way... ---*/

    su2double* buffer = new su2double [nMarkerMon];

    auto Allreduce_inplace = [buffer](int size, su2double* x) {
      SU2_MPI::Allreduce(x, buffer, size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      for(int i=0; i<size; ++i) x[i] = buffer[i];
    };

    Allreduce_inplace(nMarkerMon, SurfaceMntCoeff.CL);
    Allreduce_inplace(nMarkerMon, SurfaceMntCoeff.CD);
    Allreduce_inplace(nMarkerMon, SurfaceMntCoeff.CSF);

    for (iMarker_Monitoring = 0; iMarker_Monitoring < nMarkerMon; iMarker_Monitoring++)
      SurfaceMntCoeff.CEff[iMarker_Monitoring] = SurfaceMntCoeff.CL[iMarker_Monitoring] / (SurfaceMntCoeff.CD[iMarker_Monitoring] + EPS);

    Allreduce_inplace(nMarkerMon, SurfaceMntCoeff.CFx);
    Allreduce_inplace(nMarkerMon, SurfaceMntCoeff.CFy);
    Allreduce_inplace(nMarkerMon, SurfaceMntCoeff.CFz);

    Allreduce_inplace(nMarkerMon, SurfaceMntCoeff.CMx);
    Allreduce_inplace(nMarkerMon, SurfaceMntCoeff.CMy);
    Allreduce_inplace(nMarkerMon, SurfaceMntCoeff.CMz);

    delete [] buffer;

  }

#endif

  /*--- Update the total coefficients (note that all the nodes have the same value) ---*/

  TotalCoeff.CD            += AllBoundMntCoeff.CD;
  TotalCoeff.CL            += AllBoundMntCoeff.CL;
  TotalCoeff.CSF           += AllBoundMntCoeff.CSF;
  TotalCoeff.CEff          = TotalCoeff.CL / (TotalCoeff.CD + EPS);
  TotalCoeff.CFx           += AllBoundMntCoeff.CFx;
  TotalCoeff.CFy           += AllBoundMntCoeff.CFy;
  TotalCoeff.CFz           += AllBoundMntCoeff.CFz;
  TotalCoeff.CMx           += AllBoundMntCoeff.CMx;
  TotalCoeff.CMy           += AllBoundMntCoeff.CMy;
  TotalCoeff.CMz           += AllBoundMntCoeff.CMz;
  TotalCoeff.CoPx          += AllBoundMntCoeff.CoPx;
  TotalCoeff.CoPy          += AllBoundMntCoeff.CoPy;
  TotalCoeff.CoPz          += AllBoundMntCoeff.CoPz;
  TotalCoeff.CT            += AllBoundMntCoeff.CT;
  TotalCoeff.CQ            += AllBoundMntCoeff.CQ;
  TotalCoeff.CMerit        = TotalCoeff.CT / (TotalCoeff.CQ + EPS);

  /*--- Update the total coefficients per surface (note that all the nodes have the same value)---*/

  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    SurfaceCoeff.CL[iMarker_Monitoring]         += SurfaceMntCoeff.CL[iMarker_Monitoring];
    SurfaceCoeff.CD[iMarker_Monitoring]         += SurfaceMntCoeff.CD[iMarker_Monitoring];
    SurfaceCoeff.CSF[iMarker_Monitoring]        += SurfaceMntCoeff.CSF[iMarker_Monitoring];
    SurfaceCoeff.CEff[iMarker_Monitoring]       += SurfaceMntCoeff.CL[iMarker_Monitoring] / (SurfaceMntCoeff.CD[iMarker_Monitoring] + EPS);
    SurfaceCoeff.CFx[iMarker_Monitoring]        += SurfaceMntCoeff.CFx[iMarker_Monitoring];
    SurfaceCoeff.CFy[iMarker_Monitoring]        += SurfaceMntCoeff.CFy[iMarker_Monitoring];
    SurfaceCoeff.CFz[iMarker_Monitoring]        += SurfaceMntCoeff.CFz[iMarker_Monitoring];
    SurfaceCoeff.CMx[iMarker_Monitoring]        += SurfaceMntCoeff.CMx[iMarker_Monitoring];
    SurfaceCoeff.CMy[iMarker_Monitoring]        += SurfaceMntCoeff.CMy[iMarker_Monitoring];
    SurfaceCoeff.CMz[iMarker_Monitoring]        += SurfaceMntCoeff.CMz[iMarker_Monitoring];
  }

}

template<class V, ENUM_REGIME FlowRegime>
void CFVMFlowSolverBase<V,FlowRegime>::Friction_Forces(const CGeometry* geometry, const CConfig* config) {

  /// TODO: Major cleanup needed.

  unsigned long iVertex, iPoint, iPointNormal;
  unsigned short Boundary, Monitoring, iMarker, iMarker_Monitoring, iDim, jDim;
  su2double Viscosity = 0.0, div_vel, WallDist[3] = {0.0},
  Area, WallShearStress, TauNormal, RefTemp, RefVel2 = 0.0, RefDensity = 0.0, GradTemperature, Density = 0.0, WallDistMod, FrictionVel,
  Mach2Vel, Mach_Motion, UnitNormal[3] = {0.0}, TauElem[3] = {0.0}, TauTangent[3] = {0.0},
  Tau[3][3] = {{0.0}}, Cp, thermal_conductivity, MaxNorm = 8.0,
  Grad_Vel[3][3] = {{0.0}}, Grad_Temp[3] = {0.0},
  delta[3][3] = {{1.0, 0.0, 0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};
  su2double AxiFactor;
  const su2double *Coord = nullptr, *Coord_Normal = nullptr, *Normal = nullptr;

  string Marker_Tag, Monitoring_Tag;

  su2double Alpha = config->GetAoA()*PI_NUMBER/180.0;
  su2double Beta = config->GetAoS()*PI_NUMBER/180.0;
  su2double RefArea = config->GetRefArea();
  su2double RefLength = config->GetRefLength();
  su2double RefHeatFlux = config->GetHeat_Flux_Ref();
  su2double Gas_Constant = config->GetGas_ConstantND();
  const su2double *Origin = nullptr;

  if (config->GetnMarker_Monitoring() != 0) { Origin = config->GetRefOriginMoment(0); }

  su2double Prandtl_Lam = config->GetPrandtl_Lam();
  bool energy = config->GetEnergy_Equation();
  bool QCR = config->GetQCR();
  bool axisymmetric = config->GetAxisymmetric();

  /// TODO: Move these ifs to specialized functions.

  if (FlowRegime == COMPRESSIBLE) {
    /*--- Evaluate reference values for non-dimensionalization.
     For dynamic meshes, use the motion Mach number as a reference value
     for computing the force coefficients. Otherwise, use the freestream values,
     which is the standard convention. ---*/

    RefTemp = Temperature_Inf;
    RefDensity = Density_Inf;
    if (dynamic_grid) {
      Mach2Vel = sqrt(Gamma*Gas_Constant*RefTemp);
      Mach_Motion = config->GetMach_Motion();
      RefVel2 = (Mach_Motion*Mach2Vel)*(Mach_Motion*Mach2Vel);
    } else {
      RefVel2 = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        RefVel2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
    }
  }

  if (FlowRegime == INCOMPRESSIBLE) {
    /*--- Evaluate reference values for non-dimensionalization.
     For dimensional or non-dim based on initial values, use
     the far-field state (inf). For a custom non-dim based
     on user-provided reference values, use the ref values
     to compute the forces. ---*/

    if ((config->GetRef_Inc_NonDim() == DIMENSIONAL) ||
        (config->GetRef_Inc_NonDim() == INITIAL_VALUES)) {
      RefDensity  = Density_Inf;
      RefVel2 = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        RefVel2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
    }
    else if (config->GetRef_Inc_NonDim() == REFERENCE_VALUES) {
      RefDensity = config->GetInc_Density_Ref();
      RefVel2    = config->GetInc_Velocity_Ref()*config->GetInc_Velocity_Ref();
    }
  }

  const su2double factor = 1.0 / (0.5*RefDensity*RefArea*RefVel2);

  /*--- Variables initialization ---*/

  AllBoundViscCoeff.setZero();
  SurfaceViscCoeff.setZero();

  AllBound_HF_Visc = 0.0;  AllBound_MaxHF_Visc = 0.0;

  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    Surface_HF_Visc[iMarker_Monitoring]  = 0.0; Surface_MaxHF_Visc[iMarker_Monitoring]   = 0.0;
  }

  /*--- Loop over the Navier-Stokes markers ---*/

  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    Boundary = config->GetMarker_All_KindBC(iMarker);
    Monitoring = config->GetMarker_All_Monitoring(iMarker);

    /*--- Obtain the origin for the moment computation for a particular marker ---*/

    if (Monitoring == YES) {
      for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
        Monitoring_Tag = config->GetMarker_Monitoring_TagBound(iMarker_Monitoring);
        Marker_Tag = config->GetMarker_All_TagBound(iMarker);
        if (Marker_Tag == Monitoring_Tag)
          Origin = config->GetRefOriginMoment(iMarker_Monitoring);
      }
    }

    if ((Boundary == HEAT_FLUX) || (Boundary == ISOTHERMAL) || (Boundary == CHT_WALL_INTERFACE)) {

      /*--- Forces initialization at each Marker ---*/

      ViscCoeff.setZero(iMarker);

      HF_Visc[iMarker] = 0.0;    MaxHF_Visc[iMarker] = 0.0;

      su2double ForceViscous[MAXNDIM] = {0.0}, MomentViscous[MAXNDIM] = {0.0};
      su2double MomentX_Force[MAXNDIM] = {0.0}, MomentY_Force[MAXNDIM] = {0.0}, MomentZ_Force[MAXNDIM] = {0.0};

      /*--- Loop over the vertices to compute the forces ---*/

      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {

        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        iPointNormal = geometry->vertex[iMarker][iVertex]->GetNormal_Neighbor();

        Coord = geometry->nodes->GetCoord(iPoint);
        Coord_Normal = geometry->nodes->GetCoord(iPointNormal);

        Normal = geometry->vertex[iMarker][iVertex]->GetNormal();

        for (iDim = 0; iDim < nDim; iDim++) {
          for (jDim = 0 ; jDim < nDim; jDim++) {
            Grad_Vel[iDim][jDim] = nodes->GetGradient_Primitive(iPoint,iDim+1, jDim);
          }

          /// TODO: Move the temperature index logic to a function.

          if (FlowRegime == COMPRESSIBLE)
            Grad_Temp[iDim] = nodes->GetGradient_Primitive(iPoint,0, iDim);

          if (FlowRegime == INCOMPRESSIBLE)
            Grad_Temp[iDim] = nodes->GetGradient_Primitive(iPoint,nDim+1, iDim);
        }

        Viscosity = nodes->GetLaminarViscosity(iPoint);
        Density = nodes->GetDensity(iPoint);

        Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);

        for (iDim = 0; iDim < nDim; iDim++) {
          UnitNormal[iDim] = Normal[iDim]/Area;
        }

        /*--- Evaluate Tau ---*/

        div_vel = 0.0; for (iDim = 0; iDim < nDim; iDim++) div_vel += Grad_Vel[iDim][iDim];

        for (iDim = 0; iDim < nDim; iDim++) {
          for (jDim = 0 ; jDim < nDim; jDim++) {
            Tau[iDim][jDim] = Viscosity*(Grad_Vel[jDim][iDim] + Grad_Vel[iDim][jDim]) -
                              TWO3*Viscosity*div_vel*delta[iDim][jDim];
          }
        }

        /*--- If necessary evaluate the QCR contribution to Tau ---*/

        if (QCR) {
          su2double den_aux, c_cr1=0.3, O_ik, O_jk;
          unsigned short kDim;

          /*--- Denominator Antisymmetric normalized rotation tensor ---*/

          den_aux = 0.0;
          for (iDim = 0 ; iDim < nDim; iDim++)
            for (jDim = 0 ; jDim < nDim; jDim++)
              den_aux += Grad_Vel[iDim][jDim] * Grad_Vel[iDim][jDim];
          den_aux = sqrt(max(den_aux,1E-10));

          /*--- Adding the QCR contribution ---*/

          su2double tauQCR[MAXNDIM][MAXNDIM] = {{0.0}};

          for (iDim = 0 ; iDim < nDim; iDim++){
            for (jDim = 0 ; jDim < nDim; jDim++){
              for (kDim = 0 ; kDim < nDim; kDim++){
                O_ik = (Grad_Vel[iDim][kDim] - Grad_Vel[kDim][iDim])/ den_aux;
                O_jk = (Grad_Vel[jDim][kDim] - Grad_Vel[kDim][jDim])/ den_aux;
                tauQCR[iDim][jDim] += O_ik * Tau[jDim][kDim] + O_jk * Tau[iDim][kDim];
              }
            }
          }

          for (iDim = 0 ; iDim < nDim; iDim++)
            for (jDim = 0 ; jDim < nDim; jDim++)
              Tau[iDim][jDim] -= c_cr1 * tauQCR[iDim][jDim];
        }

        /*--- Project Tau in each surface element ---*/

        for (iDim = 0; iDim < nDim; iDim++) {
          TauElem[iDim] = 0.0;
          for (jDim = 0; jDim < nDim; jDim++) {
            TauElem[iDim] += Tau[iDim][jDim]*UnitNormal[jDim];
          }
        }

        /*--- Compute wall shear stress (using the stress tensor). Compute wall skin friction coefficient, and heat flux on the wall ---*/

        TauNormal = 0.0; for (iDim = 0; iDim < nDim; iDim++) TauNormal += TauElem[iDim] * UnitNormal[iDim];

        WallShearStress = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          TauTangent[iDim] = TauElem[iDim] - TauNormal * UnitNormal[iDim];
          CSkinFriction[iMarker][iDim][iVertex] = TauTangent[iDim] / (0.5*RefDensity*RefVel2);
          WallShearStress += TauTangent[iDim] * TauTangent[iDim];
        }
        WallShearStress = sqrt(WallShearStress);

        for (iDim = 0; iDim < nDim; iDim++) WallDist[iDim] = (Coord[iDim] - Coord_Normal[iDim]);
        WallDistMod = 0.0; for (iDim = 0; iDim < nDim; iDim++) WallDistMod += WallDist[iDim]*WallDist[iDim]; WallDistMod = sqrt(WallDistMod);

        /*--- Compute y+ and non-dimensional velocity ---*/

        FrictionVel = sqrt(fabs(WallShearStress)/Density);
        YPlus[iMarker][iVertex] = WallDistMod*FrictionVel/(Viscosity/Density);

        /*--- Compute total and maximum heat flux on the wall ---*/

        GradTemperature = 0.0;

        /// TODO: Move these ifs to specialized functions.

        if (FlowRegime == COMPRESSIBLE) {

          for (iDim = 0; iDim < nDim; iDim++)
            GradTemperature -= Grad_Temp[iDim]*UnitNormal[iDim];

          Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
          thermal_conductivity = Cp * Viscosity/Prandtl_Lam;
        }

        if (FlowRegime == INCOMPRESSIBLE) {

          if (energy)
            for (iDim = 0; iDim < nDim; iDim++)
              GradTemperature -= Grad_Temp[iDim]*UnitNormal[iDim];

          thermal_conductivity = nodes->GetThermalConductivity(iPoint);
        }

        HeatFlux[iMarker][iVertex] = -thermal_conductivity*GradTemperature*RefHeatFlux;


        /*--- Note that y+, and heat are computed at the
         halo cells (for visualization purposes), but not the forces ---*/

        if ((geometry->nodes->GetDomain(iPoint)) && (Monitoring == YES)) {

          /*--- Axisymmetric simulations ---*/

          if (axisymmetric) AxiFactor = 2.0*PI_NUMBER*geometry->nodes->GetCoord(iPoint, 1);
          else AxiFactor = 1.0;

          /*--- Force computation ---*/

          su2double Force[MAXNDIM] = {0.0}, MomentDist[MAXNDIM] = {0.0};
          for (iDim = 0; iDim < nDim; iDim++) {
            Force[iDim] = TauElem[iDim] * Area * factor * AxiFactor;
            ForceViscous[iDim] += Force[iDim];
            MomentDist[iDim] = Coord[iDim] - Origin[iDim];
          }

          /*--- Moment with respect to the reference axis ---*/

          if (iDim == 3) {
            MomentViscous[0] += (Force[2]*MomentDist[1] - Force[1]*MomentDist[2])/RefLength;
            MomentX_Force[1] += (-Force[1]*Coord[2]);
            MomentX_Force[2] += (Force[2]*Coord[1]);

            MomentViscous[1] += (Force[0]*MomentDist[2] - Force[2]*MomentDist[0])/RefLength;
            MomentY_Force[2] += (-Force[2]*Coord[0]);
            MomentY_Force[0] += (Force[0]*Coord[2]);
          }
          MomentViscous[2] += (Force[1]*MomentDist[0] - Force[0]*MomentDist[1])/RefLength;
          MomentZ_Force[0] += (-Force[0]*Coord[1]);
          MomentZ_Force[1] += (Force[1]*Coord[0]);

          HF_Visc[iMarker]          += HeatFlux[iMarker][iVertex]*Area;
          MaxHF_Visc[iMarker]       += pow(HeatFlux[iMarker][iVertex], MaxNorm);
        }
      }

      /*--- Project forces and store the non-dimensional coefficients ---*/

      if (Monitoring == YES) {
        if (nDim == 2) {
          ViscCoeff.CD[iMarker]          =  ForceViscous[0]*cos(Alpha) + ForceViscous[1]*sin(Alpha);
          ViscCoeff.CL[iMarker]          = -ForceViscous[0]*sin(Alpha) + ForceViscous[1]*cos(Alpha);
          ViscCoeff.CEff[iMarker]        = ViscCoeff.CL[iMarker] / (ViscCoeff.CD[iMarker]+EPS);
          ViscCoeff.CFx[iMarker]         = ForceViscous[0];
          ViscCoeff.CFy[iMarker]         = ForceViscous[1];
          ViscCoeff.CMz[iMarker]         = MomentViscous[2];
          ViscCoeff.CoPx[iMarker]        = MomentZ_Force[1];
          ViscCoeff.CoPy[iMarker]        = -MomentZ_Force[0];
          ViscCoeff.CT[iMarker]          = -ViscCoeff.CFx[iMarker];
          ViscCoeff.CQ[iMarker]          = -ViscCoeff.CMz[iMarker];
          ViscCoeff.CMerit[iMarker]      = ViscCoeff.CT[iMarker] / (ViscCoeff.CQ[iMarker]+EPS);
          MaxHF_Visc[iMarker]            = pow(MaxHF_Visc[iMarker], 1.0/MaxNorm);
        }
        if (nDim == 3) {
          ViscCoeff.CD[iMarker]          =  ForceViscous[0]*cos(Alpha)*cos(Beta) + ForceViscous[1]*sin(Beta) + ForceViscous[2]*sin(Alpha)*cos(Beta);
          ViscCoeff.CL[iMarker]          = -ForceViscous[0]*sin(Alpha) + ForceViscous[2]*cos(Alpha);
          ViscCoeff.CSF[iMarker]         = -ForceViscous[0]*sin(Beta)*cos(Alpha) + ForceViscous[1]*cos(Beta) - ForceViscous[2]*sin(Beta)*sin(Alpha);
          ViscCoeff.CEff[iMarker]        = ViscCoeff.CL[iMarker]/(ViscCoeff.CD[iMarker] + EPS);
          ViscCoeff.CFx[iMarker]         = ForceViscous[0];
          ViscCoeff.CFy[iMarker]         = ForceViscous[1];
          ViscCoeff.CFz[iMarker]         = ForceViscous[2];
          ViscCoeff.CMx[iMarker]         = MomentViscous[0];
          ViscCoeff.CMy[iMarker]         = MomentViscous[1];
          ViscCoeff.CMz[iMarker]         = MomentViscous[2];
          ViscCoeff.CoPx[iMarker]        = -MomentY_Force[0];
          ViscCoeff.CoPz[iMarker]        = MomentY_Force[2];
          ViscCoeff.CT[iMarker]          = -ViscCoeff.CFz[iMarker];
          ViscCoeff.CQ[iMarker]          = -ViscCoeff.CMz[iMarker];
          ViscCoeff.CMerit[iMarker]      = ViscCoeff.CT[iMarker] / (ViscCoeff.CQ[iMarker] + EPS);
          MaxHF_Visc[iMarker]            = pow(MaxHF_Visc[iMarker], 1.0/MaxNorm);
        }

        AllBoundViscCoeff.CD          += ViscCoeff.CD[iMarker];
        AllBoundViscCoeff.CL          += ViscCoeff.CL[iMarker];
        AllBoundViscCoeff.CSF         += ViscCoeff.CSF[iMarker];
        AllBoundViscCoeff.CFx         += ViscCoeff.CFx[iMarker];
        AllBoundViscCoeff.CFy         += ViscCoeff.CFy[iMarker];
        AllBoundViscCoeff.CFz         += ViscCoeff.CFz[iMarker];
        AllBoundViscCoeff.CMx         += ViscCoeff.CMx[iMarker];
        AllBoundViscCoeff.CMy         += ViscCoeff.CMy[iMarker];
        AllBoundViscCoeff.CMz         += ViscCoeff.CMz[iMarker];
        AllBoundViscCoeff.CoPx        += ViscCoeff.CoPx[iMarker];
        AllBoundViscCoeff.CoPy        += ViscCoeff.CoPy[iMarker];
        AllBoundViscCoeff.CoPz        += ViscCoeff.CoPz[iMarker];
        AllBoundViscCoeff.CT          += ViscCoeff.CT[iMarker];
        AllBoundViscCoeff.CQ          += ViscCoeff.CQ[iMarker];
        AllBound_HF_Visc              += HF_Visc[iMarker];
        AllBound_MaxHF_Visc           += pow(MaxHF_Visc[iMarker], MaxNorm);

        /*--- Compute the coefficients per surface ---*/

        for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
          Monitoring_Tag = config->GetMarker_Monitoring_TagBound(iMarker_Monitoring);
          Marker_Tag = config->GetMarker_All_TagBound(iMarker);
          if (Marker_Tag == Monitoring_Tag) {
            SurfaceViscCoeff.CL[iMarker_Monitoring]      += ViscCoeff.CL[iMarker];
            SurfaceViscCoeff.CD[iMarker_Monitoring]      += ViscCoeff.CD[iMarker];
            SurfaceViscCoeff.CSF[iMarker_Monitoring]     += ViscCoeff.CSF[iMarker];
            SurfaceViscCoeff.CEff[iMarker_Monitoring]    += ViscCoeff.CEff[iMarker];
            SurfaceViscCoeff.CFx[iMarker_Monitoring]     += ViscCoeff.CFx[iMarker];
            SurfaceViscCoeff.CFy[iMarker_Monitoring]     += ViscCoeff.CFy[iMarker];
            SurfaceViscCoeff.CFz[iMarker_Monitoring]     += ViscCoeff.CFz[iMarker];
            SurfaceViscCoeff.CMx[iMarker_Monitoring]     += ViscCoeff.CMx[iMarker];
            SurfaceViscCoeff.CMy[iMarker_Monitoring]     += ViscCoeff.CMy[iMarker];
            SurfaceViscCoeff.CMz[iMarker_Monitoring]     += ViscCoeff.CMz[iMarker];
            Surface_HF_Visc[iMarker_Monitoring]          += HF_Visc[iMarker];
            Surface_MaxHF_Visc[iMarker_Monitoring]       += pow(MaxHF_Visc[iMarker],MaxNorm);
          }
        }
      }
    }
  }

  /*--- Update some global coeffients ---*/

  AllBoundViscCoeff.CEff = AllBoundViscCoeff.CL / (AllBoundViscCoeff.CD + EPS);
  AllBoundViscCoeff.CMerit = AllBoundViscCoeff.CT / (AllBoundViscCoeff.CQ + EPS);
  AllBound_MaxHF_Visc = pow(AllBound_MaxHF_Visc, 1.0/MaxNorm);


#ifdef HAVE_MPI

  /*--- Add AllBound information using all the nodes ---*/

  if (config->GetComm_Level() == COMM_FULL) {

    auto Allreduce = [](su2double x) {
      su2double tmp = x; x = 0.0;
      SU2_MPI::Allreduce(&tmp, &x, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      return x;
    };
    AllBoundViscCoeff.CD = Allreduce(AllBoundViscCoeff.CD);
    AllBoundViscCoeff.CL = Allreduce(AllBoundViscCoeff.CL);
    AllBoundViscCoeff.CSF = Allreduce(AllBoundViscCoeff.CSF);
    AllBoundViscCoeff.CEff = AllBoundViscCoeff.CL / (AllBoundViscCoeff.CD + EPS);

    AllBoundViscCoeff.CMx = Allreduce(AllBoundViscCoeff.CMx);
    AllBoundViscCoeff.CMy = Allreduce(AllBoundViscCoeff.CMy);
    AllBoundViscCoeff.CMz = Allreduce(AllBoundViscCoeff.CMz);

    AllBoundViscCoeff.CFx = Allreduce(AllBoundViscCoeff.CFx);
    AllBoundViscCoeff.CFy = Allreduce(AllBoundViscCoeff.CFy);
    AllBoundViscCoeff.CFz = Allreduce(AllBoundViscCoeff.CFz);

    AllBoundViscCoeff.CoPx = Allreduce(AllBoundViscCoeff.CoPx);
    AllBoundViscCoeff.CoPy = Allreduce(AllBoundViscCoeff.CoPy);
    AllBoundViscCoeff.CoPz = Allreduce(AllBoundViscCoeff.CoPz);

    AllBoundViscCoeff.CT = Allreduce(AllBoundViscCoeff.CT);
    AllBoundViscCoeff.CQ = Allreduce(AllBoundViscCoeff.CQ);
    AllBoundViscCoeff.CMerit = AllBoundViscCoeff.CT / (AllBoundViscCoeff.CQ + EPS);

    AllBound_HF_Visc = Allreduce(AllBound_HF_Visc);
    AllBound_MaxHF_Visc = pow(Allreduce(pow(AllBound_MaxHF_Visc, MaxNorm)), 1.0/MaxNorm);

  }

  /*--- Add the forces on the surfaces using all the nodes ---*/

  if (config->GetComm_Level() == COMM_FULL) {

    int nMarkerMon = config->GetnMarker_Monitoring();

    /*--- Use the same buffer for all reductions. We could avoid the copy back into
     *    the original variable by swaping pointers, but it is safer this way... ---*/

    su2double* buffer = new su2double [nMarkerMon];

    auto Allreduce_inplace = [buffer](int size, su2double* x) {
      SU2_MPI::Allreduce(x, buffer, size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      for(int i=0; i<size; ++i) x[i] = buffer[i];
    };

    Allreduce_inplace(nMarkerMon, SurfaceViscCoeff.CL);
    Allreduce_inplace(nMarkerMon, SurfaceViscCoeff.CD);
    Allreduce_inplace(nMarkerMon, SurfaceViscCoeff.CSF);

    for (iMarker_Monitoring = 0; iMarker_Monitoring < nMarkerMon; iMarker_Monitoring++)
      SurfaceViscCoeff.CEff[iMarker_Monitoring] = SurfaceViscCoeff.CL[iMarker_Monitoring] /
                                                 (SurfaceViscCoeff.CD[iMarker_Monitoring] + EPS);

    Allreduce_inplace(nMarkerMon, SurfaceViscCoeff.CFx);
    Allreduce_inplace(nMarkerMon, SurfaceViscCoeff.CFy);
    Allreduce_inplace(nMarkerMon, SurfaceViscCoeff.CFz);

    Allreduce_inplace(nMarkerMon, SurfaceViscCoeff.CMx);
    Allreduce_inplace(nMarkerMon, SurfaceViscCoeff.CMy);
    Allreduce_inplace(nMarkerMon, SurfaceViscCoeff.CMz);

    Allreduce_inplace(nMarkerMon, Surface_HF_Visc);
    Allreduce_inplace(nMarkerMon, Surface_MaxHF_Visc);

    delete [] buffer;

  }

#endif

  /*--- Update the total coefficients (note that all the nodes have the same value)---*/

  TotalCoeff.CD          += AllBoundViscCoeff.CD;
  TotalCoeff.CL          += AllBoundViscCoeff.CL;
  TotalCoeff.CSF         += AllBoundViscCoeff.CSF;
  TotalCoeff.CEff         = TotalCoeff.CL / (TotalCoeff.CD + EPS);
  TotalCoeff.CFx         += AllBoundViscCoeff.CFx;
  TotalCoeff.CFy         += AllBoundViscCoeff.CFy;
  TotalCoeff.CFz         += AllBoundViscCoeff.CFz;
  TotalCoeff.CMx         += AllBoundViscCoeff.CMx;
  TotalCoeff.CMy         += AllBoundViscCoeff.CMy;
  TotalCoeff.CMz         += AllBoundViscCoeff.CMz;
  TotalCoeff.CoPx        += AllBoundViscCoeff.CoPx;
  TotalCoeff.CoPy        += AllBoundViscCoeff.CoPy;
  TotalCoeff.CoPz        += AllBoundViscCoeff.CoPz;
  TotalCoeff.CT          += AllBoundViscCoeff.CT;
  TotalCoeff.CQ          += AllBoundViscCoeff.CQ;
  TotalCoeff.CMerit       = AllBoundViscCoeff.CT / (AllBoundViscCoeff.CQ + EPS);
  Total_Heat         = AllBound_HF_Visc;
  Total_MaxHeat      = AllBound_MaxHF_Visc;

  /*--- Update the total coefficients per surface (note that all the nodes have the same value)---*/

  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    SurfaceCoeff.CL[iMarker_Monitoring]         += SurfaceViscCoeff.CL[iMarker_Monitoring];
    SurfaceCoeff.CD[iMarker_Monitoring]         += SurfaceViscCoeff.CD[iMarker_Monitoring];
    SurfaceCoeff.CSF[iMarker_Monitoring]        += SurfaceViscCoeff.CSF[iMarker_Monitoring];
    SurfaceCoeff.CEff[iMarker_Monitoring]        = SurfaceViscCoeff.CL[iMarker_Monitoring] / (SurfaceCoeff.CD[iMarker_Monitoring] + EPS);
    SurfaceCoeff.CFx[iMarker_Monitoring]        += SurfaceViscCoeff.CFx[iMarker_Monitoring];
    SurfaceCoeff.CFy[iMarker_Monitoring]        += SurfaceViscCoeff.CFy[iMarker_Monitoring];
    SurfaceCoeff.CFz[iMarker_Monitoring]        += SurfaceViscCoeff.CFz[iMarker_Monitoring];
    SurfaceCoeff.CMx[iMarker_Monitoring]        += SurfaceViscCoeff.CMx[iMarker_Monitoring];
    SurfaceCoeff.CMy[iMarker_Monitoring]        += SurfaceViscCoeff.CMy[iMarker_Monitoring];
    SurfaceCoeff.CMz[iMarker_Monitoring]        += SurfaceViscCoeff.CMz[iMarker_Monitoring];
  }

}
