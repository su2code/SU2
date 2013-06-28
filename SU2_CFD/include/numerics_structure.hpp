/*!
 * \file numerics_structure.hpp
 * \brief Headers of the main subroutines for the dumerical definition of the problem.
 *        The subroutines and functions are in the <i>numerics_structure.cpp</i>, 
 *        <i>numerics_convective.cpp</i>, <i>numerics_viscous.cpp</i>, and 
 *        <i>numerics_source.cpp</i> files.
 * \author Current Development: Stanford University.
 *         Original Structure: CADES 1.0 (2009).
 * \version 1.0.
 *
 * Stanford University Unstructured (SU2) Code
 * Copyright (C) 2012 Aerospace Design Laboratory
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <cmath>
#include <iostream>

#include "../../Common/include/config_structure.hpp"

using namespace std;

/*! 
 * \class CNumerics
 * \brief Class for defining the numerical methods.
 * \author F. Palacios.
 * \version 1.0.
 */
class CNumerics {
protected:
	unsigned short nDim,	/*!< \brief Number of dimensions. */
	nVar;					/*!< \brief Number of variables. */
	unsigned short nSpecies, 	/*!< \brief No of species present in plasma */
	nFluids;					/*!< \brief No of fluids modeled in plasma */
	double Gamma;				/*!< \brief Fluid's Gamma constant (ratio of specific heats). */
	double GammaMonatomic;	/*!< \brief Fluid's Gamma constant (ratio of specific heats) for monatomic species. */
	double GammaDiatomic;		/*!< \brief Fluid's Gamma constant (ratio of specific heats) for diatomic species. */
	double Gamma_Minus_One;		/*!< \brief Fluids's Gamma - 1.0  . */
	double Gas_Constant;				/*!< \brief Gas constant. */
	unsigned short nDiatomics, nMonatomics;

public:
	double **Flux_Tensor,	/*!< \brief Flux tensor (used for viscous and inviscid purposes. */
	*Proj_Flux_Tensor;		/*!< \brief Flux tensor projected in a direction. */
	double **tau,		/*!< \brief Viscous stress tensor. */
	**delta;			/*!< \brief Identity matrix. */
	double Laminar_Viscosity_i,	/*!< \brief Laminar viscosity at point i. */
	Laminar_Viscosity_j;		/*!< \brief Laminar viscosity at point j. */
	double Eddy_Viscosity_i,	/*!< \brief Eddy viscosity at point i. */
	Eddy_Viscosity_j;			/*!< \brief Eddy viscosity at point j. */
	double Pressure_i,	/*!< \brief Pressure at point i. */
	Pressure_j;			/*!< \brief Pressure at point j. */
	double GravityForce_i,	/*!< \brief Gravity force at point i. */
	GravityForce_j;			/*!< \brief Gravity force at point j. */
	double Density_i,	/*!< \brief Pressure at point i. */
	Density_j;			/*!< \brief Pressure at point j. */
	double DensityInc_i,	/*!< \brief Incompressible density at point i. */
	DensityInc_j;			/*!< \brief Incompressible density at point j. */
	double BetaInc2_i,	/*!< \brief Beta incompressible at point i. */
	BetaInc2_j;			/*!< \brief Beta incompressible at point j. */
	double Lambda_i,	/*!< \brief Spectral radius at point i. */
	Lambda_j;			/*!< \brief Spectral radius at point j. */
	double LambdaComb_i,	/*!< \brief Spectral radius at point i. */
	LambdaComb_j;			/*!< \brief Spectral radius at point j. */
	double SoundSpeed_i,	/*!< \brief Sound speed at point i. */
	SoundSpeed_j;			/*!< \brief Sound speed at point j. */
	double Enthalpy_i,	/*!< \brief Enthalpy at point i. */
	Enthalpy_j;			/*!< \brief Enthalpy at point j. */
	double dist_i,	/*!< \brief Distance of point i to the nearest wall. */
	dist_j;			/*!< \brief Distance of point j to the nearest wall. */
	double Temp_i,	/*!< \brief Temperature at point i. */
	Temp_j;			/*!< \brief Temperature at point j. */
	double F1;		/*!< \brief Menter's first blending function. */
	double *Und_Lapl_i, /*!< \brief Undivided laplacians at point i. */
	*Und_Lapl_j;		/*!< \brief Undivided laplacians at point j. */
	double Sensor_i,	/*!< \brief Pressure sensor at point i. */
	Sensor_j;			/*!< \brief Pressure sensor at point j. */
	double *GridVel_i,	/*!< \brief Grid velocity at point i. */
	*GridVel_j;			/*!< \brief Grid velocity at point j. */
	double *RotVel_i,	/*!< \brief Rotational velocity at point i. */
	*RotVel_j;			/*!< \brief Rotational velocity at point j. */
	double *U_i,		/*!< \brief Vector of conservative variables at point i. */
	*U_j,				/*!< \brief Vector of conservative variables at point j. */
	*U_0,				/*!< \brief Vector of conservative variables at node 0. */
	*U_1,				/*!< \brief Vector of conservative variables at node 1. */
	*U_2;				/*!< \brief Vector of conservative variables at node 2. */
	double *Psi_i,		/*!< \brief Vector of adjoint variables at point i. */
	*Psi_j;				/*!< \brief Vector of adjoint variables at point j. */
	double *DeltaU_i,	/*!< \brief Vector of linearized variables at point i. */
	*DeltaU_j;			/*!< \brief Vector of linearized variables at point j. */
	double *TurbVar_i,	/*!< \brief Vector of turbulent variables at point i. */
	*TurbVar_j;			/*!< \brief Vector of turbulent variables at point j. */
	double *LevelSetVar_i,	/*!< \brief Vector of turbulent variables at point i. */
	*LevelSetVar_j;			/*!< \brief Vector of turbulent variables at point j. */
	double *TurbPsi_i,	/*!< \brief Vector of adjoint turbulent variables at point i. */
	*TurbPsi_j;			/*!< \brief Vector of adjoint turbulent variables at point j. */
	double **ConsVar_Grad_i,	/*!< \brief Gradient of conservative variables at point i. */
	**ConsVar_Grad_j,			/*!< \brief Gradient of conservative variables at point j. */
	**ConsVar_Grad_0,			/*!< \brief Gradient of conservative variables at point 0. */
	**ConsVar_Grad_1,			/*!< \brief Gradient of conservative variables at point 1. */
	**ConsVar_Grad_2,			/*!< \brief Gradient of conservative variables at point 2. */
	**ConsVar_Grad;				/*!< \brief Gradient of conservative variables which is a scalar. */
	double **PrimVar_Grad_i,	/*!< \brief Gradient of primitive variables at point i. */
	**PrimVar_Grad_j;			/*!< \brief Gradient of primitive variables at point j. */
	double **PsiVar_Grad_i,		/*!< \brief Gradient of adjoint variables at point i. */
	**PsiVar_Grad_j;			/*!< \brief Gradient of adjoint variables at point j. */
	double **TurbVar_Grad_i,	/*!< \brief Gradient of turbulent variables at point i. */
	**TurbVar_Grad_j;			/*!< \brief Gradient of turbulent variables at point j. */
	double **LevelSetVar_Grad_i,	/*!< \brief Gradient of turbulent variables at point i. */
	**LevelSetVar_Grad_j;			/*!< \brief Gradient of turbulent variables at point j. */
	double **TurbPsi_Grad_i,	/*!< \brief Gradient of adjoint turbulent variables at point i. */
	**TurbPsi_Grad_j;			/*!< \brief Gradient of adjoint turbulent variables at point j. */
	double *AuxVar_Grad_i,		/*!< \brief Gradient of an auxiliary variable at point i. */
	*AuxVar_Grad_j;				/*!< \brief Gradient of an auxiliary variable at point i. */
	double *Coord_i,	/*!< \brief Cartesians coordinates of point i. */
	*Coord_j,			/*!< \brief Cartesians coordinates of point j. */
	*Coord_0,			/*!< \brief Cartesians coordinates of point 0 (Galerkin method, triangle). */
	*Coord_1,			/*!< \brief Cartesians coordinates of point 1 (Galerkin method, tetrahedra). */
	*Coord_2,			/*!< \brief Cartesians coordinates of point 2 (Galerkin method, triangle). */
	*Coord_3;			/*!< \brief Cartesians coordinates of point 3 (Galerkin method, tetrahedra). */
	unsigned short Neighbor_i,	/*!< \brief Number of neighbors of the point i. */
	Neighbor_j;					/*!< \brief Number of neighbors of the point j. */
	double *Normal,	/*!< \brief Normal vector, it norm is the area of the face. */
	*UnitaryNormal;		/*!< \brief Unitary normal vector. */
	double TimeStep,		/*!< \brief Time step useful in dual time method. */
	Area,				/*!< \brief Area of the face i-j. */
	Volume;				/*!< \brief Volume of the control volume around point i. */
	double Volume_n,	/*!< \brief Volume of the control volume at time n. */
	Volume_nM1,		/*!< \brief Volume of the control volume at time n-1. */
	Volume_nP1;		/*!< \brief Volume of the control volume at time n+1. */
	double *U_n,	/*!< \brief Vector of conservative variables at time n. */
	*U_nM1,		/*!< \brief Vector of conservative variables at time n-1. */
	*U_nP1;		/*!< \brief Vector of conservative variables at time n+1. */

	/*! 
	 * \brief Constructor of the class.
	 */
	CNumerics(void);

	/*! 
	 * \overload
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CNumerics(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

	/*! 
	 * \overload
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] val_nSpecies - Number of species of the problem.
	 * \param[in] val_nFluids - Number of fluids of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CNumerics(unsigned short val_nDim, unsigned short val_nVar, unsigned short val_nSpecies, unsigned short val_nFluids, CConfig *config);

	
	/*! 
	 * \overload
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] val_nSpecies - Number of species of the problem.
	 * \param[in] val_nDiatomics - Number of diatomic species of the problem.
	 * \param[in] val_nMonatomics - Number of monatomic species of the problem.
     * \param[in] config - Definition of the particular problem.	 */
	CNumerics(unsigned short val_nDim, unsigned short val_nVar, unsigned short val_nSpecies, unsigned short val_nDiatomics, unsigned short val_nMonatomics, CConfig *config);
	
	/*! 
	 * \brief Destructor of the class. 
	 */
	virtual ~CNumerics(void);

	/*! 
	 * \brief Set the solution at different times.
	 * \param[in] val_u_nM1 Conservative solution at time n-1.
	 * \param[in] val_u_n Conservative solution at time n.
	 * \param[in] val_u_nP1 Conservative solution at time n+1.
	 */
	void SetPastSol(double *val_u_nM1, double *val_u_n, double *val_u_nP1);

	/*! 
	 * \brief Set the control volume at different times.
	 * \param[in] val_volume_nM1 - Control volume at time n-1.
	 * \param[in] val_volume_n - Control volume at time n.
	 * \param[in] val_volume_nP1 - Control volume at time n+1.
	 */
	void SetPastVolume(double val_volume_nM1, double val_volume_n, double val_volume_nP1);

	/*! 
	 * \brief Set the time step.
	 * \param[in] val_timestep - Value of the time step.
	 */
	void SetTimeStep(double val_timestep);

	/*! 
	 * \brief Set the value of the conservative variables.
	 * \param[in] val_u_i - Value of the conservative variable at point i.
	 * \param[in] val_u_j - Value of the conservative variable at point j.
	 */
	void SetConservative(double *val_u_i, double *val_u_j);

	/*! 
	 * \brief Set the value of the conservative variables.
	 * \param[in] val_u_0 - Value of the conservative variable at point 0.
	 * \param[in] val_u_1 - Value of the conservative variable at point 1.
	 * \param[in] val_u_2 - Value of the conservative variable at point 2.
	 */
	void SetConservative(double *val_u_0, double *val_u_1, double *val_u_2);

	/*! 
	 * \brief Set the gradient of the conservative variables.
	 * \param[in] val_consvar_grad_i - Gradient of the conservative variable at point i.
	 * \param[in] val_consvar_grad_j - Gradient of the conservative variable at point j.
	 */
	void SetConsVarGradient(double **val_consvar_grad_i, double **val_consvar_grad_j);

	/*! 
	 * \brief Set the gradient of the conservative variables.
	 * \param[in] val_consvar_grad_0 - Gradient of the conservative variable at point 0.
	 * \param[in] val_consvar_grad_1 - Gradient of the conservative variable at point 1.
	 * \param[in] val_consvar_grad_2 - Gradient of the conservative variable at point 2.
	 */
	void SetConsVarGradient(double **val_consvar_grad_0, double **val_consvar_grad_1, double **val_consvar_grad_2);

	/*! 
	 * \brief Set the gradient of the conservative variables.
	 * \param[in] val_consvar_grad - Gradient of the conservative variable which is a scalar.
	 */
	void SetConsVarGradient(double **val_consvar_grad);

	/*! 
	 * \brief Set the gradient of the primitive variables.
	 * \param[in] val_primvar_grad_i - Gradient of the primitive variable at point i.
	 * \param[in] val_primvar_grad_j - Gradient of the primitive variable at point j.
	 */
	void SetPrimVarGradient(double **val_primvar_grad_i, double **val_primvar_grad_j);

	/*! 
	 * \brief Set the value of the adjoint variable.
	 * \param[in] val_psi_i - Value of the adjoint variable at point i.
	 * \param[in] val_psi_j - Value of the adjoint variable at point j.
	 */
	void SetAdjointVar(double *val_psi_i, double *val_psi_j);

	/*! 
	 * \brief Set the value of the linearized conservative variables.
	 * \param[in] val_deltau_i - Value of the linearized conservative variable at point i.
	 * \param[in] val_deltau_j - Value of the linearized conservative variable at point j.
	 */
	void SetLinearizedVar(double *val_deltau_i, double *val_deltau_j);

	/*! 
	 * \brief Set the gradient of the adjoint variables.
	 * \param[in] val_psivar_grad_i - Gradient of the adjoint variable at point i.
	 * \param[in] val_psivar_grad_j - Gradient of the adjoint variable at point j.
	 */
	void SetAdjointVarGradient(double **val_psivar_grad_i, double **val_psivar_grad_j);

	/*! 
	 * \brief Set the value of the turbulent variable.
	 * \param[in] val_turbvar_i - Value of the turbulent variable at point i.
	 * \param[in] val_turbvar_j - Value of the turbulent variable at point j.
	 */
	void SetTurbVar(double *val_turbvar_i, double *val_turbvar_j);

	/*! 
	 * \brief Set the gradient of the turbulent variables.
	 * \param[in] val_turbvar_grad_i - Gradient of the turbulent variable at point i.
	 * \param[in] val_turbvar_grad_j - Gradient of the turbulent variable at point j.
	 */
	void SetTurbVarGradient(double **val_turbvar_grad_i, double **val_turbvar_grad_j);

	/*! 
	 * \brief Set the value of the level set variable.
	 * \param[in] val_levelsetvar_i - Value of the level set variable at point i.
	 * \param[in] val_levelsetvar_j - Value of the level set variable at point j.
	 */
	void SetLevelSetVar(double *val_levelsetvar_i, double *val_levelsetvar_j);

	/*! 
	 * \brief Set the gradient of the level set variables.
	 * \param[in] val_levelsetvar_grad_i - Gradient of the level set variable at point i.
	 * \param[in] val_levelsetvar_grad_j - Gradient of the level set variable at point j.
	 */
	void SetLevelSetVarGradient(double **val_levelsetvar_grad_i, double **val_levelsetvar_grad_j);

	/*! 
	 * \brief Set the value of the adjoint turbulent variable.
	 * \param[in] val_turbpsivar_i - Value of the adjoint turbulent variable at point i.
	 * \param[in] val_turbpsivar_j - Value of the adjoint turbulent variable at point j.
	 */
	void SetTurbAdjointVar(double *val_turbpsivar_i, double *val_turbpsivar_j);

	/*! 
	 * \brief Set the gradient of the adjoint turbulent variables.
	 * \param[in] val_turbpsivar_grad_i - Gradient of the adjoint turbulent variable at point i.
	 * \param[in] val_turbpsivar_grad_j - Gradient of the adjoint turbulent variable at point j.
	 */
	void SetTurbAdjointGradient (double **val_turbpsivar_grad_i, double **val_turbpsivar_grad_j);

	/*! 
	 * \brief Set the value of the turbulent variable.
	 * \param[in] val_F1 - Value of the first Menter blending function.
	 */
	void SetF1blending(double val_F1);

	/*!
	 * \brief Set the gradient of the auxiliary variables.
	 * \param[in] val_auxvargrad_i - Gradient of the auxiliary variable at point i.
	 * \param[in] val_auxvargrad_j - Gradient of the auxiliary variable at point j.
	 */
	void SetAuxVarGrad(double *val_auxvargrad_i, double *val_auxvargrad_j);

	/*! 
	 * \brief Compute the primitive variables at point i -> [Temperature vel_x vel_y vel_z Pressure].
	 * \param[in] val_consvar - Conservative variables.
	 * \param[in] val_primvar - Primitive variables.
	 */
	void ConsVar2PrimVar(double *val_consvar, double *val_primvar);

	/*! 
	 * \brief Compute the gradient of the primitive variables at point i -> [Density vel_x vel_y vel_z Pressure].
	 * \param[in] val_flowsol - Conservative variables.
	 * \param[in] val_consvar_grad - Gradient of conservative variables.
	 * \param[out] val_primvar_grad - Gradient of primitive variables.
	 */
	void ConsGrad2PrimGrad(double *val_flowsol, double **val_consvar_grad, double **val_primvar_grad);

	/*! 
	 * \brief Set the laminar viscosity.
	 * \param[in] val_laminar_viscosity_i - Value of the laminar viscosity at point i.
	 * \param[in] val_laminar_viscosity_j - Value of the laminar viscosity at point j.
	 */
	void SetLaminarViscosity(double val_laminar_viscosity_i, double val_laminar_viscosity_j);

	/*! 
	 * \brief Set the eddy viscosity.
	 * \param[in] val_eddy_viscosity_i - Value of the eddy viscosity at point i.
	 * \param[in] val_eddy_viscosity_j - Value of the eddy viscosity at point j.
	 */
	void SetEddyViscosity(double val_eddy_viscosity_i, double val_eddy_viscosity_j);

	/*! 
	 * \brief Set the value of the distance from the nearest wall.
	 * \param[in] val_dist_i - Value of of the distance from point i to the nearest wall.
	 * \param[in] val_dist_j - Value of of the distance from point j to the nearest wall.
	 */
	void SetDistance(double val_dist_i, double val_dist_j);

	/*! 
	 * \brief Set coordinates of the points.
	 * \param[in] val_coord_i - Coordinates of the point i.
	 * \param[in] val_coord_j - Coordinates of the point j.
	 */
	void SetCoord(double *val_coord_i, double *val_coord_j);

	/*! 
	 * \overload
	 * \param[in] val_coord_0 - Coordinates of the point 0.
	 * \param[in] val_coord_1 - Coordinates of the point 1.
	 * \param[in] val_coord_2 - Coordinates of the point 2.
	 */
	void SetCoord(double *val_coord_0, double *val_coord_1, double *val_coord_2);

	/*! 
	 * \overload
	 * \param[in] val_coord_0 - Coordinates of the point 0.
	 * \param[in] val_coord_1 - Coordinates of the point 1.
	 * \param[in] val_coord_2 - Coordinates of the point 2.
	 * \param[in] val_coord_3 - Coordinates of the point 3.
	 */
	void SetCoord(double *val_coord_0, double *val_coord_1, double *val_coord_2, 
			double *val_coord_3);

	/*! 
	 * \brief Set the velocity of the computational grid.
	 * \param[in] val_gridvel_i - Grid velocity of the point i.
	 * \param[in] val_gridvel_j - Grid velocity of the point j.
	 */
	void SetGridVel(double *val_gridvel_i, double *val_gridvel_j);

	/*! 
	 * \brief Set the velocity of the rotational framework.
	 * \param[in] val_rotvel_i - Rotational frame velocity of the point i.
	 * \param[in] val_rotvel_j - Rotational frame velocity of the point j.
	 */
	void SetRotVel(double *val_rotvel_i, double *val_rotvel_j);

	/*! 
	 * \brief Set the value of the pressure.
	 * \param[in] val_pressure_i - Value of the pressure at point i.
	 * \param[in] val_pressure_j - Value of the pressure at point j.
	 */
	void SetPressure(double val_pressure_i, double val_pressure_j);

	/*! 
	 * \brief Set the value of the density for the incompressible solver.
	 * \param[in] val_densityinc_i - Value of the pressure at point i.
	 * \param[in] val_densityinc_j - Value of the pressure at point j.
	 */
	void SetDensityInc(double val_densityinc_i, double val_densityinc_j);

	/*! 
	 * \brief Set the value of the beta for incompressible flows.
	 * \param[in] val_betainc2_i - Value of beta for incompressible flows at point i.
	 * \param[in] val_betainc2_j - Value of beta for incompressible flows at point j.
	 */
	void SetBetaInc2(double val_betainc2_i, double val_betainc2_j);

	/*! 
	 * \brief Set the value of the sound speed.
	 * \param[in] val_soundspeed_i - Value of the sound speed at point i.
	 * \param[in] val_soundspeed_j - Value of the sound speed at point j.
	 */
	void SetSoundSpeed(double val_soundspeed_i, double val_soundspeed_j);

	/*! 
	 * \brief Set the value of the temperature.
	 * \param[in] val_temp_i - Value of the temperature at point i.
	 * \param[in] val_temp_j - Value of the temperature at point j.
	 */
	void SetTemperature(double val_temp_i, double val_temp_j);

	/*! 
	 * \brief Set the value of the enthalpy.
	 * \param[in] val_enthalpy_i - Value of the enthalpy at point i.
	 * \param[in] val_enthalpy_j - Value of the enthalpy at point j.
	 */
	void SetEnthalpy(double val_enthalpy_i, double val_enthalpy_j);

	/*! 
	 * \brief Set the value of the spectral radius.
	 * \param[in] val_lambda_i - Value of the spectral radius at point i.
	 * \param[in] val_lambda_j - Value of the spectral radius at point j.
	 */
	void SetLambda(double val_lambda_i, double val_lambda_j);
	
	/*! 
	 * \brief Set the value of the lambda for combustion problems.
	 * \param[in] val_lambdacomb_i - Value of the lambda combustion at point i.
	 * \param[in] val_lambdacomb_j - Value of the lambda combustion at point j.
	 */
	void SetLambdaComb(double val_lambdacomb_i, double val_lambdacomb_j);

	/*! 
	 * \brief Set the value of undivided laplacian.
	 * \param[in] val_und_lapl_i Undivided laplacian at point i.
	 * \param[in] val_und_lapl_j Undivided laplacian at point j.
	 */
	void SetUndivided_Laplacian(double *val_und_lapl_i, double *val_und_lapl_j);

	/*! 
	 * \brief Set the value of the pressure sensor.
	 * \param[in] val_sensor_i Pressure sensor at point i.
	 * \param[in] val_sensor_j Pressure sensor at point j.
	 */
	void SetSensor(double val_sensor_i, double val_sensor_j);

	/*! 
	 * \brief Set the number of neighbor to a point.
	 * \param[in] val_neighbor_i - Number of neighbor to point i.
	 * \param[in] val_neighbor_j - Number of neighbor to point j.
	 */
	void SetNeighbor(unsigned short val_neighbor_i, unsigned short val_neighbor_j);

	/*! 
	 * \brief Set the value of the normal vector to the face between two points.
	 * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
	 */
	void SetNormal(double *val_normal);

	/*! 
	 * \brief Set the value of the volume of the control volume.
	 * \param[in] val_volume Volume of the control volume.
	 */
	void SetVolume(double val_volume);

	/*! 
	 * \brief Get the inviscid fluxes.
	 * \param[in] val_density - Value of the density.
	 * \param[in] val_velocity - Value of the velocity.
	 * \param[in] val_pressure - Value of the pressure.
	 * \param[in] val_enthalpy - Value of the enthalpy.
	 */
	void GetInviscidFlux(double val_density, double *val_velocity, double val_pressure, double val_enthalpy);

	/*! 
	 * \brief Get the viscous fluxes.
	 * \param[in] val_primvar - Value of the primitive variables.
	 * \param[in] val_gradprimvar - Gradient of the primitive variables.
	 * \param[in] val_laminar_viscosity - Value of the laminar viscosity.
	 * \param[in] val_eddy_viscosity - Value of the eddy viscosity.
	 * \param[in] val_mach_inf - Value of the Mach number at the infinity.
	 */
	void GetViscousFlux(double *val_primvar, double **val_gradprimvar, 
			double val_laminar_viscosity, double val_eddy_viscosity, double val_mach_inf);

	/*! 
	 * \brief Compute the projected inviscid flux vector.
	 * \param[in] val_density - Pointer to the density.
	 * \param[in] val_velocity - Pointer to the velocity.
	 * \param[in] val_pressure - Pointer to the pressure.
	 * \param[in] val_enthalpy - Pointer to the enthalpy.
	 * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
	 * \param[out] val_Proj_Flux - Pointer to the projected flux.
	 */
	void GetInviscidProjFlux(double *val_density, double *val_velocity, double *val_pressure, double *val_enthalpy, 
			double *val_normal, double *val_Proj_Flux);

	/*! 
	 * \brief Compute the projected inviscid flux vector.
	 * \param[in] val_density - Pointer to the density.
	 * \param[in] val_velocity - Pointer to the velocity.
	 * \param[in] val_pressure - Pointer to the pressure.
	 * \param[in] val_enthalpy - Pointer to the enthalpy.
	 * \param[in] val_energy_vib - Pointer to the vibrational energy.
	 * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
	 * \param[out] val_Proj_Flux - Pointer to the projected flux.
	 */
	void GetInviscidProjFlux_(double *val_density, double **val_velocity,double *val_pressure, double *val_enthalpy, 
														double *val_energy_vib, double *val_normal, double *val_Proj_Flux);
		

	/*!
	 * \brief Compute the projected inviscid flux vector for incompresible simulations
	 * \param[in] val_density - Pointer to the density.
	 * \param[in] val_velocity - Pointer to the velocity.
	 * \param[in] val_pressure - Pointer to the pressure.
	 * \param[in] val_gravityforce - Pointer to the gravity force.
	 * \param[in] val_betainc2 - Value of the artificial compresibility factor.
	 * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
	 * \param[out] val_Proj_Flux - Pointer to the projected flux.
	 */
	void GetInviscidArtCompProjFlux(double *val_density, double *val_velocity, double *val_pressure, double *val_gravityforce, double *val_betainc2, 
			double *val_normal, double *val_Proj_Flux);

	/*! 
	 * \overload
	 * \brief Overloaded function for multi-species formulation (compressible flow).
	 * \brief Compute the projected inviscid flux vector.
	 * \param[in] val_density - Pointer to the density.
	 * \param[in] val_velocity - Pointer to the velocity.
	 * \param[in] val_pressure - Pointer to the pressure.
	 * \param[in] val_enthalpy - Pointer to the enthalpy.
	 * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
	 * \param[out] val_Proj_Flux - Pointer to the projected flux.
	 */
	void GetInviscidProjFlux(double *val_density, double **val_velocity, double *val_pressure, double *val_enthalpy, 
			double *val_normal, double *val_Proj_Flux);

	/*! 
	 * \brief Compute the projection of the viscous fluxes into a direction.
	 * \param[in] val_primvar - Primitive variables.
	 * \param[in] val_gradprimvar - Gradient of the primitive variables.
	 * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
	 * \param[in] val_laminar_viscosity - Laminar viscosity.
	 * \param[in] val_eddy_viscosity - Eddy viscosity.
	 */

	void GetViscousProjFlux(double *val_primvar, double **val_gradprimvar, double *val_normal, double val_laminar_viscosity,
			double val_eddy_viscosity);

	/*! 
	 * \brief Compute the projection of the viscous fluxes into a direction (artificial compresibility method).
	 * \param[in] val_primvar - Primitive variables.
	 * \param[in] val_gradprimvar - Gradient of the primitive variables.
	 * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
	 * \param[in] val_laminar_viscosity - Laminar viscosity.
	 * \param[in] val_eddy_viscosity - Eddy viscosity.
	 */

	void GetViscousArtCompProjFlux(double *val_primvar, double **val_gradprimvar, double *val_normal, double val_laminar_viscosity,
			double val_eddy_viscosity);

	/*! 
	 * \brief Compute the projection of the inviscid Jacobian matrices.
	 * \param[in] val_velocity Pointer to the velocity.
	 * \param[in] val_energy Value of the energy.
	 * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
	 * \param[in] val_scale - Scale of the projection.
	 * \param[out] val_Proj_Jac_tensor - Pointer to the projected inviscid Jacobian.
	 */
	void GetInviscidProjJac(double *val_velocity, double val_energy, double *val_normal, 
			double val_scale, double **val_Proj_Jac_tensor);

	/*! 
	 * \brief Compute the projection of the inviscid Jacobian matrices (artificial compresibility).
	 * \param[in] val_density - Value of the density.
	 * \param[in] val_velocity - Pointer to the velocity.
	 * \param[in] val_betainc2 - Value of the artificial compresibility factor.
	 * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
	 * \param[in] val_scale - Scale of the projection.
	 * \param[out] val_Proj_Jac_tensor - Pointer to the projected inviscid Jacobian.
	 */
	void GetInviscidArtCompProjJac(double val_density, double *val_velocity, double val_betainc2, double *val_normal, 
			double val_scale, double **val_Proj_Jac_tensor);

	/*! 
	 * \overload
	 * \brief Compute the projection of the inviscid Jacobian matrices.
	 * \param[in] val_velocity Pointer to the velocity.
	 * \param[in] val_energy Value of the energy.
	 * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
	 * \param[in] val_scale - Scale of the projection.
	 * \param[out] val_Proj_Jac_tensor - Pointer to the projected inviscid Jacobian.
	 */
	void GetInviscidProjJac(double **val_velocity, double *val_energy, double *val_normal, 
			double val_scale, double **val_Proj_Jac_tensor);

	
	/*! 
	 * \overload
	 * \brief Compute the projection of the inviscid Jacobian matrices.
	 * \param[in] val_velocity Pointer to the velocity.
	 * \param[in] val_energy Value of the energy.
	 * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
	 * \param[in] val_scale - Scale of the projection.
	 * \param[out] val_Proj_Jac_Tensor - Pointer to the projected inviscid Jacobian.
	 */
	void GetInviscidProjJac_(double **val_velocity, double *val_energy, double *val_normal, 
													 double val_scale, double **val_Proj_Jac_Tensor);


	/*! 
	 * \brief Compute the projection of the viscous Jacobian matrices.
	 * \param[in] val_density - Value of the density.
	 * \param[in] val_pressure - Value of the pressure.
	 * \param[in] val_laminar_viscosity - Value of the laminar viscosity.
	 * \param[in] val_eddy_viscosity - Value of the eddy viscosity.
	 * \param[in] val_dist_ij - Distance between the points.
	 * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
	 * \param[in] val_dS - Area of the face between two nodes.
	 * \param[in] val_Mean_PrimVar - Mean value of the primitive variables.
	 * \param[in] val_Proj_Visc_Flux - Pointer to the projected viscous flux.
	 * \param[out] val_Proj_Jac_Tensor_i - Pointer to the projected viscous Jacobian at point i.
	 * \param[out] val_Proj_Jac_Tensor_j - Pointer to the projected viscous Jacobian at point j.
	 */
	void GetViscousProjJacs(double val_density, double val_pressure, double val_laminar_viscosity, 
			double val_eddy_viscosity, double val_dist_ij, double *val_normal, double val_dS,
			double *val_Mean_PrimVar, double *val_Proj_Visc_Flux, double **val_Proj_Jac_Tensor_i,
			double **val_Proj_Jac_Tensor_j);

	/*! 
	 * \brief Compute the projection of the viscous Jacobian matrices.
	 * \param[in] val_density - Value of the density.
	 * \param[in] val_pressure - Value of the pressure.
	 * \param[in] val_laminar_viscosity - Value of the laminar viscosity.
	 * \param[in] val_eddy_viscosity - Value of the eddy viscosity.
	 * \param[in] val_dist_ij - Distance between the points.
	 * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
	 * \param[in] val_dS - Area of the face between two nodes.
	 * \param[in] val_Mean_PrimVar - Mean value of the primitive variables.
	 * \param[in] val_Proj_Visc_Flux - Pointer to the projected viscous flux.
	 * \param[out] val_Proj_Jac_Tensor_i - Pointer to the projected viscous Jacobian at point i.
	 * \param[out] val_Proj_Jac_Tensor_j - Pointer to the projected viscous Jacobian at point j.
	 */
	void GetViscousArtCompProjJacs(double val_density, double val_pressure, double val_laminar_viscosity, 
			double val_eddy_viscosity, double val_dist_ij, double *val_normal, double val_dS,
			double *val_Mean_PrimVar, double *val_Proj_Visc_Flux, double **val_Proj_Jac_Tensor_i,
			double **val_Proj_Jac_Tensor_j);

	/*! 
	 * \brief Computation of the matrix P, this matrix diagonalize the conservative Jacobians in 
	 *        the form $P^{-1}(A.Normal)P=Lambda$.
	 * \param[in] val_density - Value of the density.
	 * \param[in] val_velocity - Value of the velocity.
	 * \param[in] val_soundspeed - Value of the sound speed.
	 * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
	 * \param[out] val_p_tensor - Pointer to the P matrix.
	 */
	void GetPMatrix(double *val_density, double *val_velocity, double *val_soundspeed, 
			double *val_normal, double **val_p_tensor);

	/*! 
	 * \overload
	 * \brief Computation of the matrix P, this matrix diagonalize the conservative Jacobians in 
	 *        the form $P^{-1}(A.Normal)P=Lambda$.
	 * \param[in] val_density - Value of the density.
	 * \param[in] val_velocity - Value of the velocity.
	 * \param[in] val_soundspeed - Value of the sound speed.
	 * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
	 * \param[out] val_p_tensor - Pointer to the P matrix.
	 */
	void GetPMatrix(double *val_density, double **val_velocity, double *val_soundspeed, 
			double *val_normal, double **val_p_tensor);
	
	/*! 
	 * \overload
	 * \brief Computation of the matrix P, this matrix diagonalize the conservative Jacobians in 
	 *        the form $P^{-1}(A.Normal)P=Lambda$.
	 * \param[in] val_density - Value of the density.
	 * \param[in] val_velocity - Value of the velocity.
	 * \param[in] val_soundspeed - Value of the sound speed.
	 * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
	 * \param[out] val_p_tensor - Pointer to the P matrix.
	 */
	void GetPMatrix_(double *val_density, double **val_velocity, double *val_soundspeed, 
									 double *val_normal, double **val_p_tensor);
	
	/*!
	 * \brief Computation of the matrix P (artificial compresibility), this matrix diagonalize the conservative Jacobians in 
	 *        the form $P^{-1}(A.Normal)P=Lambda$.
	 * \param[in] val_density - Value of the density.
	 * \param[in] val_velocity - Value of the velocity.
	 * \param[in] val_betainv2 - Value of the compresibility factor.
	 * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
	 * \param[out] val_p_tensor - Pointer to the P matrix.
	 */
	void GetPArtCompMatrix(double *val_density, double *val_velocity, double *val_betainv2, double *val_normal, double **val_p_tensor);

	/*! 
	 * \brief Computation of the matrix P^{-1}, this matrix diagonalize the conservative Jacobians 
	 *        in the form $P^{-1}(A.Normal)P=Lambda$.
	 * \param[in] val_density - Value of the density.
	 * \param[in] val_velocity - Value of the velocity.
	 * \param[in] val_soundspeed - Value of the sound speed.
	 * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
	 * \param[out] val_invp_tensor - Pointer to inverse of the P matrix.
	 */
	void GetPMatrix_inv(double *val_density, double *val_velocity, double *val_soundspeed, 
			double *val_normal, double **val_invp_tensor);

	/*! 
	 * \overload
	 * \brief Computation of the matrix P^{-1}, this matrix diagonalize the conservative Jacobians 
	 *        in the form $P^{-1}(A.Normal)P=Lambda$.
	 * \param[in] val_density - Value of the density.
	 * \param[in] val_velocity - Value of the velocity.
	 * \param[in] val_soundspeed - Value of the sound speed.
	 * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
	 * \param[out] val_invp_tensor - Pointer to inverse of the P matrix.
	 */
	void GetPMatrix_inv(double *val_density, double **val_velocity, double *val_soundspeed, 
			double *val_normal, double **val_invp_tensor);
	
	/*! 
	 * \overload
	 * \brief Computation of the matrix P^{-1}, this matrix diagonalize the conservative Jacobians 
	 *        in the form $P^{-1}(A.Normal)P=Lambda$.
	 * \param[in] val_density - Value of the density.
	 * \param[in] val_velocity - Value of the velocity.
	 * \param[in] val_soundspeed - Value of the sound speed.
	 * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
	 * \param[out] val_invp_tensor - Pointer to inverse of the P matrix.
	 */
	void GetPMatrix_inv_(double *val_density, double **val_velocity, double *val_soundspeed, 
											 double *val_normal, double **val_invp_tensor);

	/*!
	 * \brief Computation of the matrix P^{-1} (artificial compresibility), this matrix diagonalize the conservative Jacobians 
	 *        in the form $P^{-1}(A.Normal)P=Lambda$.
	 * \param[in] val_density - Value of the density.
	 * \param[in] val_velocity - Value of the velocity.
	 * \param[in] val_betainv2 - Value of the compresibility factor.
	 * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
	 * \param[out] val_invp_tensor - Pointer to inverse of the P matrix.
	 */
	void GetPArtCompMatrix_inv(double *val_density, double *val_velocity, double *val_betainv2, double *val_normal, double **val_invp_tensor);

	/*! 
	 * \brief Computation of the projected inviscid lambda (eingenvalues).
	 * \param[in] val_velocity - Value of the velocity.
	 * \param[in] val_soundspeed - Value of the sound speed.
	 * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
	 * \param[in] val_Lambda_Vector - Pointer to Lambda matrix.
	 */
	void GetJacInviscidLambda_fabs(double *val_velocity, double val_soundspeed, 
			double *val_normal, double *val_Lambda_Vector);

	/*! 
	 * \brief Compute the numerical residual.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetResidual(double *val_residual, CConfig *config);

	/*! 
	 * \overload
	 * \param[out] val_residual_i - Pointer to the total residual at point i.
	 * \param[out] val_residual_j - Pointer to the total residual at point j.
	 */
	virtual void SetResidual(double *val_residual_i, double *val_residual_j);

	/*! 
	 * \overload
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config);

	/*! 
	 * \overload
	 * \param[out] val_resconv - Pointer to the convective residual.
	 * \param[out] val_resvisc - Pointer to the artificial viscosity residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] art_diss - If <code>TRUE</code>, the artificial dissipation is computed; otherwise <code>FALSE</code>.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetResidual(double *val_resconv, double *val_resvisc, double **val_Jacobian_i, double **val_Jacobian_j, 
			bool art_diss, CConfig *config);

	/*! 
	 * \overload
	 * \param[out] val_residual_i - Pointer to the total residual at point i.
	 * \param[out] val_residual_j - Pointer to the total viscosity residual at point j.
	 * \param[out] val_Jacobian_ii - Jacobian of the numerical method at node i (implicit computation) from node i.
	 * \param[out] val_Jacobian_ij - Jacobian of the numerical method at node i (implicit computation) from node j.
	 * \param[out] val_Jacobian_ji - Jacobian of the numerical method at node j (implicit computation) from node i.
	 * \param[out] val_Jacobian_jj - Jacobian of the numerical method at node j (implicit computation) from node j.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetResidual(double *val_residual_i, double *val_residual_j, 
			double **val_Jacobian_ii, double **val_Jacobian_ij,
			double **val_Jacobian_ji, double **val_Jacobian_jj, CConfig *config);

	/*! 
	 * \overload
	 * \param[out] val_resconv_i - Pointer to the convective residual at point i.
	 * \param[out] val_resvisc_i - Pointer to the artificial viscosity residual at point i.
	 * \param[out] val_resconv_j - Pointer to the convective residual at point j.
	 * \param[out] val_resvisc_j - Pointer to the artificial viscosity residual at point j.
	 * \param[out] val_Jacobian_ii - Jacobian of the numerical method at node i (implicit computation) from node i.
	 * \param[out] val_Jacobian_ij - Jacobian of the numerical method at node i (implicit computation) from node j.
	 * \param[out] val_Jacobian_ji - Jacobian of the numerical method at node j (implicit computation) from node i.
	 * \param[out] val_Jacobian_jj - Jacobian of the numerical method at node j (implicit computation) from node j.
	 * \param[in] local_art_diss - If 3 then 3th order dissipation, if 1, then 1st order dissipation; if 0, no dissipation.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetResidual(double *val_resconv_i, double *val_resvisc_i, double *val_resconv_j, double *val_resvisc_j, 
			double **val_Jacobian_ii, double **val_Jacobian_ij, double **val_Jacobian_ji, double **val_Jacobian_jj,
			unsigned short local_art_diss, CConfig *config);

	/*! 
	 * \overload
	 * \param[out] val_stiffmatrix_elem - Stiffness matrix for Galerkin computation.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetResidual(double **val_stiffmatrix_elem, CConfig *config);

	/*! 
	 * \overload
	 * \param[in] config - Definition of the particular problem.
	 * \param[out] val_residual - residual of the source terms
	 * \param[out] val_Jacobian_i - Jacobian of the source terms
	 */
	virtual void SetResidual(double *val_residual, double **val_Jacobian_i, CConfig *config);

};

/*! 
 * \class CUpwRoe_Flow
 * \brief Class for solving an approximate Riemann solver of Roe for the flow equations.
 * \ingroup ConvDiscr
 * \author A. Bueno (UPM) & F. Palacios (Stanford University).
 * \version 1.0.
 */
class CUpwRoe_Flow : public CNumerics {
private:
	bool implicit;
	double *Diff_U;
	double *Velocity_i, *Velocity_j, *RoeVelocity;
	double *Proj_flux_tensor_i, *Proj_flux_tensor_j;
	double *delta_wave, *delta_vel;
	double *Lambda, *Epsilon;
	double **P_Tensor, **invP_Tensor;
	double sq_vel, Proj_ModJac_Tensor_ij, Density_i, Energy_i, SoundSpeed_i, Pressure_i, Enthalpy_i, 
	Density_j, Energy_j, SoundSpeed_j, Pressure_j, Enthalpy_j, R, RoeDensity, RoeEnthalpy, RoeSoundSpeed, 
	ProjVelocity, ProjVelocity_i, ProjVelocity_j, proj_delta_vel, delta_p, delta_rho;
	unsigned short iDim, iVar, jVar, kVar;

public:

	/*! 
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CUpwRoe_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

	/*! 
	 * \brief Destructor of the class. 
	 */
	~CUpwRoe_Flow(void);

	/*! 
	 * \brief Compute the Roe's flux between two nodes i and j.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config);
};

/*! 
 * \class CUpwRoeArtComp_Flow
 * \brief Class for solving an approximate Riemann solver of Roe for the incompressible flow equations.
 * \ingroup ConvDiscr
 * \author F. Palacios.
 * \version 1.0.
 */
class CUpwRoeArtComp_Flow : public CNumerics {
private:
	bool implicit;
	bool gravity;
	double Froude;
	double *Diff_U;
	double *Velocity_i, *Velocity_j, *MeanVelocity;
	double *Proj_flux_tensor_i, *Proj_flux_tensor_j;
	double *Lambda, *Epsilon;
	double **P_Tensor, **invP_Tensor;
	double sq_vel, Proj_ModJac_Tensor_ij, Density_i, Energy_i, SoundSpeed_i, Pressure_i, Enthalpy_i, 
	Density_j, Energy_j, SoundSpeed_j, Pressure_j, Enthalpy_j, R, MeanDensity, MeanEnthalpy, MeanSoundSpeed, MeanPressure, MeanBetaInc2,
	ProjVelocity, ProjVelocity_i, ProjVelocity_j, proj_delta_vel, delta_p, delta_rho, vn;
	unsigned short iDim, iVar, jVar, kVar;

public:

	/*! 
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CUpwRoeArtComp_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

	/*! 
	 * \brief Destructor of the class. 
	 */
	~CUpwRoeArtComp_Flow(void);

	/*! 
	 * \brief Compute the Roe's flux between two nodes i and j.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config);
};

/*! 
 * \class CUpwRoe_AdjFlow
 * \brief Class for solving an approximate Riemann solver of Roe 
 *        for the adjoint flow equations.
 * \ingroup ConvDiscr
 * \author F. Palacios.
 * \version 1.0.
 */
class CUpwRoe_AdjFlow : public CNumerics {
private:
	double *Residual_Roe;
	double area, Sx, Sy, Sz, rarea, nx, ny, nz, rho_l, u_l, v_l, w_l, h_l, rho_r, 
	u_r, v_r, w_r, h_r, psi1, psi2, psi3, psi4, psi5;
	double h, u, v, w, c, psi1_l, psi2_l, psi3_l, psi4_l, psi5_l,
	psi1_r, psi2_r, psi3_r, psi4_r, psi5_r, q_l, q_r, Q_l, Q_r, vn,
	rrho_l, weight, rweight1, cc;
	double l1psi, l2psi, absQ, absQp, absQm, q2, alpha, beta_u, beta_v, beta_w, Q, l1l2p, l1l2m, eta;
	unsigned short iDim, iVar;

public:

	/*! 
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CUpwRoe_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

	/*! 
	 * \brief Destructor of the class. 
	 */
	~CUpwRoe_AdjFlow(void);

	/*! 
	 * \brief Compute the adjoint Roe's flux between two nodes i and j.
	 * \param[out] val_residual_i - Pointer to the total residual at point i.
	 * \param[out] val_residual_j - Pointer to the total residual at point j.
	 */
	void SetResidual (double *val_residual_i, double *val_residual_j);
};

/*! 
 * \class CUpwRoe_Plasma
 * \brief Class for solving an approximate Riemann solver of Roe for the plasma equations.
 * \ingroup ConvDiscr
 * \author ADL Stanford
 * \version 1.0.
 */
class CUpwRoe_Plasma : public CNumerics {
private:
	bool implicit;
	double *Diff_U;
	double **Velocity_i, **Velocity_j, **RoeVelocity;
	double *Proj_flux_tensor_i, *Proj_flux_tensor_j;
	double *delta_wave, **delta_vel;
	double *Lambda, *Epsilon;
	double *Density_i, *Energy_i, *SoundSpeed_i, *Pressure_i, *Enthalpy_i, 
	*Density_j, *Energy_j, *SoundSpeed_j, *Pressure_j, *Enthalpy_j,  *RoeDensity, *RoeEnthalpy, *RoeSoundSpeed, 
	*ProjVelocity, *ProjVelocity_i, *ProjVelocity_j;
	double **P_Tensor, **invP_Tensor;
	double sq_vel, Proj_ModJac_Tensor_ij, R,delta_p, delta_rho,proj_delta_vel; 
	unsigned short iDim, iVar, jVar, kVar, iFluids;

public:

	/*! 
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] val_nSpecies - Number of species in the problem.
	 * \param[in] val_nFluids - Number of fluids in the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CUpwRoe_Plasma(unsigned short val_nDim, unsigned short val_nVar,unsigned short val_nSpecies, unsigned short val_nFluids, CConfig *config);

	/*! 
	 * \brief Destructor of the class. 
	 */
	~CUpwRoe_Plasma(void);

	/*! 
	 * \brief Compute the Roe's flux between two nodes i and j.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config);
};

/*! 
 * \class CUpwRoe_Plasma
 * \brief Class for solving an approximate Riemann solver of Roe for the plasma equations.
 * \ingroup ConvDiscr
 * \author ADL Stanford
 * \version 1.0.
 */
class CUpwRoe_PlasmaDiatomic : public CNumerics {
private:
	bool implicit;
	double *Diff_U;
	double **Velocity_i, **Velocity_j, **RoeVelocity;
	double *Proj_flux_tensor_i, *Proj_flux_tensor_j;
	double *delta_wave, **delta_vel;
	double *Lambda, *Epsilon;
	double *Density_i, *Energy_i, *Energy_vib_i, *SoundSpeed_i, *Pressure_i, *Enthalpy_i, 
	*Density_j, *Energy_j, *Energy_vib_j, *SoundSpeed_j, *Pressure_j, *Enthalpy_j,  *RoeDensity, *RoeEnthalpy, *RoeSoundSpeed, 
	*ProjVelocity, *ProjVelocity_i, *ProjVelocity_j;
	double **P_Tensor, **invP_Tensor;
	double sq_vel, Proj_ModJac_Tensor_ij, R,delta_p, delta_rho,proj_delta_vel; 
	unsigned short iDim, iVar, jVar, kVar, iFluids;
	
public:
	
	/*! 
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] val_nSpecies - Number of species in the problem.
	 * \param[in] val_nDiatomics - Number of diatomic species in the problem.
	 * \param[in] val_nMonatomics - Number of monatomic species in the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CUpwRoe_PlasmaDiatomic(unsigned short val_nDim, unsigned short val_nVar,unsigned short val_nSpecies, unsigned short val_nDiatomics, unsigned short val_nMonatomics, CConfig *config);
	
	/*! 
	 * \brief Destructor of the class. 
	 */
	~CUpwRoe_PlasmaDiatomic(void);
	
	/*! 
	 * \brief Compute the Roe's flux between two nodes i and j.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config);
};

/*! 
 * \class CUpwAUSM_Flow
 * \brief Class for solving an approximate Riemann AUSM.
 * \ingroup ConvDiscr
 * \author F. Palacios
 * \version 1.0.
 */
class CUpwAUSM_Flow : public CNumerics {
private:
	bool implicit;
	double *Diff_U;
	double *Velocity_i, *Velocity_j, *RoeVelocity;
	double *Proj_flux_tensor_i, *Proj_flux_tensor_j;
	double *delta_wave, *delta_vel;
	double *Lambda, *Epsilon;
	double **P_Tensor, **invP_Tensor;
	double sq_vel, Proj_ModJac_Tensor_ij, Density_i, Energy_i, SoundSpeed_i, Pressure_i, Enthalpy_i, 
	Density_j, Energy_j, SoundSpeed_j, Pressure_j, Enthalpy_j, R, RoeDensity, RoeEnthalpy, RoeSoundSpeed, 
	ProjVelocity, ProjVelocity_i, ProjVelocity_j, proj_delta_vel, delta_p, delta_rho;
	unsigned short iDim, iVar, jVar, kVar;
	
public:
	
	/*! 
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CUpwAUSM_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
	
	/*! 
	 * \brief Destructor of the class. 
	 */
	~CUpwAUSM_Flow(void);
	
	/*! 
	 * \brief Compute the Roe's flux between two nodes i and j.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config);
};

/*! 
 * \class CUpwHLLC_Flow
 * \brief Class for solving an approximate Riemann AUSM.
 * \ingroup ConvDiscr
 * \author F. Palacios, based on the Joe code implementation 
 * \version 1.0.
 */
class CUpwHLLC_Flow : public CNumerics {
private:
	bool implicit;
	double *Diff_U;
	double *Velocity_i, *Velocity_j, *RoeVelocity;
	double *Proj_flux_tensor_i, *Proj_flux_tensor_j;
	double *delta_wave, *delta_vel;
	double *Lambda, *Epsilon;
	double **P_Tensor, **invP_Tensor;
	double sq_vel, Proj_ModJac_Tensor_ij, Density_i, Energy_i, SoundSpeed_i, Pressure_i, Enthalpy_i, 
	Density_j, Energy_j, SoundSpeed_j, Pressure_j, Enthalpy_j, R, RoeDensity, RoeEnthalpy, RoeSoundSpeed, 
	ProjVelocity, ProjVelocity_i, ProjVelocity_j, proj_delta_vel, delta_p, delta_rho;
	unsigned short iDim, iVar, jVar, kVar;
	
public:
	
	/*! 
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CUpwHLLC_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
	
	/*! 
	 * \brief Destructor of the class. 
	 */
	~CUpwHLLC_Flow(void);
	
	/*! 
	 * \brief Compute the Roe's flux between two nodes i and j.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config);
};

/*! 
 * \class CUpwLin_TurbSA
 * \brief Class for performing a linear upwind solver for the Spalart-Allmaras turbulence model equations.
 * \ingroup ConvDiscr
 * \author A. Bueno.
 * \version 1.0.
 */
class CUpwLin_TurbSA : public CNumerics {
private:
	double *Velocity_i;
	double *Velocity_j;

public:

	/*! 
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CUpwLin_TurbSA(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

	/*! 
	 * \brief Destructor of the class. 
	 */
	~CUpwLin_TurbSA(void);

	/*! 
	 * \brief Compute the upwind flux between two nodes i and j.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetResidual (double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config);
};

/*! 
 * \class CUpwLin_LevelSet
 * \brief Class for performing a linear upwind solver for the Level Set equations.
 * \ingroup ConvDiscr
 * \author F. Palacios.
 * \version 1.0.
 */
class CUpwLin_LevelSet : public CNumerics {
private:
	bool implicit;
	double *Velocity_i;
	double *Velocity_j;

public:

	/*! 
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CUpwLin_LevelSet(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

	/*! 
	 * \brief Destructor of the class. 
	 */
	~CUpwLin_LevelSet(void);

	/*! 
	 * \brief Compute the upwind flux between two nodes i and j.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config);
};

/*! 
 * \class CUpwLin_Combustion
 * \brief Class for performing a linear upwind solver for the combustion equations.
 * \ingroup ConvDiscr
 */
class CUpwLin_Combustion : public CNumerics {
private:
	double *RhoVel_i, *RhoVel_j;
	
public:
	
	/*! 
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CUpwLin_Combustion(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
	
	/*! 
	 * \brief Destructor of the class. 
	 */
	~CUpwLin_Combustion(void);
	
	/*! 
	 * \brief Compute the upwind flux between two nodes i and j.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetResidual (double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config);
};

/*! 
 * \class CUpwLin_TurbSST
 * \brief Class for performing a linear upwind solver for the Menter SST turbulence model equations.
 * \ingroup ConvDiscr
 * \author A. Campos.
 * \version 1.0.
 */
class CUpwLin_TurbSST : public CNumerics {
private:
	double *Velocity_i;
	double *Velocity_j;

public:

	/*!
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CUpwLin_TurbSST(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

	/*!
	 * \brief Destructor of the class.
	 */
	~CUpwLin_TurbSST(void);

	/*!
	 * \brief Compute the upwind flux between two nodes i and j.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetResidual (double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config);
};

/*!
 * \class CUpwLin_AdjTurb
 * \brief Class for performing a linear upwind solver for the adjoint turbulence equations.
 * \ingroup ConvDiscr
 * \author A. Bueno.
 * \version 1.0.
 */
class CUpwLin_AdjTurb : public CNumerics {
private:
	double *Velocity_i;

public:

	/*! 
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CUpwLin_AdjTurb(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

	/*! 
	 * \brief Destructor of the class. 
	 */
	~CUpwLin_AdjTurb(void);

	/*! 
	 * \brief Compute the adjoint upwind flux between two nodes i and j.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetResidual (double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config);
};

/*! 
 * \class CUpwSca_TurbSA
 * \brief Class for doing a scalar upwind solver for the Spalar-Allmaral turbulence model equations.
 * \ingroup ConvDiscr
 * \author A. Bueno.
 * \version 1.0.
 */
class CUpwSca_TurbSA : public CNumerics {
private:
	double *Velocity_i, *Velocity_j;
	bool implicit;
	double Density_i, Density_j,
	q_ij,
	a0, a1;
	unsigned short iDim;

public:

	/*! 
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CUpwSca_TurbSA(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

	/*! 
	 * \brief Destructor of the class. 
	 */
	~CUpwSca_TurbSA(void);

	/*! 
	 * \brief Compute the scalar upwind flux between two nodes i and j.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config);
};

/*! 
 * \class CUpwSca_TurbSST
 * \brief Class for doing a scalar upwind solver for the Menter SST turbulence model equations.
 * \ingroup ConvDiscr
 * \author A. Campos.
 * \version 1.0.
 */
class CUpwSca_TurbSST : public CNumerics {
private:
	double *Velocity_i, *Velocity_j;
	bool implicit;
	double Density_i, Density_j,
	q_ij,
	a0, a1;
	unsigned short iDim;

public:

	/*!
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CUpwSca_TurbSST(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

	/*!
	 * \brief Destructor of the class.
	 */
	~CUpwSca_TurbSST(void);

	/*!
	 * \brief Compute the scalar upwind flux between two nodes i and j.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config);
};


/*!
 * \class CUpwSca_AdjTurb
 * \brief Class for doing a scalar upwind solver for the adjoint turbulence equations.
 * \ingroup ConvDiscr
 * \author A. Bueno.
 * \version 1.0.
 */
class CUpwSca_AdjTurb : public CNumerics {
private:
	double *Velocity_i, *Velocity_j;

public:

	/*! 
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CUpwSca_AdjTurb(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

	/*! 
	 * \brief Destructor of the class. 
	 */
	~CUpwSca_AdjTurb(void);

	/*! 
	 * \param[out] val_residual_i - Pointer to the total residual at point i.
	 * \param[out] val_residual_j - Pointer to the total viscosity residual at point j.
	 * \param[out] val_Jacobian_ii - Jacobian of the numerical method at node i (implicit computation) from node i.
	 * \param[out] val_Jacobian_ij - Jacobian of the numerical method at node i (implicit computation) from node j.
	 * \param[out] val_Jacobian_ji - Jacobian of the numerical method at node j (implicit computation) from node i.
	 * \param[out] val_Jacobian_jj - Jacobian of the numerical method at node j (implicit computation) from node j.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetResidual(double *val_residual_i, double *val_residual_j, double **val_Jacobian_ii, double **val_Jacobian_ij, 
			double **val_Jacobian_ji, double **val_Jacobian_jj, CConfig *config);
};

/*! 
 * \class CCentJST_Flow
 * \brief Class for centred shceme - JST.
 * \ingroup ConvDiscr
 * \author F. Palacios.
 * \version 1.0.
 */
class CCentJST_Flow : public CNumerics {

private:
	unsigned short iDim, iVar, jVar; /*!< \brief Iteration on dimension and variables. */
	double *Diff_U, *Diff_Lapl, /*!< \brief Diference of conservative variables and undivided laplacians. */
	*Velocity_i, *Velocity_j, /*!< \brief Velocity at node 0 and 1. */
	*MeanVelocity, ProjVelocity, ProjVelocity_i, ProjVelocity_j,  /*!< \brief Mean and projected velocities. */
	Density_i, Density_j, Energy_i, Energy_j,  /*!< \brief Mean Density and energies. */
	sq_vel_i, sq_vel_j,   /*!< \brief Modulus of the velocity and the normal vector. */
	MeanDensity, MeanPressure, MeanEnthalpy, MeanEnergy, /*!< \brief Mean values of primitive variables. */
	Param_p, Param_Kappa_2, Param_Kappa_4, /*!< \brief Artificial dissipation parameters. */
	Local_Lambda_i, Local_Lambda_j, MeanLambda, /*!< \brief Local eingenvalues. */
	Phi_i, Phi_j, sc2, sc4, StretchingFactor, /*!< \brief Streching parameters. */
	*Proj_flux_tensor,  /*!< \brief Projected inviscid flux tensor. */
	Epsilon_2, Epsilon_4, cte_0, cte_1; /*!< \brief Artificial dissipation values. */
	bool implicit, /*!< \brief Implicit calculation. */
	grid_movement, /*!< \brief Modification for grid movement. */
	stretching, /*!< \brief Stretching factor. */
	rotating_frame; /*!< \brief Rotational frame. */

public:

	/*! 
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimension of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CCentJST_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

	/*! 
	 * \brief Destructor of the class. 
	 */
	~CCentJST_Flow(void);

	/*! 
	 * \brief Compute the flow residual using a JST method.
	 * \param[out] val_resconv - Pointer to the convective residual.
	 * \param[out] val_resvisc - Pointer to the artificial viscosity residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] art_diss - If <code>TRUE</code>, the artificial dissipation is computed; otherwise <code>FALSE</code>.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetResidual(double *val_resconv, double *val_resvisc, double **val_Jacobian_i, double **val_Jacobian_j, 
			bool art_diss, CConfig *config);
};

/*! 
 * \class CCentJSTArtComp_Flow
 * \brief Class for centered scheme - JST (artificial compressibility).
 * \ingroup ConvDiscr
 * \author F. Palacios.
 * \version 1.0.
 */
class CCentJSTArtComp_Flow : public CNumerics {

private:
	unsigned short iDim, iVar, jVar; /*!< \brief Iteration on dimension and variables. */
	double *Diff_U, *Diff_Lapl, /*!< \brief Diference of conservative variables and undivided laplacians. */
	*Velocity_i, *Velocity_j, /*!< \brief Velocity at node 0 and 1. */
	*MeanVelocity, ProjVelocity, ProjVelocity_i, ProjVelocity_j,  /*!< \brief Mean and projected velocities. */
	sq_vel_i, sq_vel_j,   /*!< \brief Modulus of the velocity and the normal vector. */
	MeanGravityForce, MeanDensity, MeanPressure, MeanEnthalpy, MeanEnergy, MeanBetaInc2, /*!< \brief Mean values of primitive variables. */
	Param_p, Param_Kappa_2, Param_Kappa_4, /*!< \brief Artificial dissipation parameters. */
	Local_Lambda_i, Local_Lambda_j, MeanLambda, /*!< \brief Local eingenvalues. */
	Phi_i, Phi_j, sc2, sc4, StretchingFactor, /*!< \brief Streching parameters. */
	*Proj_flux_tensor,  /*!< \brief Projected inviscid flux tensor. */
	Epsilon_2, Epsilon_4, cte_0, cte_1; /*!< \brief Artificial dissipation values. */
	bool implicit, /*!< \brief Implicit calculation. */
	grid_movement, /*!< \brief Modification for grid movement. */
	stretching, /*!< \brief Stretching factor. */
	rotating_frame, /*!< \brief Rotational frame. */
	gravity; /*!< \brief computation with gravity force. */
	double Froude; /*!< \brief Froude number. */

public:

	/*! 
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimension of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CCentJSTArtComp_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

	/*! 
	 * \brief Destructor of the class. 
	 */
	~CCentJSTArtComp_Flow(void);

	/*! 
	 * \brief Compute the flow residual using a JST method.
	 * \param[out] val_resconv - Pointer to the convective residual.
	 * \param[out] val_resvisc - Pointer to the artificial viscosity residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] art_diss - If <code>TRUE</code>, the artificial dissipation is computed; otherwise <code>FALSE</code>.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetResidual(double *val_resconv, double *val_resvisc, double **val_Jacobian_i, double **val_Jacobian_j, 
			bool art_diss, CConfig *config);
};

/*! 
 * \class CCentJST_AdjFlow
 * \brief Class for and adjoint centred scheme - JST.
 * \ingroup ConvDiscr
 * \author F. Palacios.
 * \version 1.0.
 */
class CCentJST_AdjFlow : public CNumerics {
private:
	double *Diff_Psi, *Diff_Lapl;
	double *Velocity_i, *Velocity_j;
	double *MeanPhi;
	unsigned short iDim, jDim, iVar, jVar;
	double Residual, ProjVelocity_i, ProjVelocity_j, ProjPhi, ProjPhi_Vel, sq_vel, phis1, phis2;
	double MeanPsiRho, MeanPsiE, Param_p, Param_Kappa_4, Param_Kappa_0, Local_Lambda_i, Local_Lambda_j, MeanLambda;
	double Phi_i, Phi_j, sc4, StretchingFactor, Epsilon_4, Epsilon_0;
	bool implicit, stretching, grid_movement, rotating_frame;

public:

	/*! 
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CCentJST_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

	/*! 
	 * \brief Destructor of the class. 
	 */
	~CCentJST_AdjFlow(void);

	/*! 
	 * \brief Compute the adjoint flow residual using a JST method.
	 * \param[out] val_resconv_i - Pointer to the convective residual at point i.
	 * \param[out] val_resvisc_i - Pointer to the artificial viscosity residual at point i.
	 * \param[out] val_resconv_j - Pointer to the convective residual at point j.
	 * \param[out] val_resvisc_j - Pointer to the artificial viscosity residual at point j.
	 * \param[out] val_Jacobian_ii - Jacobian of the numerical method at node i (implicit computation) from node i.
	 * \param[out] val_Jacobian_ij - Jacobian of the numerical method at node i (implicit computation) from node j.
	 * \param[out] val_Jacobian_ji - Jacobian of the numerical method at node j (implicit computation) from node i.
	 * \param[out] val_Jacobian_jj - Jacobian of the numerical method at node j (implicit computation) from node j.
	 * \param[in] local_art_diss - If 3 then 3th order dissipation, if 1, then 1st order dissipation; if 0, no dissipation.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetResidual (double *val_resconv_i, double *val_resvisc_i, double *val_resconv_j, double *val_resvisc_j, 
			double **val_Jacobian_ii, double **val_Jacobian_ij, double **val_Jacobian_ji, double **val_Jacobian_jj,
			unsigned short local_art_diss, CConfig *config);
};

/*! 
 * \class CCentJST_LinFlow
 * \brief Class for linearized centred scheme - JST.
 * \ingroup ConvDiscr
 * \author F. Palacios.
 * \version 1.0.
 */
class CCentJST_LinFlow : public CNumerics {
private:
	double *Diff_DeltaU, *Diff_Lapl;
	double *Velocity_i, *Velocity_j;
	double *MeanDeltaVel, *MeanVelocity;
	double **MeanJacobian;
	double **Jacobian_i, **Jacobian_j;
	unsigned short iDim, iVar, jVar;
	double sq_vel, Density_i, DensityEnergy_i, Energy_i, Pressure_i, Density_j, DensityEnergy_j, Energy_j, 
	Pressure_j, Param_p, Param_Kappa_4, Local_Lambda_i, Local_Lambda_j, MeanLambda, sc4, StretchingFactor, 
	Epsilon_4, MeanDeltaRho, MeanDeltaE, ProjVelocity_i, ProjVelocity_j, MeanDensity, MeanPressure, 
	MeanEnthalpy, MeanEnergy, Phi_i, Phi_j;
	bool stretching;


public:

	/*! 
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CCentJST_LinFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

	/*! 
	 * \brief Destructor of the class. 
	 */
	~CCentJST_LinFlow(void);

	/*! 
	 * \brief Compute the linearized flow residual using a JST method.
	 * \param[out] val_resconv - Pointer to the convective residual.
	 * \param[out] val_resvisc - Pointer to the artificial viscosity residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] art_diss - If <code>TRUE</code>, the artificial dissipation is computed; otherwise <code>FALSE</code>.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetResidual (double *val_resconv, double *val_resvisc, double **val_Jacobian_i, double **val_Jacobian_j, bool art_diss, CConfig *config);
};

/*! 
 * \class CCentLaxComp_Flow
 * \brief Class for computing the Lax-Friedrich centred scheme.
 * \ingroup ConvDiscr
 * \author F. Palacios.
 * \version 1.0.
 */
class CCentLaxComp_Flow : public CNumerics {
private:
	unsigned short iDim, iVar, jVar; /*!< \brief Iteration on dimension and variables. */
	double *Diff_U, /*!< \brief Difference of conservative variables. */
	*Velocity_i, *Velocity_j, /*!< \brief Velocity at node 0 and 1. */
	*MeanVelocity, ProjVelocity, ProjVelocity_i, ProjVelocity_j,  /*!< \brief Mean and projected velocities. */
	*Proj_flux_tensor,  /*!< \brief Projected inviscid flux tensor. */
	Density_i, Density_j, Energy_i, Energy_j,  /*!< \brief Mean Density and energies. */
	sq_vel_i, sq_vel_j,   /*!< \brief Modulus of the velocity and the normal vector. */
	MeanDensity, MeanPressure, MeanEnthalpy, MeanEnergy, /*!< \brief Mean values of primitive variables. */
	Param_p, Param_Kappa_0, /*!< \brief Artificial dissipation parameters. */
	Local_Lambda_i, Local_Lambda_j, MeanLambda, /*!< \brief Local eingenvalues. */
	Phi_i, Phi_j, sc0, StretchingFactor, /*!< \brief Streching parameters. */
	Epsilon_0, cte; /*!< \brief Artificial dissipation values. */
	bool implicit, /*!< \brief Implicit calculation. */
	grid_movement, /*!< \brief Modification for grid movement. */
	rotating_frame; /*!< \brief Rotational frame. */
	bool stretching;

public:

	/*! 
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimension of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CCentLaxComp_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

	/*! 
	 * \brief Destructor of the class. 
	 */
	~CCentLaxComp_Flow(void);

	/*! 
	 * \brief Compute the flow residual using a Lax method.
	 * \param[out] val_resconv - Pointer to the convective residual.
	 * \param[out] val_resvisc - Pointer to the artificial viscosity residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] art_diss - If <code>TRUE</code>, the artificial dissipation is computed; otherwise <code>FALSE</code>.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetResidual(double *val_resconv, double *val_resvisc, double **val_Jacobian_i, double **val_Jacobian_j, 
			bool art_diss, CConfig *config);
};

/*! 
 * \class CCentLaxArtComp_Flow
 * \brief Class for computing the Lax-Friedrich centered scheme (artificial compressibility).
 * \ingroup ConvDiscr
 * \author F. Palacios.
 * \version 1.0.
 */
class CCentLaxArtComp_Flow : public CNumerics {
private:
	unsigned short iDim, iVar, jVar; /*!< \brief Iteration on dimension and variables. */
	double *Diff_U, /*!< \brief Difference of conservative variables. */
	*Velocity_i, *Velocity_j, /*!< \brief Velocity at node 0 and 1. */
	*MeanVelocity, ProjVelocity, ProjVelocity_i, ProjVelocity_j,  /*!< \brief Mean and projected velocities. */
	*Proj_flux_tensor,  /*!< \brief Projected inviscid flux tensor. */
	sq_vel_i, sq_vel_j,   /*!< \brief Modulus of the velocity and the normal vector. */
	MeanGravityForce, MeanDensity, MeanPressure, MeanEnthalpy, MeanEnergy, MeanBetaInc2, /*!< \brief Mean values of primitive variables. */
	Param_p, Param_Kappa_0, /*!< \brief Artificial dissipation parameters. */
	Local_Lambda_i, Local_Lambda_j, MeanLambda, /*!< \brief Local eingenvalues. */
	Phi_i, Phi_j, sc0, StretchingFactor, /*!< \brief Streching parameters. */
	Epsilon_0, cte; /*!< \brief Artificial dissipation values. */
	bool implicit, /*!< \brief Implicit calculation. */
	grid_movement, /*!< \brief Modification for grid movement. */
	gravity, /*!< \brief Modification for for gravity force. */
	rotating_frame; /*!< \brief Rotational frame. */
	bool stretching;
	double Froude;

public:

	/*! 
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimension of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CCentLaxArtComp_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

	/*! 
	 * \brief Destructor of the class. 
	 */
	~CCentLaxArtComp_Flow(void);

	/*! 
	 * \brief Compute the flow residual using a Lax method.
	 * \param[out] val_resconv - Pointer to the convective residual.
	 * \param[out] val_resvisc - Pointer to the artificial viscosity residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] art_diss - If <code>TRUE</code>, the artificial dissipation is computed; otherwise <code>FALSE</code>.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetResidual(double *val_resconv, double *val_resvisc, double **val_Jacobian_i, double **val_Jacobian_j, 
			bool art_diss, CConfig *config);
};

/*! 
 * \class CCentLax_AdjFlow
 * \brief Class for computing the Lax-Friedrich adjoint centred scheme.
 * \ingroup ConvDiscr
 * \author F. Palacios.
 * \version 1.0.
 */
class CCentLax_AdjFlow : public CNumerics {
private:
	double *Diff_Psi;
	double *Velocity_i, *Velocity_j;
	double *MeanPhi;
	unsigned short iDim, jDim, iVar, jVar;
	double Residual, ProjVelocity_i, ProjVelocity_j, ProjPhi, ProjPhi_Vel, sq_vel, phis1, phis2, 
	MeanPsiRho, MeanPsiE, Param_p, Param_Kappa_0, Local_Lambda_i, Local_Lambda_j, MeanLambda, 
	Phi_i, Phi_j, sc2, StretchingFactor, Epsilon_0, cte_0;
	bool implicit;
	bool stretching;
  bool rotating_frame;

public:

	/*! 
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CCentLax_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

	/*! 
	 * \brief Destructor of the class. 
	 */
	~CCentLax_AdjFlow(void);

	/*! 
	 * \brief Compute the adjoint flow residual using a Lax method.
	 * \param[out] val_resconv_i - Pointer to the convective residual at point i.
	 * \param[out] val_resvisc_i - Pointer to the artificial viscosity residual at point i.
	 * \param[out] val_resconv_j - Pointer to the convective residual at point j.
	 * \param[out] val_resvisc_j - Pointer to the artificial viscosity residual at point j.
	 * \param[out] val_Jacobian_ii - Jacobian of the numerical method at node i (implicit computation) from node i.
	 * \param[out] val_Jacobian_ij - Jacobian of the numerical method at node i (implicit computation) from node j.
	 * \param[out] val_Jacobian_ji - Jacobian of the numerical method at node j (implicit computation) from node i.
	 * \param[out] val_Jacobian_jj - Jacobian of the numerical method at node j (implicit computation) from node j.
	 * \param[in] local_art_diss - If 3 then 3th order dissipation, if 1, then 1st order dissipation; if 0, no dissipation.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetResidual (double *val_resconv_i, double *val_resvisc_i, double *val_resconv_j, double *val_resvisc_j, 
			double **val_Jacobian_ii, double **val_Jacobian_ij, double **val_Jacobian_ji, double **val_Jacobian_jj,
			unsigned short local_art_diss, CConfig *config);
};

/*! 
 * \class CCentLax_LinFlow
 * \brief Class for computing the Lax-Friedrich linearized centred scheme.
 * \ingroup ConvDiscr
 * \author F. Palacios.
 * \version 1.0.
 */
class CCentLax_LinFlow : public CNumerics {
private:
	double *Diff_DeltaU;
	double *Velocity_i, *Velocity_j;
	double *MeanDeltaVel, *MeanVelocity;
	double **MeanJacobian;
	double **Jacobian_i;
	double **Jacobian_j;
	unsigned short iDim, iVar, jVar;
	double sq_vel, Density_i, DensityEnergy_i, Energy_i, Pressure_i, Density_j, 
	DensityEnergy_j, Energy_j,Pressure_j, Param_p, Param_Kappa_0, 
	Local_Lambda_i, Local_Lambda_j, MeanLambda, cte_0, StretchingFactor, 
	Epsilon_i, MeanDeltaRho, MeanDeltaE, ProjVelocity_i, ProjVelocity_j, 
	dS, MeanDensity, MeanPressure, 
	MeanEnthalpy, MeanEnergy, Phi_i, Phi_j, 
	sc2;
	bool stretching;

public:

	/*! 
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CCentLax_LinFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

	/*! 
	 * \brief Destructor of the class. 
	 */
	~CCentLax_LinFlow(void);

	/*! 
	 * \brief Compute the linearized flow residual using a Lax method.
	 * \param[out] val_resconv - Pointer to the convective residual.
	 * \param[out] val_resvisc - Pointer to the artificial viscosity residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] art_diss - If <code>TRUE</code>, the artificial dissipation is computed; otherwise <code>FALSE</code>.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetResidual(double *val_resconv, double *val_resvisc, double **val_Jacobian_i, double **val_Jacobian_j, bool art_diss, CConfig *config);
};

/*! 
 * \class CAvgGrad_Flow
 * \brief Class for computing viscous term using the average of gradients.
 * \ingroup ViscDiscr
 * \author A. Bueno, and F. Palacios.
 * \version 1.0.
 */
class CAvgGrad_Flow : public CNumerics {
private:
	unsigned short iDim, iVar; /*!< \brief Iterators in dimension an variable. */
	double *Mean_PrimVar, /*!< \brief Mean primitive variables. */
	*Velocity_i, *Velocity_j, /*!< \brief Velocity at point i and 1. */
	*Prim_Var_i, *Prim_Var_j, /*!< \brief Primitives variables at point i and 1. */
	**Mean_GradPrimVar, /*!< \brief Mean value of the gradient. */
	Mean_Laminar_Viscosity, Mean_Eddy_Viscosity, /*!< \brief Mean value of the viscosity. */
	Mean_Density,  /*!< \brief Mean density. */
	*Proj_flux_tensor, /*!< \brief Projection of the viscous fluxes. */
	dist_ij; /*!< \brief Length of the edge and face. */
	bool implicit; /*!< \brief Implicit calculus. */

public:

	/*! 
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimension of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CAvgGrad_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

	/*! 
	 * \brief Destructor of the class. 
	 */
	~CAvgGrad_Flow(void);
	/*! 
	 * \brief Compute the viscous flow residual using an average of gradients.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config);
};

/*! 
 * \class CAvgGradArtComp_Flow
 * \brief Class for computing viscous term using an average of gradients.
 * \ingroup ViscDiscr
 * \author A. Bueno, and F. Palacios.
 * \version 1.0.
 */
class CAvgGradArtComp_Flow : public CNumerics {
private:
	unsigned short iDim, iVar; /*!< \brief Iterators in dimension an variable. */
	double *Mean_PrimVar, /*!< \brief Mean primitive variables. */
	*Velocity_i, *Velocity_j, /*!< \brief Velocity at point i and 1. */
	*Prim_Var_i, *Prim_Var_j, /*!< \brief Primitives variables at point i and 1. */
	**Mean_GradPrimVar, /*!< \brief Mean value of the gradient. */
	Mean_Laminar_Viscosity, Mean_Eddy_Viscosity, /*!< \brief Mean value of the viscosity. */
	Mean_Density,  /*!< \brief Mean density. */
	*Proj_flux_tensor, /*!< \brief Projection of the viscous fluxes. */
	dist_ij; /*!< \brief Length of the edge and face. */
	bool implicit; /*!< \brief Implicit calculus. */

public:

	/*! 
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimension of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CAvgGradArtComp_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

	/*! 
	 * \brief Destructor of the class. 
	 */
	~CAvgGradArtComp_Flow(void);
	/*! 
	 * \brief Compute the viscous flow residual using an average of gradients.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config);
};

/*! 
 * \class CAvgGrad_Combustion
 * \brief Class for computing viscous term using average of gradients (Spalart-Allmaras Turbulence model).
 * \ingroup ViscDiscr
 */
class CAvgGrad_Combustion : public CNumerics {
private:
	
public:
	
	/*! 
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CAvgGrad_Combustion(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
	
	/*! 
	 * \brief Destructor of the class. 
	 */
	~CAvgGrad_Combustion(void);
	
	/*! 
	 * \brief Compute the viscous turbulence terms residual using an average of gradients.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetResidual(double *val_residual, double **Jacobian_i, double **Jacobian_j, CConfig *config);
};

/*! 
 * \class CAvgGrad_TurbSA
 * \brief Class for computing viscous term using average of gradients (Spalart-Allmaras Turbulence model).
 * \ingroup ViscDiscr
 * \author A. Bueno.
 * \version 1.0.
 */
class CAvgGrad_TurbSA : public CNumerics {
private:
	double **Mean_GradTurbVar;
	double *Proj_Mean_GradTurbVar_Kappa, *Proj_Mean_GradTurbVar_Edge;
	double *Edge_Vector;
	bool implicit;
	double sigma;
	double nu_i, nu_j, nu_e;
	double dist_ij_2;
	double proj_vector_ij;
	unsigned short iVar, iDim;
	double nu_hat_i;
	double nu_hat_j;

public:

	/*! 
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CAvgGrad_TurbSA(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

	/*! 
	 * \brief Destructor of the class. 
	 */
	~CAvgGrad_TurbSA(void);

	/*! 
	 * \brief Compute the viscous turbulence terms residual using an average of gradients.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetResidual(double *val_residual, double **Jacobian_i, double **Jacobian_j, CConfig *config);
};

/*! 
 * \class CAvgGrad_TurbSST
 * \brief Class for computing viscous term using average of gradients (Menter SST Turbulence model).
 * \ingroup ViscDiscr
 * \author A. Bueno.
 * \version 1.0.
 */
class CAvgGrad_TurbSST : public CNumerics {
private:
	double **Mean_GradTurbVar;
	double *Proj_Mean_GradTurbVar_Kappa, *Proj_Mean_GradTurbVar_Edge;
	double *Edge_Vector;
	bool implicit;
	double sigma;
	double nu_i, nu_j, nu_e;
	double dist_ij_2;
	double proj_vector_ij;
	unsigned short iVar, iDim;
	double nu_hat_i;
	double nu_hat_j;

public:

	/*!
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CAvgGrad_TurbSST(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

	/*!
	 * \brief Destructor of the class.
	 */
	~CAvgGrad_TurbSST(void);

	/*!
	 * \brief Compute the viscous turbulence terms residual using an average of gradients.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetResidual(double *val_residual, double **Jacobian_i, double **Jacobian_j, CConfig *config);
};

/*!
 * \class CAvgGrad_AdjFlow
 * \brief Class for computing the adjoint viscous terms.
 * \ingroup ViscDiscr
 * \author F. Palacios.
 * \version 1.0.
 */
class CAvgGrad_AdjFlow : public CNumerics {
private:
	double *Velocity_i;	/*!< \brief Auxiliary vector for storing the velocity of point i. */
	double *Velocity_j;	/*!< \brief Auxiliary vector for storing the velocity of point j. */
	double *Mean_Velocity;
	double *Mean_GradPsiE;	/*!< \brief Counter for dimensions of the problem. */
	double **Mean_GradPhi;	/*!< \brief Counter for dimensions of the problem. */
public:

	/*! 
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.	 
	 */	
	CAvgGrad_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

	/*! 
	 * \brief Destructor of the class.
	 */	
	~CAvgGrad_AdjFlow(void);

	/*! 
	 * \brief Residual computation.
	 * \param[out] val_residual_i - Pointer to the total residual at point i.
	 * \param[out] val_residual_j - Pointer to the total residual at point j.
	 */	
	void SetResidual(double *val_residual_i, double *val_residual_j);
};

/*! 
 * \class CAvgGradCorrected_Flow
 * \brief Class for computing viscous term using the average of gradients with a correction.
 * \ingroup ViscDiscr
 * \author A. Bueno, and F. Palacios.
 * \version 1.0.
 */
class CAvgGradCorrected_Flow : public CNumerics {
private:
	unsigned short iDim, iVar; /*!< \brief Iterators in dimension an variable. */
	double *Mean_PrimVar, /*!< \brief Mean primitive variables. */
	*Velocity_i, *Velocity_j, /*!< \brief Velocity at point i and 1. */
	*Prim_Var_i, *Prim_Var_j, /*!< \brief Primitives variables at point i and 1. */
	*Edge_Vector,  /*!< \brief Vector form point i to point j. */
	**Mean_GradPrimVar, *Proj_Mean_GradPrimVar_Edge, /*!< \brief Mean value of the gradient. */
	Mean_Laminar_Viscosity, Mean_Eddy_Viscosity, /*!< \brief Mean value of the viscosity. */
	dist_ij_2, /*!< \brief Length of the edge and face. */
	*Proj_flux_tensor; /*!< \brief Projection of the viscous fluxes. */
	bool implicit; /*!< \brief Implicit calculus. */

public:

	/*! 
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimension of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CAvgGradCorrected_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

	/*! 
	 * \brief Destructor of the class. 
	 */
	~CAvgGradCorrected_Flow(void);

	/*! 
	 * \brief Compute the viscous flow residual using an average of gradients with correction.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config);
};

/*! 
 * \class CAvgGradCorrectedArtComp_Flow
 * \brief Class for computing viscous term using an average of gradients with correction (artificial compresibility).
 * \ingroup ViscDiscr
 * \author F. Palacios.
 * \version 1.0.
 */
class CAvgGradCorrectedArtComp_Flow : public CNumerics {
private:
	unsigned short iDim, iVar; /*!< \brief Iterators in dimension an variable. */
	double *Mean_PrimVar, /*!< \brief Mean primitive variables. */
	*Velocity_i, *Velocity_j, /*!< \brief Velocity at point i and 1. */
	*Prim_Var_i, *Prim_Var_j, /*!< \brief Primitives variables at point i and 1. */
	*Edge_Vector,  /*!< \brief Vector form point i to point j. */
	**Mean_GradPrimVar, *Proj_Mean_GradPrimVar_Edge, /*!< \brief Mean value of the gradient. */
	Mean_Laminar_Viscosity, Mean_Eddy_Viscosity, /*!< \brief Mean value of the viscosity. */
	Mean_Density, /*!< \brief Mean value of the density. */
	dist_ij_2, /*!< \brief Length of the edge and face. */
	*Proj_flux_tensor; /*!< \brief Projection of the viscous fluxes. */
	bool implicit; /*!< \brief Implicit calculus. */

public:

	/*! 
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimension of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CAvgGradCorrectedArtComp_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

	/*! 
	 * \brief Destructor of the class. 
	 */
	~CAvgGradCorrectedArtComp_Flow(void);

	/*! 
	 * \brief Compute the viscous flow residual using an average of gradients with correction.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config);
};

/*! 
 * \class CAvgGradCorrected_TurbSA
 * \brief Class for computing viscous term using average of gradients with correction (Spalart-Allmaras turbulence model).
 * \ingroup ViscDiscr
 * \author A. Bueno.
 * \version 1.0.
 */
class CAvgGradCorrected_TurbSA : public CNumerics {
private:
	double **Mean_GradTurbVar;
	double *Proj_Mean_GradTurbVar_Kappa, *Proj_Mean_GradTurbVar_Edge, *Proj_Mean_GradTurbVar_Corrected;
	double *Edge_Vector;
	bool implicit;
	double sigma, nu_i, nu_j, nu_e, dist_ij_2, proj_vector_ij, nu_hat_i, nu_hat_j;
	unsigned short iVar, iDim;

public:

	/*! 
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CAvgGradCorrected_TurbSA(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

	/*! 
	 * \brief Destructor of the class. 
	 */
	~CAvgGradCorrected_TurbSA(void);

	/*! 
	 * \brief Compute the viscous turbulent residual using an average of gradients with correction.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetResidual(double *val_residual, double **Jacobian_i, double **Jacobian_j, CConfig *config);
};

/*! 
 * \class CAvgGradCorrected_TurbSST
 * \brief Class for computing viscous term using average of gradient with correction (Menter SST turbulence model).
 * \ingroup ViscDiscr
 * \author A. Bueno.
 * \version 1.0.
 */
class CAvgGradCorrected_TurbSST : public CNumerics {
private:
	double **Mean_GradTurbVar;
	double *Proj_Mean_GradTurbVar_Kappa, *Proj_Mean_GradTurbVar_Edge, *Proj_Mean_GradTurbVar_Corrected, *Edge_Vector;
	double dist_ij_2, proj_vector_ij;
	double sigma_kine_1, sigma_omega_1, sigma_kine_2, sigma_omega_2;
	bool implicit;
	unsigned short iVar, iDim;

public:

	/*!
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CAvgGradCorrected_TurbSST(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

	/*!
	 * \brief Destructor of the class.
	 */
	~CAvgGradCorrected_TurbSST(void);

	/*!
	 * \brief Compute the viscous turbulent residual using an average of gradients wtih correction.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetResidual(double *val_residual, double **Jacobian_i, double **Jacobian_j, CConfig *config);
};

/*!
 * \class CAvgGradCorrected_AdjFlow
 * \brief Class for computing the adjoint viscous terms, including correction.
 * \ingroup ViscDiscr
 * \author A. Bueno.
 * \version 1.0.
 */
class CAvgGradCorrected_AdjFlow : public CNumerics {
private:
	double *Velocity_i;	/*!< \brief Auxiliary vector for storing the velocity of point i. */
	double *Velocity_j;	/*!< \brief Auxiliary vector for storing the velocity of point j. */
	double *Mean_Velocity;
	double **Mean_GradPsiVar;	/*!< \brief Counter for dimensions of the problem. */
	double *Edge_Vector;	/*!< \brief Vector going from node 0 to node 1. */
	double *Proj_Mean_GradPsiVar_Edge;	/*!< \brief Projection of Mean_GradPsiVar onto Edge_Vector. */
	double *Mean_GradPsiE;	/*!< \brief Counter for dimensions of the problem. */
	double **Mean_GradPhi;	/*!< \brief Counter for dimensions of the problem. */

public:

	/*! 
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.	 
	 */	
	CAvgGradCorrected_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

	/*! 
	 * \brief Destructor of the class.
	 */	
	~CAvgGradCorrected_AdjFlow(void);

	/*! 
	 * \brief Compute the adjoint flow viscous residual in a non-conservative way using an average of gradients and derivative correction.
	 * \param[out] val_residual_i - Pointer to the total residual at point i.
	 * \param[out] val_residual_j - Pointer to the total residual at point j.
	 */	
	void SetResidual(double *val_residual_i, double *val_residual_j);

	/*! 
	 * \brief Compute the adjoint flow viscous residual in a non-conservative way using an average of gradients and derivative correction.
	 * \param[out] val_residual_i - Pointer to the viscous residual at point i.
	 * \param[out] val_residual_j - Pointer to the viscous residual at point j.
	 * \param[out] val_Jacobian_ii - Jacobian of the numerical method at node i (implicit computation) from node i.
	 * \param[out] val_Jacobian_ij - Jacobian of the numerical method at node i (implicit computation) from node j.
	 * \param[out] val_Jacobian_ji - Jacobian of the numerical method at node j (implicit computation) from node i.
	 * \param[out] val_Jacobian_jj - Jacobian of the numerical method at node j (implicit computation) from node j.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetResidual (double *val_residual_i, double *val_residual_j, double **val_Jacobian_ii, double **val_Jacobian_ij, 
			double **val_Jacobian_ji, double **val_Jacobian_jj, CConfig *config);
};

/*! 
 * \class CAvgGradCorrected_AdjTurb
 * \brief Class for adjoint turbulent using average of gradients with a correction.
 * \ingroup ViscDiscr
 * \author A. Bueno.
 * \version 1.0.
 */
class CAvgGradCorrected_AdjTurb : public CNumerics {
private:
	double **Mean_GradTurbPsi;
	double *Proj_Mean_GradTurbPsi_Kappa, *Proj_Mean_GradTurbPsi_Edge, *Proj_Mean_GradTurbPsi_Corrected;
	double *Edge_Vector;

public:

	/*! 
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.	 
	 */
	CAvgGradCorrected_AdjTurb(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

	/*! 
	 * \brief Destructor of the class. 
	 */
	~CAvgGradCorrected_AdjTurb(void);

	/*! 
	 * \brief Compute the adjoint turbulent residual using average of gradients and a derivative correction.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */

	void SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config);

	/*! 
	 * \overload
	 * \param[out] val_residual_i - Pointer to the total residual at point i.
	 * \param[out] val_residual_j - Pointer to the total viscosity residual at point j.
	 * \param[out] val_Jacobian_ii - Jacobian of the numerical method at node i (implicit computation) from node i.
	 * \param[out] val_Jacobian_ij - Jacobian of the numerical method at node i (implicit computation) from node j.
	 * \param[out] val_Jacobian_ji - Jacobian of the numerical method at node j (implicit computation) from node i.
	 * \param[out] val_Jacobian_jj - Jacobian of the numerical method at node j (implicit computation) from node j.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetResidual(double *val_residual_i, double *val_residual_j, double **val_Jacobian_ii, double **val_Jacobian_ij, 
			double **val_Jacobian_ji, double **val_Jacobian_jj, CConfig *config);
};

/*! 
 * \class CAvgGrad_Plasma
 * \brief Class for computing viscous term using average of gradients.
 * \ingroup ViscDiscr
 * \author ADL Stanford.
 * \version 1.0.
 */
class CAvgGrad_Plasma : public CNumerics {
private:
	unsigned short iDim, iVar; /*!< \brief Iterators in dimension an variable. */
	double *Mean_PrimVar, /*!< \brief Mean primitive variables. */
	*Velocity_i, *Velocity_j, /*!< \brief Velocity at point i and 1. */
	*Prim_Var_i, *Prim_Var_j, /*!< \brief Primitives variables at point i and 1. */
	**Mean_GradPrimVar, /*!< \brief Mean value of the gradient. */
	Mean_Laminar_Viscosity, Mean_Eddy_Viscosity, /*!< \brief Mean value of the viscosity. */
	Mean_Density,  /*!< \brief Mean density. */
	*Proj_flux_tensor, /*!< \brief Projection of the viscous fluxes. */
	dist_ij, dS; /*!< \brief Length of the edge and face. */
	bool implicit; /*!< \brief Implicit calculus. */

public:

	/*! 
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimension of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CAvgGrad_Plasma(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

	/*! 
	 * \brief Destructor of the class. 
	 */
	~CAvgGrad_Plasma(void);
	/*! 
	 * \brief Compute the viscous flow residual using an average of gradients.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config);
};

/*! 
 * \class CGalerkin_Flow
 * \brief Class for computing the stiffness matrix of the Galerkin method.
 * \ingroup ViscDiscr
 * \author F. Palacios.
 * \version 1.0.
 */
class CGalerkin_Flow : public CNumerics {
public:

	/*! 
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CGalerkin_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

	/*! 
	 * \brief Destructor of the class. 
	 */
	~CGalerkin_Flow(void);

	/*! 
	 * \brief Computing stiffness matrix of the Galerkin method.
	 * \param[out] val_stiffmatrix_elem - Stiffness matrix for Galerkin computation.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetResidual (double **val_stiffmatrix_elem, CConfig *config);
};

/*! 
 * \class CSourceNothing
 * \brief Dummy class.
 * \ingroup SourceDiscr
 * \author F. Palacios.
 * \version 1.0.
 */
class CSourceNothing : public CNumerics {
public:

	/*! 
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CSourceNothing(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

	/*! 
	 * \brief Destructor of the class. 
	 */
	~CSourceNothing(void);
};

/*! 
 * \class CSourcePieceWise_TurbSA
 * \brief Class for integrating the source terms of the Spalart-Allmaras turbulence model equation.
 * \ingroup SourceDiscr
 * \author A. Bueno.
 * \version 1.0.
 */
class CSourcePieceWise_TurbSA : public CNumerics {
private:
	double cv1_3;
	double k2;
	double cb1;
	double cw2;
	double cw3_6;
	double sigma;
	double cb2;
	double cw1;
	double DivVelocity, Vorticity;
	unsigned short iDim;
	double nu, Ji, fv1, fv2, Omega, Shat, dist_0_2, Ji_2, Ji_3;
	double r, g, g_6, glim, fw;
	double norm2_Grad;
	double dfv1, dfv2, dShat;
	double dr, dg, dfw;;
	double nu_hat_i;
	double grad_nu_hat;
	double prod_grads;
	bool implicit;

public:

	/*! 
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CSourcePieceWise_TurbSA(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

	/*! 
	 * \brief Destructor of the class. 
	 */
	~CSourcePieceWise_TurbSA(void);

	/*! 
	 * \brief Residual for source term integration.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config);
};

/*! 
 * \class CSourcePieceWise_TurbSST
 * \brief Class for integrating the source terms of the Menter SST turbulence model equations.
 * \ingroup SourceDiscr
 * \author A. Campos.
 * \version 1.0.
 */
class CSourcePieceWise_TurbSST : public CNumerics {
private:
	double gamma_1, gamma_2;
	double beta_1, beta_2;
	double sigma_omega_1, sigma_omega_2;
	double von_Karman;
	double beta_star;
	double **strain_rate, **reynolds_stress;
	double DivVelocity, Vorticity;
	double norm2_Grad;
	bool implicit;

public:

	/*!
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CSourcePieceWise_TurbSST(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

	/*!
	 * \brief Destructor of the class.
	 */
	~CSourcePieceWise_TurbSST(void);

	/*!
	 * \brief Residual for source term integration.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config);
};

/*!
 * \class CSourcePieceWise_Gravity
 * \brief Class for the source term integration of the gravity force.
 * \ingroup SourceDiscr
 * \author F. Palacios
 * \version 1.0.
 */
class CSourcePieceWise_Gravity : public CNumerics {
	double U_ref, L_ref, Froude;
	bool implicit, incompressible;
	
public:
	
	/*! 
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CSourcePieceWise_Gravity(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
	
	/*! 
	 * \brief Destructor of the class. 
	 */
	~CSourcePieceWise_Gravity(void);
	
	/*! 
	 * \brief Source term integration for the electrical potential.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetResidual(double *val_residual, CConfig *config);
};

/*!
 * \class CSourcePieceWise_Combustion
 * \brief Class for the source term integration of the combustion source term.
 * \ingroup SourceDiscr
 * \version 1.0.
 */
class CSourcePieceWise_Combustion : public CNumerics {
	double Alpha, T_i, Gas_Constant;
	bool implicit;
	
public:
	
	/*! 
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CSourcePieceWise_Combustion(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
	
	/*! 
	 * \brief Destructor of the class. 
	 */
	~CSourcePieceWise_Combustion(void);
	
	/*! 
	 * \brief Source term integration for the electrical potential.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetResidual(double *val_residual, CConfig *config);
};

/*!
 * \class CSourcePieceWise_Elec
 * \brief Class for the soruce term integration of the electrical potential equation.
 * \ingroup SourceDiscr
 * \author A. Bueno.
 * \version 1.0.
 */
class CSourcePieceWise_Elec : public CNumerics {
public:

	/*! 
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CSourcePieceWise_Elec(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

	/*! 
	 * \brief Destructor of the class. 
	 */
	~CSourcePieceWise_Elec(void);

	/*! 
	 * \brief Source term integration for the electrical potential.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetResidual(double *val_residual, CConfig *config);
};

/*! 
 * \class CSourcePieceWise_AdjFlow
 * \brief Class for source term integration in adjoint problem.
 * \ingroup SourceDiscr
 * \author F. Palacios.
 * \version 1.0.
 */
class CSourcePieceWise_AdjFlow : public CNumerics {
private:
	double *Velocity, *GradDensity, *GradInvDensity, *dPoDensity2, *alpha, *beta, *Sigma_5_vec;
	double **GradVel_o_Rho, **sigma, **Sigma_phi, **Sigma_5_Tensor, **Sigma;
public:

	/*! 
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CSourcePieceWise_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

	/*! 
	 * \brief Destructor of the class. 
	 */
	~CSourcePieceWise_AdjFlow(void);

	/*! 
	 * \brief Source term integration of the flow adjoint equation.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetResidual (double *val_residual, CConfig *config);
};

/*! 
 * \class CSourcePieceWise_AdjTurb
 * \brief Class for source term integration of the adjoint turbulent equation.
 * \ingroup SourceDiscr
 * \author A. Bueno.
 * \version 1.0.
 */
class CSourcePieceWise_AdjTurb : public CNumerics {
private:
	double **tau, *Velocity;

public:

	/*! 
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CSourcePieceWise_AdjTurb(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

	/*! 
	 * \brief Destructor of the class. 
	 */
	~CSourcePieceWise_AdjTurb(void);

	/*! 
	 * \brief Source term integration of the adjoint turbulence equation.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config);
};

/*! 
 * \class CSourcePieceWise_AdjElec
 * \brief Class for source term integration of the adjoint electric potential equation.
 * \ingroup SourceDiscr
 * \author F. Palacios.
 * \version 1.0.
 */
class CSourcePieceWise_AdjElec : public CNumerics {
public:

	/*! 
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CSourcePieceWise_AdjElec(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

	/*! 
	 * \brief Destructor of the class. 
	 */
	~CSourcePieceWise_AdjElec(void);

	/*! 
	 * \brief Source term integration of the adjoint electric potential equation.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetResidual(double *val_residual, CConfig *config);
};

/*! 
 * \class CSourcePieceWise_LinElec
 * \brief Class for source term integration of the linearized electric potential equation.
 * \ingroup SourceDiscr
 * \author F. Palacios.
 * \version 1.0.
 */
class CSourcePieceWise_LinElec : public CNumerics {
public:

	/*! 
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CSourcePieceWise_LinElec(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

	/*! 
	 * \brief Destructor of the class. 
	 */
	~CSourcePieceWise_LinElec(void);

	/*! 
	 * \brief Source term integration of the linearized electric potential equation.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetResidual(double *val_residual, CConfig *config);
};

/*! 
 * \class CSourcePieceWise_Plasma
 * \brief Class for integrating the source terms of the plasma equation.
 * \ingroup SourceDiscr
 * \author Amrita Mittal
 * \version 1.0.
 */
class CSourcePieceWise_Plasma : public CNumerics {
private:
	bool implicit;

	double r1, m1, n1, e1, r2, m2, n2, e2,r3, m3, n3, e3,T1,T2,T3, P1,P2,P3;

	double tolerance;
	double dT1_dr1,dT2_dr2,dT3_dr3,dT1_dm1,dT2_dm2,dT3_dm3,dT1_dn1,dT2_dn2,dT3_dn3,dT1_de1,dT2_de2,dT3_de3;
	double C12,C13,C23;
	double dC12_dT1,dC12_dT2,dC13_dT1,dC13_dT3,dC23_dT2,dC23_dT3;
	double nu12,nu13,nu23;
	double nu21,nu31,nu32;
	double dr1_dT1,dr2_dT2,dr3_dT3;
	double dv12_dT1,dv12_dT2,dv21_dT1,dv21_dT2;
	double dv13_dT1,dv13_dT3,dv31_dT1,dv31_dT3;
	double dv23_dT2,dv23_dT3,dv32_dT2,dv32_dT3;
	double Cv1, Cv2,Cv3;
	double gam1, gam2, gam3;
	double f13,f12,f23;
	double df13_dT1,df13_dT2,df13_dT3;
	double df23_dT1, df23_dT2,df23_dT3;
	double df12_dT1, df12_dT2,df12_dT3;
	double QT1, QT2, QT3;
	double dQT1_dT1,dQT1_dT2,dQT1_dT3;
	double dQT2_dT1,dQT2_dT2,dQT2_dT3;
	double dQT3_dT1,dQT3_dT2,dQT3_dT3;
	double dQT1_dr1, dQT1_dr2, dQT1_dr3;
	double dQT2_dr1, dQT2_dr2, dQT2_dr3;
	double dQT3_dr1, dQT3_dr2, dQT3_dr3;
	double C1, C2,C3;
	double Ck1,Ck2,Ck3;
	double eta1,eta2,eta3;
	double zeta1,zeta2,zeta3;
	double theta1,theta2,theta3;
	double phi1,phi2,phi3;
	double kf1,kf2,kf3,ke1,ke2,ke3,kb1,kb2,kb3;
	double R,dkf1_dT1,dkf2_dT2,dkf3_dT3,dke1_dT1,dke2_dT2,dke3_dT3;
	double dkb1_dT1,dkb2_dT2,dkb3_dT3;
	double dR_dr1,dR_dm1,dR_dn1,dR_de1;
	double dR_dr2,dR_dm2,dR_dn2,dR_de2;
	double dR_dr3,dR_dm3,dR_dn3,dR_de3;
	double S1,dS1_dr1,dS1_dm1,dS1_dn1,dS1_de1,dS1_dr2,dS1_dm2,dS1_dn2,dS1_de2,dS1_dr3,dS1_dm3,dS1_dn3,dS1_de3;
	double S2,dS2_dr1,dS2_dm1,dS2_dn1,dS2_de1,dS2_dr2,dS2_dm2,dS2_dn2,dS2_de2,dS2_dr3,dS2_dm3,dS2_dn3,dS2_de3;
	double S3,dS3_dr1,dS3_dm1,dS3_dn1,dS3_de1,dS3_dr2,dS3_dm2,dS3_dn2,dS3_de2,dS3_dr3,dS3_dm3,dS3_dn3,dS3_de3;
	double S4,dS4_dr1,dS4_dm1,dS4_dn1,dS4_de1,dS4_dr2,dS4_dm2,dS4_dn2,dS4_de2,dS4_dr3,dS4_dm3,dS4_dn3,dS4_de3;
	double S5,dS5_dr1,dS5_dm1,dS5_dn1,dS5_de1,dS5_dr2,dS5_dm2,dS5_dn2,dS5_de2,dS5_dr3,dS5_dm3,dS5_dn3,dS5_de3;
	double S6,dS6_dr1,dS6_dm1,dS6_dn1,dS6_de1,dS6_dr2,dS6_dm2,dS6_dn2,dS6_de2,dS6_dr3,dS6_dm3,dS6_dn3,dS6_de3;
	double S7,dS7_dr1,dS7_dm1,dS7_dn1,dS7_de1,dS7_dr2,dS7_dm2,dS7_dn2,dS7_de2,dS7_dr3,dS7_dm3,dS7_dn3,dS7_de3;
	double S8,dS8_dr1,dS8_dm1,dS8_dn1,dS8_de1,dS8_dr2,dS8_dm2,dS8_dn2,dS8_de2,dS8_dr3,dS8_dm3,dS8_dn3,dS8_de3;
	double S9,dS9_dr1,dS9_dm1,dS9_dn1,dS9_de1,dS9_dr2,dS9_dm2,dS9_dn2,dS9_de2,dS9_dr3,dS9_dm3,dS9_dn3,dS9_de3;
	double S10,dS10_dr1,dS10_dm1,dS10_dn1,dS10_de1,dS10_dr2,dS10_dm2,dS10_dn2,dS10_de2,dS10_dr3,dS10_dm3,dS10_dn3,dS10_de3;
	double S11,dS11_dr1,dS11_dm1,dS11_dn1,dS11_de1,dS11_dr2,dS11_dm2,dS11_dn2,dS11_de2,dS11_dr3,dS11_dm3,dS11_dn3,dS11_de3;
	double S12,dS12_dr1,dS12_dm1,dS12_dn1,dS12_de1,dS12_dr2,dS12_dm2,dS12_dn2,dS12_de2,dS12_dr3,dS12_dm3,dS12_dn3,dS12_de3;

	double Kb;
	double M1;
	double M3;
	double M2;
	double ec;
	double eps0;
	double Te;
	double Rc;
	double r12;
	double r13;
	double r23;
	double sigma12;
	double sigma13;
	double sigma23;

	double Ex, Ey;
	double AvgNum;

	double tol;
	double Tstart;
	double M1Avg, M2Avg, M3Avg, M1M2M3Avg3;
public:

	/*! 
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CSourcePieceWise_Plasma(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

	/*! 
	 * \brief Destructor of the class. 
	 */
	~CSourcePieceWise_Plasma(void);

	/*! 
	 * \brief Residual for source term integration.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetResidual(double *val_residual, double **val_Jacobian_i,CConfig *config);
};

/*! 
 * \class CSourcePieceWise_PlasmaDiatomic
 * \brief Class for integrating the source terms of the plasma equation.
 * \ingroup SourceDiscr
 * \author Amrita Mittal
 * \version 1.0.
 */
class CSourcePieceWise_PlasmaDiatomic : public CNumerics {
private:
	bool implicit;
	
	double r1, m1, n1, e1, r2, m2, n2, e2,r3, m3, n3, e3,T1,T2,T3, P1,P2,P3;
	
	double tolerance;
	double dT1_dr1,dT2_dr2,dT3_dr3,dT1_dm1,dT2_dm2,dT3_dm3,dT1_dn1,dT2_dn2,dT3_dn3,dT1_de1,dT2_de2,dT3_de3;
	double C12,C13,C23;
	double dC12_dT1,dC12_dT2,dC13_dT1,dC13_dT3,dC23_dT2,dC23_dT3;
	double nu12,nu13,nu23;
	double nu21,nu31,nu32;
	double dr1_dT1,dr2_dT2,dr3_dT3;
	double dv12_dT1,dv12_dT2,dv21_dT1,dv21_dT2;
	double dv13_dT1,dv13_dT3,dv31_dT1,dv31_dT3;
	double dv23_dT2,dv23_dT3,dv32_dT2,dv32_dT3;
	double Cv1, Cv2,Cv3;
	double gam1, gam2, gam3;
	double f13,f12,f23;
	double df13_dT1,df13_dT2,df13_dT3;
	double df23_dT1, df23_dT2,df23_dT3;
	double df12_dT1, df12_dT2,df12_dT3;
	double QT1, QT2, QT3;
	double dQT1_dT1,dQT1_dT2,dQT1_dT3;
	double dQT2_dT1,dQT2_dT2,dQT2_dT3;
	double dQT3_dT1,dQT3_dT2,dQT3_dT3;
	double dQT1_dr1, dQT1_dr2, dQT1_dr3;
	double dQT2_dr1, dQT2_dr2, dQT2_dr3;
	double dQT3_dr1, dQT3_dr2, dQT3_dr3;
	double C1, C2,C3;
	double Ck1,Ck2,Ck3;
	double eta1,eta2,eta3;
	double zeta1,zeta2,zeta3;
	double theta1,theta2,theta3;
	double phi1,phi2,phi3;
	double kf1,kf2,kf3,ke1,ke2,ke3,kb1,kb2,kb3;
	double R,dkf1_dT1,dkf2_dT2,dkf3_dT3,dke1_dT1,dke2_dT2,dke3_dT3;
	double dkb1_dT1,dkb2_dT2,dkb3_dT3;
	double dR_dr1,dR_dm1,dR_dn1,dR_de1;
	double dR_dr2,dR_dm2,dR_dn2,dR_de2;
	double dR_dr3,dR_dm3,dR_dn3,dR_de3;
	double S1,dS1_dr1,dS1_dm1,dS1_dn1,dS1_de1,dS1_dr2,dS1_dm2,dS1_dn2,dS1_de2,dS1_dr3,dS1_dm3,dS1_dn3,dS1_de3;
	double S2,dS2_dr1,dS2_dm1,dS2_dn1,dS2_de1,dS2_dr2,dS2_dm2,dS2_dn2,dS2_de2,dS2_dr3,dS2_dm3,dS2_dn3,dS2_de3;
	double S3,dS3_dr1,dS3_dm1,dS3_dn1,dS3_de1,dS3_dr2,dS3_dm2,dS3_dn2,dS3_de2,dS3_dr3,dS3_dm3,dS3_dn3,dS3_de3;
	double S4,dS4_dr1,dS4_dm1,dS4_dn1,dS4_de1,dS4_dr2,dS4_dm2,dS4_dn2,dS4_de2,dS4_dr3,dS4_dm3,dS4_dn3,dS4_de3;
	double S5,dS5_dr1,dS5_dm1,dS5_dn1,dS5_de1,dS5_dr2,dS5_dm2,dS5_dn2,dS5_de2,dS5_dr3,dS5_dm3,dS5_dn3,dS5_de3;
	double S6,dS6_dr1,dS6_dm1,dS6_dn1,dS6_de1,dS6_dr2,dS6_dm2,dS6_dn2,dS6_de2,dS6_dr3,dS6_dm3,dS6_dn3,dS6_de3;
	double S7,dS7_dr1,dS7_dm1,dS7_dn1,dS7_de1,dS7_dr2,dS7_dm2,dS7_dn2,dS7_de2,dS7_dr3,dS7_dm3,dS7_dn3,dS7_de3;
	double S8,dS8_dr1,dS8_dm1,dS8_dn1,dS8_de1,dS8_dr2,dS8_dm2,dS8_dn2,dS8_de2,dS8_dr3,dS8_dm3,dS8_dn3,dS8_de3;
	double S9,dS9_dr1,dS9_dm1,dS9_dn1,dS9_de1,dS9_dr2,dS9_dm2,dS9_dn2,dS9_de2,dS9_dr3,dS9_dm3,dS9_dn3,dS9_de3;
	double S10,dS10_dr1,dS10_dm1,dS10_dn1,dS10_de1,dS10_dr2,dS10_dm2,dS10_dn2,dS10_de2,dS10_dr3,dS10_dm3,dS10_dn3,dS10_de3;
	double S11,dS11_dr1,dS11_dm1,dS11_dn1,dS11_de1,dS11_dr2,dS11_dm2,dS11_dn2,dS11_de2,dS11_dr3,dS11_dm3,dS11_dn3,dS11_de3;
	double S12,dS12_dr1,dS12_dm1,dS12_dn1,dS12_de1,dS12_dr2,dS12_dm2,dS12_dn2,dS12_de2,dS12_dr3,dS12_dm3,dS12_dn3,dS12_de3;
	
	double Kb;
	double M1;
	double M3;
	double M2;
	double ec;
	double eps0;
	double Te;
	double Rc;
	double r12;
	double r13;
	double r23;
	double sigma12;
	double sigma13;
	double sigma23;
	
	double Ex, Ey;
	double AvgNum;
	
	double tol;
	double Tstart;
	double M1Avg, M2Avg, M3Avg, M1M2M3Avg3;
public:
	
	/*! 
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CSourcePieceWise_PlasmaDiatomic(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
	
	/*! 
	 * \brief Destructor of the class. 
	 */
	~CSourcePieceWise_PlasmaDiatomic(void);
	
	/*! 
	 * \brief Residual for source term integration.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetResidual(double *val_residual, double **val_Jacobian_i,CConfig *config);
};

/*! 
 * \class CSourceConservative_AdjFlow
 * \brief Class for source term integration in adjoint problem using a conservative scheme.
 * \ingroup SourceDiscr
 * \author F. Palacios.
 * \version 1.0.
 */
class CSourceConservative_AdjFlow : public CNumerics {
private:
	double *Velocity, *Residual_i, *Residual_j, *Mean_Residual;
	double **Mean_PrimVar_Grad;

public:

	/*! 
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CSourceConservative_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

	/*! 
	 * \brief Destructor of the class. 
	 */
	~CSourceConservative_AdjFlow(void);

	/*! 
	 * \brief Source term integration using a conservative scheme.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetResidual(double *val_residual, CConfig *config);
};

/*! 
 * \class CSourceConservative_AdjTurb
 * \brief Class for source term integration in adjoint turbulent problem using a conservative scheme.
 * \ingroup SourceDiscr
 * \author A. Bueno.
 * \version 1.0.
 */
class CSourceConservative_AdjTurb : public CNumerics {
public:

	/*! 
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CSourceConservative_AdjTurb(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

	/*! 
	 * \brief Destructor of the class. 
	 */
	~CSourceConservative_AdjTurb(void);

	/*! 
	 * \brief Source term integration using a conservative scheme.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config);
};

/*! 
 * \class CSourceRotationalFrame_Flow
 * \brief Class for source term for doing a rotational frame.
 * \ingroup SourceDiscr
 * \author F. Palacios.
 * \version 1.0.
 */
class CSourceRotationalFrame_Flow : public CNumerics {
public:

	/*! 
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CSourceRotationalFrame_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

	/*! 
	 * \brief Destructor of the class. 
	 */
	~CSourceRotationalFrame_Flow(void);

	/*! 
	 * \brief Residual of the rotational frame source term.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetResidual(double *val_residual, CConfig *config);
};

/*!
 * \class CSource_Template
 * \brief Dummy class.
 * \ingroup SourceDiscr
 * \author A. Lonkar.
 * \version 1.0.
 */
class CSource_Template : public CNumerics {
public:

	/*!
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config -  Name of the input config file
	 *
	 */
	CSource_Template(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);


	/*!
	 * \brief Residual for source term integration.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetResidual(double *val_residual, double **val_Jacobian_i,CConfig *config);

	/*!
	 * \brief Destructor of the class.
	 */
	~CSource_Template(void);
};

/*!
 * \class CConvectiveTemplate
 * \brief Class for setting up new method for spatial discretization of convective terms in flow Equations
 * \ingroup ConvDiscr
 * \author A. Lonkar
 * \version 1.0.
 */
class CConvective_Template : public CNumerics {
private:

	/* define private variables here */
	bool implicit;
	double *Diff_U;
	double *Velocity_i, *Velocity_j, *RoeVelocity;
	double *Proj_flux_tensor_i, *Proj_flux_tensor_j;
	double *delta_wave, *delta_vel;
	double *Lambda, *Epsilon;
	double **P_Tensor, **invP_Tensor;
	double sq_vel, Proj_ModJac_Tensor_ij, Density_i, Energy_i, SoundSpeed_i, Pressure_i, Enthalpy_i,
	Density_j, Energy_j, SoundSpeed_j, Pressure_j, Enthalpy_j, R, RoeDensity, RoeEnthalpy, RoeSoundSpeed,
	ProjVelocity, ProjVelocity_i, ProjVelocity_j, proj_delta_vel, delta_p, delta_rho;
	unsigned short iDim, iVar, jVar, kVar;

public:

	/*!
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CConvective_Template(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

	/*!
	 * \brief Destructor of the class.
	 */
	~CConvective_Template(void);

	/*!
	 * \brief Compute the Roe's flux between two nodes i and j.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config);
};

/*! 
 * \class CViscous_Template
 * \brief Class for computing viscous term using average of gradients.
 * \ingroup ViscDiscr
 * \author F. Palacios
 * \version 1.0.
 */
class CViscous_Template : public CNumerics {
private:
	
public:
	
	/*! 
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimension of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CViscous_Template(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
	
	/*! 
	 * \brief Destructor of the class. 
	 */
	~CViscous_Template(void);
	
	/*! 
	 * \brief Compute the viscous flow residual using an average of gradients.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config);
};

#include "numerics_structure.inl"
