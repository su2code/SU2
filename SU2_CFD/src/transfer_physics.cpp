/*!
 * \file transfer_structure.cpp
 * \brief Main subroutines for physics of the information transfer between zones
 * \author R. Sanchez
 * \version 4.0.1 "Cardinal"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *
 * Copyright (C) 2012-2015 SU2, the open-source CFD code.
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

#include "../include/transfer_structure.hpp"

CTransfer_FlowTraction::CTransfer_FlowTraction(void) : CTransfer() {

}

CTransfer_FlowTraction::CTransfer_FlowTraction(unsigned short val_nVar, unsigned short val_nConst, CConfig *config) : CTransfer(val_nVar, val_nConst, config) {

}

CTransfer_FlowTraction::~CTransfer_FlowTraction(void) {

}

void CTransfer_FlowTraction::GetPhysical_Constants(CSolver *flow_solution, CSolver *struct_solution,
		   	   	   	   	   	   	   	   	   	   	   CGeometry *flow_geometry, CGeometry *struct_geometry,
												   CConfig *flow_config, CConfig *struct_config){

	unsigned short iVar;

	/*--- We have to clear the traction before applying it, because we are "adding" to node and not "setting" ---*/

	for (unsigned long iPoint = 0; iPoint < struct_geometry->GetnPoint(); iPoint++){
		struct_solution->node[iPoint]->Clear_FlowTraction();
	}

  	/*--- Redimensionalize the pressure ---*/

	su2double *Velocity_ND, *Velocity_Real;
	su2double Density_ND,  Density_Real, Velocity2_Real, Velocity2_ND;

    Velocity_Real = flow_config->GetVelocity_FreeStream();
    Density_Real = flow_config->GetDensity_FreeStream();

    Velocity_ND = flow_config->GetVelocity_FreeStreamND();
    Density_ND = flow_config->GetDensity_FreeStreamND();

	Velocity2_Real = 0.0;
	Velocity2_ND = 0.0;
    for (iVar = 0; iVar < nVar; iVar++){
    	Velocity2_Real += Velocity_Real[iVar]*Velocity_Real[iVar];
    	Velocity2_ND += Velocity_ND[iVar]*Velocity_ND[iVar];
    }

    Physical_Constants[0] = Density_Real*Velocity2_Real/(Density_ND*Velocity2_ND);

	/*--- Apply a ramp to the transfer of the fluid loads ---*/

	su2double ModAmpl;
	su2double CurrentTime = struct_config->GetCurrent_DynTime();
	su2double Static_Time = struct_config->GetStatic_Time();

	bool Ramp_Load = struct_config->GetRamp_Load();
	su2double Ramp_Time = struct_config->GetRamp_Time();

    bool Sigmoid_Load = struct_config->GetSigmoid_Load();
	su2double Sigmoid_Time = struct_config->GetSigmoid_Time();
	su2double Sigmoid_K = struct_config->GetSigmoid_K();
	su2double SigAux = 0.0;

	if (CurrentTime <= Static_Time){ ModAmpl=0.0; }
	else if((CurrentTime > Static_Time) &&
			(CurrentTime <= (Static_Time + Ramp_Time)) &&
			(Ramp_Load)){
		ModAmpl = (CurrentTime-Static_Time) / Ramp_Time;
		ModAmpl = max(ModAmpl,0.0);
		ModAmpl = min(ModAmpl,1.0);
		Physical_Constants[1] = ModAmpl;
	}
	else if((CurrentTime > Static_Time) &&
			(CurrentTime <= (Static_Time + Sigmoid_Time)) &&
			(Sigmoid_Load)){
		SigAux = (CurrentTime-Static_Time) / Sigmoid_Time;
		ModAmpl = (1 / (1+exp(-1*Sigmoid_K*(SigAux - 0.5)) ) );
		ModAmpl = max(ModAmpl,0.0);
		ModAmpl = min(ModAmpl,1.0);
		Physical_Constants[1] = ModAmpl;
	}
	else{ Physical_Constants[1] = 1.0; }

}

void CTransfer_FlowTraction::GetDonor_Variable(CSolver *flow_solution, CGeometry *flow_geometry, CConfig *flow_config,
					   	   	   	   	   	   	   unsigned long Marker_Flow, unsigned long Vertex_Flow, unsigned long Point_Struct){


	unsigned short iVar, jVar;
	unsigned long Point_Flow;
	su2double *Normal_Flow;

	// Check the kind of fluid problem
	bool compressible       = (flow_config->GetKind_Regime() == COMPRESSIBLE);
	bool incompressible     = (flow_config->GetKind_Regime() == INCOMPRESSIBLE);
	bool viscous_flow       = ((flow_config->GetKind_Solver() == NAVIER_STOKES) ||
							   (flow_config->GetKind_Solver() == RANS) );

	// Parameters for the calculations
	// Pn: Pressure
	// Pinf: Pressure_infinite
	// div_vel: Velocity divergence
	// Dij: Dirac delta
	su2double Pn = 0.0, div_vel = 0.0, Dij = 0.0;
	su2double Viscosity = 0.0;
	su2double **Grad_PrimVar;
	su2double Tau[3][3] = { {0.0, 0.0, 0.0} ,
							{0.0, 0.0, 0.0} ,
							{0.0, 0.0, 0.0} } ;

	su2double Pinf = flow_solution->GetPressure_Inf();

	Point_Flow = flow_geometry->vertex[Marker_Flow][Vertex_Flow]->GetNode();
		// Get the normal at the vertex: this normal goes inside the fluid domain.
	Normal_Flow = flow_geometry->vertex[Marker_Flow][Vertex_Flow]->GetNormal();

	// Retrieve the values of pressure and viscosity
	if (incompressible){

		Pn = flow_solution->node[Point_Flow]->GetPressureInc();

		if (viscous_flow){

			Grad_PrimVar = flow_solution->node[Point_Flow]->GetGradient_Primitive();
			Viscosity = flow_solution->node[Point_Flow]->GetLaminarViscosityInc();

		}
	}
	else if (compressible){

		Pn = flow_solution->node[Point_Flow]->GetPressure();

		if (viscous_flow){

			Grad_PrimVar = flow_solution->node[Point_Flow]->GetGradient_Primitive();
			Viscosity = flow_solution->node[Point_Flow]->GetLaminarViscosity();

		}
	}

	// Calculate tn in the fluid nodes for the inviscid term --> Units of force (non-dimensional).
	for (iVar = 0; iVar < nVar; iVar++) {
		Donor_Variable[iVar] = -(Pn-Pinf)*Normal_Flow[iVar];
	}

	// Calculate tn in the fluid nodes for the viscous term

	if ((incompressible || compressible) && viscous_flow){

		// Divergence of the velocity
		div_vel = 0.0; for (iVar = 0; iVar < nVar; iVar++) div_vel += Grad_PrimVar[iVar+1][iVar];
		if (incompressible) div_vel = 0.0;

		for (iVar = 0; iVar < nVar; iVar++) {

			for (jVar = 0 ; jVar < nVar; jVar++) {
				// Dirac delta
				Dij = 0.0; if (iVar == jVar) Dij = 1.0;

				// Viscous stress
				Tau[iVar][jVar] = Viscosity*(Grad_PrimVar[jVar+1][iVar] + Grad_PrimVar[iVar+1][jVar]) -
						TWO3*Viscosity*div_vel*Dij;

				// Viscous component in the tn vector --> Units of force (non-dimensional).
				Donor_Variable[iVar] += Tau[iVar][jVar]*Normal_Flow[jVar];
			}
		}
	}

	// Redimensionalize and take into account ramp transfer of the loads
	for (iVar = 0; iVar < nVar; iVar++){
		Donor_Variable[iVar] = Donor_Variable[iVar] * Physical_Constants[0] * Physical_Constants[1];
	}

}

void CTransfer_FlowTraction::SetTarget_Variable(CSolver *fea_solution, CGeometry *fea_geometry,
												CConfig *fea_config, unsigned long Marker_Struct,
												unsigned long Vertex_Struct, unsigned long Point_Struct){

	/*--- Add to the Flow traction ---*/
	fea_solution->node[Point_Struct]->Add_FlowTraction(Target_Variable);

}



CTransfer_StructuralDisplacements::CTransfer_StructuralDisplacements(void) : CTransfer() {

}

CTransfer_StructuralDisplacements::CTransfer_StructuralDisplacements(unsigned short val_nVar, unsigned short val_nConst, CConfig *config) : CTransfer(val_nVar, val_nConst, config) {

}

CTransfer_StructuralDisplacements::~CTransfer_StructuralDisplacements(void) {

}


void CTransfer_StructuralDisplacements::GetPhysical_Constants(CSolver *struct_solution, CSolver *flow_solution,
		   	   	   	   	   	   	   	   	   	   	   CGeometry *struct_geometry, CGeometry *flow_geometry,
												   CConfig *struct_config, CConfig *flow_config){

}

void CTransfer_StructuralDisplacements::GetDonor_Variable(CSolver *struct_solution, CGeometry *struct_geometry, CConfig *struct_config,
					   	   	   	   	   	   	   	          unsigned long Marker_Struct, unsigned long Vertex_Struct, unsigned long Point_Struct){


	su2double *DisplacementDonor, *DisplacementDonor_Prev;
	unsigned short iVar;

    /*--- The displacements come from the predicted solution ---*/
    DisplacementDonor = struct_solution->node[Point_Struct]->GetSolution_Pred();

    DisplacementDonor_Prev = struct_solution->node[Point_Struct]->GetSolution_Pred_Old();

	for (iVar = 0; iVar < nVar; iVar++){
		Donor_Variable[iVar] = DisplacementDonor[iVar] - DisplacementDonor_Prev[iVar];
	}


}

void CTransfer_StructuralDisplacements::SetTarget_Variable(CSolver *flow_solution, CGeometry *flow_geometry,
														   CConfig *flow_config, unsigned long Marker_Flow,
														   unsigned long Vertex_Flow, unsigned long Point_Flow){

	su2double VarCoord[3] = {0.0, 0.0, 0.0};
	unsigned short iVar;

	for (iVar = 0; iVar < nVar; iVar++)
		VarCoord[iVar] = Target_Variable[iVar];

	flow_geometry->vertex[Marker_Flow][Vertex_Flow]->SetVarCoord(VarCoord);

}

CTransfer_StructuralDisplacements_Original::CTransfer_StructuralDisplacements_Original(void) : CTransfer() {

}

CTransfer_StructuralDisplacements_Original::CTransfer_StructuralDisplacements_Original(unsigned short val_nVar, unsigned short val_nConst, CConfig *config) : CTransfer(val_nVar, val_nConst, config) {

}

CTransfer_StructuralDisplacements_Original::~CTransfer_StructuralDisplacements_Original(void) {

}


void CTransfer_StructuralDisplacements_Original::GetPhysical_Constants(CSolver *struct_solution, CSolver *flow_solution,
		   	   	   	   	   	   	   	   	   	   	   CGeometry *struct_geometry, CGeometry *flow_geometry,
												   CConfig *struct_config, CConfig *flow_config){

}

void CTransfer_StructuralDisplacements_Original::GetDonor_Variable(CSolver *struct_solution, CGeometry *struct_geometry, CConfig *struct_config,
					   	   	   	   	   	   	   	          unsigned long Marker_Struct, unsigned long Vertex_Struct, unsigned long Point_Struct){


	su2double *Coord_Struct, *Displacement_Struct;
	unsigned short iVar;

    Coord_Struct = struct_geometry->node[Point_Struct]->GetCoord();

    /*--- The displacements come from the predicted solution ---*/
    Displacement_Struct = struct_solution->node[Point_Struct]->GetSolution_Pred();

	for (iVar = 0; iVar < nVar; iVar++){
		Donor_Variable[iVar] = Coord_Struct[iVar] + Displacement_Struct[iVar];
	}

}

void CTransfer_StructuralDisplacements_Original::SetTarget_Variable(CSolver *flow_solution, CGeometry *flow_geometry,
														   CConfig *flow_config, unsigned long Marker_Flow,
														   unsigned long Vertex_Flow, unsigned long Point_Flow){

	su2double *Coord, VarCoord[3] = {0.0, 0.0, 0.0};
	unsigned short iVar;

	Coord = flow_geometry->node[Point_Flow]->GetCoord();

	for (iVar = 0; iVar < nVar; iVar++)
		VarCoord[iVar] = Target_Variable[iVar]-Coord[iVar];

	flow_geometry->vertex[Marker_Flow][Vertex_Flow]->SetVarCoord(VarCoord);

}




CTransfer_ConservativeVars::CTransfer_ConservativeVars(void) : CTransfer() {

}

CTransfer_ConservativeVars::CTransfer_ConservativeVars(unsigned short val_nVar, unsigned short val_nConst, CConfig *config) : CTransfer(val_nVar, val_nConst, config) {

}

CTransfer_ConservativeVars::~CTransfer_ConservativeVars(void) {

}


void CTransfer_ConservativeVars::GetPhysical_Constants(CSolver *donor_solution, CSolver *target_solution,
		   	   	   	   	   	   	   	   	   	   	       CGeometry *donor_geometry, CGeometry *target_geometry,
													   CConfig *donor_config, CConfig *target_config){

}

void CTransfer_ConservativeVars::GetDonor_Variable(CSolver *donor_solution, CGeometry *donor_geometry, CConfig *donor_config,
					   	   	   	   	   	   	   	   unsigned long Marker_Donor, unsigned long Vertex_Donor, unsigned long Point_Donor){

	su2double *Solution;
	unsigned short iVar;

    /*--- Retrieve solution and set it as the donor variable ---*/
	Solution = donor_solution->node[Point_Donor]->GetSolution();

	for (iVar = 0; iVar < nVar; iVar++){
		Donor_Variable[iVar] = Solution[iVar];
	}


}

void CTransfer_ConservativeVars::SetTarget_Variable(CSolver *target_solution, CGeometry *target_geometry,
													CConfig *target_config, unsigned long Marker_Target,
													unsigned long Vertex_Target, unsigned long Point_Target){

	/*--- Set the target solution with the value of the Target Variable ---*/
	target_solution->node[Point_Target]->SetSolution(Target_Variable);

}




CTransfer_MixingPlaneInterface::CTransfer_MixingPlaneInterface(void) : CTransfer() {

}

CTransfer_MixingPlaneInterface::CTransfer_MixingPlaneInterface(unsigned short val_nVar, unsigned short val_nConst, CConfig *config){
	unsigned short iVar;
	nVar = val_nVar;

	Donor_Variable     = new su2double[nVar + 2];
	Target_Variable    = new su2double[nVar + 2];



	for (iVar = 0; iVar < nVar + 2; iVar++){
		Donor_Variable[iVar]  = 0.0;
		Target_Variable[iVar] = 0.0;
	}


}

CTransfer_MixingPlaneInterface::~CTransfer_MixingPlaneInterface(void) {
	if (Donor_Variable       != NULL) delete [] Donor_Variable;
	if (Target_Variable      != NULL) delete [] Target_Variable;
}



void CTransfer_MixingPlaneInterface::GetDonor_Variable(CSolver *donor_solution, CGeometry *donor_geometry,
								   	     	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	CConfig *donor_config, unsigned long Marker_Donor,
																														unsigned long iSpan, unsigned long rank) {

	unsigned short nDim = nVar - 2;

	Donor_Variable[0] = donor_solution->GetAverageDensity(Marker_Donor, iSpan);
	Donor_Variable[1]	= donor_solution->GetAveragePressure(Marker_Donor, iSpan);
	Donor_Variable[2] = donor_solution->GetAverageTotPressure(Marker_Donor, iSpan);
 	Donor_Variable[3] = donor_solution->GetAverageTotTemperature(Marker_Donor, iSpan);
 	Donor_Variable[4] = donor_solution->GetAverageTurboVelocity(Marker_Donor, iSpan)[0];
 	Donor_Variable[5] = donor_solution->GetAverageTurboVelocity(Marker_Donor, iSpan)[1];

 	if(nDim == 3){
 		Donor_Variable[6] = donor_solution->GetAverageTurboVelocity(Marker_Donor, iSpan)[2];
 	}

}


void CTransfer_MixingPlaneInterface::SetTarget_Variable(CSolver *target_solution, CGeometry *target_geometry,
										  CConfig *target_config, unsigned long Marker_Target,
										  unsigned long iSpan, unsigned long rank) {

	unsigned short nDim = nVar - 2;

	target_solution->SetExtAverageDensity(Marker_Target, iSpan, Target_Variable[0]);
	target_solution->SetExtAveragePressure(Marker_Target, iSpan, Target_Variable[1]);
	target_solution->SetExtAverageTotPressure(Marker_Target, iSpan, Target_Variable[2]);
  target_solution->SetExtAverageTotTemperature(Marker_Target, iSpan, Target_Variable[3]);
	target_solution->SetExtAverageTurboVelocity(Marker_Target, iSpan, 0, Target_Variable[4]);
	target_solution->SetExtAverageTurboVelocity(Marker_Target, iSpan, 1, Target_Variable[5]);


  if(nDim == 3){
  	target_solution->SetExtAverageTurboVelocity(Marker_Target, iSpan, 2, Target_Variable[6]);
  }


}

void CTransfer_MixingPlaneInterface::GetSetTurboPerformance(CSolver *donor_solution, CSolver *target_solution, unsigned short donorZone){

	/*--- Loop over the nMarker of turboperformance and get the desired values ---*/
//	TotalStaticEfficiency[iMarker_Monitoring] = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotalStaticEfficiency(iMarker_Monitoring);
//	TotalTotalEfficiency[iMarker_Monitoring]  = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotalTotalEfficiency(iMarker_Monitoring);
	target_solution->SetKineticEnergyLoss(donor_solution->GetKineticEnergyLoss(donorZone), donorZone);
	target_solution->SetTotalPressureLoss(donor_solution->GetTotalPressureLoss(donorZone), donorZone);
	target_solution->SetMassFlowIn(donor_solution->GetMassFlowIn(donorZone), donorZone);
	target_solution->SetMassFlowOut(donor_solution->GetMassFlowOut(donorZone), donorZone);
	target_solution->SetFlowAngleIn(donor_solution->GetFlowAngleIn(donorZone), donorZone);
	target_solution->SetFlowAngleOut(donor_solution->GetFlowAngleOut(donorZone), donorZone);
	target_solution->SetEulerianWork(donor_solution->GetEulerianWork(donorZone), donorZone);
	target_solution->SetTotalEnthalpyIn(donor_solution->GetTotalEnthalpyIn(donorZone), donorZone);
	target_solution->SetPressureRatio(donor_solution->GetPressureRatio(donorZone), donorZone);
	target_solution->SetPressureOut(donor_solution->GetPressureOut(donorZone), donorZone);
	target_solution->SetEnthalpyOut(donor_solution->GetEnthalpyOut(donorZone), donorZone);
	target_solution->SetMachIn(donor_solution->GetMachIn(donorZone), donorZone);
	target_solution->SetMachOut(donor_solution->GetMachOut(donorZone), donorZone);
	target_solution->SetNormalMachIn(donor_solution->GetNormalMachIn(donorZone), donorZone);
	target_solution->SetNormalMachOut(donor_solution->GetNormalMachOut(donorZone), donorZone);
	target_solution->SetVelocityOutIs(donor_solution->GetVelocityOutIs(donorZone), donorZone);
	target_solution->SetTotalPresureIn(donor_solution->GetTotalPresureIn(donorZone), donorZone);
	target_solution->SetTotalTemperatureIn(donor_solution->GetTotalTemperatureIn(donorZone), donorZone);
	target_solution->SetFlowAngleIn_BC(donor_solution->GetFlowAngleIn_BC(donorZone), donorZone);
	target_solution->SetEntropyIn(donor_solution->GetEntropyIn(donorZone), donorZone);
	target_solution->SetEntropyIn_BC(donor_solution->GetEntropyIn_BC(donorZone), donorZone);
	target_solution->SetTotalEnthalpyIn_BC(donor_solution->GetTotalEnthalpyIn_BC(donorZone), donorZone);

}

