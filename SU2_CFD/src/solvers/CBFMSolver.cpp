/*!
 * \file CBFMSolver.cpp
 * \brief Subroutines to be implemented for any new solvers
 * \author E.C. Bunschoten
 * \version 7.1.0 "Blackbird"
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


#include "../../include/solvers/CBFMSolver.hpp"
#include "../../include/variables/CBFMVariable.hpp"
#include "../../include/numerics/ReadBFMInput.hpp"
#include "../../include/numerics/BFMInterpolator.hpp"
#include "../../include/gradients/computeGradientsGreenGauss.hpp"
#include "../../include/gradients/computeGradientsLeastSquares.hpp"

CBFMSolver::CBFMSolver(void) : CSolver() { }

CBFMSolver::CBFMSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CSolver() {
    /* This function initalizes the BFM solver */

    
    rank = SU2_MPI::GetRank();
    if(rank == MASTER_NODE)
        cout << "Initiating BFM solver" << endl;
    nPoint = geometry->GetnPoint();
    nDim = geometry->GetnDim();
    nVar = N_BFM_PARAMS;
    Omega = config->GetBFM_Rotation() * PI_NUMBER / 30;
    // Defining BFM variable class
    nodes = new CBFMVariable(nPoint, nDim, N_BFM_PARAMS);
    SetBaseClassPointerToNodes();
    
    // Reading BFM geometry input file
    BFM_File_Reader = new ReadBFMInput(config, config->GetBFM_FileName());

    
    BFM_formulation = config->GetBFM_Formulation();
    if(rank == MASTER_NODE){
        switch (BFM_formulation)
        {
        case HALL:
            cout << "Body-Force Model selection: Hall" << endl;
            break;
        case THOLLET:
            cout << "Body-Force Model selection: Thollet" << endl;
            break;
        default:
            SU2_MPI::Error(string("No suitable Body-Force Model was selected "),
                    CURRENT_FUNCTION);
            break;
        }
    }
    // Filling in the blade geometry parameter names
    BFM_Parameter_Names.resize(N_BFM_PARAMS);
    BFM_Parameter_Names.at(I_AXIAL_CHORD) = "axial_coordinate";
    BFM_Parameter_Names.at(I_RADIAL_COORDINATE) = "radial_coordinate";
    BFM_Parameter_Names.at(I_BLOCKAGE_FACTOR) = "blockage_factor";
    BFM_Parameter_Names.at(I_CAMBER_NORMAL_AXIAL) = "n_ax";
    BFM_Parameter_Names.at(I_CAMBER_NORMAL_TANGENTIAL) = "n_tang";
    BFM_Parameter_Names.at(I_CAMBER_NORMAL_RADIAL) = "n_rad";
    BFM_Parameter_Names.at(I_LEADING_EDGE_AXIAL) = "ax_LE";

    // Commencing blade geometry interpolation
    if(rank == MASTER_NODE)
        cout << "Interpolating blade geometry parameters to nodes" << endl;
    Interpolator = new BFMInterpolator(BFM_File_Reader, this, geometry, config);
    Interpolator->Interpolate(BFM_File_Reader, this, geometry);
    SU2_OMP_BARRIER
    // Setting the solution as the metal blockage factor so the periodic, spatial gradient can be computed
    for(unsigned long iPoint=0; iPoint<nPoint; ++iPoint){
        nodes->SetSolution(iPoint, 0, nodes->GetAuxVar(iPoint, I_BLOCKAGE_FACTOR));
    }

    // Commencing metal blockage factor gradient computation
    const auto &solution = nodes->GetSolution();
    auto &gradient = nodes->GetGradient();
    computeGradientsGreenGauss(this, SOLUTION_GRADIENT, PERIODIC_SOL_GG, *geometry,
                             *config, solution, 0, 1, gradient);
    
    // Storing metal blockage gradient in auxilary variable gradient
    for(unsigned long iPoint=0; iPoint < nPoint; ++iPoint){
        for(unsigned short iDim=0; iDim<nDim; ++iDim){
            nodes->SetAuxVarGradient(iPoint, I_BLOCKAGE_FACTOR, iDim, nodes->GetGradient(iPoint, 0, iDim));
        }
    }

    // Interpolator class is no longer needed, so it's deleted to free up memory
    delete Interpolator;

    // Computing cylindrical projections of the node coordinates.
    ComputeCylProjections(geometry, config);
    
}

CBFMSolver::~CBFMSolver(void) { 
    delete nodes;
    delete BFM_File_Reader;
}



void CBFMSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) { 

}

void CBFMSolver::ComputeBFMSources(CSolver **solver_container, unsigned long iPoint, vector<su2double> & BFM_sources){
    // In this function, the flow source terms according to the body-force model are computed

    su2double bffac; // Body-force factor
    su2double W_ax, W_th, W_r; // Relative, cylindrical velocity components (axial, tangential, radial)
    vector<su2double*> W_cyl = {&W_ax, &W_th, &W_r}; // Vector containing the relative velocity components.
    
    // Getting node body-force factor. If 1, the BFM source terms are computed at the current node.
    bffac = nodes->GetAuxVar(iPoint, I_BODY_FORCE_FACTOR);

    // Computing the relative velocity at the current node.
    ComputeRelativeVelocity(solver_container, iPoint, W_cyl);
    
    if(bffac == 1){
        // Selection of the BFM stated in the config file.
        switch (BFM_formulation)
        {
        case HALL:
            ComputeBFMSources_Hall(solver_container, iPoint, BFM_sources, W_cyl);
            break;
        case THOLLET:
            ComputeBFMSources_Thollet(solver_container, iPoint, BFM_sources, W_cyl);
            break;
        default:
            break;
        }
    }else{
        // If the BFM is inactive at the current node, the BFM source terms are set to zero.
        for(unsigned short iDim=0; iDim<nDim+2; ++iDim){
            BFM_sources.at(iDim) = 0;
        }
    }
}

void CBFMSolver::ComputeBFMSources_Hall(CSolver **solver_container, unsigned long iPoint, vector<su2double> & BFM_sources, vector<su2double*>&W_cyl){
    /*---Halls Body-Force Model */
    /*
        This function computes the body-force source terms according to Halls BFM formulation. 
        Halls BFM generates source terms for flow deflection only and does not take blade thickness
        effects into account.
    */

    su2double F[nDim]; // body-force vector(Cartesian)
    su2double Rotation_rate = Omega; // Blade rotation rate [rad s^-1]
    su2double Nx, Nt, Nr, rotFac, blade_count; // axial, tangential, radial camber normal vector components[-], rotation factor[-], blade count[-]
    su2double W_ax, W_th, W_r, // axial, tangential, radial relative velocity [m s^-1]
     WdotN, // Dot product between camber normal vector and relative velocity vector [m s^-1]
    W_nx, W_nth, W_nr, // Normal projections of relative velocity components [m s^-1]
    W_px, W_pth, W_pr; // Parallel projections of relative velocity components [m s^-1]
    su2double W_p, W_mag; // Parallel projection relative velocity magnitude [m s^-1]
    su2double radius, // Radius [m]
    pitch, // Blade pitch [m]
    density; // Fluid density [kg m^-3]
    su2double delta; // Flow incidence angle [rad]
    su2double F_n, F_p, // Normal, parallel body-force [m s^-2]
     F_ax, F_th, F_r, // Cylindrical body-force components
     e_source; // Energy equation source term [kg m^-1 s^-3]
    
    // Getting interpolated blade geometry parameters from auxilary variable
    Nx = nodes->GetAuxVar(iPoint, I_CAMBER_NORMAL_AXIAL); 		// Camber normal component in axial direction
    Nt = nodes->GetAuxVar(iPoint, I_CAMBER_NORMAL_TANGENTIAL); 			// Camber normal component in tangential direction
    Nr = nodes->GetAuxVar(iPoint, I_CAMBER_NORMAL_RADIAL);			// Camber normal component in radial direction
    rotFac = nodes->GetAuxVar(iPoint, I_ROTATION_FACTOR);	// Rotation factor. Multiplies the body-force rotation value.
    blade_count = nodes->GetAuxVar(iPoint, I_BLADE_COUNT);

    // // Getting cylindrical, relative velocity vector
    W_ax = *W_cyl.at(0);
    W_th = *W_cyl.at(1);
    W_r = *W_cyl.at(2);
    
    // Getting local radius and computing blade pitch
    radius = nodes->GetAuxVar(iPoint, I_RADIAL_COORDINATE);
    pitch = 2 * PI_NUMBER * radius / blade_count;	// Computing pitch

    WdotN = W_ax * Nx + W_th * Nt + W_r * Nr;		// Dot product of relative velocity and camber normal vector
    W_nx = WdotN * Nx, W_nr = WdotN * Nr, W_nth = WdotN * Nt;		// Relative velocity components normal to the blade
    W_px = W_ax - W_nx, W_pr = W_r - W_nr, W_pth = W_th - W_nth;  // Relative velocity components parallel to the blade 
    
    W_p = sqrt(W_px * W_px + W_pr * W_pr + W_pth * W_pth);	// Parallel relative velocity magnitude
    // Relative velocity magnitude
    W_mag = sqrt(W_ax * W_ax + W_th * W_th + W_r * W_r);
    // Calculating the deviation angle
    delta = asin(WdotN / (W_mag + 1e-6));
    
    // Computing normal body force magnitude
    F_n = -PI_NUMBER * delta * (1 / pitch) * (1 / abs(Nt)) * W_mag * W_mag;

    // Transforming the normal and force component to cyllindrical coordinates
    density = solver_container[FLOW_SOL]->GetNodes()->GetDensity(iPoint);
    F_ax = F_n * (cos(delta) * Nx - sin(delta) * (W_px / (W_p + 1e-6)));		// Axial body-force component
    F_r = F_n * (cos(delta) * Nr - sin(delta) * (W_pr / (W_p + 1e-6)));		// Radial body-force component
    F_th = F_n * (cos(delta) * Nt - sin(delta) * (W_pth / (W_p + 1e-6)));	// Tangential body-force component
    e_source = rotFac * Rotation_rate * radius * F_th;				// Energy source term
    
    // Appending Cartesial body-forces to body-force vector
    for(int iDim=0; iDim<nDim; iDim++){
        F[iDim] = density * (F_ax*(nodes->GetAxialProjection(iPoint, iDim))
                          + F_th*(nodes->GetTangentialProjection(iPoint, iDim))
                          + F_r*(nodes->GetRadialProjection(iPoint, iDim)));
    }

    
    // Appending source term values to the body-force source term vector
    BFM_sources.at(0) = 0.0;	// empty density component

    // Appending momentum source terms
    for(int iDim = 0; iDim < nDim; iDim ++) {
        BFM_sources.at(iDim+1) = F[iDim];
    }
    // Appending energy source term
    BFM_sources.at(nDim+1) = density * e_source;

}

void CBFMSolver::ComputeBFMSources_Thollet(CSolver **solver_container, unsigned long iPoint, vector<su2double>&BFM_sources, vector<su2double*>&W_cyl){
    /*---Thollets Body-Force Model */
    /*
        This function computes the body-force source terms according to Tholletss BFM formulation. 
        Thollets BFM generates body-forces normal and parallel to the local, relative flow, enabling loss modeling.
        Additionally, compressibility effects are accounted for while computing the normal body force component and
        flow obstruction due to blade thickness is modeled through an additional set of source terms, of which the
        magnitude depends on the metal blockage factor distribution.
    */

    su2double F[nDim];  // body-force vector(Cartesian)
    su2double Rotation_rate = Omega;  // Blade rotation rate [rad s^-1]
    su2double b, // Metal blockage factor [-]
     Nx, Nt, Nr, // axial, tangential, radial camber normal vector components[-]
     rotFac, // Rotation factor [-]
     blade_count; // Blade row blade count [-]
    su2double W_ax, W_th, W_r, // axial, tangential, radial relative velocity [m s^-1]
    WdotN, // Dot product between camber normal vector and relative velocity vector [m s^-1]
    W_nx, W_nth, W_nr, // Normal projections of relative velocity components [m s^-1]
    W_px, W_pth, W_pr; // Parallel projections of relative velocity components [m s^-1]
    su2double W_p, // Parallel projected rlative velocity magnitude [m s^-1]
    W_mag, // Relative velocity magnitude [m s^-1]
    K_mach; // Normal force compressibility correction factor [-]
    su2double radius, pitch; // radius [m], blade pitch[m]
    su2double delta; // Flow incidense angle [rad]
    su2double F_n, F_p, // Normal, parallel body-force [m s^-2]
     F_ax, F_th, F_r, // Cylindrical body-force components
     e_source; // Energy equation source term [kg m^-1 s^-3]
    su2double ax, ax_le, // Axial node coordiate [m], blade leading edge coordinate[m]
     mu, // Dynamic viscosity [kg m^-1 s^-2]
     density, // Fluid density [kg m^-3]
     Re_ax, C_f; // Axial Reynolds number[-], blade friction factor [-]

    // Getting the blade geometric parameters from the auxilary variable
    b = nodes->GetAuxVar(iPoint, I_BLOCKAGE_FACTOR);			// Metal blockage factor
    Nx = nodes->GetAuxVar(iPoint, I_CAMBER_NORMAL_AXIAL); 		// Camber normal component in axial direction
    Nt = nodes->GetAuxVar(iPoint, I_CAMBER_NORMAL_TANGENTIAL); 			// Camber normal component in tangential direction
    Nr = nodes->GetAuxVar(iPoint, I_CAMBER_NORMAL_RADIAL);			// Camber normal component in radial direction
    rotFac = nodes->GetAuxVar(iPoint, I_ROTATION_FACTOR);	// Rotation factor. Multiplies the body-force rotation value.
    blade_count = nodes->GetAuxVar(iPoint, I_BLADE_COUNT);

    // Getting the relative velocity vector
    W_ax = *W_cyl.at(0);
    W_th = *W_cyl.at(1);
    W_r = *W_cyl.at(2);
    
    // Getting the local radius and computing blade pitch
    radius = nodes->GetAuxVar(iPoint, I_RADIAL_COORDINATE);
    pitch = 2 * PI_NUMBER * radius / blade_count;	// Computing pitch

    WdotN = W_ax * Nx + W_th * Nt + W_r * Nr;		// Dot product of relative velocity and camber normal vector
    W_nx = WdotN * Nx, W_nr = WdotN * Nr, W_nth = WdotN * Nt;		// Relative velocity components normal to the blade
    W_px = W_ax - W_nx, W_pr = W_r - W_nr, W_pth = W_th - W_nth;  // Relative velocity components parallel to the blade 
    
    W_p = sqrt(W_px * W_px + W_pr * W_pr + W_pth * W_pth);	// Parallel relative velocity magnitude
    // Relative velocity magnitude
    W_mag = sqrt(W_ax * W_ax + W_th * W_th + W_r * W_r);
    // Calculating the deviation angle
    delta = asin(WdotN / W_mag);
    if(isinf(delta)) cout << "infinite angle of attack" << endl;
    if(isnan(delta)) cout << "NaN angle of attack" << endl;
    // Computing the compressiblity correction factor for the normal body force
    K_mach = ComputeKMach(solver_container, iPoint, W_cyl);

    // Computing normal body force magnitude
    F_n = -PI_NUMBER * K_mach * delta * (1 / pitch) * (1 / abs(Nt)) * (1 / b) * W_mag * W_mag;

    // Getting leading edge axial coordinate and node axial coordinate
    ax_le = nodes->GetAuxVar(iPoint, I_LEADING_EDGE_AXIAL);
    ax = nodes->GetAuxVar(iPoint, I_AXIAL_COORDINATE);

    // Getting fluid density
    density = solver_container[FLOW_SOL]->GetNodes()->GetDensity(iPoint);
    
    // Computing the axial, relative Reynolds number
    // TODO: Get dynamic viscosity from fluid model
    Re_ax = abs(density * W_mag * (ax - ax_le)) / (1.716E-5);
    if(Re_ax == 0.0){
			Re_ax = (0.001 * W_mag * density) / (1.716E-5);
		}
    // Computing the blade friction factor
    // TODO: Allow for the user to set the coefficients
    C_f = 0.0592 * pow(Re_ax, -0.2);
    // Computing the parallel, loss generating body force
    F_p = -C_f * (1 / pitch) * (1 / abs(Nt)) * (1 / b) * W_mag * W_mag;

    // Transforming the normal and force component to cyllindrical coordinates
    F_ax = F_n * (cos(delta) * Nx - sin(delta) * (W_px / (W_p + 1e-6))) + F_p * W_ax / (W_mag + 1e-6);		// Axial body-force component
    F_r = F_n * (cos(delta) * Nr - sin(delta) * (W_pr / (W_p + 1e-6))) + F_p * W_r / (W_mag + 1e-6);		// Radial body-force component
    F_th = F_n * (cos(delta) * Nt - sin(delta) * (W_pth / (W_p + 1e-6)))+ F_p * W_th / (W_mag + 1e-6);	// Tangential body-force component
    e_source = rotFac * Rotation_rate * radius * F_th;				// Energy source term
    
    // Appending Cartesial body-forces to body-force vector
    for(int iDim=0; iDim<nDim; iDim++){
        F[iDim] = density * (F_ax*(nodes->GetAxialProjection(iPoint, iDim))
                          + F_th*(nodes->GetTangentialProjection(iPoint, iDim))
                          + F_r*(nodes->GetRadialProjection(iPoint, iDim)));
    }
    
    // Appending source term values to the body-force source term vector
    BFM_sources.at(0) = 0.0;	// empty density component

    // Appending momentum source terms
    for(int iDim = 0; iDim < nDim; iDim ++) {
        BFM_sources.at(iDim+1) = F[iDim];
    }
    // Appending energy source term
    BFM_sources.at(nDim+1) = density * e_source;

    // Computing source terms due to metal blockage
    ComputeBlockageSources(solver_container, iPoint, BFM_sources);
    
}

void CBFMSolver::ComputeBlockageSources(CSolver **solver_container, unsigned long iPoint, vector<su2double>&BFM_sources){
    /*
     This function computes the flow source terms due to flow obstruction caused by blade thickness effects.
     The resulting source terms are added to the BFM source terms according to Thollets BFM formulation
    */
    su2double density, // Fluid density [kg m^-3]
    energy, // Fluid energy [m^-1 s^-2]
    blockage_gradient; // Metal blockage factor gradient [m^-1]
    su2double velocity, velocity_j; // Flow velocity components [m s^-1]
    su2double b; // Metal blockage factor [-]
    su2double source_density{0}, // Density source term[kg m^-3 s^-1]
    source_energy{0}, // Energy source term
    source_momentum[nDim]; // Momentum source terms

    // Getting fluid density and energy
    density = solver_container[FLOW_SOL]->GetNodes()->GetDensity(iPoint);
    energy = solver_container[FLOW_SOL]->GetNodes()->GetEnergy(iPoint);

    // Getting interpolated metal blockage factor
    b = nodes->GetAuxVar(iPoint, I_BLOCKAGE_FACTOR);

    // Looping over dimensions to compute the divergence source terms
    for(unsigned short iDim=0; iDim<nDim; ++iDim){
        // Setting momentum source term to zero
        source_momentum[iDim] = 0;

        // Getting flow velocity
        velocity = solver_container[FLOW_SOL]->GetNodes()->GetVelocity(iPoint, iDim);
        // Getting metal blockage gradient
        blockage_gradient = nodes->GetAuxVarGradient(iPoint, I_BLOCKAGE_FACTOR, iDim);

        // Updating density and energy source terms
        source_density += density * velocity * blockage_gradient / b;
        source_energy += density * energy * velocity * blockage_gradient / b;

        // Looping over dimensions to compute momentum source terms due to metal blockage
        for(unsigned short jDim=0; jDim<nDim; jDim++){
            velocity_j = solver_container[FLOW_SOL]->GetNodes()->GetVelocity(iPoint, jDim);
            source_momentum[iDim] += density * velocity * velocity_j * blockage_gradient / b;
        }
    }

    // Subtracting metal blockage source terms from body force source terms
    BFM_sources.at(0) -= source_density;
    for(unsigned short iDim=0; iDim<nDim; ++iDim){
        BFM_sources.at(iDim + 1) -= source_momentum[iDim];
    }
    BFM_sources.at(nDim + 1) -= source_energy;


}
void CBFMSolver::ComputeCylProjections(CGeometry *geometry, CConfig *config)
{
    
	su2double *Coord, *Rotation_vector;
	su2double rot_dot_x, rot_dot_rot;
	su2double ax, radius;
	su2double axial_vector[nDim], radial_vector[nDim];
    
	for(unsigned long iPoint=0; iPoint<nPoint; iPoint++){
		// Getting node Cartesian coordinates
		Coord = geometry->nodes->GetCoord(iPoint);

		rot_dot_x = 0;	// dot product between axial unit vector and coordinate vector
		rot_dot_rot = 0;	// dot product of axial unit vector
		radius = 0;	// Radial coordinate
		ax = 0;		// Axial coordinate

		// Computing axial dot products
		for(int iDim = 0; iDim < nDim; iDim ++){
			rot_dot_x += config->GetBFM_Axis(iDim)*Coord[iDim];
			rot_dot_rot += (config->GetBFM_Axis(iDim)) * (config->GetBFM_Axis(iDim));
		}
		//

       
        
		for(int iDim = 0; iDim < nDim; iDim ++){
			axial_vector[iDim] = config->GetBFM_Axis(iDim)*Coord[iDim]; // Parallel projection of coordinate vector onto rotation axis
			radial_vector[iDim] = Coord[iDim] - axial_vector[iDim];	// Orthogonal projection of coordinate vector onto rotation axis
			ax += axial_vector[iDim]*axial_vector[iDim];	
			radius += radial_vector[iDim]*radial_vector[iDim];
            su2double rotrate = config->GetBFM_Axis(iDim);
            nodes->SetAxialProjection(iPoint, iDim, config->GetBFM_Axis(iDim));
        }
        
		// Computing absolute values for cylindrical coordinates
		radius = sqrt(radius);	// Computing radial coordinate
		ax = sqrt(ax);	// Computing axial coordinate
		if(rot_dot_x < 0.0){ax = -ax;}	// In case the coordinate parallel projection is negative, the axial coordinate is flipped

        nodes->SetAuxVar(iPoint, I_AXIAL_COORDINATE, ax);
        nodes->SetAuxVar(iPoint, I_RADIAL_COORDINATE, radius);
		// Storing cylindrical coordinates in class property
		// cyl_coordinates[iPoint][0] = ax;
		// cyl_coordinates[iPoint][1] = radius;

		// Computation of tangential and radial projection vectors
		int i_left, i_right;
		for(int iDim=0; iDim<nDim; iDim++){

			// Computing tangential unit vector components through cross product of radial and axial vectors
			i_left = (iDim + 1) % nDim;
			i_right = (iDim + 2) % nDim;
            nodes->SetTangentialProjection(iPoint, iDim, (config->GetBFM_Axis(i_left)*radial_vector[i_right] - config->GetBFM_Axis(i_right)*radial_vector[i_left])/radius);
			//Radial projection unit vector is normalized radial coordinate vector
            nodes->SetRadialProjection(iPoint, iDim, radial_vector[iDim]/radius);
		}
			
	}
}
void CBFMSolver::ComputeRelativeVelocity(CGeometry *geometry, CSolver **solver_container, CConfig *config)
{
    su2double *Coord, U_i,*Geometric_Parameters;
	su2double W_ax, W_r, W_th, rotFac;
    su2double Rotation_rate = Omega;

	for(unsigned long iPoint=0; iPoint<geometry->GetnPoint(); iPoint++){

		// Obtaining node coordinates
		Coord = geometry->nodes->GetCoord(iPoint);

		// Obtaining solution flow variables
		

		// Computing relative velocity components in axial, radial and tangential direction
		W_ax = 0;	// Relative, axial velocity
		W_r = 0;	// Relative, radial velocity
		W_th = 0;	// Relative, tangential velocity

		// Obtaining BFM geometric parameters from node
		rotFac = nodes->GetAuxVar(iPoint, I_ROTATION_FACTOR);	// Rotation factor, distinguish between rotor and stator

		// Adding to the relative velocity components through dot product of absolute velocity and respective unit vector
		for(int iDim = 0; iDim < nDim; iDim++){
            U_i = solver_container[FLOW_SOL]->GetNodes()->GetVelocity(iPoint, iDim);
            W_ax += nodes->GetAxialProjection(iPoint, iDim) * U_i;
            W_th += nodes->GetTangentialProjection(iPoint, iDim) * U_i;
            W_r += nodes->GetRadialProjection(iPoint, iDim) * U_i;
		}
		// Induced velocity due to rotation is subtracted from tangential velocity
		W_th -= rotFac * Rotation_rate * nodes->GetAuxVar(iPoint, I_RADIAL_COORDINATE);

		// Storing relative velocity vector in class data structure
        nodes->SetRelativeVelocity(iPoint, 0, W_ax);
        nodes->SetRelativeVelocity(iPoint, 1, W_th);
        nodes->SetRelativeVelocity(iPoint, 2, W_r);
		
	}
}

void CBFMSolver::ComputeRelativeVelocity(CSolver **solver_container, unsigned long iPoint, vector<su2double*> &W_cyl){
    su2double *Coord, U_i,*Geometric_Parameters;
	su2double W_ax, W_r, W_th, rotFac;
    su2double Rotation_rate = Omega;

    // Obtaining solution flow variables
    

    // Computing relative velocity components in axial, radial and tangential direction
    W_ax = 0;	// Relative, axial velocity
    W_r = 0;	// Relative, radial velocity
    W_th = 0;	// Relative, tangential velocity

    // Obtaining BFM geometric parameters from node
    rotFac = nodes->GetAuxVar(iPoint, I_ROTATION_FACTOR);	// Rotation factor, distinguish between rotor and stator

    // Adding to the relative velocity components through dot product of absolute velocity and respective unit vector
    for(int iDim = 0; iDim < nDim; iDim++){
        U_i = solver_container[FLOW_SOL]->GetNodes()->GetVelocity(iPoint, iDim);
        W_ax += nodes->GetAxialProjection(iPoint, iDim) * U_i;
        W_th += nodes->GetTangentialProjection(iPoint, iDim) * U_i;
        W_r += nodes->GetRadialProjection(iPoint, iDim) * U_i;
    }
    // Induced velocity due to rotation is subtracted from tangential velocity
    W_th -= rotFac * Rotation_rate * nodes->GetAuxVar(iPoint, I_RADIAL_COORDINATE);

    *W_cyl.at(0) = W_ax;
    *W_cyl.at(1) = W_th;
    *W_cyl.at(2) = W_r;
    // Storing relative velocity vector in class data structure
    nodes->SetRelativeVelocity(iPoint, 0, W_ax);
    nodes->SetRelativeVelocity(iPoint, 1, W_th);
    nodes->SetRelativeVelocity(iPoint, 2, W_r);
		
	
}


su2double CBFMSolver::ComputeKMach(CSolver **solver_container, unsigned long iPoint, vector<su2double*>W_cyl){
    su2double Temperature = solver_container[FLOW_SOL]->GetNodes()->GetTemperature(iPoint);
    su2double SoS = solver_container[FLOW_SOL]->GetNodes()->GetSoundSpeed(iPoint);
    su2double W_mag{0};
    for(unsigned short iDim=0; iDim<nDim; ++iDim){
        W_mag += (*W_cyl.at(iDim))*(*W_cyl.at(iDim));
    }
    W_mag = sqrt(W_mag);
    su2double M_rel = W_mag / SoS;
    if(M_rel == 1.0){
        M_rel -= 1e-6;
    }
    su2double K_prime = M_rel < 1 ? 1/sqrt(1 - M_rel * M_rel) : 2 / (PI_NUMBER * sqrt(M_rel*M_rel - 1));
    su2double K_Mach = K_prime <= 3 ? K_prime : 3;
    if (isnan(K_Mach)) K_Mach = 3;

    return K_Mach;
}