/*!
 * \file output_cgns.cpp
 * \brief Main subroutines for output solver information
 * \author T. Economon, M. Colonno
 * \version 6.2.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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

#include "../include/output_structure.hpp"

void COutput::SetCGNS_Coordinates(CConfig *config, CGeometry *geometry, unsigned short iZone) {
  
#ifdef HAVE_CGNS
  
  /*--- local CGNS variables ---*/
  int cgns_file, cgns_coord, element_dims, physical_dims, cgns_err;
  unsigned long iExtIter = config->GetExtIter();
  string base_file, buffer, elements_name;
  stringstream name, results_file;
  bool unsteady = config->GetUnsteady_Simulation();
  cgsize_t isize[3][1];
  
  /*--- Create CGNS base file name ---*/
  base_file = config->GetFlow_FileName();
  
  /*--- Add CGNS extension. ---*/
  base_file = base_file.append(".cgns");
  
  /*--- Create CGNS results file name ---*/
  if (unsteady) {
    
    buffer = config->GetFlow_FileName();
    
    results_file.str(string()); results_file << buffer;
    if (((int)iExtIter >= 0) && ((int)iExtIter < 10))      results_file << "_0000" << iExtIter;
    if (((int)iExtIter >= 10) && ((int)iExtIter < 100))    results_file << "_000" << iExtIter;
    if (((int)iExtIter >= 100) && ((int)iExtIter < 1000))    results_file << "_00" << iExtIter;
    if (((int)iExtIter >= 1000) && ((int)iExtIter < 10000))  results_file << "_0" << iExtIter;
    if ((int)iExtIter >= 10000)              results_file << iExtIter;
    results_file << ".cgns";
  }
  
  /*--- Write base file if not already done ---*/
  if (!wrote_base_file) {
    
    /*--- Write base file ---*/
    cgns_err = cg_open((char *)base_file.c_str(), CG_MODE_MODIFY, &cgns_file);
    if (cgns_err) cg_error_print();
    
    element_dims = geometry->GetnDim();    // Currently (release 2.0) only all-2D or all-3D zones permitted
    physical_dims = element_dims;
    
    isize[0][0] = (cgsize_t)nGlobal_Poin;        // vertex size
    isize[1][0] = (cgsize_t)nGlobal_Elem;        // cell size
    isize[2][0] = 0;            // boundary vertex size (zero if elements not sorted)
    
    cgns_err = cg_goto(cgns_file, cgns_base,"Zone_t", cgns_zone,"end");
    if (cgns_err) cg_error_print();

    
    /*--- write CGNS node coordinates ---*/
    cgns_err = cg_coord_write(cgns_file, cgns_base, cgns_zone, RealDouble,"CoordinateX", Coords[0], &cgns_coord);
    if (cgns_err) cg_error_print();
    cgns_err = cg_coord_write(cgns_file, cgns_base, cgns_zone, RealDouble,"CoordinateY", Coords[1], &cgns_coord);
    if (cgns_err) cg_error_print();
    if (geometry->GetnDim() == 3) {
      cgns_err = cg_coord_write(cgns_file, cgns_base, cgns_zone, RealDouble,"CoordinateZ", Coords[2], &cgns_coord);
      if (cgns_err) cg_error_print();
    }
    
    cgns_err = cg_close(cgns_file);
    if (cgns_err) cg_error_print();
    
        wrote_base_file = true;
    
  }
  
  /*--- Set up results file for this time step if necessary ---*/
  if (unsteady) {
    
    cgns_err = cg_open((char *)results_file.str().c_str(), CG_MODE_WRITE, &cgns_file);

    element_dims = geometry->GetnDim();    // Currently only all-2D or all-3D zones permitted
    physical_dims = element_dims;

    /*--- write CGNS base data (one base assumed currently) ---*/
    cgns_err = cg_base_write(cgns_file,"SU2 Base", element_dims, physical_dims, &cgns_base_results);
    if (cgns_err) cg_error_print();

    isize[0][0] = (cgsize_t)geometry->GetGlobal_nPointDomain();        // vertex size
    isize[1][0] = (cgsize_t)nGlobal_Elem;        // cell size
    isize[2][0] = 0;            // boundary vertex size (zero if elements not sorted)

    /*--- write CGNS zone data ---*/
    cgns_err = cg_zone_write(cgns_file, cgns_base_results,"SU2 Zone", isize[0],Unstructured, &cgns_zone_results);
    if (cgns_err) cg_error_print();

    cgns_err = cg_goto(cgns_file, cgns_base_results,"Zone_t", cgns_zone_results,"end");
    if (cgns_err) cg_error_print();

    /*--- Write CGNS node coordinates, if appliciable ---*/
    if (config->GetGrid_Movement()) {
      
      /*--- write CGNS node coordinates ---*/
      cgns_err = cg_coord_write(cgns_file, cgns_base_results, cgns_zone_results, RealDouble,"x", Coords[0], &cgns_coord);
      if (cgns_err) cg_error_print();
      cgns_err = cg_coord_write(cgns_file, cgns_base_results, cgns_zone_results, RealDouble,"y", Coords[1], &cgns_coord);
      if (cgns_err) cg_error_print();
      if (geometry->GetnDim() == 3) {
        cgns_err = cg_coord_write(cgns_file, cgns_base_results, cgns_zone_results, RealDouble,"z", Coords[2], &cgns_coord);
        if (cgns_err) cg_error_print();
      }
    }
    else {
      /*--- Write a CGNS link for the node coordinates ---*/
      cgns_err = cg_link_write("GridCoordinates",(char *)base_file.c_str(),"/SU2 Base/SU2 Zone/GridCoordinates");
      if (cgns_err) cg_error_print();
    }
    
    /*--- Write a CGNS link for each element type connectivity ---*/
    if (nGlobal_Tria > 0) cgns_err = cg_link_write("Triangle Elements",(char *)base_file.c_str(),"/SU2 Base/SU2 Zone/Triangle Elements");
    if (nGlobal_Quad > 0) cgns_err = cg_link_write("Quadrilateral Elements",(char *)base_file.c_str(),"/SU2 Base/SU2 Zone/Quadrilateral Elements");
    if (nGlobal_Tetr > 0) cgns_err = cg_link_write("Tetrahedral Elements",(char *)base_file.c_str(),"/SU2 Base/SU2 Zone/Tetrahedral Elements");
    if (nGlobal_Hexa > 0) cgns_err = cg_link_write("Hexahedral Elements",(char *)base_file.c_str(),"/SU2 Base/SU2 Zone/Hexahedral Elements");
    if (nGlobal_Pyra > 0) cgns_err = cg_link_write("Pyramid Elements",(char *)base_file.c_str(),"/SU2 Base/SU2 Zone/Pyramid Elements");
    if (nGlobal_Pris > 0) cgns_err = cg_link_write("Prism Elements",(char *)base_file.c_str(),"/SU2 Base/SU2 Zone/Prism Elements");
    if (nGlobal_Line > 0) cgns_err = cg_link_write("Line Elements",(char *)base_file.c_str(),"/SU2 Base/SU2 Zone/Line Elements");
    if (cgns_err) cg_error_print();

    
    /*--- Close CGNS file ---*/
    cgns_err = cg_close(cgns_file);
    if (cgns_err) cg_error_print();
    
  }


  
#else // Not built with CGNS support
  
  cout << "CGNS file requested but SU2 was built without CGNS support. No file written" << "\n"; 
  
#endif
  
}

void COutput::SetCGNS_Connectivity(CConfig *config, CGeometry *geometry, unsigned short iZone) {
  
#ifdef HAVE_CGNS
  
  /*--- local CGNS variables ---*/
  int cgns_file, element_dims, physical_dims, cgns_err;
  int cgns_section;
  unsigned long iExtIter = config->GetExtIter();
  string base_file, buffer, elements_name;
  stringstream name, results_file;
  bool unsteady = config->GetUnsteady_Simulation();
  cgsize_t isize[3][1], elem_start, elem_end;
  
  /*--- Create CGNS base file name ---*/
  base_file = config->GetFlow_FileName();
  
  /*--- Add CGNS extension. ---*/
  base_file = base_file.append(".cgns");
  
  /*--- Create CGNS results file name ---*/
  if (unsteady) {
    
    buffer = config->GetFlow_FileName();
    
    results_file.str(string()); results_file << buffer;
    if (((int)iExtIter >= 0) && ((int)iExtIter < 10))      results_file << "_0000" << iExtIter;
    if (((int)iExtIter >= 10) && ((int)iExtIter < 100))    results_file << "_000" << iExtIter;
    if (((int)iExtIter >= 100) && ((int)iExtIter < 1000))    results_file << "_00" << iExtIter;
    if (((int)iExtIter >= 1000) && ((int)iExtIter < 10000))  results_file << "_0" << iExtIter;
    if ((int)iExtIter >= 10000)              results_file << iExtIter;
    results_file << ".cgns";
  }
  
  /*--- Write base file if not already done ---*/
  if (!wrote_base_file) {
    
    /*--- Write base file ---*/
    cgns_err = cg_open((char *)base_file.c_str(), CG_MODE_WRITE, &cgns_file);
    if (cgns_err) cg_error_print();
    
    element_dims = geometry->GetnDim();    // Currently only all-2D or all-3D zones permitted
    physical_dims = element_dims;
    
    /*--- write CGNS base data (one base assumed currently) ---*/
    cgns_err = cg_base_write(cgns_file,"SU2 Base", element_dims, physical_dims, &cgns_base);
    if (cgns_err) cg_error_print();
    
    /*--- write CGNS descriptor data ---*/
    cgns_err = cg_goto(cgns_file, cgns_base,"end");
    if (cgns_err) cg_error_print();
    
    cgns_err = cg_equationset_write(physical_dims);
    if (cgns_err) cg_error_print();
    
    /*--- Write governing equations to CGNS file ---*/
    cgns_err = cg_goto(cgns_file, cgns_base,"FlowEquationSet_t",1,"end");
    if (cgns_err) cg_error_print();
    switch (config->GetKind_Solver()) {
      case EULER:
        cgns_err = cg_governing_write(Euler); break;
      case NAVIER_STOKES:
        cgns_err = cg_governing_write(NSLaminar); break;
      case RANS:
        cgns_err = cg_governing_write(NSTurbulent); break;
      default:
        break; // cgns_err = cg_governing_write(CG_UserDefined);
    }
    if (cgns_err) cg_error_print();
    
    if (unsteady) cgns_err = cg_simulation_type_write(cgns_file, cgns_base, TimeAccurate);
    else cgns_err = cg_simulation_type_write(cgns_file, cgns_base, NonTimeAccurate);
    if (cgns_err) cg_error_print();
    
    cgns_err = cg_descriptor_write("Solver Information","SU2");
    if (cgns_err) cg_error_print();
    
    isize[0][0] = (cgsize_t)geometry->GetGlobal_nPointDomain(); //;        // vertex size
    isize[1][0] = (cgsize_t)nGlobal_Elem;        // cell size
    isize[2][0] = 0;            // boundary vertex size (zero if elements not sorted)
    
    /*--- write CGNS zone data ---*/
    cgns_err = cg_zone_write(cgns_file, cgns_base,"SU2 Zone", isize[0],Unstructured, &cgns_zone);
    if (cgns_err) cg_error_print();

    cgns_err = cg_goto(cgns_file, cgns_base,"Zone_t", cgns_zone,"end");
    if (cgns_err) cg_error_print();
    
        /*--- Reference Note: CGNS element type list:
     NODE, BAR_2, BAR_3, TRI_3, TRI_6, QUAD_4, QUAD_8, QUAD_9, TETRA_4, TETRA_10, PYRA_5,
     PYRA_14, PENTA_6, PENTA_15, PENTA_18, HEXA_8, HEXA_20, HEXA_27, MIXED, PYRA_13, NGON_n, NFACE_n ---*/
    
    /*--- Write a CGNS section for each element type ---*/
    // ier = cg_section_write(int fn, int B, int Z, char *ElementSectionName, ElementType_t type,
    // cgsize_t start, cgsize_t end, int nbndry, cgsize_t *Elements, int *S);
    
    if (nGlobal_Tria > 0) {
      elem_start = 1; elem_end = (int)nGlobal_Tria;
      cgns_err = cg_section_write(cgns_file, cgns_base, cgns_zone,
                                  "Triangle Elements", TRI_3, elem_start, elem_end,
                                  0,(cgsize_t *)Conn_Tria, &cgns_section);
    }
    if (nGlobal_Quad > 0) {
      elem_start = 1; elem_end = (int)nGlobal_Quad;
      cgns_err = cg_section_write(cgns_file, cgns_base, cgns_zone,"Quadrilateral Elements", QUAD_4,
                                  elem_start, elem_end,0,(cgsize_t *)Conn_Quad, &cgns_section);
    }
    if (nGlobal_Tetr > 0) {
      elem_start = 1; elem_end = (int)nGlobal_Tetr;
      cgns_err = cg_section_write(cgns_file, cgns_base, cgns_zone,"Tetrahedral Elements", TETRA_4,
                                  elem_start, elem_end,0,(cgsize_t *)Conn_Tetr, &cgns_section);
    }
    if (nGlobal_Hexa > 0) {
      elem_start = 1; elem_end = (int)nGlobal_Hexa;
      cgns_err = cg_section_write(cgns_file, cgns_base, cgns_zone,"Hexahedral Elements", HEXA_8,
                                  elem_start, elem_end,0,(cgsize_t *)Conn_Hexa, &cgns_section);
    }
    if (nGlobal_Pyra > 0) {
      elem_start = 1; elem_end = (int)nGlobal_Pyra;
      cgns_err = cg_section_write(cgns_file, cgns_base, cgns_zone,"Pyramid Elements", PYRA_5,
                                  elem_start, elem_end,0,(cgsize_t *)Conn_Pyra, &cgns_section);
    }
    if (nGlobal_Pris > 0) {
      elem_start = 1; elem_end = (int)nGlobal_Pris;
      cgns_err = cg_section_write(cgns_file, cgns_base, cgns_zone,"Prism Elements", PENTA_6,
                                  elem_start, elem_end,0,(cgsize_t *)Conn_Pris, &cgns_section);
    }
    if (nGlobal_Line > 0) {
      elem_start = 1; elem_end = (int)nGlobal_Line;
      cgns_err = cg_section_write(cgns_file, cgns_base, cgns_zone,"Line Elements", BAR_2,
                                  elem_start, elem_end,0,(cgsize_t *)Conn_Line, &cgns_section);
    }
    if (cgns_err) cg_error_print();
    
    
    cgns_err = cg_close(cgns_file);
    if (cgns_err) cg_error_print();
    
  }
  
#else // Not built with CGNS support
  
  cout << "CGNS file requested but SU2 was built without CGNS support. No file written" << "\n"; 
  
#endif
  
}

void COutput::SetCGNS_Solution(CConfig *config, CGeometry *geometry, unsigned short iZone) {
  
#ifdef HAVE_CGNS
  
  /*--- local CGNS variables ---*/
  int cgns_file, cgns_flow, cgns_field, cgns_err;
//  int element_dims;
  unsigned long jVar, iVar, iExtIter = config->GetExtIter();
  string base_file, buffer, elements_name;
  stringstream name, results_file;
  bool unsteady = config->GetUnsteady_Simulation();
  
  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  
  /*--- Create CGNS base file name ---*/
  base_file = config->GetFlow_FileName();
  
  /*--- Add CGNS extension. ---*/
  base_file = base_file.append(".cgns");
  
  /*--- Create CGNS results file name ---*/
  if (unsteady) {
    
    buffer = config->GetFlow_FileName();
    
    results_file.str(string()); results_file << buffer;
    if (((int)iExtIter >= 0) && ((int)iExtIter < 10))      results_file << "_0000" << iExtIter;
    if (((int)iExtIter >= 10) && ((int)iExtIter < 100))    results_file << "_000" << iExtIter;
    if (((int)iExtIter >= 100) && ((int)iExtIter < 1000))    results_file << "_00" << iExtIter;
    if (((int)iExtIter >= 1000) && ((int)iExtIter < 10000))  results_file << "_0" << iExtIter;
    if ((int)iExtIter >= 10000)              results_file << iExtIter;
    results_file << ".cgns";
  }
    
    if (!unsteady) {
      
      /*--- Write base file ---*/
      cgns_err = cg_open((char *)base_file.c_str(), CG_MODE_MODIFY, &cgns_file);
      if (cgns_err) cg_error_print();
      
//      element_dims = geometry->GetnDim();    // Currently (release 2.0) only all-2D or all-3D zones permitted
      
      /*--- write CGNS descriptor data ---*/
      cgns_err = cg_goto(cgns_file, cgns_base,"end");
      if (cgns_err) cg_error_print();
      
      /*--- Create a CGNS solution node ---*/
      cgns_err = cg_sol_write(cgns_file, cgns_base, cgns_zone,"Solution", Vertex, &cgns_flow);
      if (cgns_err) cg_error_print();
      
      cgns_err = cg_goto(cgns_file, cgns_base,"Zone_t", cgns_zone,"FlowSolution_t", cgns_flow,"end");
      if (cgns_err) cg_error_print();
      
      cgns_err = cg_gridlocation_write(Vertex);
      if (cgns_err) cg_error_print();
    }

    
    //wrote_CGNS_base = true;
    else {
  
  /*--- Set up results file for this time step if necessary ---*/
    
    cgns_err = cg_open((char *)results_file.str().c_str(), CG_MODE_MODIFY, &cgns_file);
    
//    element_dims = geometry->GetnDim();    // Currently (release 2.0) only all-2D or all-3D zones permitted
//    
//    /*--- write CGNS base data (one base assumed currently) ---*/
//    cgns_err = cg_base_write(cgns_file,"SU2 Base", element_dims, physical_dims, &cgns_base);
//    if (cgns_err) cg_error_print();
        
//    /*--- write CGNS zone data ---*/
//    cgns_err = cg_zone_write(cgns_file, cgns_base,"SU2 Zone", isize[0],Unstructured, &cgns_zone);
//    if (cgns_err) cg_error_print();
    
    cgns_err = cg_goto(cgns_file, cgns_base_results,"Zone_t", cgns_zone_results,"end");
    if (cgns_err) cg_error_print();
    
    /*--- Write a CGNS solution node for this time step ---*/
    cgns_err = cg_sol_write(cgns_file, cgns_base_results, cgns_zone_results,"Solution", Vertex, &cgns_flow);
    if (cgns_err) cg_error_print();
    
    cgns_err = cg_goto(cgns_file, cgns_base_results,"Zone_t", cgns_zone_results,"FlowSolution_t", cgns_flow,"end");
    if (cgns_err) cg_error_print();
    
    cgns_err = cg_gridlocation_write(Vertex);
    if (cgns_err) cg_error_print();
    
      cgns_base = cgns_base_results;
      cgns_zone = cgns_zone_results;
  }
//  else {
//    
//    /*--- Open CGNS file for soltuion writing ---*/
//    cgns_err = cg_open((char *)base_file.c_str(), CG_MODE_MODIFY, &cgns_file);
//    cgns_base = 1; cgns_zone = 1; cgns_flow = 1;  // fix for multiple zones
//    
//  }
  
  /*  Reference Note on solution variables:
   index 0 --> (nVar_Consv-1)      = Conservative Variables
   nVar_Consv --> (2*nVar_Consv-1)    = Conservative Variable Residuals
   (2*nVar_Consv-1)+          = Additional p, M, T, laminar, eddy depending on solver used */
  
  /*--- Write conservative variables to CGNS file ---*/
  for (iVar = 0; iVar < nVar_Consv; iVar++) {
    name.str(string()); name << "Conservative Variable " << iVar+1;
    cgns_err = cg_field_write(cgns_file, cgns_base, cgns_zone, cgns_flow, RealDouble,(char *)name.str().c_str(), Data[iVar], &cgns_field);
    if (cgns_err) cg_error_print();
  }
  
  /*--- Write primitive variable residuals to CGNS file ---*/
  if (config->GetWrt_Limiters()) {
    for (jVar = 0; jVar < nVar_Consv; jVar++) {
      name.str(string()); name << "Primitive Limiter " << jVar+1;
      cgns_err = cg_field_write(cgns_file, cgns_base, cgns_zone, cgns_flow, RealDouble,(char *)name.str().c_str(), Data[iVar], &cgns_field); iVar++;
      if (cgns_err) cg_error_print();
    }
  }
  
  /*--- Write conservative variable residuals to CGNS file ---*/
  if (config->GetWrt_Residuals()) {
    for (jVar = 0; jVar < nVar_Consv; jVar++) {
      name.str(string()); name << "Conservative Residual " << jVar+1;
      cgns_err = cg_field_write(cgns_file, cgns_base, cgns_zone, cgns_flow, RealDouble,(char *)name.str().c_str(), Data[iVar], &cgns_field); iVar++;
      if (cgns_err) cg_error_print();
    }
  }
  
  /*--- Write grid velocities to CGNS file, if applicable ---*/
  if (config->GetGrid_Movement()) {
    cgns_err = cg_field_write(cgns_file, cgns_base, cgns_zone, cgns_flow, RealDouble,"Grid Velocity X", Data[iVar], &cgns_field); iVar++;
    if (cgns_err) cg_error_print();
    cgns_err = cg_field_write(cgns_file, cgns_base, cgns_zone, cgns_flow, RealDouble,"Grid Velocity Y", Data[iVar], &cgns_field); iVar++;
    if (cgns_err) cg_error_print();
    if (geometry->GetnDim() == 3) {
      cgns_err = cg_field_write(cgns_file, cgns_base, cgns_zone, cgns_flow, RealDouble,"Grid Velocity Z", Data[iVar], &cgns_field); iVar++;
      if (cgns_err) cg_error_print();
    }
  }
  
  if (compressible) {
    switch (config->GetKind_Solver()) {
        
        /*--- Write pressure and Mach data to CGNS file ---*/
      case EULER:
        cgns_err = cg_field_write(cgns_file, cgns_base, cgns_zone, cgns_flow, RealDouble,"Pressure", Data[iVar], &cgns_field); iVar++;
        if (cgns_err) cg_error_print();
        cgns_err = cg_field_write(cgns_file, cgns_base, cgns_zone, cgns_flow, RealDouble,"Mach", Data[iVar], &cgns_field); iVar++;
        if (cgns_err) cg_error_print();
        break;
        
        /*--- Write temperature and laminar viscosity to CGNS file, if applicable ---*/
      case NAVIER_STOKES:
        cgns_err = cg_field_write(cgns_file, cgns_base, cgns_zone, cgns_flow, RealDouble,"Pressure", Data[iVar], &cgns_field); iVar++;
        if (cgns_err) cg_error_print();
        cgns_err = cg_field_write(cgns_file, cgns_base, cgns_zone, cgns_flow, RealDouble,"Mach", Data[iVar], &cgns_field); iVar++;
        if (cgns_err) cg_error_print();
        cgns_err = cg_field_write(cgns_file, cgns_base, cgns_zone, cgns_flow, RealDouble,"Temperature", Data[iVar], &cgns_field); iVar++;
        if (cgns_err) cg_error_print();
        cgns_err = cg_field_write(cgns_file, cgns_base, cgns_zone, cgns_flow, RealDouble,"Viscosity", Data[iVar], &cgns_field); iVar++;
        if (cgns_err) cg_error_print();
        break;
        
        /*--- Write eddy viscosity to CGNS file, if applicable ---*/
      case RANS:
        cgns_err = cg_field_write(cgns_file, cgns_base, cgns_zone, cgns_flow, RealDouble,"Pressure", Data[iVar], &cgns_field); iVar++;
        if (cgns_err) cg_error_print();
        cgns_err = cg_field_write(cgns_file, cgns_base, cgns_zone, cgns_flow, RealDouble,"Mach", Data[iVar], &cgns_field); iVar++;
        if (cgns_err) cg_error_print();
        cgns_err = cg_field_write(cgns_file, cgns_base, cgns_zone, cgns_flow, RealDouble,"Temperature", Data[iVar], &cgns_field); iVar++;
        if (cgns_err) cg_error_print();
        cgns_err = cg_field_write(cgns_file, cgns_base, cgns_zone, cgns_flow, RealDouble,"Viscosity", Data[iVar], &cgns_field); iVar++;
        if (cgns_err) cg_error_print();
        cgns_err = cg_field_write(cgns_file, cgns_base, cgns_zone, cgns_flow, RealDouble,"Eddy Viscosity", Data[iVar], &cgns_field); iVar++;
        if (cgns_err) cg_error_print();
        break;
        
      default:
        SU2_MPI::Error("Unrecognized equation type", CURRENT_FUNCTION);
        break;
      }
  }  
  
  /*--- Close CGNS file ---*/
  cgns_err = cg_close(cgns_file);
  if (cgns_err) cg_error_print();
  
#else // Not built with CGNS support
  
  cout << "CGNS file requested but SU2 was built without CGNS support. No file written" << "\n"; 
  
#endif
  
}
