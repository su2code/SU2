/*!
 * \file CCGNSBFFileWriter.cpp
 * \brief Filewriter class for CGNS format.
 * \author E. Saetta, L. Russo, R. Tognaccini
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

#include "../../../include/output/filewriter/CCGNSFileWriter.hpp"
#include "../../../../externals/cgns/cgnslib.h"
#include "../../../../Common/include/geometry/CGeometry.hpp"
#include "../../../../Common/include/CConfig.hpp"
#include "../../../include/solvers/CSolver.hpp"
#include "../../../../Common/include/option_structure.hpp"
#include "../../../include/output/COutputLegacy.hpp"

#include <algorithm>

const string CCGNSFileWriter::fileExt = ".cgns";

CCGNSFileWriter::CCGNSFileWriter(string valFileName, CParallelDataSorter *valDataSorter) :
  CFileWriter(std::move(valFileName), valDataSorter, fileExt){}

CCGNSFileWriter::~CCGNSFileWriter(){}

void CCGNSFileWriter::Write_Data_CGNS(CConfig *config){

  if (!dataSorter->GetConnectivitySorted()){
    SU2_MPI::Error("Connectivity must be sorted.", CURRENT_FUNCTION);
  }

  #ifdef HAVE_MPI
      SU2_MPI::Barrier(MPI_COMM_WORLD);
  #endif

  startTime = SU2_MPI::Wtime();
  int iProcessor, nProcessor = size, index;
  unsigned long Buffer_Send_nVertex[1], *Buffer_Recv_nVertex = nullptr;

  unsigned long nLocalVertex = 0, MaxLocalVertex = 0;
  nLocalVertex   = dataSorter->GetnPoints();
  unsigned long nPoint_Global = dataSorter->GetnPointsGlobal();
  Buffer_Send_nVertex[0] = nLocalVertex;
  if (rank == MASTER_NODE) Buffer_Recv_nVertex = new unsigned long[nProcessor];

  SU2_MPI::Allreduce(&nLocalVertex, &MaxLocalVertex, 1,
                       MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Gather(&Buffer_Send_nVertex, 1, MPI_UNSIGNED_LONG,
                    Buffer_Recv_nVertex,  1, MPI_UNSIGNED_LONG,
                    MASTER_NODE, MPI_COMM_WORLD);

  const vector<string> fieldNames = dataSorter->GetFieldNames();
  unsigned short nVar = fieldNames.size();

  su2double *buffLocData = new su2double[MaxLocalVertex * nVar]();
  su2double *buffGlobData = nullptr;

  index = 0;
  for (int iPoint = 0; iPoint < nLocalVertex; iPoint++) {
      for(int varCount = 0; varCount < nVar; varCount++){
          buffLocData[index] = dataSorter->GetData(varCount, iPoint);
          index++;
      }
  }

  #ifdef HAVE_MPI
      SU2_MPI::Barrier(MPI_COMM_WORLD);
  #endif

  if (rank == MASTER_NODE)
      buffGlobData = new su2double[nProcessor * nVar * MaxLocalVertex]();

  SU2_MPI::Gather(buffLocData, (int)MaxLocalVertex*nVar, MPI_DOUBLE,
                  buffGlobData, (int)MaxLocalVertex*nVar, MPI_DOUBLE,
                  MASTER_NODE, MPI_COMM_WORLD);

  su2double **GlobalData;
  if (rank == MASTER_NODE){
      GlobalData = new su2double*[nVar];
      for(int i = 0; i < nVar; i++)
          GlobalData[i] = new su2double[dataSorter->GetnPointsGlobal()];

    index = 0;
    int index2 = 0;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        for(int j = 0; j < Buffer_Recv_nVertex[iProcessor]; j++){
            index  = (iProcessor*MaxLocalVertex*fieldNames.size() +
                                              j*fieldNames.size());
            for (int i = 0; i < nVar; i++)
                GlobalData[i][j+index2] = buffGlobData[index + i];
        }
        index2 += Buffer_Recv_nVertex[iProcessor];
    }
//    delete [] Buffer_Recv_nVertex;
//    delete [] buffGlobData;
  }
// delete [] buffLocData;

// Free memory: *Buffer_Recv_nVertex(solo master), *buffGlobData(solo master), *buffLocData
  #ifdef HAVE_MPI
      SU2_MPI::Barrier(MPI_COMM_WORLD);
  #endif

  if (rank == MASTER_NODE){
  #ifdef HAVE_CGNS

    /*---    Internal Variables    ---*/
    int cgns_err, cgns_file, cgns_file1, cgns_flow, cgns_field, nbndry, NormalIndex, bc_wall, bc_far, bc_nulli, coord_index;
    int nbases, base_number, cell_dim, phys_dim, nzones, zone_number, ncoords, nsections, i, S;
    int bcnum, parent_flag;
    cgsize_t isize[3][1];
    cgsize_t start, end, ElementDataSize, NormalListSize, npnts, j, k, z, p, q, NewElementDataSize, dimensions;
    cgsize_t *ConnectOffset;
    cgsize_t pnts[1][2];
    cgsize_t A[4], B[4];
    cgsize_t *Elements = NULL;
    cgsize_t *NewElements = NULL;
    PointSetType_t ptset_type;
    BCType_t bocotype;
    DataType_t datatype;
    GridLocation_t location;
    ElementType_t type;
    DataClass_t dataclass;

    unsigned short iDim = 0, nDim = dataSorter->GetnDim();
    unsigned short varStart = 2;
    if (nDim == 3) varStart++;
    char boconame;
    char basename[] = "Base";
    char zonename[] = "Zone 1";
    char coordname[32];
    char ElementSectionName[32];
    unsigned long jVar;
    unsigned short NVar, iVar, iPoint;
    stringstream name, secname, ss;
    string mesh_file, base_file, wall_string, far_string, null_string, str;
    unsigned long nLocalVertex   = dataSorter->GetnPoints();
    nElem = dataSorter->GetnElemGlobal();

    isize[0][0] = (cgsize_t)nPoint_Global;        // vertex size
    isize[1][0] = (cgsize_t)nElem;                // cell size
    isize[2][0] = 0;                              //(cgsize_t)nSurf_Elem; boundary vertex size (zero if elements not sorted)

    /*---    Open Mesh File    ---*/
    mesh_file = config->GetMesh_FileName();
    cgns_err = cg_open(mesh_file.c_str(), CG_MODE_READ, &cgns_file);
    if(cgns_err) cg_error_print();
    nbases = 1;
    nzones = 1;

    /*---      Create Solution File    ---*/
    base_file = fileName;
    cgns_err = cg_open(base_file.c_str(), CG_MODE_WRITE, &cgns_file1);
    if(cgns_err) cg_error_print;
    cell_dim = 3;
    phys_dim = 3;
    cgns_err = cg_base_write(cgns_file1, basename, cell_dim, phys_dim,  &base_number);
    if(cgns_err) cg_error_print;
    cgns_err = cg_zone_write(cgns_file1, base_number, zonename, *isize, Unstructured, &zone_number);

    /*---    Write Coordinates    ---*/
    cgns_err = cg_ncoords(cgns_file, nbases, nzones, &ncoords);
    if(cgns_err) cg_error_print();
    for (i = 1; i < ncoords + 1; i++){
        cgns_err = cg_coord_info(cgns_file, nbases, nzones, i, &datatype, coordname);
        if(cgns_err) cg_error_print();
        start = 1;
        end = isize[0][0];
        vector<double> buf(end);
        cgns_err = cg_coord_read(cgns_file, nbases, nzones, coordname,
                                 datatype, &start, &end, buf.data());
        if(cgns_err) cg_error_print();
        cgns_err = cg_coord_write(cgns_file1, base_number, zone_number, datatype, coordname,
                                  buf.data(), &coord_index);
        if(cgns_err) cg_error_print();
    }

    /*---    Write Connectivity    ---*/
    cgns_err = cg_nsections(cgns_file, nbases, nzones, &nsections);
    if(cgns_err) cg_error_print();
    for ( i = 1; i < nsections + 1; i++)    {
        cgns_err = cg_section_read(cgns_file, nbases, nzones, i,
                                   ElementSectionName, &type, &start,
                                   &end, &nbndry, &parent_flag);
        if(cgns_err) cg_error_print();
        cgns_err = cg_ElementDataSize(cgns_file, nbases, nzones, i, &ElementDataSize);
        if(cgns_err) cg_error_print();
        if(type != 20){
            Elements = new cgsize_t[ElementDataSize];
            cgns_err = cg_elements_read(cgns_file, nbases, nzones, i, Elements, NULL);
            if(cgns_err) cg_error_print();
            cgns_err = cg_section_write(cgns_file1, base_number, zone_number,
                                        ElementSectionName, type, start,
                                        end, nbndry, Elements, &S);
            if(cgns_err) cg_error_print();
        }
        else{
            Elements = new cgsize_t[ElementDataSize];
            cgns_err = cg_poly_elements_read(cgns_file, nbases, nzones, i, Elements, ConnectOffset, NULL);
            if(cgns_err) cg_error_print();
            cgns_err = cg_poly_section_write(cgns_file1, base_number, zone_number,
                                             ElementSectionName, type, start,
                                             end, nbndry, Elements, ConnectOffset, &S);
            if(cgns_err) cg_error_print();
        }
        /*---    Add Boundary Condition    ---*/
        for (j = 0; j <= config->GetnMarker_All(); j++){
            if (ElementSectionName == config->GetMarker_All_TagBound(j)){
                if (config->GetMarker_All_KindBC(j) == 26){
                    bocotype = (BCType_t)20;
                    ptset_type = (PointSetType_t)PointList;
                    cgns_err = cg_boco_write(cgns_file1, base_number, zone_number, ElementSectionName,
                                             bocotype, ptset_type, ElementDataSize, Elements, &bcnum);
                    if(cgns_err) cg_error_print();
                    location = Vertex;
                    cgns_err = cg_boco_gridlocation_write(cgns_file1, base_number, zone_number, bcnum, location);
                    if(cgns_err) cg_error_print();
                }
                else if (config->GetMarker_All_KindBC(j) == 2){
                    bocotype = (BCType_t)7;
                    ptset_type = (PointSetType_t)PointList;
                    cgns_err = cg_boco_write(cgns_file1, base_number, zone_number, ElementSectionName,
                                             bocotype, ptset_type, ElementDataSize, Elements, &bcnum);
                    if(cgns_err) cg_error_print();
                    location = Vertex;
                    cgns_err = cg_boco_gridlocation_write(cgns_file1, base_number, zone_number, bcnum, location);
                    if(cgns_err) cg_error_print();
                }
                else{
                    bocotype = (BCType_t)0;
                    ptset_type = (PointSetType_t)PointList;
                    cgns_err = cg_boco_write(cgns_file1, base_number, zone_number, ElementSectionName,
                                             bocotype, ptset_type, ElementDataSize, Elements, &bcnum);
                    if(cgns_err) cg_error_print();
                    location = Vertex;
                    cgns_err = cg_boco_gridlocation_write(cgns_file1, base_number, zone_number, bcnum, location);
                    if(cgns_err) cg_error_print();
                }
            }
        }
    }

    /*---      Close CGNS Mesh File    ---*/
    cgns_err = cg_close(cgns_file);
    if(cgns_err) cg_error_print();

    /*---    Write Solution    ---*/
    cgns_err = cg_sol_write(cgns_file1, base_number, zone_number, (char *)"Solution", Vertex, &cgns_flow);
    if(cgns_err) cg_error_print();
    cgns_err = cg_goto(cgns_file1, nbases, "Zone_t", 1, "end");
    dataclass = NormalizedByUnknownDimensional;
    cgns_err = cg_dataclass_write(dataclass);
    nbases = base_number;
    nzones = zone_number;
    if(cgns_err) cg_error_print();
    for (unsigned short iField = varStart; iField < fieldNames.size(); iField++){
        const char *fieldname = fieldNames[iField].c_str();
        cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble, fieldname, GlobalData[iField], &cgns_field);
        if(cgns_err) cg_error_print();
    }
    /*---    Close CGNS Solution File    ---*/
    cgns_err = cg_close(cgns_file1);
    if(cgns_err) cg_error_print();

  #else
    cout << "CGNS file requested but SU2 was built without CGNS support. No file written" << "\n";
  #endif
  }
  #ifdef HAVE_MPI
      SU2_MPI::Barrier(MPI_COMM_WORLD);
  #endif
/*
  if(rank==MASTER_NODE){
      for (int iVar = 0; iVar < nVar; iVar++)
          delete [] GlobalData[iVar];
          delete [] GlobalData;
  }
*/
  /*--- Compute and store the write time. ---*/

  stopTime = SU2_MPI::Wtime();
  usedTime = stopTime-startTime;
  fileSize = Determine_Filesize(fileName);

  /*--- Compute and store the bandwidth ---*/

  bandwidth = fileSize/(1.0e6)/usedTime;
}
