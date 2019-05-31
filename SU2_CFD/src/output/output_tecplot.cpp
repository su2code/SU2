/*!
 * \file output_tecplot.cpp
 * \brief Main subroutines for output solver information.
 * \author F. Palacios, T. Economon, M. Colonno
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

#include "../../include/output/COutput.hpp"

void COutput::WriteTecplotASCII_Parallel(CConfig *config, CGeometry *geometry, unsigned short val_iZone, unsigned short val_nZone, bool surf_sol) {
  
  unsigned short iVar, nDim = geometry->GetnDim();
  
  unsigned long iPoint, iElem, iNode;
  unsigned long iExtIter = config->GetExtIter();
  
  if (config->GetUnsteady_Simulation() == HARMONIC_BALANCE){
    iExtIter = config->GetiInst();
  }
  
  int iProcessor;

  ofstream Tecplot_File;
  
  string filename;
  
  if (surf_sol) filename = config->GetFilename(SurfaceFilename, ".dat");
  else filename          = config->GetFilename(VolumeFilename, ".dat");
  
  /*--- Open Tecplot ASCII file and write the header. ---*/
  
  if (rank == MASTER_NODE) {
    Tecplot_File.open(filename.c_str(), ios::out);
    Tecplot_File.precision(6);
    if (surf_sol) Tecplot_File << "TITLE = \"Visualization of the surface solution\"" << endl;
    else Tecplot_File << "TITLE = \"Visualization of the volumetric solution\"" << endl;
    
    Tecplot_File << "VARIABLES = ";
    for (iVar = 0; iVar < Variable_Names.size()-1; iVar++) {
      Tecplot_File << "\"" << Variable_Names[iVar] << "\",";
    }
    Tecplot_File << "\"" << Variable_Names[Variable_Names.size()-1] << "\"" << endl;
    
    /*--- Write the header ---*/
    
    Tecplot_File << "ZONE ";
    if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
      Tecplot_File << "STRANDID="<<SU2_TYPE::Int(iExtIter+1)<<", SOLUTIONTIME="<<config->GetDelta_UnstTime()*iExtIter<<", ";
    } else if (config->GetUnsteady_Simulation() == HARMONIC_BALANCE) {
      /*--- Compute period of oscillation & compute time interval using nTimeInstances ---*/
      su2double period = config->GetHarmonicBalance_Period();
      su2double deltaT = period/(su2double)(config->GetnTimeInstances());
      Tecplot_File << "STRANDID="<<SU2_TYPE::Int(val_iZone+1)<<", SOLUTIONTIME="<<deltaT*val_iZone<<", ";
    }
    if (nDim == 2) {
      if (surf_sol) Tecplot_File << "NODES= "<< nGlobal_Surf_Poin <<", ELEMENTS= "<< nSurf_Elem_Par <<", DATAPACKING=POINT, ZONETYPE=FELINESEG"<< endl;
      else Tecplot_File << "NODES= "<< nGlobal_Poin_Par <<", ELEMENTS= "<< nGlobal_Elem_Par <<", DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL"<< endl;
    } else {
      if (surf_sol) Tecplot_File << "NODES= "<< nGlobal_Surf_Poin <<", ELEMENTS= "<< nSurf_Elem_Par <<", DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL"<< endl;
      else Tecplot_File << "NODES= "<< nGlobal_Poin_Par <<", ELEMENTS= "<< nGlobal_Elem_Par <<", DATAPACKING=POINT, ZONETYPE=FEBRICK"<< endl;
    }

    Tecplot_File.close();
    
  }
  
#ifdef HAVE_MPI
  SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif
  
  /*--- Each processor opens the file. ---*/
  
  Tecplot_File.open(filename.c_str(), ios::out | ios::app);
  
  /*--- Write surface and volumetric solution data. ---*/
  
  for (iProcessor = 0; iProcessor < size; iProcessor++) {
    if (rank == iProcessor) {
      
        /*--- Write the node data from this proc ---*/
        
        if (surf_sol) {
            for (iPoint = 0; iPoint < nSurf_Poin_Par; iPoint++) {
        for (iVar = 0; iVar < GlobalField_Counter; iVar++)
          Tecplot_File << scientific << Parallel_Surf_Data[iVar][iPoint] << "\t";
        Tecplot_File << endl;
            
          }
        } else {
          
           for (iPoint = 0; iPoint < nParallel_Poin; iPoint++) {
          for (iVar = 0; iVar < GlobalField_Counter; iVar++)
            Tecplot_File << scientific << Parallel_Data[iVar][iPoint] << "\t";
          Tecplot_File << endl;
        }
      }
    }
    Tecplot_File.flush();
#ifdef HAVE_MPI
    SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif
  }
  
  /*--- Write connectivity data. ---*/
  
  for (iProcessor = 0; iProcessor < size; iProcessor++) {
    if (rank == iProcessor) {
      
      if (surf_sol) {
        
        for (iElem = 0; iElem < nParallel_Line; iElem++) {
          iNode = iElem*N_POINTS_LINE;
          Tecplot_File << Conn_BoundLine_Par[iNode+0] << "\t";
          Tecplot_File << Conn_BoundLine_Par[iNode+1] << "\n";
        }
        
        for (iElem = 0; iElem < nParallel_BoundTria; iElem++) {
          iNode = iElem*N_POINTS_TRIANGLE;
          Tecplot_File << Conn_BoundTria_Par[iNode+0] << "\t";
          Tecplot_File << Conn_BoundTria_Par[iNode+1] << "\t";
          Tecplot_File << Conn_BoundTria_Par[iNode+2] << "\t";
          Tecplot_File << Conn_BoundTria_Par[iNode+2] << "\n";
        }
        
        for (iElem = 0; iElem < nParallel_BoundQuad; iElem++) {
          iNode = iElem*N_POINTS_QUADRILATERAL;
          Tecplot_File << Conn_BoundQuad_Par[iNode+0] << "\t";
          Tecplot_File << Conn_BoundQuad_Par[iNode+1] << "\t";
          Tecplot_File << Conn_BoundQuad_Par[iNode+2] << "\t";
          Tecplot_File << Conn_BoundQuad_Par[iNode+3] << "\n";
        }
        
      } else {
        
      for (iElem = 0; iElem < nParallel_Tria; iElem++) {
        iNode = iElem*N_POINTS_TRIANGLE;
        Tecplot_File << Conn_Tria_Par[iNode+0] << "\t";
        Tecplot_File << Conn_Tria_Par[iNode+1] << "\t";
        Tecplot_File << Conn_Tria_Par[iNode+2] << "\t";
        Tecplot_File << Conn_Tria_Par[iNode+2] << "\n";
      }
      
      for (iElem = 0; iElem < nParallel_Quad; iElem++) {
        iNode = iElem*N_POINTS_QUADRILATERAL;
        Tecplot_File << Conn_Quad_Par[iNode+0] << "\t";
        Tecplot_File << Conn_Quad_Par[iNode+1] << "\t";
        Tecplot_File << Conn_Quad_Par[iNode+2] << "\t";
        Tecplot_File << Conn_Quad_Par[iNode+3] << "\n";
      }
      
      for (iElem = 0; iElem < nParallel_Tetr; iElem++) {
        iNode = iElem*N_POINTS_TETRAHEDRON;
        Tecplot_File << Conn_Tetr_Par[iNode+0] << "\t" << Conn_Tetr_Par[iNode+1] << "\t";
        Tecplot_File << Conn_Tetr_Par[iNode+2] << "\t" << Conn_Tetr_Par[iNode+2] << "\t";
        Tecplot_File << Conn_Tetr_Par[iNode+3] << "\t" << Conn_Tetr_Par[iNode+3] << "\t";
        Tecplot_File << Conn_Tetr_Par[iNode+3] << "\t" << Conn_Tetr_Par[iNode+3] << "\n";
      }
      
      for (iElem = 0; iElem < nParallel_Hexa; iElem++) {
        iNode = iElem*N_POINTS_HEXAHEDRON;
        Tecplot_File << Conn_Hexa_Par[iNode+0] << "\t" << Conn_Hexa_Par[iNode+1] << "\t";
        Tecplot_File << Conn_Hexa_Par[iNode+2] << "\t" << Conn_Hexa_Par[iNode+3] << "\t";
        Tecplot_File << Conn_Hexa_Par[iNode+4] << "\t" << Conn_Hexa_Par[iNode+5] << "\t";
        Tecplot_File << Conn_Hexa_Par[iNode+6] << "\t" << Conn_Hexa_Par[iNode+7] << "\n";
      }
      
      for (iElem = 0; iElem < nParallel_Pris; iElem++) {
        iNode = iElem*N_POINTS_PRISM;
        Tecplot_File << Conn_Pris_Par[iNode+0] << "\t" << Conn_Pris_Par[iNode+1] << "\t";
        Tecplot_File << Conn_Pris_Par[iNode+1] << "\t" << Conn_Pris_Par[iNode+2] << "\t";
        Tecplot_File << Conn_Pris_Par[iNode+3] << "\t" << Conn_Pris_Par[iNode+4] << "\t";
        Tecplot_File << Conn_Pris_Par[iNode+4] << "\t" << Conn_Pris_Par[iNode+5] << "\n";
      }
      
      for (iElem = 0; iElem < nParallel_Pyra; iElem++) {
        iNode = iElem*N_POINTS_PYRAMID;
        Tecplot_File << Conn_Pyra_Par[iNode+0] << "\t" << Conn_Pyra_Par[iNode+1] << "\t";
        Tecplot_File << Conn_Pyra_Par[iNode+2] << "\t" << Conn_Pyra_Par[iNode+3] << "\t";
        Tecplot_File << Conn_Pyra_Par[iNode+4] << "\t" << Conn_Pyra_Par[iNode+4] << "\t";
        Tecplot_File << Conn_Pyra_Par[iNode+4] << "\t" << Conn_Pyra_Par[iNode+4] << "\n";
      }
      }
      
    }
    Tecplot_File.flush();
#ifdef HAVE_MPI
    SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif
  }
  
  Tecplot_File.close();
  
}

#ifdef HAVE_MPI

namespace
{

/*!
 * \brief Calculate the partitioning of nodes to determine:
 * (a) For a given global node number, to which partition does it belong and what is its local node number; and
 * (b) How many nodes are held by each partition.
 */
class NodePartitioner {
public:
  /*!
   * \param[in] global_num_nodes - The total number of nodes being output
   * \param[in] num_ranks - The number of MPI ranks involved in the output
   */
  NodePartitioner(unsigned long global_num_nodes, int num_ranks)
    : m_num_ranks(num_ranks) {
    /* rank i has (1-based) global nodes m_node_range[i] + 1 through m_node_range[i + 1] */
    unsigned long nodes_per_rank = global_num_nodes / num_ranks;
    unsigned long num_extra_nodes = global_num_nodes - nodes_per_rank * num_ranks;
    m_node_range.resize(num_ranks + 1);
    m_node_range[0] = 0;
    for (int ii = 1; ii <= num_ranks; ii++) {
      m_node_range[ii] = m_node_range[ii - 1] + nodes_per_rank;
      if (num_extra_nodes > 0) {
        ++m_node_range[ii];
        --num_extra_nodes;
      }
    }
    assert(m_node_range[num_ranks] == global_num_nodes);
  }

  /*!
   * \brief Determine the MPI rank that owns a global node number and its corresponding local node number.
   * \param global_node_number[in] - The global node number; global node numbers are sequential across all MPI ranks.
   * \param owning_rank[out] - The MPI rank that owns (will output) the global node
   * \param node_number[out] - The rank-local node number for the given global node number
   */
  void GetOwningRankAndNodeNumber(unsigned long global_node_number, int &owning_rank, unsigned long &node_number)
  {
    owning_rank = static_cast<int>(global_node_number / m_node_range[1]);
    if (owning_rank >= m_num_ranks)
      owning_rank = m_num_ranks - 1;
    while(global_node_number > m_node_range[owning_rank + 1])
      ++owning_rank;
    while(global_node_number <= m_node_range[owning_rank])
      --owning_rank;
    node_number = global_node_number - m_node_range[owning_rank];
  }

  /*!
   * \brief Determine the number of nodes to be output by a particular rank
   * \param which_rank[in] - The MPI rank
   * \ret - The number of nodes that will be output by the give MPI rank.
   */
  int64_t GetRankNumNodes(int which_rank)
  {
    return static_cast<int64_t>(m_node_range[which_rank + 1] - m_node_range[which_rank]);
  }

private:
  int m_num_ranks;
  vector<unsigned long> m_node_range;
};

int64_t GetHaloNodeNumber(unsigned long global_node_number, unsigned long last_local_node, vector<unsigned long> const &halo_node_list)
{
  vector<unsigned long>::const_iterator it = lower_bound(halo_node_list.begin(), halo_node_list.end(), global_node_number);
  assert(it != halo_node_list.end());
  assert(*it == global_node_number);
  /* When C++11 is universally available, replace the following mouthful with "auto" */
  iterator_traits<vector<unsigned long>::const_iterator>::difference_type offset = distance(halo_node_list.begin(), it);
  assert(offset >= 0);
  return (int64_t)(last_local_node + offset + 1);
}

} /* namespace */

#endif /* HAVE_MPI */

void COutput::WriteTecplotBinary_Parallel(CConfig *config, CGeometry *geometry, unsigned short val_iZone, unsigned short val_nZone, bool surf_sol) {

#ifdef HAVE_TECIO
  
  /*--- Open Tecplot binary file. ---*/
  
  string filename;
  
  if (surf_sol) filename = config->GetFilename(SurfaceFilename, ".szplt");
  else filename          = config->GetFilename(VolumeFilename, ".szplt");

  
  string data_set_title = surf_sol
    ? "Visualization of the surface solution"
    : "Visualization of the volumetric solution";

  ostringstream tecplot_variable_names;
  for (size_t iVar = 0; iVar < Variable_Names.size()-1; ++iVar) {
    tecplot_variable_names << Variable_Names[iVar] << ",";
  }
  tecplot_variable_names << Variable_Names[Variable_Names.size()-1];

  void* file_handle = NULL;
  int32_t err = tecFileWriterOpen(filename.c_str(), data_set_title.c_str(), tecplot_variable_names.str().c_str(),
    FILEFORMAT_SZL, FILETYPE_FULL, (int32_t)FieldDataType_Double, NULL, &file_handle);
  if (err) cout << "Error opening Tecplot file '" << filename << "'" << endl;

#ifdef HAVE_MPI
  err = tecMPIInitialize(file_handle, MPI_COMM_WORLD, MASTER_NODE);
  if (err) cout << "Error initializing Tecplot parallel output." << endl;
#endif
  
  /*--- Define the zone(s). For 2D, and for 3D surfaces, each rank outputs a separate zone. ---*/

  int64_t num_nodes;
  int64_t num_cells;
  int32_t zone_type;
  if (surf_sol) {
    num_nodes = static_cast<int64_t>(nGlobal_Surf_Poin);
    num_cells = static_cast<int64_t>(nSurf_Elem_Par);
    if (geometry->GetnDim() == 2) 
      zone_type = ZONETYPE_FELINESEG;
    else
      zone_type = ZONETYPE_FEQUADRILATERAL;
  } else {
    num_nodes = static_cast<int64_t>(nGlobal_Poin_Par);
    num_cells = static_cast<int64_t>(nGlobal_Elem_Par);
    if (geometry->GetnDim() == 2)
      zone_type = ZONETYPE_FEQUADRILATERAL;
    else
      zone_type = ZONETYPE_FEBRICK;
  }

  bool is_unsteady = false;
  passivedouble solution_time = 0.0;
  if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
    is_unsteady = true;
    solution_time = SU2_TYPE::GetValue(config->GetDelta_UnstTime()*config->GetExtIter());
  } else if (config->GetUnsteady_Simulation() == HARMONIC_BALANCE) {
    is_unsteady = true;
    /*--- Compute period of oscillation & compute time interval using nTimeInstances ---*/
    passivedouble period = SU2_TYPE::GetValue(config->GetHarmonicBalance_Period());
    passivedouble deltaT = period/SU2_TYPE::GetValue(config->GetnTimeInstances());
    solution_time = deltaT*val_iZone;
  }

  int32_t zone;
  vector<int32_t> value_locations(GlobalField_Counter, 1); /* Nodal variables. */
  err = tecZoneCreateFE(file_handle, "Zone", zone_type, num_nodes, num_cells, NULL, NULL, &value_locations[0], NULL, 0, 0, 0, &zone);
  if (err) cout << rank << ": Error creating Tecplot zone." << endl;
  if (is_unsteady) {
    err = tecZoneSetUnsteadyOptions(file_handle, zone, solution_time, config->GetExtIter() + 1);
    if (err) cout << rank << ": Error setting Tecplot zone unsteady options." << std::endl;
  }

#ifdef HAVE_MPI

  unsigned short iVar;
  NodePartitioner node_partitioner(num_nodes, size);
  set<unsigned long> halo_nodes;
  vector<unsigned long> sorted_halo_nodes;
  vector<passivedouble> halo_var_data;
  vector<int> num_nodes_to_receive(size, 0);
  vector<int> values_to_receive_displacements(size);

  if (zone_type == ZONETYPE_FEBRICK) {

    /* We output a single, partitioned zone where each rank outputs one partition. */
    vector<int32_t> partition_owners;
    partition_owners.reserve(size);
    for (int32_t iRank = 0; iRank < size; ++iRank)
      partition_owners.push_back(iRank);
    err = tecZoneMapPartitionsToMPIRanks(file_handle, zone, size, &partition_owners[0]);
    if (err) cout << rank << ": Error assigning MPI ranks for Tecplot zone partitions." << endl;
  
    /* Gather a list of nodes we refer to but are not outputting. */

    for (unsigned long i = 0; i < nParallel_Tria * N_POINTS_TRIANGLE; ++i)
      if ((unsigned long)Conn_Tria_Par[i] <= beg_node[rank] || end_node[rank] < (unsigned long)Conn_Tria_Par[i])
        halo_nodes.insert(Conn_Tria_Par[i]);
  
    for (unsigned long i = 0; i < nParallel_Quad * N_POINTS_QUADRILATERAL; ++i)
      if ((unsigned long)Conn_Quad_Par[i] <= beg_node[rank] || end_node[rank] < (unsigned long)Conn_Quad_Par[i])
        halo_nodes.insert(Conn_Quad_Par[i]);
 
    for (unsigned long i = 0; i < nParallel_Tetr * N_POINTS_TETRAHEDRON; ++i)
      if ((unsigned long)Conn_Tetr_Par[i] <= beg_node[rank] || end_node[rank] < (unsigned long)Conn_Tetr_Par[i])
        halo_nodes.insert(Conn_Tetr_Par[i]);

    for (unsigned long i = 0; i < nParallel_Hexa * N_POINTS_HEXAHEDRON; ++i)
      if ((unsigned long)Conn_Hexa_Par[i] <= beg_node[rank] || end_node[rank] < (unsigned long)Conn_Hexa_Par[i])
        halo_nodes.insert(Conn_Hexa_Par[i]);
      
    for (unsigned long i = 0; i < nParallel_Pris * N_POINTS_PRISM; ++i)
      if ((unsigned long)Conn_Pris_Par[i] <= beg_node[rank] || end_node[rank] < (unsigned long)Conn_Pris_Par[i])
        halo_nodes.insert(Conn_Pris_Par[i]);
    
    for (unsigned long i = 0; i < nParallel_Pyra * N_POINTS_PYRAMID; ++i)
      if ((unsigned long)Conn_Pyra_Par[i] <= beg_node[rank] || end_node[rank] < (unsigned long)Conn_Pyra_Par[i])
        halo_nodes.insert(Conn_Pyra_Par[i]);

    /* Sorted list of halo nodes for this MPI rank. */
    sorted_halo_nodes.assign(halo_nodes.begin(), halo_nodes.end());
        
    /* Have to include all nodes our cells refer to or TecIO will barf, so add the halo node count to the number of local nodes. */
    int64_t partition_num_nodes = end_node[rank] - beg_node[rank] + static_cast<int64_t>(halo_nodes.size());
    int64_t partition_num_cells = nParallel_Tetr + nParallel_Hexa + nParallel_Pris + nParallel_Pyra;

    /*--- We effectively tack the halo nodes onto the end of the node list for this partition.
      TecIO will later replace them with references to nodes in neighboring partitions. */
    size_t num_halo_nodes = sorted_halo_nodes.size();
    vector<int64_t> halo_node_local_numbers(max((size_t)1, num_halo_nodes)); /* Min size 1 to avoid crashes when we access these vectors below. */
    vector<int32_t> neighbor_partitions(max((size_t)1, num_halo_nodes));
    vector<int64_t> neighbor_nodes(max((size_t)1, num_halo_nodes));
    for(int64_t i = 0; i < static_cast<int64_t>(num_halo_nodes); ++i) {
      halo_node_local_numbers[i] = end_node[rank] - beg_node[rank] + i + 1;
      int owning_rank;
      unsigned long node_number;
      node_partitioner.GetOwningRankAndNodeNumber(sorted_halo_nodes[i], owning_rank, node_number);
      neighbor_partitions[i] = owning_rank + 1; /* Partition numbers are 1-based. */
      neighbor_nodes[i] = static_cast<int64_t>(node_number);
    }
    err = tecFEPartitionCreate64(file_handle, zone, rank + 1, partition_num_nodes, partition_num_cells,
      static_cast<int64_t>(num_halo_nodes), &halo_node_local_numbers[0], &neighbor_partitions[0], &neighbor_nodes[0], 0, NULL);
    if (err) cout << rank << ": Error creating Tecplot zone partition." << endl;

    /* Gather halo node data. First, tell each rank how many nodes' worth of data we need from them. */
    for (size_t i = 0; i < num_halo_nodes; ++i)
      ++num_nodes_to_receive[neighbor_partitions[i] - 1];
    vector<int> num_nodes_to_send(size);
    SU2_MPI::Alltoall(&num_nodes_to_receive[0], 1, MPI_INT, &num_nodes_to_send[0], 1, MPI_INT, MPI_COMM_WORLD);

    /* Now send the global node numbers whose data we need,
       and receive the same from all other ranks.
       Each rank has globally consecutive node numbers,
       so we can just parcel out sorted_halo_nodes for send. */
    vector<int> nodes_to_send_displacements(size);
    vector<int> nodes_to_receive_displacements(size);
    nodes_to_send_displacements[0] = 0;
    nodes_to_receive_displacements[0] = 0;
    for(int iRank = 1; iRank < size; ++iRank) {
      nodes_to_send_displacements[iRank] = nodes_to_send_displacements[iRank - 1] + num_nodes_to_send[iRank - 1];
      nodes_to_receive_displacements[iRank] = nodes_to_receive_displacements[iRank - 1] + num_nodes_to_receive[iRank - 1];
    }
    int total_num_nodes_to_send = nodes_to_send_displacements[size - 1] + num_nodes_to_send[size - 1];
    vector<unsigned long> nodes_to_send(max(1, total_num_nodes_to_send));

    /* The terminology gets a bit confusing here. We're sending the node numbers
       (sorted_halo_nodes) whose data we need to receive, and receiving
       lists of nodes whose data we need to send. */
    if (sorted_halo_nodes.empty()) sorted_halo_nodes.resize(1); /* Avoid crash. */
    SU2_MPI::Alltoallv(&sorted_halo_nodes[0], &num_nodes_to_receive[0], &nodes_to_receive_displacements[0], MPI_UNSIGNED_LONG,
                       &nodes_to_send[0],     &num_nodes_to_send[0],    &nodes_to_send_displacements[0],    MPI_UNSIGNED_LONG,
                       MPI_COMM_WORLD);
    
    /* Now actually send and receive the data */
    vector<passivedouble> data_to_send(max(1, total_num_nodes_to_send * GlobalField_Counter));
    halo_var_data.resize(max((size_t)1, GlobalField_Counter * num_halo_nodes));
    vector<int> num_values_to_send(size);
    vector<int> values_to_send_displacements(size);
    vector<int> num_values_to_receive(size);
    size_t index = 0;
    for(int iRank = 0; iRank < size; ++iRank) {
      /* We send and receive GlobalField_Counter values per node. */
      num_values_to_send[iRank]              = num_nodes_to_send[iRank] * GlobalField_Counter;
      values_to_send_displacements[iRank]    = nodes_to_send_displacements[iRank] * GlobalField_Counter;
      num_values_to_receive[iRank]           = num_nodes_to_receive[iRank] * GlobalField_Counter;
      values_to_receive_displacements[iRank] = nodes_to_receive_displacements[iRank] * GlobalField_Counter;
      for(iVar = 0; iVar < GlobalField_Counter; ++iVar)
        for(int iNode = 0; iNode < num_nodes_to_send[iRank]; ++iNode) {
          unsigned long node_offset = nodes_to_send[nodes_to_send_displacements[iRank] + iNode] - beg_node[rank] - 1;
          data_to_send[index++] = SU2_TYPE::GetValue(Parallel_Data[iVar][node_offset]);
        }
    }
    SU2_MPI::Alltoallv(&data_to_send[0],  &num_values_to_send[0],    &values_to_send_displacements[0],    MPI_DOUBLE,
                       &halo_var_data[0], &num_values_to_receive[0], &values_to_receive_displacements[0], MPI_DOUBLE,
                       MPI_COMM_WORLD);
  }
  else {
    /* Zone will be gathered to and output by MASTER_NODE */
    int32_t partition_owner = MASTER_NODE;
    err = tecZoneMapPartitionsToMPIRanks(file_handle, zone, 1, &partition_owner);
  }

  /*--- Write surface and volumetric solution data. ---*/
  
  if (zone_type == ZONETYPE_FEBRICK) {
    std::vector<passivedouble> values_to_write(nParallel_Poin);
    for (iVar = 0; err == 0 && iVar < GlobalField_Counter; iVar++) {
      for(unsigned long i = 0; i < nParallel_Poin; ++i)
        values_to_write[i] = SU2_TYPE::GetValue(Parallel_Data[iVar][i]);
      err = tecZoneVarWriteDoubleValues(file_handle, zone, iVar + 1, rank + 1, nParallel_Poin, &values_to_write[0]);
      if (err) cout << rank << ": Error outputting Tecplot variable values." << endl;
      for (int iRank = 0; err == 0 && iRank < size; ++iRank) {
        if (num_nodes_to_receive[iRank] > 0) {
          int var_data_offset = values_to_receive_displacements[iRank] + num_nodes_to_receive[iRank] * iVar;
          err = tecZoneVarWriteDoubleValues(file_handle, zone, iVar + 1, rank + 1, static_cast<int64_t>(num_nodes_to_receive[iRank]), &halo_var_data[var_data_offset]);
          if (err) cout << rank << ": Error outputting Tecplot halo values." << endl;
        }
      }
    }
  } else {
    if (rank == MASTER_NODE) {
      vector<passivedouble> var_data;
      vector<unsigned long> num_surface_points(size);
      if (surf_sol)
        SU2_MPI::Gather(&nSurf_Poin_Par, 1, MPI_UNSIGNED_LONG, &num_surface_points[0], 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
      for(int iRank = 0; iRank < size; ++iRank) {
        int64_t rank_num_points;
        if (surf_sol)
          rank_num_points = num_surface_points[iRank];
        else
          rank_num_points = node_partitioner.GetRankNumNodes(iRank);
        if (rank_num_points > 0) {
          if (iRank == rank) { /* Output local data. */
            std::vector<passivedouble> values_to_write;
            for (iVar = 0; err == 0 && iVar < GlobalField_Counter; iVar++) {
              if (surf_sol) {
                values_to_write.resize(nSurf_Poin_Par);
                for(unsigned long i = 0; i < nSurf_Poin_Par; ++i)
                  values_to_write[i] = SU2_TYPE::GetValue(Parallel_Surf_Data[iVar][i]);
                err = tecZoneVarWriteDoubleValues(file_handle, zone, iVar + 1, 0, nSurf_Poin_Par, &values_to_write[0]);
              }
              else {
                values_to_write.resize(rank_num_points);
                for(unsigned long i = 0; i < (unsigned long)rank_num_points; ++i)
                  values_to_write[i] = SU2_TYPE::GetValue(Parallel_Data[iVar][i]);
                err = tecZoneVarWriteDoubleValues(file_handle, zone, iVar + 1, 0, rank_num_points, &values_to_write[0]);
              }
              if (err) cout << rank << ": Error outputting Tecplot variable values." << endl;
            }
          }
          else { /* Receive data from other rank. */
            var_data.resize(max((int64_t)1, GlobalField_Counter * rank_num_points));
            SU2_MPI::Recv(&var_data[0], GlobalField_Counter * rank_num_points, MPI_DOUBLE, iRank, iRank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (iVar = 0; err == 0 && iVar < GlobalField_Counter; iVar++) {
              err = tecZoneVarWriteDoubleValues(file_handle, zone, iVar + 1, 0, rank_num_points, &var_data[iVar * rank_num_points]);
              if (err) cout << rank << ": Error outputting Tecplot surface variable values." << endl;
            }
          }
        }
      }
    }
    else { /* Send data to MASTER_NODE */
      if (surf_sol)
        SU2_MPI::Gather(&nSurf_Poin_Par, 1, MPI_UNSIGNED_LONG, NULL, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
      vector<passivedouble> var_data;
      size_t var_data_size = GlobalField_Counter * (surf_sol ? nSurf_Poin_Par : nParallel_Poin);
      var_data.reserve(var_data_size);
      for (iVar = 0; err == 0 && iVar < GlobalField_Counter; iVar++)
        if (surf_sol)
          for(unsigned long i = 0; i < nSurf_Poin_Par; ++i)
            var_data.push_back(SU2_TYPE::GetValue(Parallel_Surf_Data[iVar][i]));
        else
          for(unsigned long i = 0; i < nParallel_Poin; ++i)
            var_data.push_back(SU2_TYPE::GetValue(Parallel_Data[iVar][i]));
      if (var_data.size() > 0)
        SU2_MPI::Send(&var_data[0], static_cast<int>(var_data.size()), MPI_DOUBLE, MASTER_NODE, rank, MPI_COMM_WORLD);
    }
  }

#else

  unsigned short iVar;

  vector<passivedouble> var_data;
  size_t var_data_size = nVar_Par * (surf_sol ? nSurf_Poin_Par : nParallel_Poin);
  var_data.reserve(var_data_size);
  
  if (surf_sol) {
    for (iVar = 0; err == 0 && iVar < GlobalField_Counter; iVar++) {
      for(unsigned long i = 0; i < nSurf_Poin_Par; ++i)
        var_data.push_back(SU2_TYPE::GetValue(Parallel_Surf_Data[iVar][i]));
      err = tecZoneVarWriteDoubleValues(file_handle, zone, iVar + 1, 0, nSurf_Poin_Par, &var_data[iVar * nSurf_Poin_Par]);
      if (err) cout << rank << ": Error outputting Tecplot variable value." << endl;
    }
  } else {
    for (iVar = 0; err == 0 && iVar < GlobalField_Counter; iVar++) {
      for(unsigned long i = 0; i < nParallel_Poin; ++i)
        var_data.push_back(SU2_TYPE::GetValue(Parallel_Data[iVar][i]));
      err = tecZoneVarWriteDoubleValues(file_handle, zone, iVar + 1, 0, nParallel_Poin, &var_data[iVar * nParallel_Poin]);
      if (err) cout << rank << ": Error outputting Tecplot variable value." << endl;
    }
  }

#endif /* HAVE_MPI */
  
  /*--- Write connectivity data. ---*/

  unsigned long iElem, iNode;
  
#ifdef HAVE_MPI
  if (zone_type == ZONETYPE_FEBRICK) {

    int64_t nodes[8];

    /**
     *  Each rank writes node numbers relative to the partition it is outputting (starting with node number 1).
     *  Ghost (halo) nodes identified above are numbered sequentially just beyond the end of the actual, local nodes.
     *  Note that beg_node and end_node refer to zero-based node numbering, but Conn_* contain one-based node numbers.
     */
#define MAKE_LOCAL(n) beg_node[rank] < (unsigned long)n && (unsigned long)n <= end_node[rank] \
  ? (int64_t)((unsigned long)n - beg_node[rank]) \
  : GetHaloNodeNumber(n, end_node[rank] - beg_node[rank], sorted_halo_nodes)

    for (iElem = 0; err == 0 && iElem < nParallel_Tetr; iElem++) {
      iNode = iElem*N_POINTS_TETRAHEDRON;
      nodes[0] = MAKE_LOCAL(Conn_Tetr_Par[iNode+0]);
      nodes[1] = MAKE_LOCAL(Conn_Tetr_Par[iNode+1]);
      nodes[2] = MAKE_LOCAL(Conn_Tetr_Par[iNode+2]);
      nodes[3] = nodes[2];
      nodes[4] = MAKE_LOCAL(Conn_Tetr_Par[iNode+3]);
      nodes[5] = nodes[4];
      nodes[6] = nodes[4];
      nodes[7] = nodes[4];
      err = tecZoneNodeMapWrite64(file_handle, zone, rank + 1, 1, 8, nodes);
      if (err) cout << rank << ": Error outputting Tecplot node values." << endl;
    }

    for (iElem = 0; err == 0 && iElem < nParallel_Hexa; iElem++) {
      iNode = iElem*N_POINTS_HEXAHEDRON;
      nodes[0] = MAKE_LOCAL(Conn_Hexa_Par[iNode+0]);
      nodes[1] = MAKE_LOCAL(Conn_Hexa_Par[iNode+1]);
      nodes[2] = MAKE_LOCAL(Conn_Hexa_Par[iNode+2]);
      nodes[3] = MAKE_LOCAL(Conn_Hexa_Par[iNode+3]);
      nodes[4] = MAKE_LOCAL(Conn_Hexa_Par[iNode+4]);
      nodes[5] = MAKE_LOCAL(Conn_Hexa_Par[iNode+5]);
      nodes[6] = MAKE_LOCAL(Conn_Hexa_Par[iNode+6]);
      nodes[7] = MAKE_LOCAL(Conn_Hexa_Par[iNode+7]);
      err = tecZoneNodeMapWrite64(file_handle, zone, rank + 1, 1, 8, nodes);
      if (err) cout << rank << ": Error outputting Tecplot node values." << endl;
    }
      
    for (iElem = 0; err == 0 && iElem < nParallel_Pris; iElem++) {
      iNode = iElem*N_POINTS_PRISM;
      nodes[0] = MAKE_LOCAL(Conn_Pris_Par[iNode+0]);
      nodes[1] = MAKE_LOCAL(Conn_Pris_Par[iNode+1]);
      nodes[2] = nodes[1];
      nodes[3] = MAKE_LOCAL(Conn_Pris_Par[iNode+2]);
      nodes[4] = MAKE_LOCAL(Conn_Pris_Par[iNode+3]);
      nodes[5] = MAKE_LOCAL(Conn_Pris_Par[iNode+4]);
      nodes[6] = nodes[5];
      nodes[7] = MAKE_LOCAL(Conn_Pris_Par[iNode+5]);
      err = tecZoneNodeMapWrite64(file_handle, zone, rank + 1, 1, 8, nodes);
      if (err) cout << rank << ": Error outputting Tecplot node values." << endl;
    }
    
    for (iElem = 0; err == 0 && iElem < nParallel_Pyra; iElem++) {
      iNode = iElem*N_POINTS_PYRAMID;
      nodes[0] = MAKE_LOCAL(Conn_Pyra_Par[iNode+0]);
      nodes[1] = MAKE_LOCAL(Conn_Pyra_Par[iNode+1]);
      nodes[2] = MAKE_LOCAL(Conn_Pyra_Par[iNode+2]);
      nodes[3] = MAKE_LOCAL(Conn_Pyra_Par[iNode+3]);
      nodes[4] = MAKE_LOCAL(Conn_Pyra_Par[iNode+4]);
      nodes[5] = nodes[4];
      nodes[6] = nodes[4];
      nodes[7] = nodes[4];
      err = tecZoneNodeMapWrite64(file_handle, zone, rank + 1, 1, 8, nodes);
      if (err) cout << rank << ": Error outputting Tecplot node values." << endl;
    }
  } else {
    if (rank == MASTER_NODE) {

      /* Non-hexahedral output by the master node. Output local data directly, and gather other data from the other ranks. */

      int64_t nodes[4];

      vector<unsigned long> connectivity_sizes(size);
      unsigned long unused = 0;
      SU2_MPI::Gather(&unused, 1, MPI_UNSIGNED_LONG, &connectivity_sizes[0], 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
      vector<int64_t> connectivity;
      for(int iRank = 0; iRank < size; ++iRank) {
        if (iRank == rank) {
          if (surf_sol) {
            for (iElem = 0; err == 0 && iElem < nParallel_Line; iElem++) {
              iNode = iElem*N_POINTS_LINE;
              nodes[0] = Conn_BoundLine_Par[iNode+0];
              nodes[1] = Conn_BoundLine_Par[iNode+1];
              err = tecZoneNodeMapWrite64(file_handle, zone, 0, 1, 2, nodes);
              if (err) cout << rank << ": Error outputting Tecplot node values." << endl;
            }
          
            for (iElem = 0; err == 0 && iElem < nParallel_BoundTria; iElem++) {
              iNode = iElem*N_POINTS_TRIANGLE;
              nodes[0] = Conn_BoundTria_Par[iNode+0];
              nodes[1] = Conn_BoundTria_Par[iNode+1];
              nodes[2] = Conn_BoundTria_Par[iNode+2];
              nodes[3] = Conn_BoundTria_Par[iNode+2];
              err = tecZoneNodeMapWrite64(file_handle, zone, 0, 1, 4, nodes);
              if (err) cout << rank << ": Error outputting Tecplot node values." << endl;
            }
          
            for (iElem = 0; err == 0 && iElem < nParallel_BoundQuad; iElem++) {
              iNode = iElem*N_POINTS_QUADRILATERAL;
              nodes[0] = Conn_BoundQuad_Par[iNode+0];
              nodes[1] = Conn_BoundQuad_Par[iNode+1];
              nodes[2] = Conn_BoundQuad_Par[iNode+2];
              nodes[3] = Conn_BoundQuad_Par[iNode+3];
              err = tecZoneNodeMapWrite64(file_handle, zone, 0, 1, 4, nodes);
              if (err) cout << rank << ": Error outputting Tecplot node values." << endl;
            }
          } else {
            for (iElem = 0; err == 0 && iElem < nParallel_Tria; iElem++) {
              iNode = iElem*N_POINTS_TRIANGLE;
              nodes[0] = Conn_Tria_Par[iNode+0];
              nodes[1] = Conn_Tria_Par[iNode+1];
              nodes[2] = Conn_Tria_Par[iNode+2];
              nodes[3] = Conn_Tria_Par[iNode+2];
              err = tecZoneNodeMapWrite64(file_handle, zone, 0, 1, 4, nodes);
              if (err) cout << rank << ": Error outputting Tecplot node values." << endl;
            }
        
            for (iElem = 0; err == 0 && iElem < nParallel_Quad; iElem++) {
              iNode = iElem*N_POINTS_QUADRILATERAL;
              nodes[0] = Conn_Quad_Par[iNode+0];
              nodes[1] = Conn_Quad_Par[iNode+1];
              nodes[2] = Conn_Quad_Par[iNode+2];
              nodes[3] = Conn_Quad_Par[iNode+3];
              err = tecZoneNodeMapWrite64(file_handle, zone, 0, 1, 4, nodes);
              if (err) cout << rank << ": Error outputting Tecplot node values." << endl;
            }
          }
        } else { /* Receive node map and write out. */
          connectivity.resize(max((unsigned long)1, connectivity_sizes[iRank]));
          SU2_MPI::Recv(&connectivity[0], connectivity_sizes[iRank], MPI_UNSIGNED_LONG, iRank, iRank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          err = tecZoneNodeMapWrite64(file_handle, zone, 0, 1, connectivity_sizes[iRank], &connectivity[0]);
          if (err) cout << rank << ": Error outputting Tecplot node values." << endl;
        }
      }
    } else {

      /* Non-hexahedral output by non-master node. Send what we've got to the master node. */

      unsigned long connectivity_size;
      if (surf_sol)
        connectivity_size = 2 * nParallel_Line + 4 * nParallel_BoundTria + 4 * nParallel_BoundQuad;
      else
        connectivity_size = 4 * (nParallel_Tria + nParallel_Quad);
      SU2_MPI::Gather(&connectivity_size, 1, MPI_UNSIGNED_LONG, NULL, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
      vector<int64_t> connectivity;
      connectivity.reserve(connectivity_size);
      if (surf_sol) {
        for (iElem = 0; err == 0 && iElem < nParallel_Line; iElem++) {
          iNode = iElem*N_POINTS_LINE;
          connectivity.push_back(Conn_BoundLine_Par[iNode+0]);
          connectivity.push_back(Conn_BoundLine_Par[iNode+1]);
        }
      
        for (iElem = 0; err == 0 && iElem < nParallel_BoundTria; iElem++) {
          iNode = iElem*N_POINTS_TRIANGLE;
          connectivity.push_back(Conn_BoundTria_Par[iNode+0]);
          connectivity.push_back(Conn_BoundTria_Par[iNode+1]);
          connectivity.push_back(Conn_BoundTria_Par[iNode+2]);
          connectivity.push_back(Conn_BoundTria_Par[iNode+2]);
        }
      
        for (iElem = 0; err == 0 && iElem < nParallel_BoundQuad; iElem++) {
          iNode = iElem*N_POINTS_QUADRILATERAL;
          connectivity.push_back(Conn_BoundQuad_Par[iNode+0]);
          connectivity.push_back(Conn_BoundQuad_Par[iNode+1]);
          connectivity.push_back(Conn_BoundQuad_Par[iNode+2]);
          connectivity.push_back(Conn_BoundQuad_Par[iNode+3]);
        }
      } else {
        for (iElem = 0; err == 0 && iElem < nParallel_Tria; iElem++) {
          iNode = iElem*N_POINTS_TRIANGLE;
          connectivity.push_back(Conn_Tria_Par[iNode+0]);
          connectivity.push_back(Conn_Tria_Par[iNode+1]);
          connectivity.push_back(Conn_Tria_Par[iNode+2]);
          connectivity.push_back(Conn_Tria_Par[iNode+2]);
        }
    
        for (iElem = 0; err == 0 && iElem < nParallel_Quad; iElem++) {
          iNode = iElem*N_POINTS_QUADRILATERAL;
          connectivity.push_back(Conn_Quad_Par[iNode+0]);
          connectivity.push_back(Conn_Quad_Par[iNode+1]);
          connectivity.push_back(Conn_Quad_Par[iNode+2]);
          connectivity.push_back(Conn_Quad_Par[iNode+3]);
        }
      }
      if (connectivity.empty()) connectivity.resize(1); /* Avoid crash */
      SU2_MPI::Send(&connectivity[0], connectivity_size, MPI_UNSIGNED_LONG, MASTER_NODE, rank, MPI_COMM_WORLD);
    }
  }
#else
  if (surf_sol) {

    int64_t nodes[4];

    for (iElem = 0; err == 0 && iElem < nParallel_Line; iElem++) {
      iNode = iElem*N_POINTS_LINE;
      nodes[0] = Conn_BoundLine_Par[iNode+0];
      nodes[1] = Conn_BoundLine_Par[iNode+1];
      err = tecZoneNodeMapWrite64(file_handle, zone, rank, 1, 2, nodes);
      if (err) cout << rank << ": Error outputting Tecplot node values." << endl;
    }
        
    for (iElem = 0; err == 0 && iElem < nParallel_BoundTria; iElem++) {
      iNode = iElem*N_POINTS_TRIANGLE;
      nodes[0] = Conn_BoundTria_Par[iNode+0];
      nodes[1] = Conn_BoundTria_Par[iNode+1];
      nodes[2] = Conn_BoundTria_Par[iNode+2];
      nodes[3] = Conn_BoundTria_Par[iNode+2];
      err = tecZoneNodeMapWrite64(file_handle, zone, rank, 1, 4, nodes);
      if (err) cout << rank << ": Error outputting Tecplot node values." << endl;
    }
        
    for (iElem = 0; err == 0 && iElem < nParallel_BoundQuad; iElem++) {
      iNode = iElem*N_POINTS_QUADRILATERAL;
      nodes[0] = Conn_BoundQuad_Par[iNode+0];
      nodes[1] = Conn_BoundQuad_Par[iNode+1];
      nodes[2] = Conn_BoundQuad_Par[iNode+2];
      nodes[3] = Conn_BoundQuad_Par[iNode+3];
      err = tecZoneNodeMapWrite64(file_handle, zone, rank, 1, 4, nodes);
      if (err) cout << rank << ": Error outputting Tecplot node values." << endl;
    }

  } else {

    int64_t nodes[8];

    for (iElem = 0; err == 0 && iElem < nParallel_Tria; iElem++) {
      iNode = iElem*N_POINTS_TRIANGLE;
      nodes[0] = Conn_Tria_Par[iNode+0];
      nodes[1] = Conn_Tria_Par[iNode+1];
      nodes[2] = Conn_Tria_Par[iNode+2];
      nodes[3] = Conn_Tria_Par[iNode+2];
      err = tecZoneNodeMapWrite64(file_handle, zone, rank, 1, 4, nodes);
      if (err) cout << rank << ": Error outputting Tecplot node values." << endl;
    }

    for (iElem = 0; err == 0 && iElem < nParallel_Quad; iElem++) {
      iNode = iElem*N_POINTS_QUADRILATERAL;
      nodes[0] = Conn_Quad_Par[iNode+0];
      nodes[1] = Conn_Quad_Par[iNode+1];
      nodes[2] = Conn_Quad_Par[iNode+2];
      nodes[3] = Conn_Quad_Par[iNode+3];
      err = tecZoneNodeMapWrite64(file_handle, zone, rank, 1, 4, nodes);
      if (err) cout << rank << ": Error outputting Tecplot node values." << endl;
    }

    for (iElem = 0; err == 0 && iElem < nParallel_Tetr; iElem++) {
      iNode = iElem*N_POINTS_TETRAHEDRON;
      nodes[0] = Conn_Tetr_Par[iNode+0];
      nodes[1] = Conn_Tetr_Par[iNode+1];
      nodes[2] = Conn_Tetr_Par[iNode+2];
      nodes[3] = Conn_Tetr_Par[iNode+2];
      nodes[4] = Conn_Tetr_Par[iNode+3];
      nodes[5] = Conn_Tetr_Par[iNode+3];
      nodes[6] = Conn_Tetr_Par[iNode+3];
      nodes[7] = Conn_Tetr_Par[iNode+3];
      err = tecZoneNodeMapWrite64(file_handle, zone, rank, 1, 8, nodes);
      if (err) cout << rank << ": Error outputting Tecplot node values." << endl;
    }

    for (iElem = 0; err == 0 && iElem < nParallel_Hexa; iElem++) {
      iNode = iElem*N_POINTS_HEXAHEDRON;
      nodes[0] = Conn_Hexa_Par[iNode+0];
      nodes[1] = Conn_Hexa_Par[iNode+1];
      nodes[2] = Conn_Hexa_Par[iNode+2];
      nodes[3] = Conn_Hexa_Par[iNode+3];
      nodes[4] = Conn_Hexa_Par[iNode+4];
      nodes[5] = Conn_Hexa_Par[iNode+5];
      nodes[6] = Conn_Hexa_Par[iNode+6];
      nodes[7] = Conn_Hexa_Par[iNode+7];
      err = tecZoneNodeMapWrite64(file_handle, zone, rank, 1, 8, nodes);
      if (err) cout << rank << ": Error outputting Tecplot node values." << endl;
    }
      
    for (iElem = 0; err == 0 && iElem < nParallel_Pris; iElem++) {
      iNode = iElem*N_POINTS_PRISM;
      nodes[0] = Conn_Pris_Par[iNode+0];
      nodes[1] = Conn_Pris_Par[iNode+1];
      nodes[2] = Conn_Pris_Par[iNode+1];
      nodes[3] = Conn_Pris_Par[iNode+2];
      nodes[4] = Conn_Pris_Par[iNode+3];
      nodes[5] = Conn_Pris_Par[iNode+4];
      nodes[6] = Conn_Pris_Par[iNode+4];
      nodes[7] = Conn_Pris_Par[iNode+5];
      err = tecZoneNodeMapWrite64(file_handle, zone, rank, 1, 8, nodes);
      if (err) cout << rank << ": Error outputting Tecplot node values." << endl;
    }
    
    for (iElem = 0; err == 0 && iElem < nParallel_Pyra; iElem++) {
      iNode = iElem*N_POINTS_PYRAMID;
      nodes[0] = Conn_Pyra_Par[iNode+0];
      nodes[1] = Conn_Pyra_Par[iNode+1];
      nodes[2] = Conn_Pyra_Par[iNode+2];
      nodes[3] = Conn_Pyra_Par[iNode+3];
      nodes[4] = Conn_Pyra_Par[iNode+4];
      nodes[5] = Conn_Pyra_Par[iNode+4];
      nodes[6] = Conn_Pyra_Par[iNode+4];
      nodes[7] = Conn_Pyra_Par[iNode+4];
      err = tecZoneNodeMapWrite64(file_handle, zone, rank, 1, 8, nodes);
      if (err) cout << rank << ": Error outputting Tecplot node values." << endl;
    }
      
  }

#endif
  
  err = tecFileWriterClose(&file_handle);
  if (err) cout << rank << ": Error finishing Tecplot file output." << endl;
  
#endif /* HAVE_TECIO */

}
