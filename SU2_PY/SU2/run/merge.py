## \file merge.py
#  \brief python package for merging meshes
#  \author Tom Economon, Trent Lukaczyk, Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
#  \version 2.0.4
#
# Stanford University Unstructured (SU2) Code
# Copyright (C) 2012 Aerospace Design Laboratory
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import os, sys, shutil, copy
from .. import io  as su2io
from interface import SOL as SU2_SOL

# ----------------------------------------------------------------------
#  Merge Mesh
# ----------------------------------------------------------------------

def merge( config ):
    """ info = SU2.run.merge(config)
        
        Merges mesh with:
            SU2.run.SOL()    (volume merging)
            internal scripts (surface merging)
            
        Assumptions:
            config.NUMBER_PART is set 
            Skip if config.NUMBER_PART > 1
            
        Inputs:
            config - an SU2 config
                
        Ouputs:
            info - an empty SU2 State
            
        Executes in:
            ./
    """
    
    # local copy
    konfig = copy.deepcopy(config)
    
    # check if needed
    partitions = konfig['NUMBER_PART']
    if partitions <= 1:
        return su2io.State()
    
    # special cases
    special_cases = su2io.get_specialCases(konfig)
    
    # # MERGING # #
    if 'WRT_UNSTEADY' in special_cases:
        merge_unsteady(konfig)
    else:
        merge_volume(konfig)
        merge_surface(konfig)
        
    # info out (empty)
    info = su2io.State()
    
    return info

#: merge

def merge_unsteady( config, begintime=0, endtime=None ):
    
    if not endtime:
        endtime = config.EXT_ITER
    
    # SU2_SOL handles unsteady volume merge
    merge_volume( config )
    
    # python handles unsteady surface merge
    for timestep in range( 1, endtime-begintime+1 ):
        timestep += begintime
        merge_surface( config, timestep )

    return

#: def merge_unsteady()

def merge_volume( config ):
    """ SU2.io.merge.merge_volume(config)
        general volume surface merging with SU2_SOL
    """
    
    SU2_SOL( config )
    
    return

#: merge_volume( config )
    

def merge_surface( config, timestep = None ):
    """ SU2.io.merge.merge_surface(config,timestep=None)
        tecplot surface file merging
    """
  
    if timestep: timestep = '%05d' % (timestep-1)

    # pull config params
    math_problem               = config.MATH_PROBLEM
    surface_flow_filename      = config.SURFACE_FLOW_FILENAME
    surface_adjoint_filename   = config.SURFACE_ADJ_FILENAME
    mesh_filename              = config.MESH_FILENAME
    output_format              = config.OUTPUT_FORMAT
    partitions                 = config.NUMBER_PART

    assert output_format == 'TECPLOT' , 'can only merge tecplot surfaces'
    
    # -------------------------------------------------------------------
    #  Surface File Merging
    # -------------------------------------------------------------------
    
    # Write the surface .plt file
    npoint = 0; nelem = 0; addpoint = [ 0 ]
    for domain in range(int(partitions)):
        if not timestep:
            if math_problem == "DIRECT": input_file = open("%s_%s.plt" % (surface_flow_filename, domain+1))
            elif math_problem == "ADJOINT": input_file = open("%s_%s.plt" % (surface_adjoint_filename, domain+1))
        else:
            if math_problem == "DIRECT": input_file = open("%s_%s_%s.plt" % (surface_flow_filename, domain+1, timestep))
            elif math_problem == "ADJOINT": input_file = open("%s_%s_%s.plt" % (surface_adjoint_filename, domain+1, timestep))         
        for line in input_file:
            if "VARIABLES" in line:
                variables_list = line
            elif "ZONE" in line:
                if not timestep:
                    local_point = int(line.replace("ZONE"," ").split(",")[0].split("=")[1])
                    npoint = npoint + local_point
                    if domain == 0: addpoint.append(local_point)
                    else: addpoint.append(local_point+addpoint[domain])
                    nelem = nelem + int(line.replace("ZONE"," ").split(",")[1].split("=")[1])
                    datapacking = line.replace("ZONE"," ").split(",")[2].split("=")[1]
                    kindzone = line.replace("ZONE"," ").split(",")[3].split("=")[1]
                else:
                    global_ID = int(line.replace("ZONE"," ").split(",")[0].split("=")[1])
                    global_time = line.replace("ZONE"," ").split(",")[1].split("=")[1]
                    local_point = int(line.replace("ZONE"," ").split(",")[2].split("=")[1])
                    npoint = npoint + local_point
                    if domain == 0: addpoint.append(local_point)
                    else: addpoint.append(local_point+addpoint[domain])
                    nelem = nelem + int(line.replace("ZONE"," ").split(",")[3].split("=")[1])
                    datapacking = line.replace("ZONE"," ").split(",")[4].split("=")[1]
                    kindzone = line.replace("ZONE"," ").split(",")[5].split("=")[1]
        input_file.close()
    
    # Open the output file
    if not timestep:
        if math_problem == "DIRECT": output_file = open("%s.plt" % (surface_flow_filename) ,"w")
        elif math_problem == "ADJOINT": output_file = open("%s.plt" % (surface_adjoint_filename) ,"w")
    else:
        if math_problem == "DIRECT": output_file = open("%s_%s.plt" % (surface_flow_filename, timestep) ,"w")
        elif math_problem == "ADJOINT": output_file = open("%s_%s.plt" % (surface_adjoint_filename, timestep) ,"w")
    
    # Write the header
    output_file.write("TITLE = \"Visualization of the surface grid\"\n")
    
    # Write the coordinates and solution
    output_file.write("%s" % variables_list)
    if not timestep:
        output_file.write("ZONE NODES=%s , ELEMENTS=%s, DATAPACKING=%s, ZONETYPE=%s" % (npoint, nelem, datapacking, kindzone))
    else:
        output_file.write("ZONE STRANDID=%s SOLUTIONTIME=%s NODES=%s , ELEMENTS=%s, DATAPACKING=%s, ZONETYPE=%s" % (global_ID, global_time, npoint, nelem, datapacking, kindzone))
    for domain in range(int(partitions)):
        if not timestep:
            if math_problem == "DIRECT": input_file = open("%s_%s.plt" % (surface_flow_filename, domain+1))
            elif math_problem == "ADJOINT": input_file = open("%s_%s.plt" % (surface_adjoint_filename, domain+1))
        else:
            if math_problem == "DIRECT": input_file = open("%s_%s_%s.plt" % (surface_flow_filename, domain+1, timestep))
            elif math_problem == "ADJOINT": input_file = open("%s_%s_%s.plt" % (surface_adjoint_filename, domain+1, timestep))  
        while 1:
            line = input_file.readline()
            if "ZONE" in line:
                if not timestep:
                    ipoint = int(line.replace("ZONE"," ").split(",")[0].split("=")[1])
                else:
                    ipoint = int(line.replace("ZONE"," ").split(",")[2].split("=")[1])
                for list_points in range(int(ipoint)):
                    new_line = input_file.readline()
                    output_file.write("%s" %  new_line)
            if not line: break
        input_file.close()
    
    # Write the elements
    for domain in range(int(partitions)):
        if not timestep:
            if math_problem == "DIRECT": input_file = open("%s_%s.plt" % (surface_flow_filename, domain+1))
            elif math_problem == "ADJOINT": input_file = open("%s_%s.plt" % (surface_adjoint_filename, domain+1))
        else:
            if math_problem == "DIRECT": input_file = open("%s_%s_%s.plt" % (surface_flow_filename, domain+1, timestep))
            elif math_problem == "ADJOINT": input_file = open("%s_%s_%s.plt" % (surface_adjoint_filename, domain+1, timestep))
        while 1:
            line = input_file.readline()
            if "ZONE" in line:
                if not timestep:
                    ipoint = int(line.replace("ZONE"," ").split(",")[0].split("=")[1])
                    ielem = int(line.replace("ZONE"," ").split(",")[1].split("=")[1])
                else:
                    ipoint = int(line.replace("ZONE"," ").split(",")[2].split("=")[1])
                    ielem = int(line.replace("ZONE"," ").split(",")[3].split("=")[1])
                for list_points in range(int(ipoint)):
                    new_line = input_file.readline()
                for list_elems in range(int(ielem)):
                    elem_info = input_file.readline().replace("\t"," ").split(" ")
                    cont = 0
                    for iter in elem_info:
                        elem_info[cont] = ("%s" % (int(iter)+addpoint[domain]))
                        cont = cont+1;
                    output_file.write("%s\n" %  ("\t".join(elem_info)))
            if not line: break
        input_file.close()
    
    for domain in range(int(partitions)):
        if not timestep:
            if math_problem == "DIRECT": os.remove("%s_%s.plt" % (surface_flow_filename, domain+1))
            elif math_problem == "ADJOINT": os.remove("%s_%s.plt" % (surface_adjoint_filename, domain+1))
        else:
            if math_problem == "DIRECT":  os.remove("%s_%s_%s.plt" % (surface_flow_filename, domain+1, timestep))
            elif math_problem == "ADJOINT":  os.remove("%s_%s_%s.plt" % (surface_adjoint_filename, domain+1, timestep))
        
    return

#: def merge_surface()


