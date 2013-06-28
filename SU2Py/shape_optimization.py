#!/usr/bin/python

## \file shape_optimization.py
#  \brief Python script for doing a the shape optimization.
#  \author Current Development: Stanford University.
#          Original Structure: CADES 1.0 (2009).
#  \version 1.1.
#
#  /usr/bin/python
#  /usr/local/epd-6.3-2-rh5-x86_64/bin/python (oscar)
#  /usr/local/bin/python2.7 (gauss)
#  /Library/Frameworks/Python.framework/Versions/2.6/bin/python (laptop)
#  /usr/local/Python-2.6.7/python (desktop)
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


from optparse import OptionParser
from scipy.optimize import fmin_bfgs
from scipy.optimize import fmin_slsqp
import time, os, shutil, subprocess, numpy

counter = 0

parser=OptionParser()
parser.add_option("-f", "--file", dest="filename",
                  help="read config from FILE", metavar="FILE")
parser.add_option("-p", "--partitions", dest="partitions", default=1,
                  help="number of PARTITIONS", metavar="PARTITIONS")
parser.add_option("-g", "--gradient", dest="gradient", default="Adjoint",
                  help="Method for computing the GRADIENT (Adjoint or FinDiff)", metavar="GRADIENT")

(options, args)=parser.parse_args()

if int(options.partitions) == 1:
    parallel = False
else : parallel = True

output_format = "TECPLOT"
free_surface = "NO"
equivalent_area = "NO"
rotating_frame = "NO"
    
# Read design variable information
for line in open(options.filename):
    if "DEFINITION_DV=" in line:
        CompleteDVInfo = line.split("=")[1]
        UnitaryDVInfo = CompleteDVInfo.split(";")
        KindDVInfo= [ [ 0 ] for i in range(len(UnitaryDVInfo)) ]
        MarkerInfo= [ [ 0 ] for i in range(len(UnitaryDVInfo)) ]
        DesignVar= [ [ 0 ] for i in range(len(UnitaryDVInfo)) ]
        for i in range(len(UnitaryDVInfo)) :
            CleanUnitaryDVInfo = UnitaryDVInfo[i].strip().strip("(").strip(")")
            GeneralInfo = CleanUnitaryDVInfo.split("|")
            KindDVInfo[i] = GeneralInfo[0].strip().split(",")
            MarkerInfo[i] = GeneralInfo[1].strip().split(",")
            DesignVar[i] = GeneralInfo[2].strip().split(",")
        for i in range(len(UnitaryDVInfo)) :
            KindDVInfo[i][0] = KindDVInfo[i][0].replace('101','MACH_NUMBER')
            KindDVInfo[i][0] = KindDVInfo[i][0].replace('102','AOA')
            KindDVInfo[i][0] = KindDVInfo[i][0].replace('10','FFD_ROTATION')
            KindDVInfo[i][0] = KindDVInfo[i][0].replace('11','FFD_CAMBER')
            KindDVInfo[i][0] = KindDVInfo[i][0].replace('12','FFD_THICKNESS')
            KindDVInfo[i][0] = KindDVInfo[i][0].replace('13','FFD_VOLUME')
            KindDVInfo[i][0] = KindDVInfo[i][0].replace('1','HICKS_HENNE')
            KindDVInfo[i][0] = KindDVInfo[i][0].replace('4','NACA_4DIGITS')
            KindDVInfo[i][0] = KindDVInfo[i][0].replace('5','DISPLACEMENT')
            KindDVInfo[i][0] = KindDVInfo[i][0].replace('6','ROTATION')
            KindDVInfo[i][0] = KindDVInfo[i][0].replace('7','FFD_CONTROL_POINT')
            KindDVInfo[i][0] = KindDVInfo[i][0].replace('8','FFD_DIHEDRAL_ANGLE')
            KindDVInfo[i][0] = KindDVInfo[i][0].replace('9','FFD_TWIST_ANGLE')
            KindDVInfo[i][1] = float(KindDVInfo[i][1])
            for j in range(len(DesignVar[i])) :
                DesignVar[i][j] = float(DesignVar[i][j])
            for j in range(len(MarkerInfo[i])) :
                MarkerInfo[i][j] = MarkerInfo[i][j]
    elif "OUTPUT_FORMAT=" in line:
        output_format = line.split("=")[1].strip()
    elif "FREE_SURFACE=" in line:
        free_surface = line.split("=")[1].strip()
    elif "EQUIV_AREA=" in line:
        equivalent_area = line.split("=")[1].strip()
    elif "ROTATING_FRAME=" in line:
        rotating_frame = line.split("=")[1].strip()  

nVar = len(DesignVar)
x_old = numpy.zeros(nVar)
x_new = numpy.zeros(nVar)
x_incr = numpy.zeros(nVar)

if (output_format == "PARAVIEW"):
    History_file = "history_opt_" + options.filename.replace(".cfg",".csv")
elif (output_format == "TECPLOT"):
    History_file = "history_opt_" + options.filename.replace(".cfg",".plt")
    
Config_MDC_file = "config_opt_MDC_" + options.filename
Config_CFD_file = "config_opt_CFD_" + options.filename
Config_GRAD_file = "config_opt_GRAD_" + options.filename
Mesh_MDC_file = "mesh_opt_MDC_" + options.filename.replace(".cfg",".su2")
if (output_format == "PARAVIEW"): optimization_file = "optimization_" + options.filename.replace(".cfg",".csv")
elif (output_format == "TECPLOT"): optimization_file = "optimization_" + options.filename.replace(".cfg",".plt")

# Create optimization history file
Opt_file = open(optimization_file,"w")

if (output_format == "TECPLOT"):
    Opt_file.write( "TITLE = \"SU2 Optimization\"\n" )
    Opt_file.write( "VARIABLES = ")
    
Opt_file.write( "\"Iteration\",\"CLift\",\"CDrag\",\"CSideForce\",\"CMx\",\"CMy\",\"CMz\",\"CFx\",\"CFy\",\"CFz\",\"CEff\"" )
if free_surface == "YES": Opt_file.write( ",\"CFreeSurface\"\n" )
elif rotating_frame == "YES": Opt_file.write( ",\"CMerit\",\"CT\",\"CQ\"\n" )
elif equivalent_area == "YES": Opt_file.write( ",\"CEquivArea\",\"CNearField\"\n" )
else: Opt_file.write( "\n" )

Opt_file.close()

###########################################################################
# Function for animating optimization files in Tecplot
def add_tecplot_iter(filename,iter):
  
  # Open the solution file for this design iteration
  # and find the "ZONE ..." line
  New_file = open("temp_file.plt","w")
  for line in open(filename):
    if "ZONE" in line:
      aug_line = line.replace('ZONE','ZONE STRANDID=' + str(iter) + ', SOLUTIONTIME=' + str(iter) + ',',1)
      New_file.write(aug_line)
    else:
      New_file.write(line)
  New_file.close()
  
  # Overwrite old solution with new solution file
  os.rename("temp_file.plt",filename)

###########################################################################
# Function for dpdating the computational grid
def update_grid(x):

# Copy all the design variables 
    global x_new, x_old
    for iDesignVar in range(len(DesignVar)):
        x_new[iDesignVar] = x[iDesignVar]
        
# Change the parameters of the design variables
    output_file = open(Config_MDC_file,"w")
    for line in open(options.filename):
        if "DV_KIND=" in line:
            for iDesignVar in range(len(DesignVar)):
                if iDesignVar == 0: output_file.write( "DV_KIND= ")
                output_file.write("%s" % KindDVInfo[iDesignVar][0])
                if iDesignVar != len(DesignVar)-1: output_file.write( " ")
                else: output_file.write( "\n")
        elif "DV_MARKER=" in line:
            iDesignVar = 0
            if len(MarkerInfo[iDesignVar]) == 1:
                output_file.write( "DV_MARKER= ( %s ) \n" % MarkerInfo[iDesignVar][0])
            else:
                for j in range(len(MarkerInfo[iDesignVar])):
                    if j == 0:  output_file.write( "DV_MARKER= ( %s, " % MarkerInfo[iDesignVar][j])
                    if j > 0 and j < len(MarkerInfo[iDesignVar])-1: output_file.write( "%s, " % MarkerInfo[iDesignVar][j])
                    if j == len(MarkerInfo[iDesignVar])-1:  output_file.write( "%s ) \n" % MarkerInfo[iDesignVar][j])               
        elif "DV_PARAM=" in line:
            for iDesignVar in range(len(DesignVar)):
                if iDesignVar == 0: output_file.write( "DV_PARAM= ")
                for j in range(len(DesignVar[iDesignVar])):
                    if j == 0:  output_file.write( "( %s, " % DesignVar[iDesignVar][j])
                    if j > 0 and j < len(DesignVar[iDesignVar])-1: output_file.write( "%s, " % DesignVar[iDesignVar][j])
                    if j == len(DesignVar[iDesignVar])-1:  output_file.write( "%s )" % DesignVar[iDesignVar][j])
                if iDesignVar != len(DesignVar)-1: output_file.write("; ")
                else: output_file.write( "\n")
        elif "DV_VALUE_OLD=" in line:
            for iDesignVar in range(len(DesignVar)):
                if iDesignVar == 0: output_file.write( "DV_VALUE_OLD= ")
                step = float(x_old[iDesignVar])
                output_file.write("%.15f" % step)
                if iDesignVar != len(DesignVar)-1: output_file.write( " ")
                else: output_file.write( "\n")
        elif "DV_VALUE_NEW=" in line:
            for iDesignVar in range(len(DesignVar)):
                if iDesignVar == 0: output_file.write( "DV_VALUE_NEW= ")
                step = float(x_new[iDesignVar]) + 1E-12
                output_file.write("%.15f" % step)
                if iDesignVar != len(DesignVar)-1: output_file.write( " ")
                else: output_file.write( "\n")
        elif "MESH_FILENAME=" in line:
            output_file.write("MESH_FILENAME= %s \n" % Mesh_MDC_file)
        elif "MESH_OUT_FILENAME=" in line:
            output_file.write("MESH_OUT_FILENAME= %s \n" % Mesh_MDC_file)
        else:
            output_file.write(line)
    output_file.close()

# Apply grid update
    os.system ("SU2_MDC " + Config_MDC_file)

# Copy solution
    for iDesignVar in range(len(DesignVar)):
        x_old[iDesignVar] = x_new[iDesignVar]

# Remove configuration files                        
    os.remove(Config_MDC_file)

###########################################################################
# Function for creating the output files
def file_output(x):

    free_surface = "NO"
    equivalent_area = "NO"
    rotating_frame = "NO"
        
    output_format = "TECPLOT"

# Read the optimization problem
    for line in open(options.filename):
        if "OBJFUNC=" in line:
            objfunc = line.split("=")[1].strip()
        elif "OBJFUNC_SCALE=" in line:
            scale = float(line.split("=")[1].strip())
        elif "FREE_SURFACE=" in line:
            free_surface = line.split("=")[1].strip()
        elif "EQUIV_AREA=" in line:
            equivalent_area = line.split("=")[1].strip()
        elif "ROTATING_FRAME=" in line:
            rotating_frame = line.split("=")[1].strip()  
        elif "OUTPUT_FORMAT=" in line:
            output_format = line.split("=")[1].strip()

# Change grid file for CFD Computation over the deformed grid
    output_file = open(Config_CFD_file,"w")
    for line in open(options.filename):
        if "MATH_PROBLEM=" in line:
            output_file.write("MATH_PROBLEM= DIRECT \n")
        elif "MESH_FILENAME=" in line:
            output_file.write("MESH_FILENAME= %s \n" % Mesh_MDC_file)
        elif "CONV_FILENAME=" in line:
            output_file.write("CONV_FILENAME= %s \n" % History_file.replace(".csv","").replace(".plt",""))
        else: 
            output_file.write(line)
    output_file.close()
    
# CFD computation to evaluate the objective function
    global counter    
    counter = counter + 1

    shutil.copyfile(Mesh_MDC_file,"mesh_ShapeIter_"+str(counter)+".su2")

    if (output_format == "PARAVIEW"):
        os.rename("surface_flow.vtk","surface_flow_ShapeIter_"+str(counter)+".vtk")
        os.rename("flow.vtk","flow_ShapeIter_"+str(counter)+".vtk")
    elif (output_format == "TECPLOT"):
        os.rename("surface_flow.plt","surface_flow_ShapeIter_"+str(counter)+".plt")
        add_tecplot_iter("surface_flow_ShapeIter_"+str(counter)+".plt",counter)
        os.rename("flow.plt","flow_ShapeIter_"+str(counter)+".plt")
        add_tecplot_iter("flow_ShapeIter_"+str(counter)+".plt",counter)
         
    if KindDVInfo[0][0] == "FFD_CONTROL_POINT" or KindDVInfo[0][0] == "FFD_DIHEDRAL_ANGLE" or KindDVInfo[0][0] == "FFD_TWIST_ANGLE" or KindDVInfo[0][0] == "FFD_ROTATION" or KindDVInfo[0][0] == "FFD_CAMBER" or KindDVInfo[0][0] == "FFD_THICKNESS" or KindDVInfo[0][0] == "FFD_VOLUME":
        if (output_format == "PARAVIEW"):
            os.rename("deformed_chunk.vtk","FFDBox_ShapeIter_"+str(counter)+".vtk")
        if (output_format == "TECPLOT"):
            os.rename("deformed_chunk.plt","FFDBox_ShapeIter_"+str(counter)+".plt")
            add_tecplot_iter("FFDBox_ShapeIter_"+str(counter)+".plt",counter)        
        
    if equivalent_area == "YES":
        if (output_format == "PARAVIEW"):
            os.rename("EquivArea.csv","EquivArea_ShapeIter_"+str(counter)+".csv")
        if (output_format == "TECPLOT"):
            os.rename("EquivArea.plt","EquivArea_ShapeIter_"+str(counter)+".plt")
            add_tecplot_iter("EquivArea_ShapeIter_"+str(counter)+".plt",counter)

    if free_surface == "YES":
        if (output_format == "PARAVIEW"):
            os.rename("LevelSet.vtk","LevelSet_ShapeIter_"+str(counter)+".vtk")
        if (output_format == "TECPLOT"):
            os.rename("LevelSet.plt","LevelSet_ShapeIter_"+str(counter)+".plt")
            add_tecplot_iter("LevelSet_ShapeIter_"+str(counter)+".plt",counter)

# Read the Objective function    
    for line in open(History_file): 
        if "ERROR" in line: error = 0.0
        
    CLift_Value = float(line.split(",")[1])
    CDrag_Value = float (line.split(",")[2])
    CSideForce_Value = float (line.split(",")[3])
    CMx_Value = float (line.split(",")[4])
    CMy_Value = float (line.split(",")[5])
    CMz_Value = float (line.split(",")[6])
    CFx_Value = float (line.split(",")[7])
    CFy_Value = float (line.split(",")[8])
    CFz_Value = float (line.split(",")[9])
    CEff_Value = float (line.split(",")[10])
    if free_surface == "YES":
        CFreeSurface_Value = float (line.split(",")[11])
    if rotating_frame == "YES":
        CMerit_Value = float (line.split(",")[11])
        CT_Value = float (line.split(",")[12])
        CQ_Value = float (line.split(",")[13])
    if equivalent_area == "YES":
        CEquivArea_Value = float (line.split(",")[11])
        CNearField_Value = float (line.split(",")[12])

# Remove configuration and mesh files            
    os.remove(Config_CFD_file)
    os.remove(History_file)

# Output optimization
    Opt_file = open(optimization_file,"a")

    Opt_file.write( "%d, \t %.10f, \t %.10f, \t %.10f, \t %.10f, \t %.10f, \t %.10f, \t %.10f, \t %.10f, \t %.10f, \t %.10f" % (counter, CLift_Value, CDrag_Value, CSideForce_Value, CMx_Value, CMy_Value, CMz_Value, CFx_Value, CFy_Value, CFz_Value, CEff_Value) )
    
    if rotating_frame == "YES": Opt_file.write( ", \t %.10f, \t %.10f, \t %.10f\n" % (CMerit_Value, CT_Value, CQ_Value) )
    elif equivalent_area == "YES": Opt_file.write( ", \t %.10f, \t %.10f\n" % (CEquivArea_Value, CNearField_Value) )
    elif free_surface == "YES": Opt_file.write( ", \t %.10f\n" % (CFreeSurface_Value) )
    else: Opt_file.write( "\n" )
        
    Opt_file.close()

###########################################################################
# Function for computing the objective function
def f(x):

    free_surface = "NO"
    equivalent_area = "NO"
    rotating_frame = "NO"
        
    output_format = "TECPLOT"

# Read the optimization problem
    for line in open(options.filename):
        if "OBJFUNC=" in line:
            objfunc = line.split("=")[1].strip()
        elif "OBJFUNC_SCALE=" in line:
            scale = float(line.split("=")[1].strip())
        elif "FREE_SURFACE=" in line:
            free_surface = line.split("=")[1].strip()
        elif "EQUIV_AREA=" in line:
            equivalent_area = line.split("=")[1].strip()
        elif "ROTATING_FRAME=" in line:
            rotating_frame = line.split("=")[1].strip()  
        elif "OUTPUT_FORMAT=" in line:
            output_format = line.split("=")[1].strip()
    
# Update the computationa grid to the design variable state
    update_grid(x)

# Change grid file for CFD Computation over the deformed grid
    output_file = open(Config_CFD_file,"w")
    for line in open(options.filename):
        if "MATH_PROBLEM=" in line:
            output_file.write("MATH_PROBLEM= DIRECT \n")
        elif "MESH_FILENAME=" in line:
            output_file.write("MESH_FILENAME= %s \n" % Mesh_MDC_file)
        elif "CONV_FILENAME=" in line:
            output_file.write("CONV_FILENAME= %s \n" % History_file.replace(".csv","").replace(".plt",""))
        else: 
            output_file.write(line)
    output_file.close()
    
# CFD computation to evaluate the objective function
    global counter    
    counter = counter + 1
    if parallel: os.system ("parallel_computation.py -p %s -f %s -o True" % (int(options.partitions), Config_CFD_file))
    else: os.system ("SU2_CFD " + Config_CFD_file)

    shutil.copyfile(Mesh_MDC_file,"mesh_ShapeIter_"+str(counter)+".su2")

    if (output_format == "PARAVIEW"):
        os.rename("surface_flow.vtk","surface_flow_ShapeIter_"+str(counter)+".vtk")
        os.rename("flow.vtk","flow_ShapeIter_"+str(counter)+".vtk")
    elif (output_format == "TECPLOT"):
        os.rename("surface_flow.plt","surface_flow_ShapeIter_"+str(counter)+".plt")
        add_tecplot_iter("surface_flow_ShapeIter_"+str(counter)+".plt",counter)
        os.rename("flow.plt","flow_ShapeIter_"+str(counter)+".plt")
        add_tecplot_iter("flow_ShapeIter_"+str(counter)+".plt",counter)
         
    if KindDVInfo[0][0] == "FFD_CONTROL_POINT" or KindDVInfo[0][0] == "FFD_DIHEDRAL_ANGLE" or KindDVInfo[0][0] == "FFD_TWIST_ANGLE" or KindDVInfo[0][0] == "FFD_ROTATION" or KindDVInfo[0][0] == "FFD_CAMBER" or KindDVInfo[0][0] == "FFD_THICKNESS" or KindDVInfo[0][0] == "FFD_VOLUME":
        if (output_format == "PARAVIEW"):
            os.rename("deformed_chunk.vtk","FFDBox_ShapeIter_"+str(counter)+".vtk")
        if (output_format == "TECPLOT"):
            os.rename("deformed_chunk.plt","FFDBox_ShapeIter_"+str(counter)+".plt")
            add_tecplot_iter("FFDBox_ShapeIter_"+str(counter)+".plt",counter)        
        
    if equivalent_area == "YES":
        if (output_format == "PARAVIEW"):
            os.rename("EquivArea.csv","EquivArea_ShapeIter_"+str(counter)+".csv")
        if (output_format == "TECPLOT"):
            os.rename("EquivArea.plt","EquivArea_ShapeIter_"+str(counter)+".plt")
            add_tecplot_iter("EquivArea_ShapeIter_"+str(counter)+".plt",counter)

    if free_surface == "YES":
        if (output_format == "PARAVIEW"):
            os.rename("LevelSet.vtk","LevelSet_ShapeIter_"+str(counter)+".vtk")
        if (output_format == "TECPLOT"):
            os.rename("LevelSet.plt","LevelSet_ShapeIter_"+str(counter)+".plt")
            add_tecplot_iter("LevelSet_ShapeIter_"+str(counter)+".plt",counter)

# Read the Objective function    
    for line in open(History_file): 
        if "ERROR" in line: error = 0.0
        
    CLift_Value = float(line.split(",")[1])
    CDrag_Value = float (line.split(",")[2])
    CSideForce_Value = float (line.split(",")[3])
    CMx_Value = float (line.split(",")[4])
    CMy_Value = float (line.split(",")[5])
    CMz_Value = float (line.split(",")[6])
    CFx_Value = float (line.split(",")[7])
    CFy_Value = float (line.split(",")[8])
    CFz_Value = float (line.split(",")[9])
    CEff_Value = float (line.split(",")[10])
    if free_surface == "YES":
        CFreeSurface_Value = float (line.split(",")[11])
    if rotating_frame == "YES":
        CMerit_Value = float (line.split(",")[11])
        CT_Value = float (line.split(",")[12])
        CQ_Value = float (line.split(",")[13])
    if equivalent_area == "YES":
        CEquivArea_Value = float (line.split(",")[11])
        CNearField_Value = float (line.split(",")[12])

# Remove configuration and mesh files            
    os.remove(Config_CFD_file)
    os.remove(History_file)

# Output optimization
    Opt_file = open(optimization_file,"a")

    Opt_file.write( "%d, \t %.10f, \t %.10f, \t %.10f, \t %.10f, \t %.10f, \t %.10f, \t %.10f, \t %.10f, \t %.10f, \t %.10f" % (counter, CLift_Value, CDrag_Value, CSideForce_Value, CMx_Value, CMy_Value, CMz_Value, CFx_Value, CFy_Value, CFz_Value, CEff_Value) )
    
    if rotating_frame == "YES": Opt_file.write( ", \t %.10f, \t %.10f, \t %.10f\n" % (CMerit_Value, CT_Value, CQ_Value) )
    elif equivalent_area == "YES": Opt_file.write( ", \t %.10f, \t %.10f\n" % (CEquivArea_Value, CNearField_Value) )
    elif free_surface == "YES": Opt_file.write( ", \t %.10f\n" % (CFreeSurface_Value) )
    else: Opt_file.write( "\n" )
        
    Opt_file.close()

# There is a change in the sign in case we need to maximize a functional
    sign = 1.0;
    if objfunc == "EFFICIENCY": sign = -1.0;
    if objfunc == "LIFT": sign = -1.0;
    if objfunc == "THRUST": sign = -1.0;
    if objfunc == "FIGURE_OF_MERIT": sign = -1.0;
    
    if objfunc == "DRAG":  return sign*CDrag_Value*scale
    if objfunc == "LIFT":  return sign*CLift_Value*scale
    if objfunc == "SIDEFORCE":  return sign*CSideForce_Value*scale
    if objfunc == "EFFICIENCY": return sign*CEff_Value*scale
    if objfunc == "EQUIVALENT_AREA": return sign*CEquivArea_Value*scale
    if objfunc == "NEARFIELD_PRESSURE": return sign*CNearField_Value*scale
    if objfunc == "FORCE_X": return sign*CFx_Value*scale
    if objfunc == "FORCE_Y": return sign*CFy_Value*scale
    if objfunc == "FORCE_Z": return sign*CFz_Value*scale
    if objfunc == "FIGURE_OF_MERIT": return sign*CMerit_Value*scale
    if objfunc == "THRUST": return sign*CT_Value*scale
    if objfunc == "TORQUE": return sign*CQ_Value*scale
    if objfunc == "FREESURFACE": return sign*CFreeSurface_Value*scale

###########################################################################
# Function for computing the derivative of the objective function
def df(x):

    free_surface = "NO"
    equivalent_area = "NO"
    rotating_frame = "NO"
    output_format = "TECPLOT"

# Update the computationa grid to the design variable state
    update_grid(x)

# Read the optimization problem
    for line in open(options.filename):
        if "OBJFUNC=" in line:
            objfunc = line.split("=")[1].strip()
        elif "OBJFUNC_SCALE=" in line:
            scale = float(line.split("=")[1].strip())
        elif "FREE_SURFACE=" in line:
            free_surface = line.split("=")[1].strip()
        elif "OUTPUT_FORMAT=" in line:
            output_format = line.split("=")[1].strip()

# Change grid file for Gradient Computation over the deformed grid
    output_file = open(Config_GRAD_file,"w")
    for line in open(options.filename):
        if "CADJ_OBJFUNC=" in line:
            output_file.write("CADJ_OBJFUNC= %s \n" % objfunc)
        elif "MESH_FILENAME=" in line:
            output_file.write("MESH_FILENAME= %s \n" % Mesh_MDC_file)
        elif "CONV_FILENAME=" in line:
            output_file.write("CONV_FILENAME= %s \n" % History_file.replace(".csv","").replace(".plt",""))
        else:
            output_file.write(line)
    output_file.close()
    
# Compute Gradient using continuous adjoint strategy
    if (options.gradient == "Adjoint"):
        if parallel: os.system ("continuous_adjoint.py -p %s -f %s" % (int(options.partitions), Config_GRAD_file))
        else: os.system ("continuous_adjoint.py -f " + Config_GRAD_file)
    elif (options.gradient == "FinDiff"):
        if parallel: os.system ("finite_differences.py -p %s -f %s" % (int(options.partitions), Config_GRAD_file))
        else: os.system ("finite_differences.py -f " + Config_GRAD_file)

# Read the objective function gradient
    r = numpy.zeros(nVar)

    if (options.gradient == "Adjoint"):
        if (output_format == "PARAVIEW"):
            cadj_grad = open("cont_adj_" + Config_GRAD_file.replace(".cfg",".csv"),"r")
        elif (output_format == "TECPLOT"):
            cadj_grad = open("cont_adj_" + Config_GRAD_file.replace(".cfg",".plt"),"r")
    elif (options.gradient == "FinDiff"):
        if (output_format == "PARAVIEW"):
            cadj_grad = open("fin_diff_" + Config_GRAD_file.replace(".cfg",".csv"),"r")
        elif (output_format == "TECPLOT"):
            cadj_grad = open("fin_diff_" + Config_GRAD_file.replace(".cfg",".plt"),"r")
        
# There is a change in the sign in case we need to maximize a functional
    sign = 1.0;
    if objfunc == "EFFICIENCY": sign = -1.0;
    if objfunc == "LIFT": sign = -1.0;
    if objfunc == "THRUST": sign = -1.0;
    if objfunc == "FIGURE_OF_MERIT": sign = -1.0;

    line = cadj_grad.readline()
    iVar = 0
    for line in cadj_grad:
        if options.gradient == "Adjoint":
            r[iVar] = sign*float(line.split(",")[1])*scale
        elif options.gradient == "FinDiff":
            if objfunc == "LIFT": r[iVar] = sign*float(line.split(",")[1])*scale
            elif objfunc == "DRAG": r[iVar] = sign*float(line.split(",")[2])*scale
            elif objfunc == "SIDEFORCE": r[iVar] = sign*float(line.split(",")[3])*scale
            elif objfunc == "MOMENT_X": r[iVar] = sign*float(line.split(",")[4])*scale
            elif objfunc == "MOMENT_Y": r[iVar] = sign*float(line.split(",")[5])*scale
            elif objfunc == "MOMENT_Z": r[iVar] = sign*float(line.split(",")[6])*scale
            elif objfunc == "FORCE_X": r[iVar] = sign*float(line.split(",")[7])*scale
            elif objfunc == "FORCE_Y": r[iVar] = sign*float(line.split(",")[8])*scale
            elif objfunc == "FORCE_Z": r[iVar] = sign*float(line.split(",")[9])*scale
            elif objfunc == "EFFICIENCY": r[iVar] = sign*float(line.split(",")[10])*scale
            elif objfunc == "EQUIVALENT_AREA": r[iVar] = sign*float(line.split(",")[11])*scale
            elif objfunc == "NEARFIELD_PRESSURE": r[iVar] = sign*float(line.split(",")[12])*scale
            elif objfunc == "FIGURE_OF_MERIT": r[iVar] = sign*float(line.split(",")[11])*scale
            elif objfunc == "THRUST": r[iVar] = sign*float(line.split(",")[12])*scale
            elif objfunc == "TORQUE": r[iVar] = sign*float(line.split(",")[13])*scale
            elif objfunc == "FREESURFACE": r[iVar] = sign*float(line.split(",")[11])*scale
            
        iVar = iVar + 1
        
    cadj_grad.close()

# Remove configuration file        
    os.remove(Config_GRAD_file)

    return r

###########################################################################
# Function for computing the equality constraint
def eqcons(x):

    free_surface = "NO"
    equivalent_area = "NO"
    rotating_frame = "NO"

# Read the optimization problem
    for line in open(options.filename):
        if "CONST_EQ=" in line:
            CompleteEqConst = line.split("=")[1]
            UnitaryEqConst = CompleteEqConst.split(",")
        elif "CONST_EQ_SCALE=" in line:
            CompleteEqConst_Scale = line.split("=")[1]
            UnitaryEqConst_Scale = CompleteEqConst_Scale.split(",")
        elif "CONST_EQ_VALUE=" in line:
            CompleteEqConst_Value = line.split("=")[1]
            UnitaryEqConst_Value = CompleteEqConst_Value.split(",")
        elif "FREE_SURFACE=" in line:
            free_surface = line.split("=")[1].strip()
        elif "EQUIV_AREA=" in line:
            equivalent_area = line.split("=")[1].strip()
        elif "ROTATING_FRAME=" in line:
            rotating_frame = line.split("=")[1].strip()  

    if UnitaryEqConst[0].strip() != "NONE" :

# Update the computationa grid to the design variable state
        update_grid(x)

# Change grid file for CFD Computation over the deformed grid
        output_file = open(Config_CFD_file,"w")
        for line in open(options.filename):
            if "MATH_PROBLEM=" in line:
                output_file.write("MATH_PROBLEM= DIRECT \n")
            elif "MESH_FILENAME=" in line:
                output_file.write("MESH_FILENAME= %s \n" % Mesh_MDC_file)
            elif "CONV_FILENAME=" in line:
                output_file.write("CONV_FILENAME= %s \n" % History_file.replace(".csv","").replace(".plt",""))
            else: 
                output_file.write(line)
        output_file.close()
    
# CFD computation to evaluate Drag and Lift Coefficient
        if parallel: os.system ("parallel_computation.py -p %s -f %s -o False" % (int(options.partitions), Config_CFD_file))
        else: os.system ("SU2_CFD " + Config_CFD_file)

# Read the Drag, Lift and Pitching Coefficient    
        for line in open(History_file): 
            if "ERROR" in line: error = 0.0

        CLift_Value = float(line.split(",")[1])
        CDrag_Value = float (line.split(",")[2])
        CSideForce_Value = float (line.split(",")[3])
        CMx_Value = float (line.split(",")[4])
        CMy_Value = float (line.split(",")[5])
        CMz_Value = float (line.split(",")[6])
        CFx_Value = float (line.split(",")[7])
        CFy_Value = float (line.split(",")[8])
        CFz_Value = float (line.split(",")[9])
        CEff_Value = float (line.split(",")[10])
        if equivalent_area == "YES":
            CEquivArea_Value = float (line.split(",")[11])
            CNearField_Value = float (line.split(",")[12])
        if rotating_frame == "YES":
            CMerit_Value = float (line.split(",")[11])
            CT_Value = float (line.split(",")[12])
            CQ_Value = float (line.split(",")[13])
        if free_surface == "YES":
            CFreeSurface_Value = float (line.split(",")[11])

# Remove configuration and mesh files            
        os.remove(Config_CFD_file)
        os.remove(History_file)

# Return the value of the inequality constraint
        r = numpy.zeros(len(UnitaryEqConst))
        
        for i in range(len(UnitaryEqConst)) :
            if UnitaryEqConst[i].strip() == "DRAG": Value = CDrag_Value
            elif UnitaryEqConst[i].strip() == "LIFT": Value = CLift_Value
            elif UnitaryEqConst[i].strip() == "SIDEFORCE": Value = CSideForce_Value
            elif UnitaryEqConst[i].strip() == "PRESSURE": Value = CPress_Value
            elif UnitaryEqConst[i].strip() == "MOMENT_X": Value = CMx_Value
            elif UnitaryEqConst[i].strip() == "MOMENT_Y": Value = CMy_Value
            elif UnitaryEqConst[i].strip() == "MOMENT_Z": Value = CMz_Value
            elif UnitaryEqConst[i].strip() == "FORCE_X": Value = CFx_Value
            elif UnitaryEqConst[i].strip() == "FORCE_Y": Value = CFy_Value
            elif UnitaryEqConst[i].strip() == "FORCE_Z": Value = CFz_Value
            elif UnitaryEqConst[i].strip() == "EFFICIENCY": Value = CEff_Value
            elif UnitaryEqConst[i].strip() == "EQUIVALENT_AREA": Value = CEquivArea_Value
            elif UnitaryEqConst[i].strip() == "NEARFIELD_PRESSURE": Value = CNearField_Value
            elif UnitaryEqConst[i].strip() == "FIGURE_OF_MERIT": Value = CMerit_Value
            elif UnitaryEqConst[i].strip() == "THRUST": Value = CT_Value
            elif UnitaryEqConst[i].strip() == "TORQUE": Value = CQ_Value
            elif UnitaryEqConst[i].strip() == "FREESURFACE": Value = CFreeSurface_Value
            
            r[i] = (Value-float(UnitaryEqConst_Value[i]))

    else:
        r = numpy.zeros(0)

    return r

###########################################################################
# Function for computing the derivative of the equality constraint
def deqcons(x):

    free_surface = "NO"
    equivalent_area = "NO"
    rotating_frame = "NO"
    output_format = "TECPLOT"

# Read the optimization problem
    for line in open(options.filename):
        if "CONST_EQ=" in line:
            CompleteEqConst = line.split("=")[1]
            UnitaryEqConst = CompleteEqConst.split(",")
        elif "CONST_EQ_SCALE=" in line:
            CompleteEqConst_Scale = line.split("=")[1]
            UnitaryEqConst_Scale = CompleteEqConst_Scale.split(",")
        elif "CONST_EQ_VALUE=" in line:
            CompleteEqConst_Value = line.split("=")[1]
            UnitaryEqConst_Value = CompleteEqConst_Value.split(",")
        elif "FREE_SURFACE=" in line:
            free_surface = line.split("=")[1].strip()
        elif "EQUIV_AREA=" in line:
            equivalent_area = line.split("=")[1].strip()
        elif "ROTATING_FRAME=" in line:
            rotating_frame = line.split("=")[1].strip()  
        elif "OUTPUT_FORMAT=" in line:
            output_format = line.split("=")[1].strip()

    # Read the conts function gradient
    r = numpy.zeros([len(UnitaryEqConst), nVar])
    
    for i in range(len(UnitaryEqConst)) :
        if UnitaryEqConst[i].strip() != "NONE" :
            # Update the computationa grid to the design variable state
            update_grid(x)

            # Change grid file for Gradient Computation over the deformed grid
            output_file = open(Config_GRAD_file,"w")
            for line in open(options.filename):
                if "CADJ_OBJFUNC=" in line:
                    output_file.write("CADJ_OBJFUNC= %s \n" % UnitaryEqConst[i])
                elif "MESH_FILENAME=" in line:
                    output_file.write("MESH_FILENAME= %s \n" % Mesh_MDC_file)
                elif "CONV_FILENAME=" in line:
                    output_file.write("CONV_FILENAME= %s \n" % History_file.replace(".csv","").replace(".plt",""))
                else:
                    output_file.write(line)
            output_file.close()
    
            # Compute Gradient using continuous adjoint strategy
            if parallel: os.system ("continuous_adjoint.py -p %s -f %s" % (int(options.partitions), Config_GRAD_file))
            else: os.system ("continuous_adjoint.py -f " + Config_GRAD_file)

            if (options.gradient == "Adjoint"):
                if (output_format == "PARAVIEW"):
                    cadj_grad = open("cont_adj_" + Config_GRAD_file.replace(".cfg",".csv"),"r")
                elif (output_format == "TECPLOT"):
                    cadj_grad = open("cont_adj_" + Config_GRAD_file.replace(".cfg",".plt"),"r")
            elif (options.gradient == "FinDiff"):
                if (output_format == "PARAVIEW"):
                    cadj_grad = open("fin_diff_" + Config_GRAD_file.replace(".cfg",".csv"),"r")
                elif (output_format == "TECPLOT"):
                    cadj_grad = open("fin_diff_" + Config_GRAD_file.replace(".cfg",".plt"),"r")        
                
            line = cadj_grad.readline()
            iVar = 0
            for line in cadj_grad:
                if (options.gradient == "Adjoint"):
                    r[i][iVar] = float(line.split(",")[1])*float(UnitaryEqConst_Scale[i])
                elif (options.gradient == "FinDiff"):
                    if UnitaryEqConst[i] == "LIFT": r[i][iVar] = float(line.split(",")[1])*float(UnitaryEqConst_Scale[i])
                    elif UnitaryEqConst[i] == "DRAG": r[i][iVar] = float(line.split(",")[2])*float(UnitaryEqConst_Scale[i])
                    elif UnitaryEqConst[i] == "SIDEFORCE": r[i][iVar] = float(line.split(",")[3])*float(UnitaryEqConst_Scale[i])
                    elif UnitaryEqConst[i] == "MOMENT_X": r[i][iVar] = float(line.split(",")[4])*float(UnitaryEqConst_Scale[i])
                    elif UnitaryEqConst[i] == "MOMENT_Y": r[i][iVar] = float(line.split(",")[5])*float(UnitaryEqConst_Scale[i])
                    elif UnitaryEqConst[i] == "MOMENT_Z": r[i][iVar] = float(line.split(",")[6])*float(UnitaryEqConst_Scale[i])
                    elif UnitaryEqConst[i] == "FORCE_X": r[i][iVar] = float(line.split(",")[7])*float(UnitaryEqConst_Scale[i])
                    elif UnitaryEqConst[i] == "FORCE_Y": r[i][iVar] = float(line.split(",")[8])*float(UnitaryEqConst_Scale[i])
                    elif UnitaryEqConst[i] == "FORCE_Z": r[i][iVar] = float(line.split(",")[9])*float(UnitaryEqConst_Scale[i])
                    elif UnitaryEqConst[i] == "EFFICIENCY": r[i][iVar] = float(line.split(",")[10])*float(UnitaryEqConst_Scale[i])
                    elif UnitaryEqConst[i] == "EQUIVALENT_AREA": r[i][iVar] = float(line.split(",")[11])*float(UnitaryEqConst_Scale[i])
                    elif UnitaryEqConst[i] == "NEARFIELD_PRESSURE": r[i][iVar] = float(line.split(",")[12])*float(UnitaryEqConst_Scale[i])
                    elif UnitaryEqConst[i] == "FIGURE_OF_MERIT": r[i][iVar] = float(line.split(",")[11])*float(UnitaryEqConst_Scale[i])
                    elif UnitaryEqConst[i] == "THRUST": r[i][iVar] = float(line.split(",")[12])*float(UnitaryEqConst_Scale[i])
                    elif UnitaryEqConst[i] == "TORQUE": r[i][iVar] = float(line.split(",")[13])*float(UnitaryEqConst_Scale[i])
                    elif UnitaryEqConst[i] == "FREESURFACE": r[i][iVar] = float(line.split(",")[11])*float(UnitaryEqConst_Scale[i])
                    
                iVar = iVar + 1
            cadj_grad.close()

            # Remove configuration file        
            os.remove(Config_GRAD_file)
        else:
            r = numpy.zeros([0, nVar])
        
    return r

###########################################################################
# Function for computing the inequality constraint
def ieqcons(x):

    free_surface = "NO"
    equivalent_area = "NO"
    rotating_frame = "NO"

# Read the optimization problem
    for line in open(options.filename):
        if "CONST_IEQ=" in line:
            CompleteIeqConst = line.split("=")[1]
            UnitaryIeqConst = CompleteIeqConst.split(",")
        elif "CONST_IEQ_SCALE=" in line:
            CompleteIeqConst_Scale = line.split("=")[1]
            UnitaryIeqConst_Scale = CompleteIeqConst_Scale.split(",")
        elif "CONST_IEQ_VALUE=" in line:
            CompleteIeqConst_Value = line.split("=")[1]
            UnitaryIeqConst_Value = CompleteIeqConst_Value.split(",")
        elif "CONST_IEQ_SIGN=" in line:
            CompleteIeqConst_Sign = line.split("=")[1]
            UnitaryIeqConst_Sign = CompleteIeqConst_Sign.split(",")
        elif "FREE_SURFACE=" in line:
            free_surface = line.split("=")[1].strip()
        elif "EQUIV_AREA=" in line:
            equivalent_area = line.split("=")[1].strip()
        elif "ROTATING_FRAME=" in line:
            rotating_frame = line.split("=")[1].strip()  

    if UnitaryIeqConst[0].strip() != "NONE" :

# Update the computationa grid to the design variable state
        update_grid(x)

# Change grid file for CFD Computation over the deformed grid
        output_file = open(Config_CFD_file,"w")
        for line in open(options.filename):
            if "MATH_PROBLEM=" in line:
                output_file.write("MATH_PROBLEM= DIRECT \n")
            elif "MESH_FILENAME=" in line:
                output_file.write("MESH_FILENAME= %s \n" % Mesh_MDC_file)
            elif "CONV_FILENAME=" in line:
                output_file.write("CONV_FILENAME= %s \n" % History_file.replace(".csv","").replace(".plt",""))
            else: 
                output_file.write(line)
        output_file.close()
    
# CFD computation to evaluate Drag and Lift Coefficient
        if parallel: os.system ("parallel_computation.py -p %s -f %s -o False" % (int(options.partitions), Config_CFD_file))
        else: os.system ("SU2_CFD " + Config_CFD_file)

# Read the Drag, Lift and Pitching Coefficient    
        for line in open(History_file): 
            if "ERROR" in line: error = 0.0
    
        CLift_Value = float(line.split(",")[1])
        CDrag_Value = float (line.split(",")[2])
        CSideForce_Value = float (line.split(",")[3])
        CMx_Value = float (line.split(",")[4])
        CMy_Value = float (line.split(",")[5])
        CMz_Value = float (line.split(",")[6])
        CFx_Value = float (line.split(",")[7])
        CFy_Value = float (line.split(",")[8])
        CFz_Value = float (line.split(",")[9])
        CEff_Value = float (line.split(",")[10])
        if free_surface == "YES":
            CFreeSurface_Value = float (line.split(",")[11])
        if equivalent_area == "YES":
            CEquivArea_Value = float (line.split(",")[11])
            CNearField_Value = float (line.split(",")[12])
        if rotating_frame == "YES":
            CMerit_Value = float (line.split(",")[11])
            CT_Value = float (line.split(",")[12])
            CQ_Value = float (line.split(",")[13])

# Remove configuration and mesh files            
        os.remove(Config_CFD_file)
        os.remove(History_file)

# Return the value of the inequality constraint
        r = numpy.zeros(len(UnitaryIeqConst))
        
        for i in range(len(UnitaryIeqConst)) :
            if UnitaryIeqConst[i].strip() == "DRAG": Value = CDrag_Value
            elif UnitaryIeqConst[i].strip() == "LIFT": Value = CLift_Value
            elif UnitaryIeqConst[i].strip() == "SIDEFORCE": Value = CSideForce_Value
            elif UnitaryIeqConst[i].strip() == "MOMENT_X": Value = CMx_Value
            elif UnitaryIeqConst[i].strip() == "MOMENT_Y": Value = CMy_Value
            elif UnitaryIeqConst[i].strip() == "MOMENT_Z": Value = CMz_Value
            elif UnitaryIeqConst[i].strip() == "FORCE_X": Value = CFx_Value
            elif UnitaryIeqConst[i].strip() == "FORCE_Y": Value = CFy_Value
            elif UnitaryIeqConst[i].strip() == "FORCE_Z": Value = CFz_Value
            elif UnitaryIeqConst[i].strip() == "EFFICIENCY": Value = CEff_Value
            elif UnitaryIeqConst[i].strip() == "EQUIVALENT_AREA": Value = CEquivArea_Value
            elif UnitaryIeqConst[i].strip() == "NEARFIELD_PRESSURE": Value = CNearField_Value
            elif UnitaryIeqConst[i].strip() == "FIGURE_OF_MERIT": Value = CMerit_Value
            elif UnitaryIeqConst[i].strip() == "THRUST": Value = CT_Value
            elif UnitaryIeqConst[i].strip() == "TORQUE": Value = CQ_Value
            elif UnitaryIeqConst[i].strip() == "FREESURFACE": Value = CFreeSurface_Value
            
            if UnitaryIeqConst_Sign[i].strip() == "GREATER": r[i] = (Value-float(UnitaryIeqConst_Value[i]))*float(UnitaryIeqConst_Scale[i])
            if UnitaryIeqConst_Sign[i].strip() == "LESS": r[i] = (float(UnitaryIeqConst_Value[i])-Value)*float(UnitaryIeqConst_Scale[i])

    else:
        r = numpy.zeros(0)

    return r

###########################################################################
# Function for computing the derivative of the inequality constraint
def dieqcons(x):

    free_surface = "NO"
    equivalent_area = "NO"
    rotating_frame = "NO"
    output_format = "TECPLOT"

# Read the optimization problem
    for line in open(options.filename):
        if "CONST_IEQ=" in line:
            CompleteIeqConst = line.split("=")[1]
            UnitaryIeqConst = CompleteIeqConst.split(",")
        elif "CONST_IEQ_SCALE=" in line:
            CompleteIeqConst_Scale = line.split("=")[1]
            UnitaryIeqConst_Scale = CompleteIeqConst_Scale.split(",")
        elif "CONST_IEQ_VALUE=" in line:
            CompleteIeqConst_Value = line.split("=")[1]
            UnitaryIeqConst_Value = CompleteIeqConst_Value.split(",")
        elif "CONST_IEQ_SIGN=" in line:
            CompleteIeqConst_Sign = line.split("=")[1]
            UnitaryIeqConst_Sign = CompleteIeqConst_Sign.split(",")
        elif "FREE_SURFACE=" in line:
            free_surface = line.split("=")[1].strip()
        elif "EQUIV_AREA=" in line:
            equivalent_area = line.split("=")[1].strip()
        elif "ROTATING_FRAME=" in line:
            rotating_frame = line.split("=")[1].strip()  
        elif "OUTPUT_FORMAT=" in line:
            output_format = line.split("=")[1].strip()

    # Read the conts function gradient
    r = numpy.zeros([len(UnitaryIeqConst), nVar])
    
    for i in range(len(UnitaryIeqConst)) :
        if UnitaryIeqConst[i].strip() != "NONE" :
            # Update the computationa grid to the design variable state
            update_grid(x)

            # Change grid file for Gradient Computation over the deformed grid
            output_file = open(Config_GRAD_file,"w")
            for line in open(options.filename):
                if "CADJ_OBJFUNC=" in line:
                    output_file.write("CADJ_OBJFUNC= %s \n" % UnitaryIeqConst[i])
                elif "MESH_FILENAME=" in line:
                    output_file.write("MESH_FILENAME= %s \n" % Mesh_MDC_file)
                elif "CONV_FILENAME=" in line:
                    output_file.write("CONV_FILENAME= %s \n" % History_file.replace(".csv","").replace(".plt",""))
                else:
                    output_file.write(line)
            output_file.close()

            # Compute Gradient using continuous adjoint strategy
            if (options.gradient == "Adjoint"):
                if parallel: os.system ("continuous_adjoint.py -p %s -f %s" % (int(options.partitions), Config_GRAD_file))
                else: os.system ("continuous_adjoint.py -f " + Config_GRAD_file)

            if (options.gradient == "Adjoint"):
                if (output_format == "PARAVIEW"):
                    cadj_grad = open("cont_adj_" + Config_GRAD_file.replace(".cfg",".csv"),"r")
                elif (output_format == "TECPLOT"):
                    cadj_grad = open("cont_adj_" + Config_GRAD_file.replace(".cfg",".plt"),"r")
            elif (options.gradient == "FinDiff"):
                if (output_format == "PARAVIEW"):
                    cadj_grad = open("fin_diff_" + Config_GRAD_file.replace(".cfg",".csv"),"r")
                elif (output_format == "TECPLOT"):
                    cadj_grad = open("fin_diff_" + Config_GRAD_file.replace(".cfg",".plt"),"r")

            sign = 1.0
            if UnitaryIeqConst_Sign[i].strip() == "GREATER": sign = 1.0
            if UnitaryIeqConst_Sign[i].strip() == "LESS": sign = -1.0
                
            line = cadj_grad.readline()
            iVar = 0
            for line in cadj_grad:
                if (options.gradient == "Adjoint"):
                    r[i][iVar] = sign*float(line.split(",")[1])*float(UnitaryIeqConst_Scale[i])
                elif (options.gradient == "FinDiff"):
                    if UnitaryIeqConst[i] == "LIFT": r[iVar] = sign*float(line.split(",")[1])*float(UnitaryIeqConst_Scale[i])
                    elif UnitaryIeqConst[i] == "DRAG": r[iVar] = sign*float(line.split(",")[2])*float(UnitaryIeqConst_Scale[i])
                    elif UnitaryIeqConst[i] == "SIDEFORCE": r[iVar] = sign*float(line.split(",")[3])*float(UnitaryIeqConst_Scale[i])
                    elif UnitaryIeqConst[i] == "MOMENT_X": r[iVar] = sign*float(line.split(",")[4])*float(UnitaryIeqConst_Scale[i])
                    elif UnitaryIeqConst[i] == "MOMENT_Y": r[iVar] = sign*float(line.split(",")[5])*float(UnitaryIeqConst_Scale[i])
                    elif UnitaryIeqConst[i] == "MOMENT_Z": r[iVar] = sign*float(line.split(",")[6])*float(UnitaryIeqConst_Scale[i])
                    elif UnitaryIeqConst[i] == "FORCE_X": r[iVar] = sign*float(line.split(",")[7])*float(UnitaryIeqConst_Scale[i])
                    elif UnitaryIeqConst[i] == "FORCE_Y": r[iVar] = sign*float(line.split(",")[8])*float(UnitaryIeqConst_Scale[i])
                    elif UnitaryIeqConst[i] == "FORCE_Z": r[iVar] = sign*float(line.split(",")[9])*float(UnitaryIeqConst_Scale[i])
                    elif UnitaryIeqConst[i] == "EFFICIENCY": r[iVar] = sign*float(line.split(",")[10])*float(UnitaryIeqConst_Scale[i])
                    elif UnitaryIeqConst[i] == "EQUIVALENT_AREA": r[iVar] = sign*float(line.split(",")[11])*float(UnitaryIeqConst_Scale[i])
                    elif UnitaryIeqConst[i] == "NEARFIELD_PRESSURE": r[iVar] = sign*float(line.split(",")[12])*float(UnitaryIeqConst_Scale[i])
                    elif UnitaryIeqConst[i] == "FIGURE_OF_MERIT": r[iVar] = sign*float(line.split(",")[11])*float(UnitaryIeqConst_Scale[i])
                    elif UnitaryIeqConst[i] == "THRUST": r[iVar] = sign*float(line.split(",")[12])*float(UnitaryIeqConst_Scale[i])
                    elif UnitaryIeqConst[i] == "TORQUE": r[iVar] = sign*float(line.split(",")[13])*float(UnitaryIeqConst_Scale[i])
                    elif UnitaryIeqConst[i] == "FREESURFACE": r[iVar] = sign*float(line.split(",")[11])*float(UnitaryIeqConst_Scale[i])
                iVar = iVar + 1
            cadj_grad.close()

            # Remove configuration file        
            os.remove(Config_GRAD_file)

        else:
            r = numpy.zeros([0, nVar])
        
    return r


# Create a copy of the original grid
for line in open(options.filename):
    if "MESH_FILENAME=" in line:
        original_grid = line.split("=")[1].strip()

output_file = open(Mesh_MDC_file,"w")
for line in open(original_grid):
    output_file.write(line)
output_file.close()

x0 = numpy.zeros(nVar)
    
fmin_slsqp(f, x0, f_eqcons=eqcons, f_ieqcons=ieqcons, bounds=[], fprime=df, fprime_eqcons=deqcons, fprime_ieqcons=dieqcons, args=(), iter=100, acc=1e-12, iprint=2, full_output=0, epsilon=1.0e-10)
