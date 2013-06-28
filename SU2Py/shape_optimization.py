#!/usr/bin/python

## \file shape_optimization.py
#  \brief Python script for doing a the shape optimization.
#  \author Current Development: Stanford University.
#          Original Structure: CADES 1.0 (2009).
#  \version 1.0
#  /usr/local/epd-6.3-2-rh5-x86_64/bin/python
#  /usr/bin/python
#  /Library/Frameworks/Python.framework/Versions/2.6/bin/python
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
import time, os, subprocess, numpy

mesh_adapt = False
counter = 0
output_format = "PARAVIEW"

parser=OptionParser()
parser.add_option("-f", "--file", dest="filename",
                  help="read config from FILE", metavar="FILE")
parser.add_option("-p", "--partitions", dest="partitions", default=1,
                  help="number of PARTITIONS", metavar="PARTITIONS")
parser.add_option("-g", "--gradient_scale", dest="gradient_scale", default=0.0001,
                  help="value of the of GRADIENT_SCALE", metavar="GRADIENT_SCALE")
parser.add_option("-c", "--constraints", dest="constraints", default="False",
                  help="optimization with CONSTRAINTS", metavar="CONSTRAINTS")

(options, args)=parser.parse_args()

gradient_scale =  float(options.gradient_scale)
if int(options.partitions) == 1:
    parallel = False
else : parallel = True

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
    elif "OUTPUT_FORMAT" in line:
                output_format = line.split("=")[1].strip()

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
if (output_format == "PARAVIEW"):
    optimization_file = "optimization_" + options.filename.replace(".cfg",".csv")
elif (output_format == "TECPLOT"):
    optimization_file = "optimization_" + options.filename.replace(".cfg",".plt")

# Create optimization history file
Opt_file = open(optimization_file,"w")
if (output_format == "PARAVIEW"):
    Opt_file.write( "\"Iteration\",\"CLift\",\"CDrag\",\"CSideForce\",\"CPress\",\"CMx\",\"CMy\",\"CMz\",\"CEff\",\"CEquivArea\",\"CNearField\"\n" )
elif (output_format == "TECPLOT"):
    Opt_file.write( "TITLE = \"SU2 Optimization\"\n" )
    Opt_file.write( "VARIABLES = \"Iteration\",\"CLift\",\"CDrag\",\"CSideForce\",\"CPress\",\"CMx\",\"CMy\",\"CMz\",\"CEff\",\"CEquivArea\",\"CNearField\"\n" )
    Opt_file.write( "ZONE T= \"Optimization history\"\n" )

Opt_file.close()

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
    os.system ("./SU2_MDC " + Config_MDC_file)

# Copy solution
    for iDesignVar in range(len(DesignVar)):
        x_old[iDesignVar] = x_new[iDesignVar]

#  Adapt once the computational grid (OF is DRAG)
    if mesh_adapt:
        os.system ("./mesh_adaptation.py -o True -f " + Config_MDC_file)

# Remove configuration files                        
    os.remove(Config_MDC_file)

def f(x):
    print ("-------> SUBROTUTINE f(x) <---------") 

# Read the optimization problem
    for line in open(options.filename):
        if "OBJFUNC=" in line:
            objfunc = line.split("=")[1].strip()
            
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
        elif "OUTPUT_FORMAT" in line:
            output_format = line.split("=")[1].strip()
            output_file.write(line)
        else: 
            output_file.write(line)
    output_file.close()
    
# CFD computation to evaluate Drag and Lift Coefficient
    global counter    
    counter = counter + 1
    if parallel: os.system ("./parallel_computation.py -p %s -f %s -o True" % (int(options.partitions), Config_CFD_file))
    else: os.system ("./SU2_CFD " + Config_CFD_file)
    
    if KindDVInfo[0][0] == "HICKS_HENNE" or KindDVInfo[0][0] == "HICKS_HENNE_NORMAL" or KindDVInfo[0][0] == "HICKS_HENNE_SHOCK" or KindDVInfo[0][0] == "NACA_4DIGITS" or KindDVInfo[0][0] == "DISPLACEMENT" or KindDVInfo[0][0] == "ROTATION": 
        os.rename("surface_flow.csv","surface_flow_ShapeIter_"+str(counter)+".csv")
    if KindDVInfo[0][0] == "FFD_CONTROL_POINT" or KindDVInfo[0][0] == "FFD_DIHEDRAL_ANGLE" or KindDVInfo[0][0] == "FFD_TWIST_ANGLE" or KindDVInfo[0][0] == "FFD_ROTATION" or KindDVInfo[0][0] == "FFD_CAMBER" or KindDVInfo[0][0] == "FFD_THICKNESS" or KindDVInfo[0][0] == "FFD_VOLUME":
        if (output_format == "PARAVIEW"):
            os.rename("deformed_chunk.vtk","FFDBox_ShapeIter_"+str(counter)+".vtk")
        if (output_format == "TECPLOT"):
            os.rename("deformed_chunk.plt","FFDBox_ShapeIter_"+str(counter)+".plt")
 #       os.rename("deformed_0chunk.vtk","chunk0_ShapeIter_"+str(counter)+".vtk")
 #       os.rename("deformed_1chunk.vtk","chunk1_ShapeIter_"+str(counter)+".vtk")
 #       os.rename("deformed_2chunk.vtk","chunk2_ShapeIter_"+str(counter)+".vtk")

        

    if (output_format == "PARAVIEW"):
        os.rename("surface_flow.vtk","surface_flow_ShapeIter_"+str(counter)+".vtk")
        os.rename("flow.vtk","flow_ShapeIter_"+str(counter)+".vtk")
    elif (output_format == "TECPLOT"):
        os.rename("surface_flow.plt","surface_flow_ShapeIter_"+str(counter)+".plt")
        os.rename("flow.plt","flow_ShapeIter_"+str(counter)+".plt")
        
    if objfunc == "EQUIVALENT_AREA":
        os.rename("equiv_area.csv","equiv_area_ShapeIter_"+str(counter)+".csv")

# Read the Drag, Lift and Pitching Coefficient    
    for line in open(History_file): 
        if "ERROR" in line: error = 0.0
    CLift_Value = float(line.split(",")[1])
    CDrag_Value = float (line.split(",")[2])
    CSideForce_Value = float (line.split(",")[3])
    CPress_Value = float (line.split(",")[4])
    CMx_Value = float (line.split(",")[5])
    CMy_Value = float (line.split(",")[6])
    CMz_Value = float (line.split(",")[7])
    CEff_Value = float (line.split(",")[8])
    CEquivArea_Value = float (line.split(",")[9])
    CNearField_Value = float (line.split(",")[10])

# Remove configuration and mesh files            
    os.remove(Config_CFD_file)
    os.remove(History_file)

# Output optimization
    Opt_file = open(optimization_file,"a")
    Opt_file.write( "%d, \t %.10f, \t %.10f, \t %.10f, \t %.10f, \t %.10f, \t %.10f, \t %.10f, \t %.10f, \t %.10f, \t %.10f\n" % (counter, CLift_Value, CDrag_Value, CSideForce_Value, CPress_Value, CMx_Value, CMy_Value, CMz_Value, CEff_Value, CEquivArea_Value, CNearField_Value) )
    Opt_file.close()

    if objfunc == "DRAG":
        return CDrag_Value*gradient_scale

    if objfunc == "LIFT":
        return CLift_Value*gradient_scale

    if objfunc == "SIDEFORCE":
        return CSideForce_Value*gradient_scale
    
    if objfunc == "EFFICIENCY":
        return CEff_Value*gradient_scale

    if objfunc == "EQUIVALENT_AREA":
        return CEquivArea_Value*gradient_scale

    if objfunc == "NEARFIELD_PRESSURE":
        return CNearField_Value*gradient_scale

def eqcons(x):
    print ("-------> SUBROTUTINE eqcons(x) <---------") 

    r = numpy.zeros(1)
    r[0] = 0.0
    return r

def deqcons(x):
    print ("-------> SUBROTUTINE deqcons(x) <---------") 

    r = numpy.zeros(nVar)
    for iVar in range(nVar):
        r[iVar] = 1E-6
    return r

def ieqcons(x):
    print ("-------> SUBROTUTINE ieqcons(x) <---------") 

# Read the optimization problem
    for line in open(options.filename):
        if "CONSTRAINT=" in line:
            constraint = line.split("=")[1].strip()
        elif "MIN_LIFT=" in line:
            min_lift = float(line.split("=")[1].strip())
        elif "MIN_PITCH=" in line:
            min_pitch = float(line.split("=")[1].strip())
        elif "MAX_PITCH=" in line:
            max_pitch = float(line.split("=")[1].strip())

    if constraint != "NONE" :
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
        if parallel: os.system ("./parallel_computation.py -p %s -f %s -o False" % (int(options.partitions), Config_CFD_file))
        else: os.system ("./SU2_CFD " + Config_CFD_file)

# Read the Drag, Lift and Pitching Coefficient    
        for line in open(History_file): 
            if "ERROR" in line: error = 0.0
        CLift_Value = float(line.split(",")[1])
        CDrag_Value = float (line.split(",")[2])
        CSideForce_Value = float (line.split(",")[3])
        CPress_Value = float (line.split(",")[4])
        CMx_Value = float (line.split(",")[5])
        CMy_Value = float (line.split(",")[6])
        CMz_Value = float (line.split(",")[7])
        CEff_Value = float (line.split(",")[8])
        CEquivArea_Value = float (line.split(",")[9])
        CNearField_Value = float (line.split(",")[10])

# Remove configuration and mesh files            
        os.remove(Config_CFD_file)
        os.remove(History_file)

# Return the objective function
        r = numpy.zeros(1)
        r[0] = (CLift_Value-min_lift)*gradient_scale

    else:
        r = numpy.zeros(1)
        r[0] = +1.0

    return r

def dieqcons(x):
    print ("-------> SUBROTUTINE dieqcons(x) <---------") 

# Read the optimization problem
    for line in open(options.filename):
        if "CONSTRAINT=" in line:
            constraint = line.split("=")[1].strip()
        elif "MIN_LIFT=" in line:
            min_lift = float(line.split("=")[1].strip())
        elif "MIN_PITCH=" in line:
            min_pitch = float(line.split("=")[1].strip())
        elif "MAX_PITCH=" in line:
            max_pitch = float(line.split("=")[1].strip())

    if constraint != "NONE" :
# Update the computationa grid to the design variable state
        update_grid(x)

# Change grid file for Gradient Computation over the deformed grid
        output_file = open(Config_GRAD_file,"w")
        for line in open(options.filename):
            if "CADJ_OBJFUNC=" in line:
                output_file.write("CADJ_OBJFUNC= %s \n" % constraint)
            elif "MESH_FILENAME=" in line:
                output_file.write("MESH_FILENAME= %s \n" % Mesh_MDC_file)
            elif "CONV_FILENAME=" in line:
                output_file.write("CONV_FILENAME= %s \n" % History_file.replace(".csv","").replace(".plt",""))
            elif "OUTPUT_FORMAT" in line:
                output_format = line.split("=")[1].strip()
                output_file.write(line)
            else:
                output_file.write(line)
        output_file.close()
    
# Compute Gradient using continuous adjoint strategy
        if parallel: os.system ("./continuous_adjoint.py -p %s -f %s" % (int(options.partitions), Config_GRAD_file))
        else: os.system ("./continuous_adjoint.py -f " + Config_GRAD_file)

# Read the objective function gradient
        r = numpy.zeros(nVar)
        if (output_format == "PARAVIEW"):
            cadj_grad = open("cont_adj_" + Config_GRAD_file.replace(".cfg",".csv"),"r")
        elif (output_format == "TECPLOT"):
            cadj_grad = open("cont_adj_" + Config_GRAD_file.replace(".cfg",".plt"),"r")

        line = cadj_grad.readline()
        iVar = 0
        for line in cadj_grad:
            if KindDVInfo[0][0] == "HICKS_HENNE": r[iVar] = float(line.split(",")[3])*gradient_scale
            if KindDVInfo[0][0] == "HICKS_HENNE_NORMAL": r[iVar] = float(line.split(",")[3])*gradient_scale
            if KindDVInfo[0][0] == "HICKS_HENNE_SHOCK": r[iVar] = float(line.split(",")[3])*gradient_scale
            if KindDVInfo[0][0] == "NACA_4DIGITS": r[iVar] = float(line.split(",")[4])*gradient_scale
            if KindDVInfo[0][0] == "DISPLACEMENT": r[iVar] = float(line.split(",")[4])*gradient_scale
            if KindDVInfo[0][0] == "ROTATION": r[iVar] = float(line.split(",")[7])*gradient_scale
            if KindDVInfo[0][0] == "FFD_CONTROL_POINT": r[iVar] = float(line.split(",")[8])*gradient_scale
            if KindDVInfo[0][0] == "FFD_DIHEDRAL_ANGLE": r[iVar] = float(line.split(",")[8])*gradient_scale
            if KindDVInfo[0][0] == "FFD_TWIST_ANGLE": r[iVar] = float(line.split(",")[4])*gradient_scale
            if KindDVInfo[0][0] == "FFD_ROTATION": r[iVar] = float(line.split(",")[8])*gradient_scale
            if KindDVInfo[0][0] == "FFD_CAMBER": r[iVar] = float(line.split(",")[4])*gradient_scale
            if KindDVInfo[0][0] == "FFD_THICKNESS": r[iVar] = float(line.split(",")[4])*gradient_scale
            if KindDVInfo[0][0] == "FFD_VOLUME": r[iVar] = float(line.split(",")[4])*gradient_scale
            iVar = iVar + 1
        cadj_grad.close()

# Remove configuration file        
        os.remove(Config_GRAD_file)

    else:
        r = numpy.zeros(nVar)
        for iVar in range(nVar):
            r[iVar] = 1E-6

    return r


def df(x):
    print ("-------> SUBROTUTINE df(x) <---------") 

# Update the computationa grid to the design variable state
    update_grid(x)

# Read the optimization problem
    for line in open(options.filename):
        if "OBJFUNC=" in line:
            objfunc = line.split("=")[1].strip()

# Change grid file for Gradient Computation over the deformed grid
    output_file = open(Config_GRAD_file,"w")
    for line in open(options.filename):
        if "CADJ_OBJFUNC=" in line:
            output_file.write("CADJ_OBJFUNC= %s \n" % objfunc)
        elif "MESH_FILENAME=" in line:
            output_file.write("MESH_FILENAME= %s \n" % Mesh_MDC_file)
        elif "CONV_FILENAME=" in line:
            output_file.write("CONV_FILENAME= %s \n" % History_file.replace(".csv","").replace(".plt",""))
        elif "OUTPUT_FORMAT" in line:
            output_format = line.split("=")[1].strip()
            output_file.write(line)
        else:
            output_file.write(line)
    output_file.close()
    
# Compute Gradient using continuous adjoint strategy
    if parallel: os.system ("./continuous_adjoint.py -p %s -f %s" % (int(options.partitions), Config_GRAD_file))
    else: os.system ("./continuous_adjoint.py -f " + Config_GRAD_file)

# Read the objective function gradient
    r = numpy.zeros(nVar)
    if (output_format == "PARAVIEW"):
        cadj_grad = open("cont_adj_" + Config_GRAD_file.replace(".cfg",".csv"),"r")
    elif (output_format == "TECPLOT"):
        cadj_grad = open("cont_adj_" + Config_GRAD_file.replace(".cfg",".plt"),"r")
            
    line = cadj_grad.readline()
    iVar = 0
    for line in cadj_grad:
        if KindDVInfo[0][0] == "HICKS_HENNE": r[iVar] = float(line.split(",")[3])*gradient_scale
        if KindDVInfo[0][0] == "HICKS_HENNE_NORMAL": r[iVar] = float(line.split(",")[3])*gradient_scale
        if KindDVInfo[0][0] == "HICKS_HENNE_SHOCK": r[iVar] = float(line.split(",")[3])*gradient_scale
        if KindDVInfo[0][0] == "NACA_4DIGITS": r[iVar] = float(line.split(",")[4])*gradient_scale
        if KindDVInfo[0][0] == "DISPLACEMENT": r[iVar] = float(line.split(",")[4])*gradient_scale
        if KindDVInfo[0][0] == "ROTATION": r[iVar] = float(line.split(",")[7])*gradient_scale
        if KindDVInfo[0][0] == "FFD_CONTROL_POINT": r[iVar] = float(line.split(",")[8])*gradient_scale
        if KindDVInfo[0][0] == "FFD_DIHEDRAL_ANGLE": r[iVar] = float(line.split(",")[8])*gradient_scale
        if KindDVInfo[0][0] == "FFD_TWIST_ANGLE": r[iVar] = float(line.split(",")[4])*gradient_scale
        if KindDVInfo[0][0] == "FFD_ROTATION": r[iVar] = float(line.split(",")[8])*gradient_scale
        if KindDVInfo[0][0] == "FFD_CAMBER": r[iVar] = float(line.split(",")[4])*gradient_scale
        if KindDVInfo[0][0] == "FFD_THICKNESS": r[iVar] = float(line.split(",")[4])*gradient_scale
        if KindDVInfo[0][0] == "FFD_VOLUME": r[iVar] = float(line.split(",")[4])*gradient_scale
        iVar = iVar + 1

        
    cadj_grad.close()

# Remove configuration file        
    os.remove(Config_GRAD_file)
    
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

if options.constraints != "False":
    fmin_slsqp(f, x0, fprime=df, args=(), f_ieqcons=ieqcons, fprime_ieqcons=dieqcons, f_eqcons=eqcons, fprime_eqcons=deqcons, iprint=2, full_output=True, iter=100, acc=1E-10)
else : 
    fmin_bfgs(f, x0, fprime=df, args=(), gtol=1e-06, norm=numpy.inf, epsilon=None, maxiter=20, full_output=True, disp=True, retall=True, callback=None)

