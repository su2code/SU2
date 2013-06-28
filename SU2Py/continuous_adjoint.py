#!/usr/bin/python

## \file continuous_adjoint.py
#  \brief Python script for doing the continuous adjoint computation using the SU2 suite.
#  \author Current Development: Stanford University.
#          Original Structure: CADES 1.0 (2009).
#  \version 1.0.
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


import os, time, sys
from optparse import OptionParser

compute_adj = True
compute_flow = True

parser=OptionParser()
parser.add_option("-f", "--file", dest="filename",
                  help="read config from FILE", metavar="FILE")
parser.add_option("-p", "--partitions", dest="partitions", default=1,
                  help="number of PARTITIONS", metavar="PARTITIONS")
parser.add_option("-s", "--step", dest="step", default=1E-4,
                  help="finite difference STEP", metavar="STEP")

(options, args)=parser.parse_args()

if int(options.partitions) == 1:
    parallel = False
else : parallel = True

FinDiff_Step = 1E-4
FinDiff_Step = float(options.step)
output_format = "PARAVIEW"

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
        # Be careful with the remplace option... two digitd must be analized before one
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
    elif "OUTPUT_FORMAT" in line:
        output_format = line.split("=")[1].strip()

if (output_format == "PARAVIEW"):
    History_file = "history_" + options.filename.replace(".cfg",".csv")
    gradient_file = "cont_adj_" + options.filename.replace(".cfg",".csv")
elif (output_format == "TECPLOT"):
    History_file = "history_" + options.filename.replace(".cfg",".plt")
    gradient_file = "cont_adj_" + options.filename.replace(".cfg",".plt")
objfunc_grad_file = "objfunc_grad_adj_" + options.filename.replace(".cfg",".dat")
Config_GPC_file = "config_GPC_" + options.filename
Config_CFD_file = "config_CFD_" + options.filename

# Create gradient file, for the header it allways take the kind of the first design variable
Grad_file = open(gradient_file,"w")
if (output_format == "TECPLOT"): Grad_file.write( "VARIABLES=" )
if KindDVInfo[0][0] == "HICKS_HENNE": Grad_file.write( "\"FinDiff_Step\",\"Up/Down\",\"Loc_Max\",\"Gradient\"\n" )
if KindDVInfo[0][0] == "NACA_4DIGITS": Grad_file.write( "\"FinDiff_Step\",\"1st_digit\",\"2nd_digit\",\"3rd&4th_digits\",\"Gradient\"\n" )
if KindDVInfo[0][0] == "DISPLACEMENT": Grad_file.write( "\"FinDiff_Step\",\"x_Disp\",\"y_Disp\",\"z_Disp\",\"Gradient\"\n" )
if KindDVInfo[0][0] == "ROTATION": Grad_file.write( "\"FinDiff_Step\",\"x_Orig\",\"y_Orig\",\"z_Orig\",\"x_End\",\"y_End\",\"z_End\",\"Gradient\"\n" )
if KindDVInfo[0][0] == "FFD_CONTROL_POINT": Grad_file.write( "\"FinDiff_Step\",\"Chunk\",\"xIndex\",\"yIndex\",\"zIndex\",\"xAxis\",\"yAxis\",\"zAxis\",\"Gradient\"\n" )
if KindDVInfo[0][0] == "FFD_DIHEDRAL_ANGLE": Grad_file.write( "\"FinDiff_Step\",\"Chunk\",\"x_Orig\",\"y_Orig\",\"z_Orig\",\"x_End\",\"y_End\",\"z_End\",\"Gradient\"\n" )
if KindDVInfo[0][0] == "FFD_TWIST_ANGLE": Grad_file.write( "\"FinDiff_Step\",\"Chunk\",\"x_Orig\",\"y_Orig\",\"z_Orig\",\"x_End\",\"y_End\",\"z_End\",\"Gradient\"\n" )
if KindDVInfo[0][0] == "FFD_ROTATION": Grad_file.write( "\"FinDiff_Step\",\"Chunk\",\"x_Orig\",\"y_Orig\",\"z_Orig\",\"x_End\",\"y_End\",\"z_End\",\"Gradient\"\n" )
if KindDVInfo[0][0] == "FFD_CAMBER": Grad_file.write( "\"FinDiff_Step\",\"Chunk\",\"xIndex\",\"yIndex\",\"Gradient\"\n" )
if KindDVInfo[0][0] == "FFD_THICKNESS": Grad_file.write( "\"FinDiff_Step\",\"Chunk\",\"xIndex\",\"yIndex\",\"Gradient\"\n" )
if KindDVInfo[0][0] == "FFD_VOLUME": Grad_file.write( "\"FinDiff_Step\",\"Chunk\",\"xIndex\",\"yIndex\",\"Gradient\"\n" )
if KindDVInfo[0][0] == "MACH_NUMBER": Grad_file.write( "\"FinDiff_Step\",\"Gradient\"\n" )
if KindDVInfo[0][0] == "AOA": Grad_file.write( "\"FinDiff_Step\",\"Gradient\"\n" )
Grad_file.close()

# Change the parameters to do a flow simulation
output_file = open(Config_CFD_file,"w")
for line in open(options.filename):
    if "MATH_PROBLEM=" in line:
        output_file.write("MATH_PROBLEM= DIRECT \n")
    elif "RESTART_FLOW_FILENAME=" in line:
        restart_flow_file = line.split("=")[1].strip()
        output_file.write(line)
    else:
        output_file.write(line)
output_file.close()

if compute_flow:
    if parallel: os.system ("./parallel_computation.py -p %s -f %s -o False" % (int(options.partitions), Config_CFD_file))
    else: os.system ("./SU2_CFD " + Config_CFD_file)

# Change the parameters to do an adjoint simulation
output_file = open(Config_CFD_file,"w")
for line in open(options.filename):
    if "MATH_PROBLEM=" in line:
        output_file.write("MATH_PROBLEM= ADJOINT \n")
    elif "SOLUTION_FLOW_FILENAME=" in line:
        output_file.write("SOLUTION_FLOW_FILENAME= %s \n" % restart_flow_file)
    elif "CONV_FILENAME=" in line:
        if (output_format == "PARAVIEW"): output_file.write("CONV_FILENAME= %s \n" % History_file.replace(".csv",""))
        elif (output_format == "TECPLOT"): output_file.write("CONV_FILENAME= %s \n" % History_file.replace(".plt",""))
    else:
        output_file.write(line)
output_file.close()

if compute_adj:
    if parallel: os.system ("./parallel_computation.py -p %s -f %s -d False -o False" % (int(options.partitions), Config_CFD_file))
    else: os.system ("./SU2_CFD " + Config_CFD_file)

# Change the parameters of the design variables
for iDesignVar in range(len(DesignVar)):
    if KindDVInfo[iDesignVar][0] != "MACH_NUMBER" and  KindDVInfo[iDesignVar][0] != "AOA":
        output_file = open(Config_GPC_file,"w")
        for line in open(options.filename):
            if "DV_KIND=" in line:
                output_file.write("DV_KIND= "+ KindDVInfo[iDesignVar][0] + "\n")
            elif "DV_MARKER=" in line:
                if len(MarkerInfo[iDesignVar]) == 1:
                   output_file.write( "DV_MARKER= ( %s ) \n" % MarkerInfo[iDesignVar][0])
                else:
                    for j in range(len(MarkerInfo[iDesignVar])):
                        if j == 0:  output_file.write( "DV_MARKER= ( %s, " % MarkerInfo[iDesignVar][j])
                        if j > 0 and j < len(MarkerInfo[iDesignVar])-1: output_file.write( "%s, " % MarkerInfo[iDesignVar][j])
                        if j == len(MarkerInfo[iDesignVar])-1:  output_file.write( "%s ) \n" % MarkerInfo[iDesignVar][j])
            elif "DV_PARAM=" in line:
                for j in range(len(DesignVar[iDesignVar])):
                    if j == 0:  output_file.write( "DV_PARAM= ( %s, " % DesignVar[iDesignVar][j])
                    if j > 0 and j < len(DesignVar[iDesignVar])-1: output_file.write( "%s, " % DesignVar[iDesignVar][j])
                    if j == len(DesignVar[iDesignVar])-1:  output_file.write( "%s ) \n" % DesignVar[iDesignVar][j])
            elif "DV_VALUE_OLD=" in line:
                output_file.write("DV_VALUE_OLD= 0.0 \n")
            elif "DV_VALUE_NEW=" in line:
                output_file.write("DV_VALUE_NEW= %s \n" % FinDiff_Step)
            elif "GRAD_OBJFUNC_FILENAME=" in line:
                output_file.write("GRAD_OBJFUNC_FILENAME= %s \n" % objfunc_grad_file)                   
            else:
                output_file.write(line)
        output_file.close()

# Compute gradient with continuous adjoint
        os.system ("./SU2_GPC " + Config_GPC_file)

# Read gradient from file
        grad = open(objfunc_grad_file,"r")
        line = grad.readline()
        line = grad.readline()        
        Gradient = float(line)
        grad.close()

# Output gradient
        Grad_file = open(gradient_file,"a")
        for j in range(len(DesignVar[iDesignVar])):
            if j == 0:  Grad_file.write( "%f, %s, " % (FinDiff_Step, DesignVar[iDesignVar][j]))
            if j > 0 and j < len(DesignVar[iDesignVar])-1: Grad_file.write( "%s, " % DesignVar[iDesignVar][j])
            if j == len(DesignVar[iDesignVar])-1:  Grad_file.write( "%s, %.10f\n" % (DesignVar[iDesignVar][j], Gradient))
        Grad_file.close()
        
# Remove configuration and mesh files            
        os.remove(Config_GPC_file)
        os.remove(objfunc_grad_file)

    else:
# The design variable does not imply surface movement
        for line in open(History_file): 
            if "ERROR" in line: error = 0.0
        Gradient = float(line.split(",")[2])

# Output gradient
        Grad_file = open(gradient_file,"a")
        Grad_file.write( "%f, %.10f, " % (FinDiff_Step, Gradient))
        Grad_file.close()


os.remove(Config_CFD_file)
os.remove(History_file)
