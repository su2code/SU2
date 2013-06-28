#!/usr/bin/python

## \file finite_differences.py
#  \brief Python script for doing the finite differences computation using the SU2 suite.
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


import os, time
from optparse import OptionParser

parser = OptionParser()
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

val_FDStep = [1E-4]
val_FDStep[0] = float(options.step)

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
            KindDVInfo[i][0] = KindDVInfo[i][0].replace('1','HICKS_HENNE')
            KindDVInfo[i][0] = KindDVInfo[i][0].replace('4','NACA_4DIGITS')
            KindDVInfo[i][0] = KindDVInfo[i][0].replace('5','DISPLACEMENT')
            KindDVInfo[i][0] = KindDVInfo[i][0].replace('6','ROTATION')
            KindDVInfo[i][0] = KindDVInfo[i][0].replace('7','FFD_CONTROL_POINT')
            KindDVInfo[i][0] = KindDVInfo[i][0].replace('8','FFD_DIHEDRAL_ANGLE')
            KindDVInfo[i][0] = KindDVInfo[i][0].replace('9','FFD_TWIST_ANGLE')
            KindDVInfo[i][0] = KindDVInfo[i][0].replace('10','FFD_ROTATION')
            KindDVInfo[i][0] = KindDVInfo[i][0].replace('11','FFD_CAMBER')
            KindDVInfo[i][0] = KindDVInfo[i][0].replace('12','FFD_THICKNESS')
            KindDVInfo[i][0] = KindDVInfo[i][0].replace('101','MACH_NUMBER')
            KindDVInfo[i][0] = KindDVInfo[i][0].replace('102','AOA')
            KindDVInfo[i][1] = float(KindDVInfo[i][1])
            for j in range(len(DesignVar[i])) :
                DesignVar[i][j] = float(DesignVar[i][j])
    elif "OUTPUT_FORMAT" in line:
        output_format = line.split("=")[1].strip()

if (output_format == "PARAVIEW"):
    History_file = "history_" + options.filename.replace(".cfg",".csv")
    gradient_file = "fin_diff_" + options.filename.replace(".cfg",".csv")
elif (output_format == "TECPLOT"):
    History_file = "history_" + options.filename.replace(".cfg",".plt")
    gradient_file = "fin_diff_" + options.filename.replace(".cfg",".plt")
Config_MDC_file = "config_MDC_" + options.filename
Config_CFD_file = "config_CFD_" + options.filename
Mesh_MDC_file = "mesh_MDC_" + options.filename.replace(".cfg",".su2")

# Create gradient file
Grad_file = open(gradient_file,"w")
if (output_format == "TECPLOT"): Grad_file.write( "VARIABLES=" )
if KindDVInfo[0][0] == "HICKS_HENNE": Grad_file.write("\"FinDiff_Step\",\"Up/Down\",\"Loc_Max\" \
,\"Grad_CLift\",\"Grad_CDrag\",\"Grad_CSideForce\",\"Grad_CPress\",\"Grad_CMx\",\"Grad_CMy\",\"Grad_CMz\",\"Grad_CEff\",\"Grad_CEquivArea\",\"Grad_CNearField\" \
,\"CLift_Value\",\"CDrag_Value\",\"CSideForce_Value\",\"CPress_Value\",\"CMx_Value\",\"CMy_Value\",\"CMz_Value\",\"CEff_Value\",\"CEquivArea_Value\",\"CNearField_Value\" \
,\"Ini_CLift\",\"Ini_CDrag\",\"Ini_CSideForce\",\"Ini_CPress\",\"Ini_CMx\",\"Ini_CMy\",\"Ini_CMz\",\"Ini_CEff\",\"Ini_CEquivArea\",\"Ini_CNearField\"\n" )
if KindDVInfo[0][0] == "NACA_4DIGITS": Grad_file.write("\"FinDiff_Step\",\"1st_digit\",\"2nd_digit\",\"3rd&4th_digits\" \
,\"Grad_CLift\",\"Grad_CDrag\",\"Grad_CSideForce\",\"Grad_CPress\",\"Grad_CMx\",\"Grad_CMy\",\"Grad_CMz\",\"Grad_CEff\",\"Grad_CEquivArea\",\"Grad_CNearField\" \
,\"CLift_Value\",\"CDrag_Value\",\"CSideForce_Value\",\"CPress_Value\",\"CMx_Value\",\"CMy_Value\",\"CMz_Value\",\"CEff_Value\",\"CEquivArea_Value\",\"CNearField_Value\" \
,\"Ini_CLift\",\"Ini_CDrag\",\"Ini_CSideForce\",\"Ini_CPress\",\"Ini_CMx\",\"Ini_CMy\",\"Ini_CMz\",\"Ini_CEff\",\"Ini_CEquivArea\",\"Ini_CNearField\"\n" )
if KindDVInfo[0][0] == "DISPLACEMENT": Grad_file.write("\"FinDiff_Step\",\"x_Disp\",\"y_Disp\",\"z_Disp\" \
,\"Grad_CLift\",\"Grad_CDrag\",\"Grad_CSideForce\",\"Grad_CPress\",\"Grad_CMx\",\"Grad_CMy\",\"Grad_CMz\",\"Grad_CEff\",\"Grad_CEquivArea\",\"Grad_CNearField\" \
,\"CLift_Value\",\"CDrag_Value\",\"CSideForce_Value\",\"CPress_Value\",\"CMx_Value\",\"CMy_Value\",\"CMz_Value\",\"CEff_Value\",\"CEquivArea_Value\",\"CNearField_Value\" \
,\"Ini_CLift\",\"Ini_CDrag\",\"Ini_CSideForce\",\"Ini_CPress\",\"Ini_CMx\",\"Ini_CMy\",\"Ini_CMz\",\"Ini_CEff\",\"Ini_CEquivArea\",\"Ini_CNearField\"\n" )
if KindDVInfo[0][0] == "ROTATION": Grad_file.write("\"FinDiff_Step\",\"x_Orig\",\"y_Orig\",\"z_Orig\",\"x_End\",\"y_End\",\"z_End\" \
,\"Grad_CLift\",\"Grad_CDrag\",\"Grad_CSideForce\",\"Grad_CPress\",\"Grad_CMx\",\"Grad_CMy\",\"Grad_CMz\",\"Grad_CEff\",\"Grad_CEquivArea\",\"Grad_CNearField\" \
,\"CLift_Value\",\"CDrag_Value\",\"CSideForce_Value\",\"CPress_Value\",\"CMx_Value\",\"CMy_Value\",\"CMz_Value\",\"CEff_Value\",\"CEquivArea_Value\",\"CNearField_Value\" \
,\"Ini_CLift\",\"Ini_CDrag\",\"Ini_CSideForce\",\"Ini_CPress\",\"Ini_CMx\",\"Ini_CMy\",\"Ini_CMz\",\"Ini_CEff\",\"Ini_CEquivArea\",\"Ini_CNearField\"\n" )
if KindDVInfo[0][0] == "FFD_CONTROL_POINT": Grad_file.write("\"FinDiff_Step\",\"Chunk\",\"xIndex\",\"yIndex\",\"zIndex\",\"xAxis\",\"yAxis\",\"zAxis\" \
,\"Grad_CLift\",\"Grad_CDrag\",\"Grad_CSideForce\",\"Grad_CPress\",\"Grad_CMx\",\"Grad_CMy\",\"Grad_CMz\",\"Grad_CEff\",\"Grad_CEquivArea\",\"Grad_CNearField\" \
,\"CLift_Value\",\"CDrag_Value\",\"CSideForce_Value\",\"CPress_Value\",\"CMx_Value\",\"CMy_Value\",\"CMz_Value\",\"CEff_Value\",\"CEquivArea_Value\",\"CNearField_Value\" \
,\"Ini_CLift\",\"Ini_CDrag\",\"Ini_CSideForce\",\"Ini_CPress\",\"Ini_CMx\",\"Ini_CMy\",\"Ini_CMz\",\"Ini_CEff\",\"Ini_CEquivArea\",\"Ini_CNearField\"\n" )
if KindDVInfo[0][0] == "FFD_DIHEDRAL_ANGLE": Grad_file.write("\"FinDiff_Step\",\"Chunk\",\"x_Orig\",\"y_Orig\",\"z_Orig\",\"x_End\",\"y_End\",\"z_End\" \
,\"Grad_CLift\",\"Grad_CDrag\",\"Grad_CSideForce\",\"Grad_CPress\",\"Grad_CMx\",\"Grad_CMy\",\"Grad_CMz\",\"Grad_CEff\",\"Grad_CEquivArea\",\"Grad_CNearField\" \
,\"CLift_Value\",\"CDrag_Value\",\"CSideForce_Value\",\"CPress_Value\",\"CMx_Value\",\"CMy_Value\",\"CMz_Value\",\"CEff_Value\",\"CEquivArea_Value\",\"CNearField_Value\" \
,\"Ini_CLift\",\"Ini_CDrag\",\"Ini_CSideForce\",\"Ini_CPress\",\"Ini_CMx\",\"Ini_CMy\",\"Ini_CMz\",\"Ini_CEff\",\"Ini_CEquivArea\",\"Ini_CNearField\"\n" )
if KindDVInfo[0][0] == "FFD_TWIST_ANGLE": Grad_file.write("\"FinDiff_Step\",\"Chunk\",\"x_Orig\",\"y_Orig\",\"z_Orig\",\"x_End\",\"y_End\",\"z_End\" \
,\"Grad_CLift\",\"Grad_CDrag\",\"Grad_CSideForce\",\"Grad_CPress\",\"Grad_CMx\",\"Grad_CMy\",\"Grad_CMz\",\"Grad_CEff\",\"Grad_CEquivArea\",\"Grad_CNearField\" \
,\"CLift_Value\",\"CDrag_Value\",\"CSideForce_Value\",\"CPress_Value\",\"CMx_Value\",\"CMy_Value\",\"CMz_Value\",\"CEff_Value\",\"CNearField_Value\" \
,\"Ini_CLift\",\"Ini_CDrag\",\"Ini_CSideForce\",\"Ini_CPress\",\"Ini_CMx\",\"Ini_CMy\",\"Ini_CMz\",\"Ini_CEff\",\"Ini_CEquivArea\",\"Ini_CNearField\"\n" )
if KindDVInfo[0][0] == "FFD_ROTATION": Grad_file.write("\"FinDiff_Step\",\"Chunk\",\"x_Orig\",\"y_Orig\",\"z_Orig\",\"x_End\",\"y_End\",\"z_End\" \
,\"Grad_CLift\",\"Grad_CDrag\",\"Grad_CSideForce\",\"Grad_CPress\",\"Grad_CMx\",\"Grad_CMy\",\"Grad_CMz\",\"Grad_CEff\",\"Grad_CEquivArea\",\"Grad_CNearField\" \
,\"CLift_Value\",\"CDrag_Value\",\"CSideForce_Value\",\"CPress_Value\",\"CMx_Value\",\"CMy_Value\",\"CMz_Value\",\"CEff_Value\",\"CEquivArea_Value\",\"CNearField_Value\" \
,\"Ini_CLift\",\"Ini_CDrag\",\"Ini_CSideForce\",\"Ini_CPress\",\"Ini_CMx\",\"Ini_CMy\",\"Ini_CMz\",\"Ini_CEff\",\"Ini_CEquivArea\",\"Ini_CNearField\"\n" )
if KindDVInfo[0][0] == "FFD_CAMBER": Grad_file.write("\"FinDiff_Step\",\"Chunk\",\"xIndex\",\"yIndex\" \
,\"Grad_CLift\",\"Grad_CDrag\",\"Grad_CSideForce\",\"Grad_CPress\",\"Grad_CMx\",\"Grad_CMy\",\"Grad_CMz\",\"Grad_CEff\",\"Grad_CEquivArea\",\"Grad_CNearField\" \
,\"CLift_Value\",\"CDrag_Value\",\"CSideForce_Value\",\"CPress_Value\",\"CMx_Value\",\"CMy_Value\",\"CMz_Value\",\"CEff_Value\",\"CEquivArea_Value\",\"CNearField_Value\" \
,\"Ini_CLift\",\"Ini_CDrag\",\"Ini_CSideForce\",\"Ini_CPress\",\"Ini_CMx\",\"Ini_CMy\",\"Ini_CMz\",\"Ini_CEff\",\"Ini_CEquivArea\",\"Ini_CNearField\"\n" )
if KindDVInfo[0][0] == "FFD_THICKNESS": Grad_file.write("\"FinDiff_Step\",\"Chunk\",\"xIndex\",\"yIndex\" \
,\"Grad_CLift\",\"Grad_CDrag\",\"Grad_CSideForce\",\"Grad_CPress\",\"Grad_CMx\",\"Grad_CMy\",\"Grad_CMz\",\"Grad_CEff\",\"Grad_CEquivArea\",\"Grad_CNearField\" \
,\"CLift_Value\",\"CDrag_Value\",\"CSideForce_Value\",\"CPress_Value\",\"CMx_Value\",\"CMy_Value\",\"CMz_Value\",\"CEff_Value\",\"CEquivArea_Value\",\"CNearField_Value\" \
,\"Ini_CLift\",\"Ini_CDrag\",\"Ini_CSideForce\",\"Ini_CPress\",\"Ini_CMx\",\"Ini_CMy\",\"Ini_CMz\",\"Ini_CEff\",\"Ini_CEquivArea\",\"Ini_CNearField\"\n" )
if KindDVInfo[0][0] == "MACH_NUMBER": Grad_file.write("\"FinDiff_Step\" \
,\"Grad_CLift\",\"Grad_CDrag\",\"Grad_CSideForce\",\"Grad_CPress\",\"Grad_CMx\",\"Grad_CMy\",\"Grad_CMz\",\"Grad_CEff\",\"Grad_CEquivArea\",\"Grad_CNearField\" \
,\"CLift_Value\",\"CDrag_Value\",\"CSideForce_Value\",\"CPress_Value\",\"CMx_Value\",\"CMy_Value\",\"CMz_Value\",\"CEff_Value\",\"CEquivArea_Value\",\"CNearField_Value\" \
,\"Ini_CLift\",\"Ini_CDrag\",\"Ini_CSideForce\",\"Ini_CPress\",\"Ini_CMx\",\"Ini_CMy\",\"Ini_CMz\",\"Ini_CEff\",\"Ini_CEquivArea\",\"Ini_CNearField\"\n" )
if KindDVInfo[0][0] == "AOA": Grad_file.write("\"FinDiff_Step\" \
,\"Grad_CLift\",\"Grad_CDrag\",\"Grad_CSideForce\",\"Grad_CPress\",\"Grad_CMx\",\"Grad_CMy\",\"Grad_CMz\",\"Grad_CEff\",\"Grad_CEquivArea\",\"Grad_CNearField\" \
,\"CLift_Value\",\"CDrag_Value\",\"CSideForce_Value\",\"CPress_Value\",\"CMx_Value\",\"CMy_Value\",\"CMz_Value\",\"CEff_Value\",\"CNearField_Value\" \
,\"Ini_CLift\",\"Ini_CDrag\",\"Ini_CSideForce\",\"Ini_CPress\",\"Ini_CMx\",\"Ini_CMy\",\"Ini_CMz\",\"Ini_CEff\",\"Ini_CEquivArea\",\"Ini_CNearField\"\n" )
Grad_file.close()

# Initial computation to evaluate adiensional coefficient
output_file = open(Config_CFD_file,"w")
for line in open(options.filename):
    if "CONV_FILENAME=" in line:
        if (output_format == "PARAVIEW"): output_file.write("CONV_FILENAME= %s \n" % History_file.replace(".csv",""))
        elif (output_format == "TECPLOT"): output_file.write("CONV_FILENAME= %s \n" % History_file.replace(".plt",""))
    elif "MATH_PROBLEM=" in line:
        output_file.write("MATH_PROBLEM= DIRECT \n")
    else:
        output_file.write(line)
output_file.close()
if parallel: os.system ("./parallel_computation.py -p %s -f %s" % (int(options.partitions), Config_CFD_file))
else: os.system ("./SU2_CFD " + Config_CFD_file)

for line in open(History_file):
    if "ERROR" in line: error = 0.0
Ini_CLift = float (line.split(",")[1])
Ini_CDrag = float (line.split(",")[2])
Ini_CSideForce = float (line.split(",")[3])
Ini_CPress = float (line.split(",")[4])
Ini_CMx = float (line.split(",")[5])
Ini_CMy = float (line.split(",")[6])
Ini_CMz = float (line.split(",")[7])
Ini_CEff = float (line.split(",")[8])
Ini_CEquivArea = float (line.split(",")[9])
Ini_CNearField = float (line.split(",")[10])

for iDesignVar in range(len(DesignVar)):
    for FDStep in val_FDStep:
# Change the parameters of the design variables
        output_file = open(Config_MDC_file,"w")
        for line in open(options.filename):
            if "DV_KIND=" in line:
                if KindDVInfo[iDesignVar][0] == "MACH_NUMBER": output_file.write("DV_KIND= NO_DEFORMATION \n")
                elif KindDVInfo[iDesignVar][0] == "AOA": output_file.write("DV_KIND= NO_DEFORMATION \n")
                else: output_file.write("DV_KIND= "+ KindDVInfo[iDesignVar][0] + "\n")
            elif "DV_MARKER=" in line:
                for j in range(len(MarkerInfo[iDesignVar])):
                    if j == 0:  output_file.write( "DV_MARKER= ( %s, " % MarkerInfo[iDesignVar][j])
                    if j > 0 and j < len(MarkerInfo[iDesignVar])-1: output_file.write( "%s, " % MarkerInfo[iDesignVar][j])
                    if j == len(MarkerInfo[iDesignVar])-1:  output_file.write( "%s ) \n" % MarkerInfo[iDesignVar][j])
            elif "DV_PARAM=" in line:
                for j in range(len(DesignVar[iDesignVar])):
                    if j == 0:  output_file.write( "DV_PARAM= ( %s, " % DesignVar[iDesignVar][j])
                    if j > 0 and j < len(DesignVar[iDesignVar])-1: output_file.write( "%s, " % DesignVar[iDesignVar][j])
                    if j == len(DesignVar[iDesignVar])-1:  output_file.write( "%s ) \n" % DesignVar[iDesignVar][j])
            elif "MACH_NUMBER=" in line:
                if KindDVInfo[iDesignVar][0] == "MACH_NUMBER": mach_number = float(line.split("=")[1].strip())
                output_file.write(line)
            elif "AoA=" in line:
                if KindDVInfo[iDesignVar][0] == "AoA": aoa = float(line.split("=")[1].strip())
                output_file.write(line)
            elif "DV_VALUE_OLD=" in line:
                output_file.write("DV_VALUE_OLD= 0.0 \n")
            elif "DV_VALUE_NEW=" in line:
                output_file.write("DV_VALUE_NEW= %s \n" % FDStep)
            elif "MESH_OUT_FILENAME=" in line:
                output_file.write("MESH_OUT_FILENAME= %s \n" % Mesh_MDC_file)                       
            else:
                output_file.write(line)
        output_file.close()
            
# Apply grid deformation
        os.system ("./SU2_MDC " + Config_MDC_file)

# Change grid file for CFD Computation over the deformed grid
        output_file = open(Config_CFD_file,"w")
        for line in open(Config_MDC_file):
            if "MESH_FILENAME=" in line:
                output_file.write("MESH_FILENAME= %s \n" % Mesh_MDC_file)
            elif "MATH_PROBLEM=" in line:
                output_file.write("MATH_PROBLEM= DIRECT \n")
            elif "MACH_NUMBER=" in line:
                if KindDVInfo[iDesignVar][0] == "MACH_NUMBER": output_file.write("MACH_NUMBER= %f \n" % (mach_number + float(FDStep)) )
                output_file.write(line)
            elif "AOA=" in line:
                if KindDVInfo[iDesignVar][0] == "AOA": output_file.write("AOA= %f \n" % (aoa + float(FDStep)) )
                output_file.write(line)
            elif "CONV_FILENAME=" in line:
                if (output_format == "PARAVIEW"): output_file.write("CONV_FILENAME= %s \n" % History_file.replace(".csv",""))
                elif (output_format == "TECPLOT"): output_file.write("CONV_FILENAME= %s \n" % History_file.replace(".plt","")) 
            else:
                output_file.write(line)
        output_file.close()

# CFD computation to evaluate adimensional coefficient
        if parallel: os.system ("./parallel_computation.py -p %s -f %s" % (int(options.partitions), Config_CFD_file))
        else: os.system ("./SU2_CFD " + Config_CFD_file)

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

        Grad_CLift = (CLift_Value-Ini_CLift) / float(FDStep)
        Grad_CDrag = (CDrag_Value-Ini_CDrag) / float(FDStep)
        Grad_CSideForce = (CSideForce_Value-Ini_CSideForce) / float(FDStep)
        Grad_CPress = (CPress_Value-Ini_CPress) / float(FDStep)
        Grad_CMx = (CMx_Value-Ini_CMx) / float(FDStep)
        Grad_CMy = (CMy_Value-Ini_CMy) / float(FDStep)
        Grad_CMz = (CMz_Value-Ini_CMz) / float(FDStep)
        Grad_CEff = (CEff_Value-Ini_CEff) / float(FDStep)
        Grad_CEquivArea = (CEquivArea_Value-Ini_CEquivArea) / float(FDStep)
        Grad_CNearField = (CNearField_Value-Ini_CNearField) / float(FDStep)

# Output gradient
        Grad_file = open(gradient_file,"a")
        if len(DesignVar[iDesignVar]) == 1:
            Grad_file.write( "%.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f\n" % (FDStep, Grad_CLift, Grad_CDrag,  Grad_CSideForce,  Grad_CPress, Grad_CMx, Grad_CMy, Grad_CMz, Grad_CEff, Grad_CEquivArea, Grad_CNearField, CLift_Value, CDrag_Value, CSideForce_Value, CPress_Value, CMx_Value, CMy_Value, CMz_Value, CEff_Value, CEquivArea_Value, CNearField_Value, Ini_CLift, Ini_CDrag, Ini_CSideForce, Ini_CPress, Ini_CMx, Ini_CMy, Ini_CMz, Ini_CEff, Ini_CEquivArea, Ini_CNearField))
        else:
            for j in range(len(DesignVar[iDesignVar])):
                if j == 0:  Grad_file.write( "%.10f, %s, " % (FDStep, DesignVar[iDesignVar][j]))
                if j > 0 and j < len(DesignVar[iDesignVar])-1: Grad_file.write( "%s, " % DesignVar[iDesignVar][j])
                if j == len(DesignVar[iDesignVar])-1:  Grad_file.write( "%s, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f\n" % (DesignVar[iDesignVar][j], Grad_CLift, Grad_CDrag,  Grad_CSideForce,  Grad_CPress, Grad_CMx, Grad_CMy, Grad_CMz, Grad_CEff, Grad_CEquivArea, Grad_CNearField, CLift_Value, CDrag_Value, CSideForce_Value, CPress_Value, CMx_Value, CMy_Value, CMz_Value, CEff_Value, CEquivArea_Value, CNearField_Value, Ini_CLift, Ini_CDrag, Ini_CSideForce, Ini_CPress, Ini_CMx, Ini_CMy, Ini_CMz, Ini_CEff, Ini_CEquivArea, Ini_CNearField))
        Grad_file.close()
            
# Remove configuration and mesh files            
        os.remove(Config_CFD_file)
        os.remove(Config_MDC_file)
        os.remove(Mesh_MDC_file)
        os.remove(History_file)

