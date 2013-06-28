#!/usr/bin/python

## \file libSU2.py
#  \brief Support Functions for SU2 Python Scripts
#  \author Current Development: Stanford University.
#          Original Structure: CADES 1.0 (2009).
#  \version 1.1.
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



import os, time, sys, shutil, glob
import numpy



# -------------------------------------------------------------------
# Get SU2 Configuration Parameters
# -------------------------------------------------------------------

def Get_ConfigParams(filename):
      
    # initialize output dictionary
    data_dict = {}   
    
    # process each line
    for line in open(filename):   
        
        # remove line returns
        line = line.strip('\r\n')
        # make sure it has useful data
        if (not "=" in line) or (line[0] == '%'):
            continue
        # split across equals sign
        line = line.split("=")
        this_param = line[0].strip()
        this_value = line[1].strip()
        
        # handle parameter types
        for case in switch(this_param):
            
            # space delimited lists
            if case("DV_KIND"):
                data_dict[this_param] = this_value.split(" ")
                break
            
            # semicolon delimited lists of dvs
            if case("DV_PARAM"):
                # remove white space
                info_General = ''.join(this_value.split())
                # split by semicolon
                info_General = info_General.split(';')
                # build list of dv params, convert string to float
                dv_Parameters = []
                for this_dvParam in info_General:
                    this_dvParam = this_dvParam.strip('()')
                    this_dvParam = [ float(x) for x in this_dvParam.split(",") ]   
                    dv_Parameters = dv_Parameters + [this_dvParam]
                data_dict[this_param] = dv_Parameters
                break            

            # unitary design variable definition
            if case("DEFINITION_DV"):
                # remove white space
                this_value = ''.join(this_value.split())                
                # split into unitary definitions
                info_Unitary = this_value.split(";")
                # process each Design Variable
                dv_Kind       = []
                dv_Scale      = []
                dv_Markers    = []
                dv_Parameters = []
                for this_General in info_Unitary:
                    # split each unitary definition into one general definition
                    info_General = this_General.strip("()").split("|") # check for needed strip()?
                    # split information for dv Kinds
                    info_Kind    = info_General[0].split(",")
                    # pull processed dv values
                    this_dvKind       = get_dvKind( float( info_Kind[0] ) )     
                    this_dvScale      = float( info_Kind[1] )
                    this_dvMarkers    = info_General[1].split(",")
                    if this_dvKind=='MACH_NUMBER' or this_dvKind=='AOA':
                        this_dvParameters = []
                    else:
                        this_dvParameters = [ float(x) for x in info_General[2].split(",") ]                    
                    # add to lists
                    dv_Kind       = dv_Kind       + [this_dvKind]
                    dv_Scale      = dv_Scale      + [this_dvScale]
                    dv_Markers    = dv_Markers    + [this_dvMarkers]
                    dv_Parameters = dv_Parameters + [this_dvParameters]
                # store in a dictionary
                dv_Definitions = { 'Kind'       : dv_Kind       ,
                                   'Scale'      : dv_Scale      ,
                                   'Markers'    : dv_Markers    ,
                                   'Parameters' : dv_Parameters }
                # save to output dictionary
                data_dict[this_param] = dv_Definitions
                break
            
            # float parameters
            if case("MACH_NUMBER")            : pass
            if case("AoA")                    :             
                data_dict[this_param] = float(this_value)
                break                 
            
            # otherwise
            # string parameters
            if case():
                data_dict[this_param] = this_value
                break              
            
            #: if case DEFINITION_DV
                        
        #: for case
        
    #: for line
            
    return data_dict
    
#: def Get_ConfigParams()



# -------------------------------------------------------------------
# Set SU2 Configuration Parameters
# -------------------------------------------------------------------

def Set_ConfigParams(filename,param_dict):
    
    temp_filename = "temp.cfg"
    shutil.copy(filename,temp_filename)
    output_file = open(filename,"w")
    
    for raw_line in open(temp_filename):
        # remove line returns
        line = raw_line.strip('\r\n')
        
        # make sure it has useful data
        if not "=" in line:
            output_file.write(raw_line)
            continue
        
        # split across equals sign
        line = line.split("=")
        this_param = line[0].strip()
        old_value  = line[1].strip()
        
        # skip if parameter unwanted
        if not param_dict.has_key(this_param):
            output_file.write(raw_line)
            continue
        
        # start writing parameter
        new_value = param_dict[this_param] 
        output_file.write(this_param + "= ")
        
        # handle parameter types
        for case in switch(this_param):  
              
            # space delimited list of floats
            if case("DV_VALUE_NEW") : pass
            if case("DV_VALUE_OLD") :
                n_lists = len(new_value)
                for i_value in range(n_lists):
                    output_file.write("%s" % new_value[i_value])
                    if i_value+1 < n_lists:
                        output_file.write(", ")               
                break
            
            # space delimited list of strings
            if case("DV_KIND") :
                n_lists = len(new_value)
                for i_value in range(n_lists):
                    output_file.write(new_value[i_value])
                    if i_value+1 < n_lists:
                        output_file.write(", ")               
                break            
            
            # comma delimited list of strings
            if case("DV_MARKER") : 
                output_file.write("( ")
                n_lists = len(new_value)
                for i_value in range(n_lists):
                    output_file.write(new_value[i_value])
                    if i_value+1 < n_lists:
                        output_file.write(", ")
                output_file.write(" )") 
                break                
            
            # semicolon delimited lists of comma delimited lists
            if case("DV_PARAM") :
                for i_value in range(len(new_value)):
                    output_file.write("( ")
                    this_list = new_value[i_value]
                    n_lists = len(new_value[i_value])
                    for j_value in range(n_lists):
                        output_file.write("%.15f" % this_list[j_value])
                        if j_value+1 < n_lists:
                            output_file.write(", ")   
                    output_file.write(") ")
                    if i_value+1 < len(new_value):
                        output_file.write("; ")            
                break

            # float parameters
            if case("MACH_NUMBER")            : pass
            if case("AoA")                    :             
                output_file.write("%.15f" % new_value)
                break       
            
            # default, assume string, integer or unformatted float 
            if case():
                output_file.write('%s' % new_value)
                break                         
                
        #: for case
        
        # next line
        output_file.write("\n")        
        
    #: for each line
    
    output_file.close()
    os.remove( temp_filename )
    
#: def Set_ConfigParams()



# -------------------------------------------------------------------
# Read SU2_GPC Gradient Values
# -------------------------------------------------------------------

def get_GradientVals( Grad_filename , scale = 1.0):
        
    # open file and skip first line
    gradfile = open(Grad_filename)
    gradfile.readline()
    
    # read values
    grad_vals = []
    for line in gradfile:
        line = line.strip()
        if len(line) == 0:
            break
        grad_vals.append(float(line) * scale)
    #: for each line
    
    return grad_vals

#: def get_GradientVals()



# -------------------------------------------------------------------
# Get Objective Function Sign
# -------------------------------------------------------------------

def get_ObjFunSign( ObjFun_name ):
    
    # flip sign for maximization problems
    if ObjFun_name == "LIFT"            : return -1.0
    if ObjFun_name == "EFFICIENCY"      : return -1.0
    if ObjFun_name == "THRUST"          : return -1.0
    if ObjFun_name == "FIGURE_OF_MERIT" : return -1.0
    
    # otherwise
    return 1.0

#: def get_ObjFunSign()



# -------------------------------------------------------------------
# Read Objective Function Values from History File
# -------------------------------------------------------------------

def get_ObjFunVals( History_filename , special_cases, scale=1.0):
    
    # skip to last line of history file
    for line in open(History_filename): 
        if "ERROR" in line: error = 0.0
    # split by comma
    line = line.split(",")
        
    # store desired value
    ObjFun_values = {}
    ObjFun_values["LIFT"]               = float (line[1])  * scale
    ObjFun_values["DRAG"]               = float (line[2])  * scale
    ObjFun_values["SIDEFORCE"]          = float (line[3])  * scale
    ObjFun_values["MOMENT_X"]           = float (line[4])  * scale
    ObjFun_values["MOMENT_Y"]           = float (line[5])  * scale
    ObjFun_values["MOMENT_Z"]           = float (line[6])  * scale
    ObjFun_values["FORCE_X"]            = float (line[7])  * scale
    ObjFun_values["FORCE_Y"]            = float (line[8])  * scale
    ObjFun_values["FORCE_Z"]            = float (line[9])  * scale
    ObjFun_values["EFFICIENCY"]         = float (line[10]) * scale
    
    for key in special_cases: 
        if key == "FREE_SURFACE"   : 
            ObjFun_values["FREE_SURFACE"]       = float (line[11])  * scale
        if key == "ROTATING_FRAME" : 
            ObjFun_values["FIGURE_OF_MERIT"]    = float (line[11]) * scale    
            ObjFun_values["THRUST"]             = float (line[12]) * scale
            ObjFun_values["TORQUE"]             = float (line[13]) * scale
        if key == "EQUIV_AREA"     : 
            ObjFun_values["EQUIVALENT_AREA"]    = float (line[11]) * scale 
            ObjFun_values["NEARFIELD_PRESSURE"] = float (line[12]) * scale           
    
    return ObjFun_values
    
#: def get_ObjFunVals()



# -------------------------------------------------------------------
# Get Adjoint Filename Prefix
# -------------------------------------------------------------------

def get_AdjointPrefix(cadj_objfunc):
    
    # return desired prefix
    if cadj_objfunc == "DRAG"               : return "cd" 
    if cadj_objfunc == "LIFT"               : return "cl"
    if cadj_objfunc == "SIDEFORCE"          : return "csf"
    if cadj_objfunc == "MOMENT_X"           : return "cmx"
    if cadj_objfunc == "MOMENT_Y"           : return "cmy"
    if cadj_objfunc == "MOMENT_Z"           : return "cmz"
    if cadj_objfunc == "FORCE_X"            : return "cfx"
    if cadj_objfunc == "FORCE_Y"            : return "cfy"
    if cadj_objfunc == "FORCE_Z"            : return "cfz"
    if cadj_objfunc == "EFFICIENCY"         : return "eff"
    if cadj_objfunc == "EQUIVALENT_AREA"    : return "ea"
    if cadj_objfunc == "NEARFIELD_PRESSURE" : return "nfp"
    if cadj_objfunc == "THRUST"             : return "ct"
    if cadj_objfunc == "TORQUE"             : return "cq"
    if cadj_objfunc == "FIGURE_OF_MERIT"    : return "merit"    
    if cadj_objfunc == "FREESURFACE"        : return "fs"  
    
    # otherwise...
    raise Exception('Unrecognized objective function')
    
#: def get_AdjointPrefix()
    
  
  
# -------------------------------------------------------------------
# Get Design Variable Kind Name
# -------------------------------------------------------------------

def get_dvKind(kindID):
    
    if kindID == 1   : return "HICKS_HENNE" 
    if kindID == 4   : return "NACA_4DIGITS"
    if kindID == 5   : return "DISPLACEMENT"
    if kindID == 6   : return "ROTATION"
    if kindID == 7   : return "FFD_CONTROL_POINT"
    if kindID == 8   : return "FFD_DIHEDRAL_ANGLE"
    if kindID == 9   : return "FFD_TWIST_ANGLE"
    if kindID == 10  : return "FFD_ROTATION"    
    if kindID == 11  : return "FFD_CAMBER"
    if kindID == 12  : return "FFD_THICKNESS"
    if kindID == 13  : return "FFD_VOLUME"
    if kindID == 101 : return "MACH_NUMBER"
    if kindID == 102 : return "AOA"     
    
    # otherwise...
    raise Exception('Unrecognized Design Variable Kind')    

#: def get_dvKind()
  
  
    
# -------------------------------------------------------------------
# Get Gradient File Header
# -------------------------------------------------------------------

def get_GradFileFormat(grad_type,plot_format,kindID,special_cases=[]):
    
    # start header, build a list of strings and join at the end
    header       = []
    write_format = []
    
    # handle plot formating
    if   plot_format == 'TECPLOT': 
        header.append('VARIABLES=')
    elif plot_format == 'PARAVIEW':
        pass
    else: raise Exception('output plot format not recognized')
    
    # CASE: CONTINUOUS ADJOINT
    if grad_type == 'CONTINUOUS_ADJOINT':
        header.append(r'"iVar","Gradient","FinDiff_Step"')
        write_format.append(r'%d, %.10f, %f')
        
    # CASE: FINITE DIFFERENCE    
    elif grad_type == 'FINITE_DIFFERENCE':
        header.append(r'"iVar","Grad_CLift","Grad_CDrag","Grad_CLDRatio","Grad_CSideForce","Grad_CMx","Grad_CMy","Grad_CMz","Grad_CFx","Grad_CFy","Grad_CFz"')
        write_format.append(r'%d, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f')
        
        for key in special_cases: 
            if key == "FREE_SURFACE"   : 
                header.append(r',"Grad_CFreeSurface"')
                write_format.append(", %.10f ")
            if key == "ROTATING_FRAME" : 
                header.append(r',"Grad_CMerit","Grad_CT","Grad_CQ"')
                write_format.append(", %.10f, %.10f, %.10f")
            if key == "EQUIV_AREA"     : 
                header.append(r',"Grad_CEquivArea","Grad_CNearField"') 
                write_format.append(", %.10f, %.10f")        
    # otherwise...
    else: raise Exception('Unrecognized Gradient Type')          
        
    # DESIGN VARIABLE PARAMETERS
    if   kindID == "HICKS_HENNE"        : 
        header.append(r',"Up/Down","Loc_Max"')
        write_format.append(r', %s, %s')
    elif kindID == "NACA_4DIGITS"       : 
        header.append(r',"1st_digit","2nd_digit","3rd&4th_digits"')
        write_format.append(r', %s, %s, %s')
    elif kindID == "DISPLACEMENT"       : 
        header.append(r',"x_Disp","y_Disp","z_Disp"')
        write_format.append(r', %s, %s, %s')
    elif kindID == "ROTATION"           : 
        header.append(r',"x_Orig","y_Orig","z_Orig","x_End","y_End","z_End"')
        write_format.append(r', %s, %s, %s, %s, %s, %s')
    elif kindID == "FFD_CONTROL_POINT"  : 
        header.append(r',"Chunk","xIndex","yIndex","zIndex","xAxis","yAxis","zAxis"')
        write_format.append(r', %s, %s, %s, %s, %s, %s, %s')
    elif kindID == "FFD_DIHEDRAL_ANGLE" : 
        header.append(r',"Chunk","x_Orig","y_Orig","z_Orig","x_End","y_End","z_End"')
        write_format.append(r', %s, %s, %s, %s, %s, %s, %s')
    elif kindID == "FFD_TWIST_ANGLE"    : 
        header.append(r',"Chunk","x_Orig","y_Orig","z_Orig","x_End","y_End","z_End"')
        write_format.append(r', %s, %s, %s, %s, %s, %s, %s')
    elif kindID == "FFD_ROTATION"       : 
        header.append(r',"Chunk","x_Orig","y_Orig","z_Orig","x_End","y_End","z_End"')
        write_format.append(r', %s, %s, %s, %s, %s, %s, %s')
    elif kindID == "FFD_CAMBER"         : 
        header.append(r',"Chunk","xIndex","yIndex"')
        write_format.append(r', %s, %s, %s')
    elif kindID == "FFD_THICKNESS"      : 
        header.append(r',"Chunk","xIndex","yIndex"')
        write_format.append(r', %s, %s, %s')
    elif kindID == "FFD_VOLUME"         : 
        header.append(r',"Chunk","xIndex","yIndex"')
        write_format.append(r', %s, %s, %s')
    elif kindID == "MACH_NUMBER"        : pass
    elif kindID == "AOA"                : pass
    
    # otherwise...
    else: raise Exception('Unrecognized Design Variable Kind') 
    
    # FINISH FORMAT
    if grad_type == 'FINITE_DIFFERENCE':    
        header.append(r',"FinDiff_Step"'+'\n')  
        write_format.append(r', %.10f'+'\n')
    elif grad_type == 'CONTINUOUS_ADJOINT':
        header.append('\n')  
        write_format.append('\n')
        
    header       = ''.join(header)
    write_format = ''.join(write_format)
    
    return [header,write_format]
        
#: def get_GradFileFormat()
    
  
  
# -------------------------------------------------------------------
# Get Extension Name
# -------------------------------------------------------------------

def get_ExtensionName(output_format):
    
    if (output_format == "PARAVIEW") : return ".csv"
    if (output_format == "TECPLOT")  : return ".plt"        
    if (output_format == "SOLUTION") : return ".dat"  
    if (output_format == "RESTART")  : return ".dat"  
    if (output_format == "CONFIG")   : return ".cfg"  

    # otherwise
    raise Exception("Output Format Unknown")

#: def get_ExtensionName()
   
        

# -------------------------------------------------------------------
# Switch Class
# -------------------------------------------------------------------  

class switch(object):
# because if elif elif is not readable
# source: Brian Beck, PSF License, ActiveState Code
#         http://code.activestate.com/recipes/410692/
    
    def __init__(self, value):
        self.value = value
        self.fall = False

    def __iter__(self):
        """Return the match method once, then stop"""
        yield self.match
        raise StopIteration
    
    def match(self, *args):
        """Indicate whether or not to enter a case suite"""
        if self.fall or not args:
            return True
        elif self.value in args: # changed for v1.5, see below
            self.fall = True
            return True
        else:
            return False
        
#: class switch()


