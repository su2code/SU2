#!/usr/bin/env python 

## \file libSU2.py
#  \brief Support Functions for SU2 Python Scripts
#  \author Trent Lukaczyk, Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
#  \version 2.0.2
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

import os, time, sys, pickle, errno, copy
import shutil, glob
import numpy

try:
    import scipy.io
    scipy_loaded = True
except ImportError:
    scipy_loaded = False

# -------------------------------------------------------------------
#  Get SU2 Configuration Parameters
# -------------------------------------------------------------------

def Get_ConfigParams(filename):
      
    # initialize output dictionary
    data_dict = {}   
    
    input_file = open(filename)
    
    # process each line
    while 1:
        # read the line
        line = input_file.readline()
        if not line:
            break
        
        # remove line returns
        line = line.strip('\r\n')
        # make sure it has useful data
        if (not "=" in line) or (line[0] == '%'):
            continue
        # split across equals sign
        line = line.split("=")
        this_param = line[0].strip()
        this_value = line[1].strip()
        
        assert not data_dict.has_key(this_param) , ('Config file has multiple specifications of %s' % this_param )
        for case in switch(this_param):
            
            # comma delimited lists of strings with or without paren's
            if case("TASKS")          : pass
            if case("GRADIENTS")      : pass
            if case("DV_KIND")        : pass
            if case("CONST_EQ")       : pass
            if case("CONST_IEQ")      : pass
            if case("CONST_IEQ_SIGN") : 
                # remove white space
                this_value = ''.join(this_value.split())   
                # split by comma
                data_dict[this_param] = this_value.split(",")
                break
            
            # semicolon delimited lists of comma delimited lists of floats
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
            
            # comma delimited lists of floats
            if case("DV_VALUE_OLD")    : pass
            if case("DV_VALUE_NEW")    : pass
            if case("CONST_EQ_SCALE")  : pass
            if case("CONST_EQ_VALUE")  : pass
            if case("CONST_IEQ_SCALE") : pass
            if case("CONST_IEQ_VALUE") :             
                # remove white space
                this_value = ''.join(this_value.split())                
                # split by comma, map to float, store in dictionary
                data_dict[this_param] = map(float,this_value.split(","))
                break              

            # float parameters
            if case("MACH_NUMBER")            : pass
            if case("AoA")                    : pass
            if case("OBJFUNC_SCALE")          : pass
            if case("FIN_DIFF_STEP")          : pass
            if case("WRT_SOL_FREQ")           :
                data_dict[this_param] = float(this_value)
                break                 
            
            # int parameters
            if case("NUMBER_PART")            : pass
            if case("AVAILABLE_PROC")         : pass
            if case("EXT_ITER")               : pass
            if case("TIME_INSTANCES")         : pass
            if case("ADAPT_CYCLES")           :
                data_dict[this_param] = int(this_value)
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
                dv_Definitions = { 'KIND'   : dv_Kind       ,
                                   'SCALE'  : dv_Scale      ,
                                   'MARKER' : dv_Markers    ,
                                   'PARAM'  : dv_Parameters }
                # save to output dictionary
                data_dict[this_param] = dv_Definitions
                break  
            
            # unitary objective definition
            if case('OPT_OBJFUNC'):
                # remove white space
                this_value = ''.join(this_value.split())                
                # split by scale
                this_value = this_value.split("*")
                this_def = {'OBJECTIVE':this_value[0]}
                this_scale = 1.0
                if len(this_value) > 1:
                    this_scale = float( this_value[1] )
                this_def['SCALE'] = this_scale
                # save to output dictionary
                data_dict[this_param] = this_def
                break
            
            # unitary constraint definition
            if case('OPT_CONSTR'):
                # remove white space
                this_value = ''.join(this_value.split())                    
                # check for none case
                if this_value == 'NONE':
                    data_dict[this_param] = {'EQUALITY':{}, 'INEQUALITY':{}}
                    break                    
                # split definitions
                this_value = this_value.split(';')
                this_def = {}
                for this_con in this_value:
                    if not this_con: continue # if no definition
                    # defaults
                    this_obj = 'NONE'
                    this_sgn = '='
                    this_scl = 1.0
                    this_val = 0.0
                    # split scale if present
                    this_con = this_con.split('*')
                    if len(this_con) > 1:
                        this_scl = float( this_con[1] )
                    this_con = this_con[0]
                    # find sign
                    for this_sgn in ['<','>','=']:
                        if this_sgn in this_con: break
                    # split sign, store objective and value
                    this_con = this_con.strip('()').split(this_sgn)
                    assert len(this_con) == 2 , 'incorrect constraint definition'
                    this_obj = this_con[0]
                    this_val = float( this_con[1] )
                    # store in dictionary
                    this_def[this_obj] = { 'SIGN'  : this_sgn ,
                                           'VALUE' : this_val ,
                                           'SCALE' : this_scl  }
                #: for each constraint definition
                # sort constraints by type
                this_sort = { 'EQUALITY'   : {} ,
                              'INEQUALITY' : {}  }
                for key,value in this_def.iteritems():
                    if value['SIGN'] == '=':
                        this_sort['EQUALITY'][key]   = value
                    else:
                        this_sort['INEQUALITY'][key] = value
                #: for each definition                
                # save to output dictionary
                data_dict[this_param] = this_sort
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
#  Set SU2 Configuration Parameters
# -------------------------------------------------------------------

def Set_ConfigParams(filename,param_dict):
    
    temp_filename = "temp.cfg"
    shutil.copy(filename,temp_filename)
    output_file = open(filename,"w")

    # break pointers
    param_dict = copy.deepcopy(param_dict)
    
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
              
            # comma delimited list of floats
            if case("DV_VALUE_NEW") : pass
            if case("DV_VALUE_OLD") :
                n_lists = len(new_value)
                for i_value in range(n_lists):
                    output_file.write("%s" % new_value[i_value])
                    if i_value+1 < n_lists:
                        output_file.write(", ")               
                break
            
            # comma delimited list of strings no paren's
            if case("DV_KIND")            : pass
            if case("TASKS")              : pass
            if case("GRADIENTS")          :            
                if not isinstance(new_value,list):
                    new_value = [ new_value ]
                n_lists = len(new_value)
                for i_value in range(n_lists):
                    output_file.write(new_value[i_value])
                    if i_value+1 < n_lists:
                        output_file.write(", ")               
                break            
            
            # comma delimited list of strings inside paren's
            if case("DV_MARKER") : 
                if not isinstance(new_value,list):
                    new_value = [ new_value ]                
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
                assert isinstance(new_value,list) , 'incorrect specification of DV_PARAM'
                if not isinstance(new_value[0],list): new_value = [ new_value ]
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
            
            # int parameters
            if case("NUMBER_PART")            : pass
            if case("EXT_ITER")               :
                output_file.write("%i" % new_value)
                break                
            
            # default, assume string, integer or unformatted float 
            if case():
                output_file.write('%s' % new_value)
                break                         
                
        #: for case
        
        # remove from param dictionary
        del param_dict[this_param]
        
        # next line
        output_file.write("\n")        
        
    #: for each line
    
    # check that all params were used
    for this_param in param_dict.keys():
        if not this_param in ['JOB_NUMBER']:
            print ( 'Warning: Parameter %s not found in config file and was not written' % (this_param) )
        
    output_file.close()
    os.remove( temp_filename )
    
#: def Set_ConfigParams()



# -------------------------------------------------------------------
#  Read SU2_GPC Gradient Values
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
#  Read All Data from History File
# -------------------------------------------------------------------

def read_History( History_filename ):
    
    extension = os.path.splitext( History_filename )[1]
    
    # open history file
    history_file = open(History_filename)
    
    # skip first line
    line = history_file.readline()

    # process header
    if extension == '.plt':
        line = history_file.readline()
        line = line.split("=")[1].strip()
    line = line.split(",")
    Variables = [ x.strip('"') for x in line ]
    n_Vars = len(Variables)
    
    # header name to config file name map
    map_dict = get_HeaderMap()
        
    # map header names
    for i_Var in range(n_Vars):
        this_variable = Variables[i_Var]
        if map_dict.has_key(this_variable):
            Variables[i_Var] = map_dict[this_variable]
    
    # start history data dictionary
    history_data = dict.fromkeys(Variables,[])
    
    # skip another line
    line = history_file.readline()
    
    # read all data rows
    while 1:
        # read line
        line = history_file.readline()
        if not line:
            break
        
        # split line
        line_data = line.strip().split(',')
        line_data = [ float(x.strip()) for x in line_data ]  
        
        # store to dictionary
        for i_Var in range(n_Vars):
            this_variable = Variables[i_Var] # can't assume dictionary keys are in order
            history_data[this_variable] = history_data[this_variable] + [ line_data[i_Var] ]
    
    # done
    history_file.close()              
    return history_data
    
#: def read_History()



# -------------------------------------------------------------------
#  Define Dictionary Map for Header Names
# -------------------------------------------------------------------

def get_HeaderMap():
    
    # header name to config file name map
    map_dict = { "Iteration"       : "ITERATION"          ,
                 "CLift"           : "LIFT"               ,
                 "CDrag"           : "DRAG"               ,
                 "CSideForce"      : "SIDEFORCE"          ,
                 "CMx"             : "MOMENT_X"           ,
                 "CMy"             : "MOMENT_Y"           ,
                 "CMz"             : "MOMENT_Z"           ,
                 "CFx"             : "FORCE_X"            ,
                 "CFy"             : "FORCE_Y"            ,
                 "CFz"             : "FORCE_Z"            ,
                 "CL/CD"           : "EFFICIENCY"         ,
                 "CEff"            : "EFFICIENCY"         ,
                 "CFreeSurface"    : "FREE_SURFACE"       ,
                 "CMerit"          : "FIGURE_OF_MERIT"    ,
                 "CQ"              : "TORQUE"             ,
                 "CT"              : "THRUST"             ,
                 "CEquivArea"      : "EQUIVALENT_AREA"    ,
                 "CNearFieldOF"    : "NEARFIELD_PRESSURE" ,
                 "CWave"           : "NOISE"              ,
                 "Time(min)"       : "TIME"                }    
    
    return map_dict

#: def get_HeaderMap()



# -------------------------------------------------------------------
#  Read Objective Function Values from History File
# -------------------------------------------------------------------

def get_ObjFunVals( History_filename , special_cases=[] ):
    
    # read the history data
    history_data = read_History(History_filename)
    
    # list of objective functions to pull
    objfun_names = [ "LIFT"               ,
                     "DRAG"               ,
                     "SIDEFORCE"          ,
                     "MOMENT_X"           ,
                     "MOMENT_Y"           ,
                     "MOMENT_Z"           ,
                     "FORCE_X"            ,
                     "FORCE_Y"            ,
                     "FORCE_Z"            ,
                     "EFFICIENCY"         ,
                     "FREESURFACE"        ,
                     "FIGURE_OF_MERIT"    ,
                     "TORQUE"             ,
                     "THRUST"             ,
                     "EQUIVALENT_AREA"    ,
                     "NEARFIELD_PRESSURE" ,
                     "NOISE"               ]

    # pull only objective functions
    ObjFun_values = {}
    for this_objfun in objfun_names:
        if history_data.has_key(this_objfun):
            ObjFun_values[this_objfun] = history_data[this_objfun] 
    
    # for unsteady cases, average time-accurate objective function values
    if 'WRT_UNSTEADY' in special_cases:
        for key,value in ObjFun_values.iteritems():
            ObjFun_values[key] = sum(value)/len(value)
    
    # otherwise, keep only last value
    else:
        for key,value in ObjFun_values.iteritems():
            ObjFun_values[key] = value[-1]
                    
    return ObjFun_values
    
#: def get_ObjFunVals()



# -------------------------------------------------------------------
#  Get Objective Function Sign
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
#  Get Constraint Sign
# -------------------------------------------------------------------

def get_ConSign( sign ):
    
    sign_map = { '>' :  1.0 ,
                 '<' : -1.0  }
    assert not sign=='=' , 'Sign "=" not valid'
    
    return sign_map[sign]

#: def get_ConSign()


# -------------------------------------------------------------------
#  Get Adjoint Filename Prefix
# -------------------------------------------------------------------

def get_AdjointPrefix(adj_objfunc):
    
    # adjoint name map
    name_map = { "DRAG"               : "cd"    ,
                 "LIFT"               : "cl"    ,
                 "SIDEFORCE"          : "csf"   ,
                 "MOMENT_X"           : "cmx"   ,
                 "MOMENT_Y"           : "cmy"   ,
                 "MOMENT_Z"           : "cmz"   ,
                 "FORCE_X"            : "cfx"   ,
                 "FORCE_Y"            : "cfy"   ,
                 "FORCE_Z"            : "cfz"   ,
                 "EFFICIENCY"         : "eff"   ,
                 "EQUIVALENT_AREA"    : "ea"    ,
                 "NEARFIELD_PRESSURE" : "nfp"   ,
                 "THRUST"             : "ct"    ,
                 "TORQUE"             : "cq"    ,
                 "FIGURE_OF_MERIT"    : "merit" ,  
                 "FREESURFACE"        : "fs"    ,
                 "NOISE"              : "fwh"    }
    
    # if none or false, return map
    if not adj_objfunc:
        return name_map
    
    # return desired objective function prefix
    elif name_map.has_key(adj_objfunc):
        return name_map[adj_objfunc]
    
    # otherwise...
    else:
        raise Exception('Unrecognized objective function')
    
#: def get_AdjointPrefix()
    
# -------------------------------------------------------------------
#  Add a Prefix
# -------------------------------------------------------------------
# might more properly be called a suffix...
def add_Prefix(base_name,prefix):
    
    base_name = os.path.splitext(base_name)    
    prefixed_name = base_name[0] + '_' + prefix + base_name[1]
    
    return prefixed_name
    
#: def add_Prefix()
  
  
# -------------------------------------------------------------------
#  Get Design Variable Kind Name
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
#  Get Gradient File Header
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
    
    # Case: continuous adjoint
    if grad_type == 'CONTINUOUS_ADJOINT':
        header.append(r'"iVar","Gradient","FinDiff_Step"')
        write_format.append(r'%4d, %.10f, %f')
        
    # Case: finite difference  
    elif grad_type == 'FINITE_DIFFERENCE':
        header.append(r'"iVar","Grad_CLift","Grad_CDrag","Grad_CLDRatio","Grad_CSideForce","Grad_CMx","Grad_CMy","Grad_CMz","Grad_CFx","Grad_CFy","Grad_CFz"')
        write_format.append(r'%4d, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f')
        
        for key in special_cases: 
            if key == "FREE_SURFACE"   : 
                header.append(r',"Grad_CFreeSurface"')
                write_format.append(", %.10f ")
            if key == "ROTATING_FRAME" : 
                header.append(r',"Grad_CMerit","Grad_CT","Grad_CQ"')
                write_format.append(", %.10f, %.10f, %.10f")
            if key == "EQUIV_AREA"     : 
                header.append(r',"Grad_CEquivArea","Grad_CNearFieldOF"') 
                write_format.append(", %.10f, %.10f")  
            if key == "AEROACOUSTIC_EULER"   : 
                header.append(r',"Grad_CWave"')
                write_format.append(", %.10f ")
    # otherwise...
    else: raise Exception('Unrecognized Gradient Type')          
        
    # design variable parameters
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
    
    # finite difference step
    if grad_type == 'FINITE_DIFFERENCE':    
        header.append(r',"FinDiff_Step"')  
        write_format.append(r', %.10f')
    
    # finish format
    header.append('\n')  
    write_format.append('\n')
        
    header       = ''.join(header)
    write_format = ''.join(write_format)
    
    return [header,write_format]
        
#: def get_GradFileFormat()
    
    
    
# -------------------------------------------------------------------
#  Get Optimization File Header
# -------------------------------------------------------------------    
    
def get_OptFileFormat(plot_format,special_cases=[]):
    
    # start header, build a list of strings and join at the end
    header_list   = []
    header_format = ''
    write_format  = []
    
    # handle plot formating
    if   plot_format == 'TECPLOT': 
        header_format = header_format + 'VARIABLES='
    elif plot_format == 'PARAVIEW':
        pass
    else: raise Exception('output plot format not recognized')

    # start header
    header_list.extend(["Iteration","CLift","CDrag","CSideForce","CMx","CMy","CMz","CFx","CFy","CFz","CEff"])
    write_format.append(r'%4d, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f')
        
    # special cases
    for key in special_cases: 
        if key == "FREE_SURFACE"   : 
            header_list.extend(["CFreeSurface"])
            write_format.append(r', %.10f ')
        if key == "ROTATING_FRAME" : 
            header_list.extend(["CMerit","CT","CQ"])
            write_format.append(r', %.10f, %.10f, %.10f')
        if key == "EQUIV_AREA"     : 
            header_list.extend(["CEquivArea","CNearFieldOF"]) 
            write_format.append(r', %.10f, %.10f')  
        if key == "AEROACOUSTIC_EULER"   : 
            header_list.extend(["CWave"])
            write_format.append(r', %.10f ')
    
    # finish formats   
    header_format = (header_format) + ('"') + ('","').join(header_list) + ('"') + (' \n')
    write_format  = ''.join(write_format)  + ' \n'
            
    # build list of objective function names
    header_vars = []
    map_dict = get_HeaderMap()
    for variable in header_list:
        assert map_dict.has_key(variable) , 'unrecognized header variable'
        header_vars.append(map_dict[variable])
    
    # done
    return [header_format,header_vars,write_format]
        
#: def get_OptFileFormat()
  
  
  
# -------------------------------------------------------------------
#  Get Extension Name
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
#  Check Special Case
# -------------------------------------------------------------------
#   returns a list of special physical problems that were
#   specified in the config file, and set to 'yes'
def get_SpecialCases(data_dict):
    
    all_special_cases = [ 'FREE_SURFACE'       ,
                          'ROTATING_FRAME'     ,
                          'EQUIV_AREA'         ,
                          'WRT_UNSTEADY'       ,
                          'AEROACOUSTIC_EULER'  ]
    
    special_cases = []
    for key in all_special_cases:
        if data_dict.has_key(key) and data_dict[key] == 'YES':
            special_cases.append(key)
        if data_dict.has_key('PHYSICAL_PROBLEM') and data_dict['PHYSICAL_PROBLEM'] == key:
            special_cases.append(key)
  
    # no support for more than one special case (except for noise)
    if len(special_cases) > 1 and 'AEROACOUSTIC_EULER' not in special_cases:
        error_str = 'Currently cannot support ' + ' and '.join(special_cases) + ' at once'
        raise Exception(error_str)   
    
    if (data_dict['WRT_SOL_FREQ'] != 1) and ('WRT_UNSTEADY' in special_cases):
        raise Exception('Must set WRT_SOL_FREQ= 1 for WRT_UNSTEADY= YES')
  
    # Special case for time-spectral
    if data_dict.has_key('UNSTEADY_SIMULATION') and data_dict['UNSTEADY_SIMULATION'] == 'TIME_SPECTRAL':
        special_cases.append('TIME_SPECTRAL')

    return special_cases

#: def get_SpecialCases()



# -------------------------------------------------------------------
#  Load a Dictionary of Data
# -------------------------------------------------------------------

def load_data( file_name, var_names=None ,
               file_format = 'infer'     ,
               core_name   = 'python_data'    ):
    """ loads dictionary of data from python pickle or matlab struct
        default looks for variable 'python_data' in file_name
        file_format = pickle, will return any python object
        file_format = matlab, will return strings or float lists and 
        requires scipy.io.loadmat
        file_format = infer (default), will infer format from extention 
        ('.mat','.pkl')
    """
    
    # process file format
    if file_format == 'infer':
        if   os.path.splitext(file_name)[1] == '.mat':
            file_format = 'matlab'
        elif os.path.splitext(file_name)[1] == '.pkl':
            file_format = 'pickle'
    assert file_format in ['matlab','pickle'] , 'unsupported file format'
        
    
    # LOAD MATLAB
    if file_format == 'matlab' and scipy_loaded:
        input_data = scipy.io.loadmat( file_name        = file_name ,
                                       squeeze_me       = False     ,
                                       chars_as_strings = True      ,
                                       struct_as_record = True       )
        # pull core variable
        assert input_data.has_key(core_name) , 'core data not found'
        input_data = input_data[core_name]
        
        # convert recarray to dictionary
        input_data = rec2dict(input_data)
        
    # LOAD PICKLE
    elif file_format == 'pickle':
        input_data = load_pickle(file_name)
        # pull core variable
        assert input_data.has_key(core_name) , 'core data not found'
        input_data = input_data[core_name]
        
    #: if file_format
    
    # load specified varname into dictionary
    if var_names != None:
        # check for one item name array
        if isinstance(var_names,str):
            var_names = [var_names,]     
        for key in input_data.keys():
            if not key in var_names:
                del input_data[key]
        #: for key
    #: if var_names
    
    return input_data

#: def load_data()



# -------------------------------------------------------------------
#  Save a Dictionary of Data
# -------------------------------------------------------------------


def save_data( file_name, data_dict, append=False ,
               file_format = 'infer'             ,
               core_name='python_data'            ):
    """ saves data file from matlab 5 and later 
        default saves variable 'python_data' in file_name
        will save nested dictionaries into nested matlab structures
        cannot save classes and modules
        uses scipy.io.loadmat
    """        
    
    # process file format
    if file_format == 'infer':
        if   os.path.splitext(file_name)[1] == '.mat':
            file_format = 'matlab'
        elif os.path.splitext(file_name)[1] == '.pkl':
            file_format = 'pickle'
    assert file_format in ['matlab','pickle'] , 'unsupported file format'
        
    # if appending needed 
    # TODO: don't overwrite other core_names
    if append == True and os.path.exists(file_name):
        # load old data
        data_dict_old = load_data( file_name   = file_name   ,
                                   var_names   = None        ,
                                   file_format = file_format ,
                                   core_name   = core_name    )
        # check for keys not in new data
        for key,value in data_dict_old.iteritems():
            if not data_dict.has_key(key):
                data_dict[key] = value
        #: for each dict item
    #: if append
    
    # save to core name
    data_dict = {core_name : data_dict}
    
    # SAVE MATLAB
    if file_format == 'matlab':
        # bunch it
        data_dict = bunch(data_dict)
        # save it
        scipy.io.savemat( file_name = file_name ,
                          mdict     = data_dict,
                          format    = '5',        # matlab 5 .mat format
                          oned_as   = 'column' )
    elif file_format == 'pickle':
        # save it
        save_pickle(file_name,data_dict)
        
    #: if file_format
    
    return

#: def save_data()     
     


# -------------------------------------------------------------------
#  Load Pickle
# -------------------------------------------------------------------

def load_pickle(file_name):
    """ assumes first entry is a list of all following data names
        returns dictionary of data
    """
    pkl_file = open(file_name,'rb')
    names = safe_unpickle.loadf(pkl_file)
    data_dict = dict.fromkeys(names,[])
    for key in names:
        data_dict[key] = safe_unpickle.loadf(pkl_file)
    pkl_file.close()
    return data_dict

#: def load_pickle()



# -------------------------------------------------------------------
#  Save Pickle
# -------------------------------------------------------------------

def save_pickle(file_name,data_dict):
    """ saves first entry as a list of all following data names
        saves dictionary of data
    """
    pkl_file = open(file_name,'wb')
    names = data_dict.keys()
    pickle.dump(names,pkl_file)
    for key in names:
        pickle.dump(data_dict[key],pkl_file)
    pkl_file.close()

#: def save_pickle()



# -------------------------------------------------------------------
#  Safe UnPickle
# -------------------------------------------------------------------

class safe_unpickle(pickle.Unpickler):
    ''' adds some safety to unpickling
        checks that only supported classes are loaded 
        original source from http://nadiana.com/python-pickle-insecure#comment-144 
    '''
    
    # modules and classes considered safe
    PICKLE_SAFE = {
        'copy_reg'        : set(['_reconstructor']), 
        '__builtin__'     : set(['object'])          ,
        'tasks_general'   : set(['General_Task'])    , # SU2 Specific
        'tasks_project'   : set(['Project','Job'])   , 
        'tasks_su2'       : set(['Decomp','Deform','Direct','Cont_Adjoint',
                                 'Multiple_Cont_Adjoint','Finite_Diff','Adapt']) ,
        'numpy'           : set(['dtype','ndarray']) ,
        'numpy.core.multiarray' : set(['scalar','_reconstruct']) ,
    }
    
    # check for save module/class
    def find_class(self, module, name):
        if not module in self.PICKLE_SAFE:
            raise pickle.UnpicklingError(
                'Attempting to unpickle unsafe module %s' % module
            )
        __import__(module)
        mod = sys.modules[module]
        if not name in self.PICKLE_SAFE[module]:
            raise pickle.UnpicklingError(
                'Attempting to unpickle unsafe class %s' % name
            )
        klass = getattr(mod, name)
        return klass
 
    # extend the load() and loads() methods
    @classmethod
    def loadf(self, pickle_file): # loads a file like pickle.load()
        return self(pickle_file).load()
    @classmethod
    def loads(self, pickle_string): #loads a string like pickle.loads()
        return self(StringIO.StringIO(pickle_string)).load()    



# -------------------------------------------------------------------
#  Flatten a List
# -------------------------------------------------------------------

def flatten_list(input_list): 
    ''' flatten an irregular list of lists of any depth
    '''  
    output_list = []
    for value in input_list:
        if isinstance(value,list):
            output_list.extend( flatten_list(value) )
        else:
            output_list.append(value)
    return output_list  

#: def flatten_list()



# -------------------------------------------------------------------
#  Append Lists in a Nested Dictionary
# -------------------------------------------------------------------

def append_nestdict(base_dict,add_dict):
    
    # ensure add_dict and base dict don't point to same object
    add_dict = copy.deepcopy(add_dict)
    
    # append add_dict keys
    for key in add_dict.keys():
        
        # ensure base_dict key exists and is a list
        if not base_dict.has_key(key):
            if isinstance( add_dict[key] , dict ):
                base_dict[key] = {}
            else:
                base_dict[key] = []
        elif not (    isinstance( base_dict[key] , list ) 
                   or isinstance( base_dict[key] , dict ) ):
            assert not isinstance( add_dict[key] , dict ) , 'base[key] is not a dictionary while add[key] is'
            base_dict[key] = [base_dict[key]]
        
        # append list or telescope
        if isinstance( base_dict[key] , dict ):
            append_nestdict(base_dict[key],add_dict[key]) # telescope
        else:
            base_dict[key].append(add_dict[key])
    
    #: for add_dict[key]
       
    # base_dict will be updated through its pointer
    return 

#: def append_nestdict()
   
            
# -------------------------------------------------------------------
#  Convert Record Array to Dictionary
# -------------------------------------------------------------------
 
def rec2dict(array_in):
    """ converts numpy record array to dictionary of lists 
        assumes array comes from scipy.io.loadmat, with
        squeeze_me = False and struct_as_record = True
    """
    
    assert isinstance(array_in,numpy.ndarray) , 'input must be a numpy record array'
    
    # make sure it's not an object array
    if array_in.dtype == numpy.dtype('object'):
        array_in = array_in.tolist()
    
    # get record keys/names
    keys = array_in.dtype.names
    
    # start output dictionary
    dataout = dict.fromkeys(keys,[])
    
    for key in keys:
        
        # squeeze_me option puts all items in a two-dim array
        value = array_in[key].tolist()[0][0]
        
        # convert string
        if isinstance(value[0],unicode):
            value = str(value[0])
            
        # convert array
        elif isinstance(value,numpy.ndarray):
            # check for another struct level
            if value.dtype.names == None:
                value = value.tolist()
            # telescoping
            else:
                value = rec2dict(value)
        
        # store value         
        dataout[key] = value        
        
    return dataout

#: def rec2dict()



# -------------------------------------------------------------------
#  Switch Class
# -------------------------------------------------------------------  
# source: Brian Beck, PSF License, ActiveState Code
#         http://code.activestate.com/recipes/410692/

class switch(object):
    
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
        elif self.value in args: 
            self.fall = True
            return True
        else:
            return False
        
#: class switch()



# -------------------------------------------------------------------
#  Bunch Class
# ------------------------------------------------------------------- 

class bunch:
    """ replicates dictionary functionality with class dot structure
        note not all dictionary methods have been implemented
    """
    
    def __init__(self, d):
        for k, v in d.items():
            if isinstance(v, dict):
                if len(v): v = bunch(v)
                else:      v = []
            self.__dict__[k] = v
            
    def __dict__(self):
        return self.__dict__
    
    # items
    def keys(self):
        return self.__dict__.keys()
    def values(self):
        return self.__dict__.values()
    def items(self):
        return self.__dict__.items()
    
    # dictionary get/set/etc
    def __getitem__(self,k):
        return self.__dict__[k]
    def __setitem__(self,k,v):
        self.__dict__[k] = v   
    def __delitem__(self,k):
        del self.__dict__[k]
    def __str__(self):
        print_format = '%s: %s'
        state = []
        for k,v in self.__dict__.items():
            if isinstance(v,bunch):
                v = '%i-item bunch' % len(v.items())
            state.append(print_format % (k,v) )
        return '\n'.join(state)

#: class bunch



# -------------------------------------------------------------------
#  Output Redirection 
# -------------------------------------------------------------------
# original source: original source: http://stackoverflow.com/questions/6796492/python-temporarily-redirect-stdout-stderr
class redirect_output(object):
    ''' Temporarily redirects sys.stdout and sys.stderr when used in
        a with contextmanager
    '''
    def __init__(self, stdout=None, stderr=None):
        
        _newout = False
        _newerr = False
        
        if isinstance(stdout,str):
            stdout = open(stdout,'a')
            _newout = True            
        if isinstance(stderr,str):
            stderr = open(stderr,'a')
            _newerr = True                   
                
        self._stdout = stdout or sys.stdout
        self._stderr = stderr or sys.stderr
        self._newout = _newout
        self._newerr = _newerr

    def __enter__(self):
        self.old_stdout, self.old_stderr = sys.stdout, sys.stderr
        self.old_stdout.flush(); self.old_stderr.flush()
        sys.stdout, sys.stderr = self._stdout, self._stderr

    def __exit__(self, exc_type, exc_value, traceback):
        self._stdout.flush(); self._stderr.flush()
        sys.stdout = self.old_stdout
        sys.stderr = self.old_stderr
        
        if self._newout:
            self._stdout.close()
        if self._newerr:
            self._stderr.close()           

#: class redirect_output()

# -------------------------------------------------------------------
#  File Lock Class
# -------------------------------------------------------------------  
# source: Evan Fosmark, BSD license
#         http://www.evanfosmark.com/2009/01/cross-platform-file-locking-support-in-python/
class FileLock(object):
    """ A file locking mechanism that has context-manager support so 
        you can use it in a with statement. This should be relatively cross
        compatible as it doesn't rely on msvcrt or fcntl for the locking.
    """
 
    def __init__(self, file_name, timeout=10, delay=.05):
        """ Prepare the file locker. Specify the file to lock and optionally
            the maximum timeout and the delay between each attempt to lock.
        """
        self.is_locked = False
        self.lockfile = os.path.join(os.getcwd(), "%s.lock" % file_name)
        self.file_name = file_name
        self.timeout = timeout
        self.delay = delay
 
 
    def acquire(self):
        """ Acquire the lock, if possible. If the lock is in use, it check again
            every `wait` seconds. It does this until it either gets the lock or
            exceeds `timeout` number of seconds, in which case it throws 
            an exception.
        """
        start_time = time.time()
        while True:
            try:
                self.fd = os.open(self.lockfile, os.O_CREAT|os.O_EXCL|os.O_RDWR)
                break;
            except OSError as e:
                if e.errno != errno.EEXIST:
                    raise 
                if (time.time() - start_time) >= self.timeout:
                    raise FileLockException("Timeout occured.")
                time.sleep(self.delay)
        self.is_locked = True
 
 
    def release(self):
        """ Get rid of the lock by deleting the lockfile. 
            When working in a `with` statement, this gets automatically 
            called at the end.
        """
        if self.is_locked:
            os.close(self.fd)
            os.unlink(self.lockfile)
            self.is_locked = False
 
 
    def __enter__(self):
        """ Activated when used in the with statement. 
            Should automatically acquire a lock to be used in the with block.
        """
        if not self.is_locked:
            self.acquire()
        return self
 
 
    def __exit__(self, type, value, traceback):
        """ Activated at the end of the with statement.
            It automatically releases the lock if it isn't locked.
        """
        if self.is_locked:
            self.release()
 
 
    def __del__(self):
        """ Make sure that the FileLock instance doesn't leave a lockfile
            lying around.
        """
        self.release()
        
class FileLockException(Exception):
    pass

#: class filelock
