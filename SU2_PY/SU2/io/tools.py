#!/usr/bin/env python 

## \file tools.py
#  \brief file i/o functions
#  \author T. Lukaczyk, F. Palacios
#  \version 4.1.2 "Cardinal"
#
# SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
#                      Dr. Thomas D. Economon (economon@stanford.edu).
#
# SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
#                 Prof. Piero Colonna's group at Delft University of Technology.
#                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
#                 Prof. Rafael Palacios' group at Imperial College London.
#
# Copyright (C) 2012-2016 SU2, the open-source CFD code.
#
# SU2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# SU2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SU2. If not, see <http://www.gnu.org/licenses/>.

# -------------------------------------------------------------------
#  Imports
# -------------------------------------------------------------------

import os, time, sys, pickle, errno, copy
import shutil, glob
from SU2.util import ordered_bunch

# -------------------------------------------------------------------
#  Read SU2_DOT Gradient Values
# -------------------------------------------------------------------

def read_gradients( Grad_filename , scale = 1.0):
    """ reads the raw gradients from the gradient file
        returns a list of floats
    """
        
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

#: def read_gradients()


# -------------------------------------------------------------------
#  Read All Data from a Plot File
# -------------------------------------------------------------------

def read_plot( filename ):
    """ reads a plot file
        returns an ordered bunch with the headers for keys
        and a list of each header's floats for values.
    """
    
    extension = os.path.splitext( filename )[1]
    
    # open history file
    plot_file = open(filename)
    
    # title?
    line = plot_file.readline()
    if line.startswith('TITLE'):
        title = line.split('=')[1] .strip() # not used right now
        line = plot_file.readline()

    # process header
    if '=' in line:
        line = line.split("=")[1].strip()
    line = line.split(",")
    Variables = [ x.strip('" ') for x in line ]
    n_Vars = len(Variables)
    
    # initialize plot data dictionary
    plot_data = ordered_bunch.fromkeys(Variables)
    # must default each value to avoid pointer problems
    for key in plot_data.keys(): plot_data[key] = [] 
    
    # zone list
    zones = []
        
    # read all data rows
    while 1:
        # read line
        line = plot_file.readline()
        if not line:
            break
        
        #zone?
        if line.startswith('ZONE'):
            zone = line.split('=')[1].strip('" ')
            zones.append(zone)
            continue
        
        # split line
        line_data = line.strip().split(',')
        line_data = [ float(x.strip()) for x in line_data ]  
        
        # store to dictionary
        for i_Var in range(n_Vars):
            this_variable = Variables[i_Var] 
            plot_data[this_variable] = plot_data[this_variable] + [ line_data[i_Var] ]
    
    #: for each line

    # check for number of zones
    if len(zones) > 1:
        raise IOError , 'multiple zones not supported'
    
    # done
    plot_file.close()              
    return plot_data


# -------------------------------------------------------------------
#  Read All Data from History File
# -------------------------------------------------------------------

def read_history( History_filename ):
    """ reads a history file
        returns an ordered bunch with the history file headers for keys
        and a list of each header's floats for values.
        if header is an optimization objective, its name is mapped to 
        the optimization name.
        Iter and Time(min) headers are mapped to ITERATION and TIME
        respectively.
    """
    
    # read plot file
    plot_data = read_plot( History_filename )
    
    # initialize history data dictionary
    history_data = ordered_bunch()    
    
    # header name to config file name map
    map_dict = get_headerMap()    
    
    # map header names
    for key in plot_data.keys():
        if map_dict.has_key(key):
            var = map_dict[key]
        else:
            var = key
        history_data[var] = plot_data[key]
    
    return history_data
    
#: def read_history()



# -------------------------------------------------------------------
#  Define Dictionary Map for Header Names
# -------------------------------------------------------------------

def get_headerMap():
    """ returns a dictionary that maps history file header names
        to optimization problem function names
    """
    # header name to config file name map
    map_dict = { "Iteration"       : "ITERATION"               ,
                 "CLift"           : "LIFT"                    ,
                 "CDrag"           : "DRAG"                    ,
                 "CSideForce"      : "SIDEFORCE"               ,
                 "Cp_Diff"         : "INVERSE_DESIGN_PRESSURE" ,
                 "HeatFlux_Diff"   : "INVERSE_DESIGN_HEATFLUX" ,
                 "HeatFlux_Total"  : "TOTAL_HEATFLUX"          ,
                 "HeatFlux_Maximum": "MAXIMUM_HEATFLUX"        ,
                 "CMx"             : "MOMENT_X"                ,
                 "CMy"             : "MOMENT_Y"                ,
                 "CMz"             : "MOMENT_Z"                ,
                 "CFx"             : "FORCE_X"                 ,
                 "CFy"             : "FORCE_Y"                 ,
                 "CFz"             : "FORCE_Z"                 ,
                 "CL/CD"           : "EFFICIENCY"              ,
                 "CEff"            : "EFFICIENCY"              ,
                 "CFreeSurface"    : "FREE_SURFACE"            ,
                 "CMerit"          : "FIGURE_OF_MERIT"         ,
                 "CQ"              : "TORQUE"                  ,
                 "CT"              : "THRUST"                  ,
                 "CEquivArea"      : "EQUIVALENT_AREA"         ,
                 "CNearFieldOF"    : "NEARFIELD_PRESSURE"      ,
                 "Avg_TotalPress"  : "AVG_TOTAL_PRESSURE"      ,
                 "FluxAvg_Pressure": "AVG_OUTLET_PRESSURE"     ,
                 "FluxAvg_Density" : "FLUXAVG_OUTLET_DENSITY"  ,
                 "FluxAvg_Velocity": "FLUXAVG_OUTLET_VELOCITY" ,
                 "Avg_Mach"        : "AVG_OUTLET_MACH"         ,
                 "Avg_Temperature" : "AVG_OUTLET_TEMPERATURE"  ,
                 "MassFlowRate"    : "MASS_FLOW_RATE"          ,
                 "Time(min)"       : "TIME"                    ,
                 "D(CLift)"        : "D_LIFT"                  ,
                 "D(CDrag)"        : "D_DRAG"                  ,
                 "D(CSideForce)"   : "D_SIDEFORCE"             ,
                 "D(CMx)"          : "D_MOMENT_X"              ,
                 "D(CMy)"          : "D_MOMENT_Y"              ,
                 "D(CMz)"          : "D_MOMENT_Z"              ,
                 "D(CFx)"          : "D_FORCE_X"               ,
                 "D(CFy)"          : "D_FORCE_Y"               ,
                 "D(CFz)"          : "D_FORCE_Z"               ,
                 "D(CL/CD)"        : "D_EFFICIENCY"}
    
    return map_dict

#: def get_headerMap()


# -------------------------------------------------------------------
#  Optimizer Function Names
# -------------------------------------------------------------------

# Aerodynamic Optimizer Function Names
optnames_aero = [ "LIFT"                    ,
                  "DRAG"                    ,
                  "SIDEFORCE"               ,
                  "MOMENT_X"                ,
                  "MOMENT_Y"                ,
                  "MOMENT_Z"                ,
                  "FORCE_X"                 ,
                  "FORCE_Y"                 ,
                  "FORCE_Z"                 ,
                  "EFFICIENCY"              ,
                  "FREE_SURFACE"            ,
                  "FIGURE_OF_MERIT"         ,
                  "TORQUE"                  ,
                  "THRUST"                  ,
                  "AVG_TOTAL_PRESSURE"      ,
                  "AVG_OUTLET_PRESSURE"     ,
                  "AVG_OUTLET_DENSITY"      ,
                  "AVG_OUTLET_VELOCITY"     ,
                  "MASS_FLOW_RATE"          ,
                  "OUTFLOW_GENERALIZED"     ,
                  "EQUIVALENT_AREA"         ,
                  "NEARFIELD_PRESSURE"      ,
                  "INVERSE_DESIGN_PRESSURE" ,
                  "INVERSE_DESIGN_HEATFLUX" ,
                  "TOTAL_HEATFLUX"          ,
                  "MAXIMUM_HEATFLUX"        ]
#: optnames_aero

optnames_stab = [ "D_LIFT_D_ALPHA"               ,
                  "D_DRAG_D_ALPHA"               ,
                  "D_SIDEFORCE_D_ALPHA"          ,
                  "D_MOMENT_X_D_ALPHA"           ,
                  "D_MOMENT_Y_D_ALPHA"           ,
                  "D_MOMENT_Z_D_ALPHA"           ,
                ]

# Geometric Optimizer Function Names
optnames_geo = [ "MAX_THICKNESS"      ,
                 "1/4_THICKNESS"      ,
                 "1/3_THICKNESS"      ,
                 "1/2_THICKNESS"      ,
                 "2/3_THICKNESS"      ,
                 "3/4_THICKNESS"      ,
                 "AREA"               ,
                 "AOA"                ,
                 "CHORD"              ,
                 "MAX_THICKNESS_SEC1" ,
                 "MAX_THICKNESS_SEC2" ,
                 "MAX_THICKNESS_SEC3" ,
                 "MAX_THICKNESS_SEC4" ,
                 "MAX_THICKNESS_SEC5" ,
                 "1/4_THICKNESS_SEC1" ,
                 "1/4_THICKNESS_SEC2" ,
                 "1/4_THICKNESS_SEC3" ,
                 "1/4_THICKNESS_SEC4" ,
                 "1/4_THICKNESS_SEC5" ,
                 "1/3_THICKNESS_SEC1" ,
                 "1/3_THICKNESS_SEC2" ,
                 "1/3_THICKNESS_SEC3" ,
                 "1/3_THICKNESS_SEC4" ,
                 "1/3_THICKNESS_SEC5" ,
                 "1/2_THICKNESS_SEC1" ,
                 "1/2_THICKNESS_SEC2" ,
                 "1/2_THICKNESS_SEC3" ,
                 "1/2_THICKNESS_SEC4" ,
                 "1/2_THICKNESS_SEC5" ,
                 "2/3_THICKNESS_SEC1" ,
                 "2/3_THICKNESS_SEC2" ,
                 "2/3_THICKNESS_SEC3" ,
                 "2/3_THICKNESS_SEC4" ,
                 "2/3_THICKNESS_SEC5" ,
                 "3/4_THICKNESS_SEC1" ,
                 "3/4_THICKNESS_SEC2" ,
                 "3/4_THICKNESS_SEC3" ,
                 "3/4_THICKNESS_SEC4" ,
                 "3/4_THICKNESS_SEC5" ,
                 "AREA_SEC1"          ,
                 "AREA_SEC2"          ,
                 "AREA_SEC3"          ,
                 "AREA_SEC4"          ,
                 "AREA_SEC5"          ,
                 "AOA_SEC1"           ,
                 "AOA_SEC2"           ,
                 "AOA_SEC3"           ,
                 "AOA_SEC4"           ,
                 "AOA_SEC5"           ,
                 "CHORD_SEC1"         ,
                 "CHORD_SEC2"         ,
                 "CHORD_SEC3"         ,
                 "CHORD_SEC4"         ,
                 "CHORD_SEC5"         ,
                 "VOLUME"              ]
#: optnames_geo

grad_names_directdiff = ["D_LIFT"                  ,
                         "D_DRAG"                  ,
                         "D_SIDEFORCE"             ,
                         "D_MOMENT_X"              ,
                         "D_MOMENT_Y"              ,
                         "D_MOMENT_Z"              ,
                         "D_FORCE_X"               ,
                         "D_FORCE_Y"               ,
                         "D_FORCE_Z"               ,
                         "D_EFFICIENCY"]

grad_names_map = { "LIFT"      : "D_LIFT"           ,
                   "DRAG"      : "D_DRAG"           ,
                   "SIDEFORCE" : "D_SIDEFORCE" ,
                   "MOMENT_X"  : "D_MOMENT_X"   ,
                   "MOMENT_Y"  : "D_MOMENT_Y"   ,
                   "MOMENT_Z"  : "D_MOMENT_Z"   ,
                   "FORCE_X"   : "D_FORCE_X"     ,
                   "FORCE_Y"   : "D_FORCE_Y"     ,
                   "FORCE_Z"   : "D_FORCE_Z"     ,
                   "EFFICIENCY" : "D_EFFICIENCY"}
# -------------------------------------------------------------------
#  Read Aerodynamic Function Values from History File
# -------------------------------------------------------------------

def read_aerodynamics( History_filename , special_cases=[], final_avg=0 ):
    """ values = read_aerodynamics(historyname, special_cases=[])
        read aerodynamic function values from history file
        
        Outputs:
            dictionary with function keys and thier values
            if special cases has 'UNSTEADY_SIMULATION', returns time averaged data
            otherwise returns final value from history file
    """
    
    # read the history data
    history_data = read_history(History_filename)
    
    # list of functions to pull
    func_names = optnames_aero + grad_names_directdiff

    # pull only these functions
    Func_Values = ordered_bunch()
    for this_objfun in func_names:
        if history_data.has_key(this_objfun):
            Func_Values[this_objfun] = history_data[this_objfun] 
    
    # for unsteady cases, average time-accurate objective function values
    if 'UNSTEADY_SIMULATION' in special_cases:
        for key,value in Func_Values.iteritems():
            Func_Values[key] = sum(value)/len(value)
         
    # average the final iterations   
    elif final_avg:
        for key,value in Func_Values.iteritems():
            # only the last few iterations
            i_fin = min([final_avg,len(value)])
            value = value[-i_fin:]
            Func_Values[key] = sum(value)/len(value)
    
    # otherwise, keep only last value
    else:
        for key,value in Func_Values.iteritems():
            Func_Values[key] = value[-1]
                    
    return Func_Values
    
#: def read_aerodynamics()



# -------------------------------------------------------------------
#  Get Objective Function Sign
# -------------------------------------------------------------------

def get_objectiveSign( ObjFun_name ):
    """ returns -1 for maximization problems:
            LIFT
            EFFICIENCY
            THRUST
            FIGURE_OF_MERIT
            MASS_FLOW_RATE
            AVG_OUTLET_PRESSURE
            AVG_TOTAL_PRESSURE
        returns +1 otherwise
    """
    
    # flip sign for maximization problems
    if ObjFun_name == "LIFT"            : return -1.0
    if ObjFun_name == "EFFICIENCY"      : return -1.0
    if ObjFun_name == "THRUST"          : return -1.0
    if ObjFun_name == "FIGURE_OF_MERIT" : return -1.0
    if ObjFun_name == "MASS_FLOW_RATE" : return -1.0
    if ObjFun_name == "AVG_TOTAL_PRESSURE" : return -1.0
    if ObjFun_name == "AVG_OUTLET_PRESSURE" : return -1.0
    
    # otherwise
    return 1.0

#: def get_objectiveSign()


# -------------------------------------------------------------------
#  Get Constraint Sign
# -------------------------------------------------------------------

def get_constraintSign( sign ):
    """ gets +/-1 given a constraint sign < or > respectively
        inequality constraint is posed as c(x) < 0
    """
    sign_map = { '>' : -1.0 ,
                 '<' : +1.0  }
    assert not sign=='=' , 'Sign "=" not valid'
    
    return sign_map[sign]

#: def get_constraintSign()


# -------------------------------------------------------------------
#  Get Adjoint Filename Suffix
# -------------------------------------------------------------------

def get_adjointSuffix(objective_function=None):
    """ gets the adjoint suffix given an objective function """
    
    # adjoint name map
    name_map = { "DRAG"                    : "cd"        ,
                 "LIFT"                    : "cl"        ,
                 "SIDEFORCE"               : "csf"       ,
                 "MOMENT_X"                : "cmx"       ,
                 "MOMENT_Y"                : "cmy"       ,
                 "MOMENT_Z"                : "cmz"       ,
                 "FORCE_X"                 : "cfx"       ,
                 "FORCE_Y"                 : "cfy"       ,
                 "FORCE_Z"                 : "cfz"       ,
                 "EFFICIENCY"              : "eff"       ,
                 "INVERSE_DESIGN_PRESSURE" : "invpress"  ,
                 "INVERSE_DESIGN_HEAT"     : "invheat"   ,
                 "MAXIMUM_HEATFLUX"        : "maxheat"   ,
                 "TOTAL_HEATFLUX"          : "totheat"   ,
                 "EQUIVALENT_AREA"         : "ea"        ,
                 "NEARFIELD_PRESSURE"      : "nfp"       ,
                 "THRUST"                  : "ct"        ,
                 "TORQUE"                  : "cq"        ,
                 "FIGURE_OF_MERIT"         : "merit"     ,
                 "AVG_TOTAL_PRESSURE"      : "pt"        ,
                 "AVG_OUTLET_PRESSURE"     : "pe"        ,
                 "MASS_FLOW_RATE"          : "mfr"       ,
                 "OUTFLOW_GENERALIZED"       : "chn"       ,
                 "FREE_SURFACE"            : "fs"       }
    
    # if none or false, return map
    if not objective_function:
        return name_map
    
    # return desired objective function suffix
    elif name_map.has_key(objective_function):
        return name_map[objective_function]
    
    # otherwise...
    else:
        raise Exception('Unrecognized adjoint function name')
    
#: def get_adjointSuffix()
    
# -------------------------------------------------------------------
#  Add a Suffix
# -------------------------------------------------------------------

def add_suffix(base_name,suffix):
    """ suffix_name = add_suffix(base_name,suffix)
        adds suffix to a filename, accounting for file type extension
        example:
            base_name   = 'input.txt'
            suffix      = 'new'
            suffix_name = 'input_new.txt'
    """
    
    base_name = os.path.splitext(base_name)    
    suffix_name = base_name[0] + '_' + suffix + base_name[1]
    
    return suffix_name
    
#: def add_suffix()



# -------------------------------------------------------------------
#  Get Design Variable ID Map
# -------------------------------------------------------------------

def get_dvMap():
    """ get dictionary that maps design variable 
        kind id number to name """
    dv_map = { 1   : "HICKS_HENNE"           ,
               2   : "COSINE_BUMP"           ,
               3   : "SPHERICAL"             ,
               4   : "NACA_4DIGITS"          ,
               5   : "DISPLACEMENT"          ,
               6   : "ROTATION"              ,
               7   : "FFD_CONTROL_POINT"     ,
               8   : "FFD_DIHEDRAL_ANGLE"    ,
               9   : "FFD_TWIST_ANGLE"       ,
               10  : "FFD_ROTATION"          ,
               11  : "FFD_CAMBER"            ,
               12  : "FFD_THICKNESS"         ,
               14  : "FOURIER"               ,
               15  : "FFD_CONTROL_POINT_2D"  ,
               16  : "FFD_CAMBER_2D"         ,
               17  : "FFD_THICKNESS_2D"      ,
               19  : "CUSTOM"                ,
               101 : "MACH_NUMBER"           ,
               102 : "AOA"                    }
    
    return dv_map

#: def get_dvMap()

# -------------------------------------------------------------------
#  Get Design Variable Kind Name from ID
# -------------------------------------------------------------------
def get_dvKind( kindID ):
    """ get design variable kind name from id number """
    dv_map = get_dvMap()
    try: 
        return dv_map[ kindID ]
    except KeyError: 
        raise Exception('Unrecognized Design Variable ID')
# def get_dvKind()

# -------------------------------------------------------------------
#  Get Design Variable Kind ID from Name
# -------------------------------------------------------------------
def get_dvID( kindName ):
    """ get design variable kind id number from name """
    dv_map = get_dvMap()
    id_map = dict((v,k) for (k,v) in dv_map.iteritems())
    try: 
        return id_map[ kindName ]
    except KeyError: 
        raise Exception('Unrecognized Design Variable Name: %s' , kindName)
#: def get_dvID()
  
  
    
# -------------------------------------------------------------------
#  Get Gradient File Header
# -------------------------------------------------------------------

def get_gradFileFormat(grad_type,plot_format,kindID,special_cases=[]):
    
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
        header.append(r'"iVar","Grad_CLift","Grad_CDrag","Grad_CLDRatio","Grad_CSideForce","Grad_CMx","Grad_CMy","Grad_CMz","Grad_CFx","Grad_CFy","Grad_CFz","Grad_HeatFlux_Total","Grad_HeatFlux_Maximum"')
        write_format.append(r'%4d, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f')
        
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
            if key == "1D_OUTPUT"     :
                header.append(r',"Grad_Avg_TotalPress","Grad_Avg_Mach","Grad_Avg_Temperature","Grad_MassFlowRate","Grad_FluxAvg_Pressure","Grad_FluxAvg_Density","Grad_FluxAvg_Velocity","Grad_FluxAvg_Enthalpy"')
                write_format.append(", %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f")
            if key == "INV_DESIGN_CP"     :
                header.append(r',"Grad_Cp_Diff"')
                write_format.append(", %.10f")
            if key == "INV_DESIGN_HEATFLUX"     :
                header.append(r',"Grad_HeatFlux_Diff"')
                write_format.append(", %.10f")
            if key =="OUTFLOW_GENERALIZED"    :
                header.append(r',"Grad_Chain_Rule"')
                write_format.append(", %.10f")

    # otherwise...
    else: raise Exception('Unrecognized Gradient Type')          
        
    # design variable parameters
    if kindID == "FFD_CONTROL_POINT_2D"  :
        header.append(r',"FFD_Box_ID","xIndex","yIndex","xAxis","yAxis"')
        write_format.append(r', %s, %s, %s, %s, %s')
    elif kindID == "FFD_CAMBER_2D"         :
        header.append(r',"FFD_Box_ID","xIndex"')
        write_format.append(r', %s, %s')
    elif kindID == "FFD_THICKNESS_2D"      :
        header.append(r',"FFD_Box_ID","xIndex"')
        write_format.append(r', %s, %s')
    elif kindID == "HICKS_HENNE"        :
        header.append(r',"Up/Down","Loc_Max"')
        write_format.append(r', %s, %s')
    elif kindID == "GAUSS_BUMP"       :
        header.append(r',"Up/Down","Loc_Max","Size_Bump"')
        write_format.append(r', %s, %s, %s')
    elif kindID == "FAIRING"       :
        header.append(r',"ControlPoint_Index","Theta_Disp","R_Disp"')
        write_format.append(r', %s, %s, %s')
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
        header.append(r',"FFD_Box_ID","xIndex","yIndex","zIndex","xAxis","yAxis","zAxis"')
        write_format.append(r', %s, %s, %s, %s, %s, %s, %s')
    elif kindID == "FFD_DIHEDRAL_ANGLE" : 
        header.append(r',"FFD_Box_ID","x_Orig","y_Orig","z_Orig","x_End","y_End","z_End"')
        write_format.append(r', %s, %s, %s, %s, %s, %s, %s')
    elif kindID == "FFD_TWIST_ANGLE"    : 
        header.append(r',"FFD_Box_ID","x_Orig","y_Orig","z_Orig","x_End","y_End","z_End"')
        write_format.append(r', %s, %s, %s, %s, %s, %s, %s')
    elif kindID == "FFD_ROTATION"       : 
        header.append(r',"FFD_Box_ID","x_Orig","y_Orig","z_Orig","x_End","y_End","z_End"')
        write_format.append(r', %s, %s, %s, %s, %s, %s, %s')
    elif kindID == "FFD_CAMBER"         : 
        header.append(r',"FFD_Box_ID","xIndex","yIndex"')
        write_format.append(r', %s, %s, %s')
    elif kindID == "FFD_THICKNESS"      : 
        header.append(r',"FFD_Box_ID","xIndex","yIndex"')
        write_format.append(r', %s, %s, %s')
    elif kindID == "MACH_NUMBER"        : pass
    elif kindID == "AOA"                : pass
    elif kindID == "CUSTOM"             : pass
    
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
        
#: def get_gradFileFormat()
    
    
    
# -------------------------------------------------------------------
#  Get Optimization File Header
# -------------------------------------------------------------------    
    
def get_optFileFormat(plot_format,special_cases=None):
    
    if special_cases is None: special_cases = []
    
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
    header_list.extend(["Iteration","CLift","CDrag","CSideForce","CMx","CMy","CMz","CFx","CFy","CFz","CEff","HeatFlux_Total","HeatFlux_Maximum"])
    write_format.append(r'%4d, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f')
        
    # special cases
    for key in special_cases: 
        if key == "FREE_SURFACE" :
            header_list.extend(["CFreeSurface"])
            write_format.append(r', %.10f ')
        if key == "ROTATING_FRAME" : 
            header_list.extend(["CMerit","CT","CQ"])
            write_format.append(r', %.10f, %.10f, %.10f')
        if key == "EQUIV_AREA"     : 
            header_list.extend(["CEquivArea","CNearFieldOF"]) 
            write_format.append(r', %.10f, %.10f')
        if key == "1D_OUTPUT":
            header_list.extend(["Avg_TotalPress","Avg_Mach","Avg_Temperature","MassFlowRate","FluxAvg_Pressure","FluxAvg_Density","FluxAvg_Velocity","FluxAvg_Enthalpy"])
            write_format.append(r', %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f')
        if key == "INV_DESIGN_CP"     :
            header_list.extend(["Cp_Diff"])
            write_format.append(r', %.10f')
        if key == "INV_DESIGN_HEATFLUX"     :
            header_list.extend(["HeatFlux_Diff"])
            write_format.append(r', %.10f')
        if key =="OUTFLOW_GENERALIZED"    :
            header_list.exted(["Chain_Rule"])
            write_format.append([r", %.10f"])

    # finish formats
    header_format = (header_format) + ('"') + ('","').join(header_list) + ('"') + (' \n')
    write_format  = ''.join(write_format)  + ' \n'
            
    # build list of objective function names
    header_vars = []
    map_dict = get_headerMap()
    for variable in header_list:
        assert map_dict.has_key(variable) , 'unrecognized header variable'
        header_vars.append(map_dict[variable])
    
    # done
    return [header_format,header_vars,write_format]
        
#: def get_optFileFormat()
  
  
  
# -------------------------------------------------------------------
#  Get Extension Name
# -------------------------------------------------------------------

def get_extension(output_format):
  
    if (output_format == "PARAVIEW")        : return ".csv"
    if (output_format == "TECPLOT")         : return ".dat"
    if (output_format == "TECPLOT_BINARY")  : return ".plt"
    if (output_format == "SOLUTION")        : return ".dat"  
    if (output_format == "RESTART")         : return ".dat"  
    if (output_format == "CONFIG")          : return ".cfg"  

    # otherwise
    raise Exception("Output Format Unknown")

#: def get_extension()



# -------------------------------------------------------------------
#  Check Special Case
# -------------------------------------------------------------------
def get_specialCases(config):
    """ returns a list of special physical problems that were
        specified in the config file, and set to 'yes'
    """
    
    all_special_cases = [ 'FREE_SURFACE'                     ,
                          'ROTATING_FRAME'                   ,
                          'EQUIV_AREA'                       ,
                          '1D_OUTPUT'                        ,
                          'INV_DESIGN_CP'                    ,
                          'INV_DESIGN_HEATFLUX'              ,
                          'OUTFLOW_GENERALIZED'                ]
    
    special_cases = []
    for key in all_special_cases:
        if config.has_key(key) and config[key] == 'YES':
            special_cases.append(key)
        if config.has_key('PHYSICAL_PROBLEM') and config['PHYSICAL_PROBLEM'] == key:
            special_cases.append(key)
            
    if config.get('UNSTEADY_SIMULATION','NO') != 'NO':
        special_cases.append('UNSTEADY_SIMULATION')
     
    # no support for more than one special case
    if len(special_cases) > 1:
        error_str = 'Currently cannot support ' + ' and '.join(special_cases) + ' at once'
        raise Exception(error_str)   
    
    if (config['WRT_SOL_FREQ'] != 1) and ('WRT_UNSTEADY' in special_cases):
        raise Exception('Must set WRT_SOL_FREQ= 1 for WRT_UNSTEADY= YES')
  
    # Special case for time-spectral
    if config.has_key('UNSTEADY_SIMULATION') and config['UNSTEADY_SIMULATION'] == 'TIME_SPECTRAL':
        special_cases.append('TIME_SPECTRAL')

    # Special case for rotating frame
    if config.has_key('GRID_MOVEMENT_KIND') and config['GRID_MOVEMENT_KIND'] == 'ROTATING_FRAME':
        special_cases.append('ROTATING_FRAME')
        
    return special_cases

#: def get_specialCases()

# -------------------------------------------------------------------
#  Check Fluid Structure Interaction
# -------------------------------------------------------------------
def get_multizone(config):
    """ returns a list of special physical problems that were
        specified in the config file, and set to 'yes'
    """
    
    all_multizone_problems = ['FLUID_STRUCTURE_INTERACTION']
    
    multizone = []
    for key in all_multizone_problems:
        if config.has_key('PHYSICAL_PROBLEM') and config['PHYSICAL_PROBLEM'] == key:
            multizone.append(key)
            
    return multizone

#: def get_multizone()


def next_folder(folder_format,num_format='%03d'):
    """ folder = next_folder(folder_format,num_format='%03d')
        finds the next folder with given format
        
        Inputs:
            folder_format - folder name with wild card (*) to mark expansion
            num_format    - %d formating to expand the wild card with
            
        Outputs:
            folder - a folder with the next index number inserted in 
            the wild card, first index is 1
    """
    
    assert '*' in folder_format , 'wildcard (*) missing in folder_format name'
    
    folders = glob.glob(folder_format)
    split   = folder_format.split('*')
    folder  = folder_format.replace('*',num_format)
    
    if folders:
        # find folder number, could be done with regex...
        max_folder = max(folders)
        if split[0]:
            max_folder = max_folder.split(split[0])[1]
        if split[1]:
            max_folder = max_folder.rsplit(split[1])[0]            
        
        # last folder number
        max_i = int(max_folder)
        
        # increment folder number
        folder = folder % (max_i+1)
    else:
        # first folder, number 1
        folder = folder % 1

    return folder


def expand_part(name,config):
    names = [name]
    return names

def expand_time(name,config):
    if 'UNSTEADY_SIMULATION' in get_specialCases(config):
        n_time = config['UNST_ADJOINT_ITER']
        name_pat = add_suffix(name,'%05d')
        names = [name_pat%i for i in range(n_time)]
    else:
        names = [name]
    return names
        
def make_link(src,dst):
    """ make_link(src,dst)
        makes a relative link
        Inputs:
            src - source file
            dst - destination to place link
        
        Windows links currently unsupported, will copy file instead
    """
    
    assert os.path.exists(src) , 'source file does not exist \n%s' % src
    
    if os.name == 'nt':
        # can't make a link in windows, need to look for other options
        if os.path.exists(dst): os.remove(dst)
        shutil.copy(src,dst)
    
    else:
        # find real file, incase source itself is a link
        src = os.path.realpath(src) 
        
        # normalize paths
        src = os.path.normpath(src)
        dst = os.path.normpath(dst)        
        
        # check for self referencing
        if src == dst: return        
        
        # find relative folder path
        srcfolder = os.path.join( os.path.split(src)[0] ) + '/'
        dstfolder = os.path.join( os.path.split(dst)[0] ) + '/'
        srcfolder = os.path.relpath(srcfolder,dstfolder)
        src = os.path.join( srcfolder, os.path.split(src)[1] )
        
        # make unix link
        if os.path.exists(dst): os.remove(dst)
        os.symlink(src,dst)
    
def restart2solution(config,state={}):
    """ restart2solution(config,state={})
        moves restart file to solution file, 
        optionally updates state
        direct or adjoint is read from config
        adjoint objective is read from config
    """

    # direct solution
    if config.MATH_PROBLEM == 'DIRECT':
        restart  = config.RESTART_FLOW_FILENAME
        solution = config.SOLUTION_FLOW_FILENAME        
        # expand unsteady time
        restarts  = expand_time(restart,config)
        solutions = expand_time(solution,config)
        # move
        for res,sol in zip(restarts,solutions):
            shutil.move( res , sol )
        # update state
        if state: state.FILES.DIRECT = solution
        
    # adjoint solution
    elif any([config.MATH_PROBLEM == 'CONTINUOUS_ADJOINT', config.MATH_PROBLEM == 'DISCRETE_ADJOINT']):
        restart  = config.RESTART_ADJ_FILENAME
        solution = config.SOLUTION_ADJ_FILENAME           
        # add suffix
        func_name = config.OBJECTIVE_FUNCTION
        suffix    = get_adjointSuffix(func_name)
        restart   = add_suffix(restart,suffix)
        solution  = add_suffix(solution,suffix)
        # expand unsteady time
        restarts  = expand_time(restart,config)
        solutions = expand_time(solution,config)        
        # move
        for res,sol in zip(restarts,solutions):
            shutil.move( res , sol )
        # udpate state
        ADJ_NAME = 'ADJOINT_' + func_name
        if state: state.FILES[ADJ_NAME] = solution
        
    else:
        raise Exception, 'unknown math problem'

