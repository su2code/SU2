#!/usr/bin/env python

## \file tools.py
#  \brief file i/o functions
#  \author T. Lukaczyk, F. Palacios
#  \version 6.2.0 "Falcon"
#
# The current SU2 release has been coordinated by the
# SU2 International Developers Society <www.su2devsociety.org>
# with selected contributions from the open-source community.
#
# The main research teams contributing to the current release are:
#  - Prof. Juan J. Alonso's group at Stanford University.
#  - Prof. Piero Colonna's group at Delft University of Technology.
#  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
#  - Prof. Rafael Palacios' group at Imperial College London.
#  - Prof. Vincent Terrapon's group at the University of Liege.
#  - Prof. Edwin van der Weide's group at the University of Twente.
#  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
#
# Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
#                      Tim Albring, and the SU2 contributors.
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

import os
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
        raise IOError('multiple zones not supported')
    
    # done
    plot_file.close()              
    return plot_data


# -------------------------------------------------------------------
#  Read All Data from History File
# -------------------------------------------------------------------

def read_history( History_filename, nZones = 1):
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
    map_dict = get_headerMap(nZones)    
    
    # map header names
    for key in plot_data.keys():
        if key in map_dict:
            var = map_dict[key]
        else:
            var = key
        history_data[var] = plot_data[key]
    
    return history_data
    
#: def read_history()



# -------------------------------------------------------------------
#  Define Dictionary Map for Header Names
# -------------------------------------------------------------------

def get_headerMap(nZones = 1):

    """ returns a dictionary that maps history file header names
        to optimization problem function names
    """
    # header name to config file name map
    history_header_map = { "Iteration"       : "ITERATION"               ,
                 "CL"              : "LIFT"                    ,
                 "CD"              : "DRAG"                    ,
                 "CSF"             : "SIDEFORCE"               ,
                 "Cp_Diff"         : "INVERSE_DESIGN_PRESSURE" ,
                 "HeatFlux_Diff"   : "INVERSE_DESIGN_HEATFLUX" ,
                 "HeatFlux_Total"  : "TOTAL_HEATFLUX"          ,
                 "HeatFlux_Maximum": "MAXIMUM_HEATFLUX"        ,
                 "Temperature_Total": "TOTAL_TEMPERATURE"        ,
                 "CMx"             : "MOMENT_X"                ,
                 "CMy"             : "MOMENT_Y"                ,
                 "CMz"             : "MOMENT_Z"                ,
                 "CFx"             : "FORCE_X"                 ,
                 "CFy"             : "FORCE_Y"                 ,
                 "CFz"             : "FORCE_Z"                 ,
                 "CL/CD"           : "EFFICIENCY"              ,
                 "AoA"             : "AOA"                     ,
                 "Custom_ObjFunc"  : "CUSTOM_OBJFUNC"          ,
                 "CMerit"          : "FIGURE_OF_MERIT"         ,
                 "Buffet_Metric"   : "BUFFET"                  ,
                 "CQ"              : "TORQUE"                  ,
                 "CT"              : "THRUST"                  ,
                 "CEquivArea"      : "EQUIVALENT_AREA"         ,
                 "CNearFieldOF"    : "NEARFIELD_PRESSURE"      ,
                 "Avg_TotalPress"  : "SURFACE_TOTAL_PRESSURE"  ,
                 "Avg_Press"       : "SURFACE_STATIC_PRESSURE" ,
                 "Avg_MassFlow"    : "SURFACE_MASSFLOW"        ,
                 "Avg_Mach"        : "SURFACE_MACH"            ,
                 "Uniformity"                : "SURFACE_UNIFORMITY"            ,
                 "Secondary_Strength"        : "SURFACE_SECONDARY"            ,
                 "Momentum_Distortion"       : "SURFACE_MOM_DISTORTION"            ,
                 "Secondary_Over_Uniformity" : "SURFACE_SECOND_OVER_UNIFORM"            ,
                 "Pressure_Drop"        : "SURFACE_PRESSURE_DROP"            ,
                 "ComboObj"        : "COMBO"                   ,
                 "Time(min)"       : "TIME"                    ,
                 'Time(min)"\n'    : "TIME"                    , # TDE hack for paraview
                 "D(CL)"           : "D_LIFT"                  ,
                 "D(CD)"           : "D_DRAG"                  ,
                 "D(CSF)"          : "D_SIDEFORCE"             ,
                 "D(CMx)"          : "D_MOMENT_X"              ,
                 "D(CMy)"          : "D_MOMENT_Y"              ,
                 "D(CMz)"          : "D_MOMENT_Z"              ,
                 "D(CFx)"          : "D_FORCE_X"               ,
                 "D(CFy)"          : "D_FORCE_Y"               ,
                 "D(CFz)"          : "D_FORCE_Z"               ,
                 "D(CL/CD)"        : "D_EFFICIENCY"            ,
                 "D(Custom_ObjFunc)" : "D_CUSTOM_OBJFUNC"      ,
                 "D(HeatFlux_Total)" : "D_HEAT"                ,
                 "D(HeatFlux_Maximum)" : "D_HEAT_MAX"          ,
                 "TotalPressureLoss_1"     : "TOTAL_PRESSURE_LOSS"    ,
                 "KineticEnergyLoss_1"     : "KINETIC_ENERGY_LOSS"    ,
                 "EntropyGen_" + str(getTurboPerfIndex(nZones)) : "ENTROPY_GENERATION"     ,                   
                 "FlowAngleOut_1"          : "FLOW_ANGLE_OUT"         ,
                 "FlowAngleIn_1"           : "FLOW_ANGLE_IN"          ,
                 "MassFlowIn_1"            : "MASS_FLOW_IN"           ,
                 "MassFlowOut_1"           : "MASS_FLOW_OUT"          ,
                 "PressureRatio_1"         : "PRESSURE_RATIO"         ,
                 "TotalEfficiency_" + str(getTurboPerfIndex(nZones))  : "TOTAL_EFFICIENCY"       ,
                 "TotalStaticEfficiency_3" : "TOTAL_STATIC_EFFICIENCY",
                 "D(TotalPressureLoss_0)"  : "D_TOTAL_PRESSURE_LOSS"  ,
                 "D(TotalEfficiency_0)"       : "D_TOTAL_EFFICIENCY"       ,
                 "D(TotalPressureLoss_0)"     : "D_TOTAL_PRESSURE_LOSS"    ,
                 "D(KineticEnergyLoss_0)"     : "D_KINETIC_ENERGY_LOSS"    ,
                 "D(TotalStaticEfficiency_0)" : "D_TOTAL_STATIC_EFFICIENCY",
                 "D(FlowAngleOut_0)"          : "D_FLOW_ANGLE_OUT"         ,
                 "D(FlowAngleIn_0)"           : "D_FLOW_ANGLE_IN"          ,
                 "D(MassFlowIn_0)"            : "D_MASS_FLOW_IN"           ,
                 "D(MassFlowOut_0)"           : "D_MASS_FLOW_OUT"          ,
                 "D(PressureRatio_0)"         : "D_PRESSURE_RATIO"         ,
                 "D(EnthalpyOut_0)"           : "D_ENTHALPY_OUT"           ,
                 "D(TotalEnthalpy_0)"         : "D_TOTAL_ENTHALPY_OUT"     ,
                 "D(Uniformity)"                : "D_SURFACE_UNIFORMITY"            ,
                 "D(Secondary_Strength)"        : "D_SURFACE_SECONDARY"             ,
                 "D(Momentum_Distortion)"       : "D_SURFACE_MOM_DISTORTION"        ,
                 "D(Secondary_Over_Uniformity)" : "D_SURFACE_SECOND_OVER_UNIFORM"   ,
                 "D(Pressure_Drop)"             : "D_SURFACE_PRESSURE_DROP"         }
 
    return history_header_map        

def getTurboPerfIndex(nZones = 1):

  if int(nZones) > 1:
    index = int(nZones) + int(int(nZones)/2.0) + 1
  else: 
    index = 1
  return index


#: def get_headerMap()


# -------------------------------------------------------------------
#  Optimizer Function Names
# -------------------------------------------------------------------

# Aerodynamic Optimizer Function Names

optnames_aero = [ "LIFT"                        ,
                  "DRAG"                        ,
                  "SIDEFORCE"                   ,
                  "MOMENT_X"                    ,
                  "MOMENT_Y"                    ,
                  "MOMENT_Z"                    ,
                  "FORCE_X"                     ,
                  "FORCE_Y"                     ,
                  "FORCE_Z"                     ,
                  "EFFICIENCY"                  ,
                  "FIGURE_OF_MERIT"             ,
                  "BUFFET"                      ,
                  "TORQUE"                      ,
                  "THRUST"                      ,
                  "SURFACE_TOTAL_PRESSURE"      ,
                  "SURFACE_STATIC_PRESSURE"     ,
                  "SURFACE_MASSFLOW"            ,
                  "SURFACE_MACH"                ,
                  "SURFACE_UNIFORMITY"          ,
                  "SURFACE_SECONDARY"           ,
                  "SURFACE_MOM_DISTORTION"      ,
                  "SURFACE_SECOND_OVER_UNIFORM" ,
                  "SURFACE_PRESSURE_DROP"       ,
                  "EQUIVALENT_AREA"             ,
                  "NEARFIELD_PRESSURE"          ,
                  "INVERSE_DESIGN_PRESSURE"     ,
                  "INVERSE_DESIGN_HEATFLUX"     ,
                  "TOTAL_HEATFLUX"              ,
                  "MAXIMUM_HEATFLUX"            ,
                  "CUSTOM_OBJFUNC"              ,
                  "COMBO"]

# Turbo performance optimizer Function Names
optnames_turbo = ["TOTAL_PRESSURE_LOSS"     ,
                  "KINETIC_ENERGY_LOSS"     ,
                  "ENTROPY_GENERATION"      ,
                  "EULERIAN_WORK"           ,
                  "FLOW_ANGLE_IN"           ,
                  "FLOW_ANGLE_OUT"          ,
                  "MASS_FLOW_IN"            ,
                  "MASS_FLOW_OUT"           ,
                  "PRESSURE_RATIO"          ,
                  "TOTAL_EFFICIENCY"        ,
                  "TOTAL_STATIC_EFFICIENCY" ,
                 ]
#: optnames_aero

optnames_stab = [ "D_LIFT_D_ALPHA"               ,
                  "D_DRAG_D_ALPHA"               ,
                  "D_SIDEFORCE_D_ALPHA"          ,
                  "D_MOMENT_X_D_ALPHA"           ,
                  "D_MOMENT_Y_D_ALPHA"           ,
                  "D_MOMENT_Z_D_ALPHA"           ,
                ]

#: Multipoint Optimizer Function Names

# optnames_multi = ['{}_{}'.format('MULTIPOINT', a) for a in optnames_aero]

optnames_multi = [ "MULTIPOINT_LIFT"               ,
                   "MULTIPOINT_DRAG"               ,
                   "MULTIPOINT_SIDEFORCE"          ,
                   "MULTIPOINT_MOMENT_X"           ,
                   "MULTIPOINT_MOMENT_Y"           ,
                   "MULTIPOINT_MOMENT_Z"           ,
                   "MULTIPOINT_CUSTOM_OBJFUNC"]

# Geometric Optimizer Function Names
optnames_geo = [ "AIRFOIL_AREA"                   ,
                 "AIRFOIL_THICKNESS"              ,
                 "AIRFOIL_CHORD"                  ,
                 "AIRFOIL_LE_RADIUS"              ,
                 "AIRFOIL_TOC"                    ,
                 "AIRFOIL_ALPHA"                  ,
                 "FUSELAGE_VOLUME"        ,
                 "FUSELAGE_WETTED_AREA"   ,
                 "FUSELAGE_MIN_WIDTH"     ,
                 "FUSELAGE_MAX_WIDTH"     ,
                 "FUSELAGE_MIN_WATERLINE_WIDTH"  ,
                 "FUSELAGE_MAX_WATERLINE_WIDTH"  ,
                 "FUSELAGE_MIN_HEIGHT"    ,
                 "FUSELAGE_MAX_HEIGHT"    ,
                 "FUSELAGE_MAX_CURVATURE" ,
                 "WING_VOLUME"            ,
                 "WING_MIN_THICKNESS" ,
                 "WING_MAX_THICKNESS" ,
                 "WING_MIN_CHORD"         ,
                 "WING_MAX_CHORD"         ,
                 "WING_MIN_LE_RADIUS"     ,
                 "WING_MAX_LE_RADIUS"     ,
                 "WING_MIN_TOC"           ,
                 "WING_MAX_TOC"           ,
                 "WING_OBJFUN_MIN_TOC"    ,
                 "WING_MAX_TWIST"         ,
                 "WING_MAX_CURVATURE"     ,
                 "WING_MAX_DIHEDRAL"           ,
                 "NACELLE_VOLUME"            ,
                 "NACELLE_MIN_THICKNESS"     ,
                 "NACELLE_MAX_THICKNESS"     ,
                 "NACELLE_MIN_CHORD"         ,
                 "NACELLE_MAX_CHORD"         ,
                 "NACELLE_MIN_LE_RADIUS"     ,
                 "NACELLE_MAX_LE_RADIUS"     ,
                 "NACELLE_MIN_TOC"           ,
                 "NACELLE_MAX_TOC"           ,
                 "NACELLE_OBJFUN_MIN_TOC"    ,
                 "NACELLE_MAX_TWIST"         ]
                 
PerStation = []
for i in range(20):
    PerStation.append("STATION" + str(i) + "_AREA")
    PerStation.append("STATION" + str(i) + "_LENGTH")
    PerStation.append("STATION" + str(i) + "_WIDTH")
    PerStation.append("STATION" + str(i) + "_WATERLINE_WIDTH")
    PerStation.append("STATION" + str(i) + "_HEIGHT")
    PerStation.append("STATION" + str(i) + "_THICKNESS")
    PerStation.append("STATION" + str(i) + "_CHORD")
    PerStation.append("STATION" + str(i) + "_LE_RADIUS")
    PerStation.append("STATION" + str(i) + "_TOC")
    PerStation.append("STATION" + str(i) + "_TWIST")

optnames_geo.extend(PerStation)
                 
#: optnames_geo

grad_names_directdiff = ["D_LIFT",
                         "D_DRAG",
                         "D_SIDEFORCE",
                         "D_MOMENT_X",
                         "D_MOMENT_Y",
                         "D_MOMENT_Z",
                         "D_FORCE_X",
                         "D_FORCE_Y",
                         "D_FORCE_Z",
                         "D_EFFICIENCY",
                         "D_CUSTOM_OBJFUNC",
                         "D_HEAT",
                         "D_MAX_HEAT",
                         "D_TOTAL_PRESSURE_LOSS",
                         "D_TOTAL_EFFICIENCY",
                         "D_TOTAL_PRESSURE_LOSS",
                         "D_KINETIC_ENERGY_LOSS",
                         "D_TOTAL_STATIC_EFFICIENCY",
                         "D_FLOW_ANGLE_OUT",
                         "D_FLOW_ANGLE_IN",
                         "D_MASS_FLOW_IN",
                         "D_MASS_FLOW_OUT",
                         "D_PRESSURE_RATIO",
                         "D_ENTHALPY_OUT",
                         "D_TOTAL_ENTHALPY_OUT",
                         "D_SURFACE_UNIFORMITY",
                         "D_SURFACE_SECONDARY",
                         "D_SURFACE_MOM_DISTORTION",
                         "D_SURFACE_SECOND_OVER_UNIFORM",
                         "D_SURFACE_PRESSURE_DROP"]

grad_names_map = ordered_bunch()
grad_names_map.MASS_FLOW_IN = "D_MASS_FLOW_IN"
grad_names_map.MOMENT_Z = "D_MOMENT_Z"
grad_names_map.FLOW_ANGLE_OUT = "D_FLOW_ANGLE_OUT"
grad_names_map.MASS_FLOW_OUT = "D_MASS_FLOW_OUT"
grad_names_map.FLOW_ANGLE_IN = "D_FLOW_ANGLE_IN"
grad_names_map.FORCE_Z = "D_FORCE_Z"
grad_names_map.FORCE_Y = "D_FORCE_Y"
grad_names_map.FORCE_X = "D_FORCE_X"
grad_names_map.TOTAL_EFFICIENCY = "D_TOTAL_EFFICIENCY"
grad_names_map.TOTAL_STATIC_EFFICIENCY = "D_TOTAL_STATIC_EFFICIENCY"
grad_names_map.PRESSURE_RATIO = "D_PRESSURE_RATIO"
grad_names_map.EFFICIENCY = "D_EFFICIENCY"
grad_names_map.DRAG = "D_DRAG"
grad_names_map.LIFT = "D_LIFT"
grad_names_map.TOTAL_ENTHALPY_OUT = "D_TOTAL_ENTHALPY_OUT"
grad_names_map.TOTAL_PRESSURE_LOSS = "D_TOTAL_PRESSURE_LOSS"
grad_names_map.MOMENT_Y = "D_MOMENT_Y"
grad_names_map.MOMENT_X="D_MOMENT_X"
grad_names_map.SIDEFORCE = "D_SIDEFORCE"
grad_names_map.ENTHALPY_OUT = "D_ENTHALPY_OUT"
grad_names_map.KINETIC_ENERGY_LOSS = "D_KINETIC_ENERGY_LOSS"
grad_names_map.CUSTOM_OBJFUNC = "D_CUSTOM_OBJFUNC"
grad_names_map.HEAT = "D_HEAT"
grad_names_map.MAX_HEAT = "D_MAX_HEAT"
grad_names_map.SURFACE_UNIFORMITY = "D_SURFACE_UNIFORMITY"
grad_names_map.SURFACE_SECONDARY = "D_SURFACE_SECONDARY"
grad_names_map.SURFACE_MOM_DISTORTION = "D_SURFACE_MOM_DISTORTION"
grad_names_map.SURFACE_SECOND_OVER_UNIFORM = "D_SURFACE_SECOND_OVER_UNIFORM"
grad_names_map.SURFACE_PRESSURE_DROP = "D_SURFACE_PRESSURE_DROP"

# per-surface functions
per_surface_map = {"LIFT"       :   "CL" ,
                  "DRAG"        :   "CD" ,
                  "SIDEFORCE"   :   "CSF"  ,
                  "MOMENT_X"    :   "CMx"   ,
                  "MOMENT_Y"    :   "CMy"   ,
                  "MOMENT_Z"    :   "CMz"   ,
                  "FORCE_X"     :   "CFx"   ,
                  "FORCE_Y"     :   "CFy"   ,
                  "FORCE_Z"     :   "CFz"   ,
                  "EFFICIENCY"  :   "CL/CD" }

# -------------------------------------------------------------------
#  Include per-surface output from History File
# ------------------------------------------------------------------- 
def update_persurface(config, state):
    # Update the header map (checking to make sure entries are not duplicated)
    header_map = get_headerMap()
    for base in per_surface_map:
        base2 = per_surface_map[base]
        for marker in config['MARKER_MONITORING']:
            if not (base2+'_'+marker) in header_map:
                header_map[base2+'_'+marker] = base2+'_'+marker
    # Update the function values in state to include the per-surface quantities
    if 'DIRECT' in state['HISTORY']:
        for base in per_surface_map:
            base2 = per_surface_map[base]
            for marker in config['MARKER_MONITORING']:
                if (base2+'_'+marker) in state['HISTORY']['DIRECT']:
                    state['FUNCTIONS'][base2+'_'+marker] = state['HISTORY']['DIRECT'][base2+'_'+marker][-1]
                    
# -------------------------------------------------------------------
#  Read Aerodynamic Function Values from History File
# -------------------------------------------------------------------

def read_aerodynamics( History_filename , nZones = 1, special_cases=[], final_avg=0 ):
    """ values = read_aerodynamics(historyname, special_cases=[])
        read aerodynamic function values from history file
        
        Outputs:
            dictionary with function keys and thier values
            if special cases has 'UNSTEADY_SIMULATION', returns time averaged data
            otherwise returns final value from history file
    """
    
    # read the history data
    history_data = read_history(History_filename, nZones)
    
    # list of functions to pull
    func_names = optnames_aero + grad_names_directdiff + optnames_turbo

    # pull only these functions
    Func_Values = ordered_bunch()
    for this_objfun in func_names:
        if this_objfun in history_data:
            Func_Values[this_objfun] = history_data[this_objfun] 
    
    # for unsteady cases, average time-accurate objective function values
    if 'UNSTEADY_SIMULATION' in special_cases and not final_avg:
        for key,value in Func_Values.items():
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
            SURFACE_TOTAL_PRESSURE
            SURFACE_STATIC_PRESSURE
            SURFACE_MASSFLOW
            SURFACE_MACH
            TOTAL_STATIC_EFFICIENCY
        returns +1 otherwise
    """
    
    # flip sign for maximization problems
    if ObjFun_name == "LIFT"            : return -1.0
    if ObjFun_name == "EFFICIENCY"      : return -1.0
    if ObjFun_name == "THRUST"          : return -1.0
    if ObjFun_name == "FIGURE_OF_MERIT" : return -1.0
    if ObjFun_name == "SURFACE_TOTAL_PRESSURE"  : return -1.0
    if ObjFun_name == "SURFACE_STATIC_PRESSURE" : return -1.0
    if ObjFun_name == "SURFACE_MASSFLOW"        : return -1.0
    if ObjFun_name == "SURFACE_MACH"            : return -1.0
    if ObjFun_name == "TOTAL_STATIC_EFFICIENCY" :return -1.0
    
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
    name_map = { "DRAG"                        : "cd"        ,
                 "LIFT"                        : "cl"        ,
                 "SIDEFORCE"                   : "csf"       ,
                 "MOMENT_X"                    : "cmx"       ,
                 "MOMENT_Y"                    : "cmy"       ,
                 "MOMENT_Z"                    : "cmz"       ,
                 "FORCE_X"                     : "cfx"       ,
                 "FORCE_Y"                     : "cfy"       ,
                 "FORCE_Z"                     : "cfz"       ,
                 "EFFICIENCY"                  : "eff"       ,
                 "INVERSE_DESIGN_PRESSURE"     : "invpress"  ,
                 "INVERSE_DESIGN_HEAT"         : "invheat"   ,
                 "MAXIMUM_HEATFLUX"            : "maxheat"   ,
                 "TOTAL_HEATFLUX"              : "totheat"   ,
                 "EQUIVALENT_AREA"             : "ea"        ,
                 "NEARFIELD_PRESSURE"          : "nfp"       ,
                 "THRUST"                      : "ct"        ,
                 "TORQUE"                      : "cq"        ,
                 "FIGURE_OF_MERIT"             : "merit"     ,
                 "BUFFET"                      : "buffet"    ,
                 "SURFACE_TOTAL_PRESSURE"      : "pt"        ,
                 "SURFACE_STATIC_PRESSURE"     : "pe"        ,
                 "SURFACE_MASSFLOW"            : "mfr"       ,
                 "SURFACE_MACH"                : "mach"      ,
                 "SURFACE_UNIFORMITY"          : "uniform"   ,
                 "SURFACE_SECONDARY"           : "second"    ,
                 "SURFACE_MOM_DISTORTION"      : "distort"   ,
                 "SURFACE_SECOND_OVER_UNIFORM" : "sou"       ,
                 "SURFACE_PRESSURE_DROP"       : "dp"        ,
                 "CUSTOM_OBJFUNC"              : "custom"    ,
                 "KINETIC_ENERGY_LOSS"         : "ke"        ,
                 "TOTAL_PRESSURE_LOSS"         : "pl"        ,
                 "ENTROPY_GENERATION"          : "entg"      ,
                 "EULERIAN_WORK"               : "ew"        ,
                 "FLOW_ANGLE_OUT"              : "fao"       ,
                 "FLOW_ANGLE_IN"               : "fai"       ,
                 "MASS_FLOW_OUT"               : "mfo"       ,
                 "MASS_FLOW_IN"                : "mfi"       ,
                 "TOTAL_EFFICIENCY"            : "teff"      ,
                 "TOTAL_STATIC_EFFICIENCY"     : "tseff"     ,
                 "COMBO"                       : "combo"}
    
    # if none or false, return map
    if not objective_function:
        return name_map
    else:
        # remove white space
        objective = ''.join(objective_function.split())
        objective = objective.split(",")
        nObj = len(objective)
        if (nObj>1):
            return "combo"
        if objective[0] in name_map:
            return name_map[objective[0]]
    
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
               2   : "SURFACE_BUMP"          ,
               4   : "NACA_4DIGITS"          ,
               5   : "TRANSLATION"           ,
               6   : "ROTATION"              ,
               7   : "FFD_CONTROL_POINT"     ,
               8   : "FFD_DIHEDRAL_ANGLE"    ,
               9   : "FFD_TWIST_ANGLE"       ,
               10  : "FFD_ROTATION"          ,
               11  : "FFD_CAMBER"            ,
               12  : "FFD_THICKNESS"         ,
               19  : "FFD_TWIST"             ,
               22  : "FFD_NACELLE"           ,
               23  : "FFD_GULL"              ,
               25  : "FFD_ROTATION"          ,
               15  : "FFD_CONTROL_POINT_2D"  ,
               16  : "FFD_CAMBER_2D"         ,
               17  : "FFD_THICKNESS_2D"      ,
               20  : "FFD_TWIST_2D"          ,
               50  : "CUSTOM"                ,
               51  : "CST"                   ,
               101 : "ANGLE_OF_ATTACK"       ,
               102 : "FFD_ANGLE_OF_ATTACK"                    }
    
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
    id_map = dict((v,k) for (k,v) in dv_map.items())
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
    if (plot_format == 'TECPLOT') or (plot_format == 'TECPLOT_BINARY'): 
        header.append('VARIABLES=')
    elif (plot_format == 'PARAVIEW') or (plot_format == 'PARAVIEW_BINARY'):
        pass
    else: raise Exception('output plot format not recognized')
    
    # Case: continuous adjoint
    if grad_type == 'CONTINUOUS_ADJOINT':
        header.append(r'"iVar","Gradient","FinDiff_Step"')
        write_format.append(r'%4d, %.10f, %f')
        
    # Case: finite difference  
    elif grad_type == 'FINITE_DIFFERENCE':
        header.append(r'"iVar","Grad_CL","Grad_CD","Grad_CSF","Grad_CMx","Grad_CMy","Grad_CMz","Grad_CFx","Grad_CFy","Grad_CFz","Grad_CL/CD","Grad_Custom_ObjFunc","Grad_HeatFlux_Total","Grad_HeatFlux_Maximum","Grad_Temperature_Total"')
        write_format.append(r'%4d, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f')

        for key in special_cases: 
            if key == "ROTATING_FRAME" : 
                header.append(r',"Grad_CMerit","Grad_CT","Grad_CQ"')
                write_format.append(", %.10f, %.10f, %.10f")
            if key == "EQUIV_AREA"     : 
                header.append(r',"Grad_CEquivArea","Grad_CNearFieldOF"') 
                write_format.append(", %.10f, %.10f")
            if key == "ENGINE"     :
                header.append(r',"Grad_AeroCDrag","Grad_SolidCDrag","Grad_Radial_Distortion","Grad_Circumferential_Distortion"')
                write_format.append(", %.10f, %.10f, %.10f, %.10f")
            if key == "1D_OUTPUT"     :
                header.append(r',"Grad_Avg_TotalPress","Grad_Avg_Mach","Grad_Avg_Temperature","Grad_MassFlowRate","Grad_Avg_Pressure","Grad_Avg_Density","Grad_Avg_Velocity","Grad_Avg_Enthalpy"')
                write_format.append(", %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f")
            if key == "INV_DESIGN_CP"     :
                header.append(r',"Grad_Cp_Diff"')
                write_format.append(", %.10f")
            if key == "INV_DESIGN_HEATFLUX"     :
                header.append(r',"Grad_HeatFlux_Diff"')
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
    elif kindID == "SURFACE_BUMP"        :
        header.append(r',"Loc_Start","Loc_End","Loc_Max"')
        write_format.append(r', %s, %s, %s')
    elif kindID == "CST"        :
        header.append(r',"Up/Down","Kulfan number", "Total Kulfan numbers"')
        write_format.append(r', %s, %s', '%s')
    elif kindID == "FAIRING"       :
        header.append(r',"ControlPoint_Index","Theta_Disp","R_Disp"')
        write_format.append(r', %s, %s, %s')
    elif kindID == "NACA_4DIGITS"       :
        header.append(r',"1st_digit","2nd_digit","3rd&4th_digits"')
        write_format.append(r', %s, %s, %s')
    elif kindID == "TRANSLATION"       : 
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
    elif kindID == "ANGLE_OF_ATTACK"      : pass
    elif kindID == "FFD_ANGLE_OF_ATTACK"  : pass
    
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
    
def get_optFileFormat(plot_format,special_cases=None, nZones = 1):
    
    if special_cases is None: special_cases = []
    
    # start header, build a list of strings and join at the end
    header_list   = []
    header_format = ''
    write_format  = []
    
    # handle plot formating
    if (plot_format == 'TECPLOT') or (plot_format == 'TECPLOT_BINARY'): 
        header_format = header_format + 'VARIABLES='
    elif (plot_format == 'PARAVIEW') or (plot_format == 'PARAVIEW_BINARY'):
        pass
    else: raise Exception('output plot format not recognized')

    # start header
    header_list.extend(["Iteration","CL","CD","CSF","CMx","CMy","CMz","CFx","CFy","CFz","CL/CD","Custom_ObjFunc","HeatFlux_Total","HeatFlux_Maximum","Temperature_Total"])
    write_format.append(r'%4d, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f')
        
    # special cases
    for key in special_cases: 
        if key == "ROTATING_FRAME" : 
            header_list.extend(["CMerit","CT","CQ"])
            write_format.append(r', %.10f, %.10f, %.10f')
        if key == "EQUIV_AREA"     : 
            header_list.extend(["CEquivArea","CNearFieldOF"]) 
            write_format.append(r', %.10f, %.10f')
        if key == "ENGINE"     :
            header_list.extend(["AeroCDrag","SolidCDrag","Radial_Distortion","Circumferential_Distortion"])
            write_format.append(r', %.10f, %.10f, %.10f, %.10f')
        if key == "1D_OUTPUT":
            header_list.extend(["AreaAvg_TotalPress","AreaAvg_Mach","AreaAvg_Temperature","MassFlowRate","Avg_Pressure","Avg_Density","Avg_Velocity","Avg_Enthalpy"])
            write_format.append(r', %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f')
        if key == "INV_DESIGN_CP"     :
            header_list.extend(["Cp_Diff"])
            write_format.append(r', %.10f')
        if key == "INV_DESIGN_HEATFLUX"     :
            header_list.extend(["HeatFlux_Diff"])
            write_format.append(r', %.10f')

    # finish formats
    header_format = (header_format) + ('"') + ('","').join(header_list) + ('"') + (' \n')
    write_format  = ''.join(write_format)  + ' \n'
            
    # build list of objective function names
    header_vars = []
    map_dict = get_headerMap(nZones)
    for variable in header_list:
        assert variable in map_dict, 'unrecognized header variable'
        header_vars.append(map_dict[variable])
    
    # done
    return [header_format,header_vars,write_format]
        
#: def get_optFileFormat()
  
  
  
# -------------------------------------------------------------------
#  Get Extension Name
# -------------------------------------------------------------------

def get_extension(output_format):
  
    if (output_format == "PARAVIEW")        : return ".csv"
    if (output_format == "PARAVIEW_BINARY") : return ".csv"
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
    
    all_special_cases = [ 'ROTATING_FRAME'                   ,
                          'EQUIV_AREA'                       ,
                          '1D_OUTPUT'                        ,
                          'INV_DESIGN_CP'                    ,
                          'INV_DESIGN_HEATFLUX'              ]
    
    special_cases = []
    for key in all_special_cases:
        if key in config and config[key] == 'YES':
            special_cases.append(key)
        if 'PHYSICAL_PROBLEM' in config and config['PHYSICAL_PROBLEM'] == key:
            special_cases.append(key)
            
    if config.get('UNSTEADY_SIMULATION','NO') != 'NO':
        special_cases.append('UNSTEADY_SIMULATION')
     
    # no support for more than one special case
    if len(special_cases) > 1:
        error_str = 'Currently cannot support ' + ' and '.join(special_cases) + ' at once'
        raise Exception(error_str)   
    
    if (config['WRT_SOL_FREQ'] != 1) and ('WRT_UNSTEADY' in special_cases):
        raise Exception('Must set WRT_SOL_FREQ= 1 for WRT_UNSTEADY= YES')
  
    # Special case for harmonic balance
    if 'UNSTEADY_SIMULATION' in config and config['UNSTEADY_SIMULATION'] == 'HARMONIC_BALANCE':
        special_cases.append('HARMONIC_BALANCE')

    # Special case for rotating frame
    if 'GRID_MOVEMENT_KIND' in config and config['GRID_MOVEMENT_KIND'] == 'ROTATING_FRAME':
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
        if 'PHYSICAL_PROBLEM' in config and config['PHYSICAL_PROBLEM'] == key:
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
        if not isinstance(name, list):
            name_pat = add_suffix(name,'%05d')
            names = [name_pat%i for i in range(n_time)]
        else:
            for n in range(len(name)):
                name_pat = add_suffix(name[n], '%05d')
                names    = [name_pat%i for i in range(n_time)]
    else:
        if not isinstance(name, list):
            names = [name]
        else:
            names = name
    return names

def expand_zones(name, config):
    if int(config.NZONES) > 1:
        if not isinstance(name, list):
            name_pat = add_suffix(name,'%d')
            names = [name_pat%i for i in range(int(config.NZONES))]
        else:
            for n in range(len(name)):
                name_pat[i] = add_suffix(name, '%d')
                names[i]    = [name_pat%i for i in range(int(config.NZONES))]
    else:
        if not isinstance(name, list):
            names = [name]
        else:
            names = name
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
        
        # expand zones
        restarts  = expand_zones(restart,config)
        solutions = expand_zones(solution,config)
        # expand unsteady time
        restarts  = expand_time(restarts,config)
        solutions = expand_time(solutions,config)
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
        # expand zones
        restarts  = expand_zones(restart,config)
        solutions = expand_zones(solution,config)
        # expand unsteady time
        restarts  = expand_time(restarts,config)
        solutions = expand_time(solutions,config)
        # move
        for res,sol in zip(restarts,solutions):
            shutil.move( res , sol )
        # udpate state
        if "," in func_name:
            func_name="COMBO"
        ADJ_NAME = 'ADJOINT_' + func_name
        if state: state.FILES[ADJ_NAME] = solution
        
    else:
        raise Exception('unknown math problem')

