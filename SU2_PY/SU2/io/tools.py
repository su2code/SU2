#!/usr/bin/env python

## \file tools.py
#  \brief file i/o functions
#  \author T. Lukaczyk, F. Palacios
#  \version 8.0.0 "Harrier"
#
# SU2 Project Website: https://su2code.github.io
#
# The SU2 Project is maintained by the SU2 Foundation
# (http://su2foundation.org)
#
# Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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
from .historyMap import history_header_map as historyOutFields

# -------------------------------------------------------------------
#  Read SU2_DOT Gradient Values
# -------------------------------------------------------------------


def read_gradients(Grad_filename, scale=1.0):
    """reads the raw gradients from the gradient file
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


def read_plot(filename):
    """reads a plot file
    returns an ordered bunch with the headers for keys
    and a list of each header's floats for values.
    """

    extension = os.path.splitext(filename)[1]

    # open history file
    plot_file = open(filename)

    # title?
    line = plot_file.readline()
    if line.startswith("TITLE"):
        title = line.split("=")[1].strip()  # not used right now
        line = plot_file.readline()

    if line.startswith("VARIABLES"):
        line = plot_file.readline()

    line = line.split(",")
    Variables = [x.strip().strip('"') for x in line]
    n_Vars = len(Variables)

    # initialize plot data dictionary
    plot_data = ordered_bunch.fromkeys(Variables)
    # must default each value to avoid pointer problems
    for key in plot_data.keys():
        plot_data[key] = []

    # zone list
    zones = []

    # read all data rows
    while 1:
        # read line
        line = plot_file.readline()
        if not line:
            break

        # zone?
        if line.startswith("ZONE"):
            zone = line.split("=")[1].strip('" ')
            zones.append(zone)
            continue

        # split line
        line_data = line.strip().split(",")
        line_data = [float(x.strip()) for x in line_data]

        # store to dictionary
        for i_Var in range(n_Vars):
            this_variable = Variables[i_Var]
            plot_data[this_variable] = plot_data[this_variable] + [line_data[i_Var]]

    #: for each line

    # check for number of zones
    if len(zones) > 1:
        raise IOError("multiple zones not supported")

    # done
    plot_file.close()
    return plot_data


# -------------------------------------------------------------------
#  Read All Data from History File
# -------------------------------------------------------------------


def read_history(History_filename, nZones=1):
    """reads a history file
    returns an ordered bunch with the history file headers for keys
    and a list of each header's floats for values.
    if header is an optimization objective, its name is mapped to
    the optimization name.
    Iter and Time(min) headers are mapped to ITERATION and TIME
    respectively.
    """

    # read plot file
    plot_data = read_plot(History_filename)

    # initialize history data dictionary
    history_data = ordered_bunch()

    # map header names
    for key in plot_data.keys():
        var = key
        for field in historyOutFields:

            if key == historyOutFields[field]["HEADER"] and nZones == 1:
                var = field

            if key.split("[")[0] == historyOutFields[field]["HEADER"] and nZones > 1:
                var = field + "[" + key.split("[")[1]

        history_data[var] = plot_data[key]

    return history_data


#: def read_history()


# -------------------------------------------------------------------
#  Define Dictionary Map for Header Names
# -------------------------------------------------------------------


def get_headerMap(nZones=1):

    headerMap = dict()
    for outputField in historyOutFields:
        headerMap[outputField] = historyOutFields[outputField]["HEADER"]

    return headerMap


def getTurboPerfIndex(nZones=1):

    if int(nZones) > 1:
        index = int(nZones) + int(int(nZones) / 2.0) + 1
    else:
        index = 1
    return index


#: def get_headerMap()


# -------------------------------------------------------------------
#  Optimizer Function Names
# -------------------------------------------------------------------

#: optnames_stab

optnames_stab = [
    "D_LIFT_D_ALPHA",
    "D_DRAG_D_ALPHA",
    "D_SIDEFORCE_D_ALPHA",
    "D_MOMENT_X_D_ALPHA",
    "D_MOMENT_Y_D_ALPHA",
    "D_MOMENT_Z_D_ALPHA",
]

#: Multipoint Optimizer Function Names

# optnames_multi = ['{}_{}'.format('MULTIPOINT', a) for a in optnames_aero]

optnames_multi = [
    "MULTIPOINT_LIFT",
    "MULTIPOINT_DRAG",
    "MULTIPOINT_SIDEFORCE",
    "MULTIPOINT_MOMENT_X",
    "MULTIPOINT_MOMENT_Y",
    "MULTIPOINT_MOMENT_Z",
    "MULTIPOINT_CUSTOM_OBJFUNC",
]

# Geometric Optimizer Function Names
optnames_geo = [
    "AIRFOIL_AREA",
    "AIRFOIL_THICKNESS",
    "AIRFOIL_CHORD",
    "AIRFOIL_LE_RADIUS",
    "AIRFOIL_TOC",
    "AIRFOIL_ALPHA",
    "FUSELAGE_VOLUME",
    "FUSELAGE_WETTED_AREA",
    "FUSELAGE_MIN_WIDTH",
    "FUSELAGE_MAX_WIDTH",
    "FUSELAGE_MIN_WATERLINE_WIDTH",
    "FUSELAGE_MAX_WATERLINE_WIDTH",
    "FUSELAGE_MIN_HEIGHT",
    "FUSELAGE_MAX_HEIGHT",
    "FUSELAGE_MAX_CURVATURE",
    "WING_VOLUME",
    "WING_MIN_THICKNESS",
    "WING_MAX_THICKNESS",
    "WING_MIN_CHORD",
    "WING_MAX_CHORD",
    "WING_MIN_LE_RADIUS",
    "WING_MAX_LE_RADIUS",
    "WING_MIN_TOC",
    "WING_MAX_TOC",
    "WING_OBJFUN_MIN_TOC",
    "WING_MAX_TWIST",
    "WING_MAX_CURVATURE",
    "WING_MAX_DIHEDRAL",
    "NACELLE_VOLUME",
    "NACELLE_MIN_THICKNESS",
    "NACELLE_MAX_THICKNESS",
    "NACELLE_MIN_CHORD",
    "NACELLE_MAX_CHORD",
    "NACELLE_MIN_LE_RADIUS",
    "NACELLE_MAX_LE_RADIUS",
    "NACELLE_MIN_TOC",
    "NACELLE_MAX_TOC",
    "NACELLE_OBJFUN_MIN_TOC",
    "NACELLE_MAX_TWIST",
]

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

# per-surface functions
per_surface_map = {
    "LIFT": "CL",
    "DRAG": "CD",
    "SIDEFORCE": "CSF",
    "MOMENT_X": "CMx",
    "MOMENT_Y": "CMy",
    "MOMENT_Z": "CMz",
    "FORCE_X": "CFx",
    "FORCE_Y": "CFy",
    "FORCE_Z": "CFz",
    "EFFICIENCY": "CL/CD",
}

# -------------------------------------------------------------------
#  Include per-surface output from History File
# -------------------------------------------------------------------
def update_persurface(config, state):
    # Update the header map (checking to make sure entries are not duplicated)
    header_map = get_headerMap()
    for base in per_surface_map:
        base2 = per_surface_map[base]
        for marker in config["MARKER_MONITORING"]:
            if not (base2 + "_" + marker) in header_map:
                header_map[base2 + "_" + marker] = base2 + "_" + marker
    # Update the function values in state to include the per-surface quantities
    if "DIRECT" in state["HISTORY"]:
        for base in per_surface_map:
            base2 = per_surface_map[base]
            for marker in config["MARKER_MONITORING"]:
                if (base2 + "_" + marker) in state["HISTORY"]["DIRECT"]:
                    state["FUNCTIONS"][base2 + "_" + marker] = state["HISTORY"][
                        "DIRECT"
                    ][base2 + "_" + marker][-1]


# -------------------------------------------------------------------
#  Read Aerodynamic Function Values from History File
# -------------------------------------------------------------------


def read_aerodynamics(
    History_filename, nZones=1, special_cases=[], final_avg=0, wnd_fct="SQUARE"
):
    """values = read_aerodynamics(historyname, special_cases=[])
    read aerodynamic function values from history file

    Outputs:
        dictionary with function keys and thier values
        if special cases has 'TIME_MARCHING', returns time averaged data
        otherwise returns final value from history file
    """

    # read the history data
    history_data = read_history(History_filename, nZones)

    # pull only these functions
    Func_Values = ordered_bunch()
    for this_objfun in historyOutFields:
        if nZones == 1:
            if this_objfun in history_data:
                if (
                    historyOutFields[this_objfun]["TYPE"] == "COEFFICIENT"
                    or historyOutFields[this_objfun]["TYPE"] == "D_COEFFICIENT"
                ):
                    Func_Values[this_objfun] = history_data[this_objfun]
        else:
            for iZone in range(nZones):
                if this_objfun + "[" + str(iZone) + "]" in history_data:
                    if (
                        historyOutFields[this_objfun]["TYPE"] == "COEFFICIENT"
                        or historyOutFields[this_objfun]["TYPE"] == "D_COEFFICIENT"
                    ):
                        Func_Values[
                            this_objfun + "[" + str(iZone) + "]"
                        ] = history_data[this_objfun + "[" + str(iZone) + "]"]

    if "TIME_MARCHING" in special_cases:
        # for unsteady cases, average time-accurate objective function values
        for key, value in Func_Values.items():
            if historyOutFields[key]["TYPE"] == "COEFFICIENT":
                if not history_data.get("TAVG_" + key):
                    raise KeyError(
                        "Key "
                        + historyOutFields["TAVG_" + key]["HEADER"]
                        + " was not found in history output."
                    )
                Func_Values[key] = history_data["TAVG_" + key][-1]
            elif historyOutFields[key]["TYPE"] == "D_COEFFICIENT":
                if not history_data.get("TAVG_" + key):
                    raise KeyError(
                        "Key "
                        + historyOutFields["TAVG_" + key]["HEADER"]
                        + " was not found in history output."
                    )
                Func_Values[key] = history_data["TAVG_" + key][-1]
    else:
        # in steady cases take only last value.
        for key, value in Func_Values.iteritems():
            if not history_data.get(key):
                raise KeyError(
                    "Key "
                    + historyOutFields[key]["HEADER"]
                    + " was not found in history output."
                )
            Func_Values[key] = value[-1]

    return Func_Values


#: def read_aerodynamics()

# -------------------------------------------------------------------
#  Get Objective Function Sign
# -------------------------------------------------------------------


def get_objectiveSign(ObjFun_name):
    """returns -1 for maximization problems:
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
    if ObjFun_name == "LIFT":
        return -1.0
    if ObjFun_name == "EFFICIENCY":
        return -1.0
    if ObjFun_name == "THRUST":
        return -1.0
    if ObjFun_name == "FIGURE_OF_MERIT":
        return -1.0
    if ObjFun_name == "SURFACE_TOTAL_PRESSURE":
        return -1.0
    if ObjFun_name == "SURFACE_STATIC_PRESSURE":
        return -1.0
    if ObjFun_name == "SURFACE_MASSFLOW":
        return -1.0
    if ObjFun_name == "SURFACE_MACH":
        return -1.0
    if ObjFun_name == "TOTAL_STATIC_EFFICIENCY":
        return -1.0

    # otherwise
    return 1.0


#: def get_objectiveSign()


# -------------------------------------------------------------------
#  Get Constraint Sign
# -------------------------------------------------------------------


def get_constraintSign(sign):
    """gets +/-1 given a constraint sign < or > respectively
    inequality constraint is posed as c(x) < 0
    """
    sign_map = {">": -1.0, "<": +1.0}
    assert not sign == "=", 'Sign "=" not valid'

    return sign_map[sign]


#: def get_constraintSign()


# -------------------------------------------------------------------
#  Get Adjoint Filename Suffix
# -------------------------------------------------------------------


def get_adjointSuffix(objective_function=None):
    """gets the adjoint suffix given an objective function"""

    # adjoint name map
    name_map = {
        "DRAG": "cd",
        "LIFT": "cl",
        "SIDEFORCE": "csf",
        "MOMENT_X": "cmx",
        "MOMENT_Y": "cmy",
        "MOMENT_Z": "cmz",
        "FORCE_X": "cfx",
        "FORCE_Y": "cfy",
        "FORCE_Z": "cfz",
        "EFFICIENCY": "eff",
        "INVERSE_DESIGN_PRESSURE": "invpress",
        "INVERSE_DESIGN_HEAT": "invheat",
        "MAXIMUM_HEATFLUX": "maxheat",
        "TOTAL_HEATFLUX": "totheat",
        "EQUIVALENT_AREA": "ea",
        "NEARFIELD_PRESSURE": "nfp",
        "THRUST": "ct",
        "TORQUE": "cq",
        "FIGURE_OF_MERIT": "merit",
        "BUFFET": "buffet",
        "SURFACE_TOTAL_PRESSURE": "pt",
        "SURFACE_STATIC_PRESSURE": "pe",
        "SURFACE_MASSFLOW": "mfr",
        "SURFACE_MACH": "mach",
        "SURFACE_UNIFORMITY": "uniform",
        "SURFACE_SECONDARY": "second",
        "SURFACE_MOM_DISTORTION": "distort",
        "SURFACE_SECOND_OVER_UNIFORM": "sou",
        "SURFACE_PRESSURE_DROP": "dp",
        "CUSTOM_OBJFUNC": "custom",
        "KINETIC_ENERGY_LOSS": "ke",
        "TOTAL_PRESSURE_LOSS": "pl",
        "ENTROPY_GENERATION": "entg",
        "EULERIAN_WORK": "ew",
        "FLOW_ANGLE_OUT": "fao",
        "FLOW_ANGLE_IN": "fai",
        "MASS_FLOW_OUT": "mfo",
        "MASS_FLOW_IN": "mfi",
        "TOTAL_EFFICIENCY": "teff",
        "TOTAL_STATIC_EFFICIENCY": "tseff",
        "COMBO": "combo",
    }

    # if none or false, return map
    if not objective_function:
        return name_map
    else:
        # remove white space
        objective = "".join(objective_function.split())
        objective = objective.split(",")
        nObj = len(objective)
        if nObj > 1:
            return "combo"
        if objective[0] in name_map:
            return name_map[objective[0]]

        # otherwise...
        else:
            raise Exception("Unrecognized adjoint function name")


#: def get_adjointSuffix()

# -------------------------------------------------------------------
#  Add a Suffix
# -------------------------------------------------------------------


def add_suffix(base_name, suffix):
    """suffix_name = add_suffix(base_name,suffix)
    adds suffix to a filename, accounting for file type extension
    example:
        base_name   = 'input.txt'
        suffix      = 'new'
        suffix_name = 'input_new.txt'
    """
    if isinstance(base_name, list):
        suffix_name = []
        for name in base_name:
            name_split = os.path.splitext(name)
            suffix_name.append(name_split[0] + "_" + suffix + name_split[1])
    else:
        base_name = os.path.splitext(base_name)
        suffix_name = base_name[0] + "_" + suffix + base_name[1]

    return suffix_name


#: def add_suffix()


# -------------------------------------------------------------------
#  Get Design Variable ID Map
# -------------------------------------------------------------------


def get_dvMap():
    """get dictionary that maps design variable
    kind id number to name"""
    dv_map = {
        0: "NO_DEFORMATION",
        1: "TRANSLATION",
        2: "ROTATION",
        3: "SCALE",
        10: "FFD_SETTING",
        11: "FFD_CONTROL_POINT",
        12: "FFD_NACELLE",
        13: "FFD_GULL",
        14: "FFD_CAMBER",
        15: "FFD_TWIST",
        16: "FFD_THICKNESS",
        18: "FFD_ROTATION",
        19: "FFD_CONTROL_POINT_2D",
        20: "FFD_CAMBER_2D",
        21: "FFD_THICKNESS_2D",
        23: "FFD_CONTROL_SURFACE",
        24: "FFD_ANGLE_OF_ATTACK",
        30: "HICKS_HENNE",
        31: "PARABOLIC",
        32: "NACA_4DIGITS",
        33: "AIRFOIL",
        34: "CST",
        35: "SURFACE_BUMP",
        36: "SURFACE_FILE",
        40: "DV_EFIELD",
        41: "DV_YOUNG",
        42: "DV_POISSON",
        43: "DV_RHO",
        44: "DV_RHO_DL",
        50: "TRANSLATE_GRID",
        51: "ROTATE_GRID",
        52: "SCALE_GRID",
        101: "ANGLE_OF_ATTACK",
    }

    return dv_map


#: def get_dvMap()

# -------------------------------------------------------------------
#  Get Design Variable Kind Name from ID
# -------------------------------------------------------------------
def get_dvKind(kindID):
    """get design variable kind name from id number"""
    dv_map = get_dvMap()
    try:
        return dv_map[kindID]
    except KeyError:
        raise Exception("Unrecognized Design Variable ID")


# def get_dvKind()

# -------------------------------------------------------------------
#  Get Design Variable Kind ID from Name
# -------------------------------------------------------------------
def get_dvID(kindName):
    """get design variable kind id number from name"""
    dv_map = get_dvMap()
    id_map = dict((v, k) for (k, v) in dv_map.items())
    try:
        return id_map[kindName]
    except KeyError:
        raise Exception("Unrecognized Design Variable Name: %s", kindName)


#: def get_dvID()


# -------------------------------------------------------------------
#  Get Gradient File Header
# -------------------------------------------------------------------


def get_gradFileFormat(grad_type, plot_format, kindID, special_cases=[]):

    # start header, build a list of strings and join at the end
    header = []
    write_format = []

    # handle plot formating
    if plot_format == "TECPLOT":
        header.append("VARIABLES=")
    elif plot_format == "CSV":
        pass
    else:
        raise Exception("output plot format not recognized")

    # Case: continuous adjoint
    if grad_type == "CONTINUOUS_ADJOINT":
        header.append(r'"iVar","Gradient","FinDiff_Step"')
        write_format.append(r"%4d, %.10f, %f")

    # Case: finite difference
    elif grad_type == "FINITE_DIFFERENCE":
        header.append(
            r'"iVar","Grad_CL","Grad_CD","Grad_CSF","Grad_CMx","Grad_CMy","Grad_CMz","Grad_CFx","Grad_CFy","Grad_CFz","Grad_CL/CD","Grad_Custom_ObjFunc","Grad_HeatFlux_Total","Grad_HeatFlux_Maximum","Grad_Temperature_Total"'
        )
        write_format.append(
            r"%4d, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f"
        )

        for key in special_cases:
            if key == "ROTATING_FRAME":
                header.append(r',"Grad_CMerit","Grad_CT","Grad_CQ"')
                write_format.append(", %.10f, %.10f, %.10f")
            if key == "EQUIV_AREA":
                header.append(r',"Grad_CEquivArea","Grad_CNearFieldOF"')
                write_format.append(", %.10f, %.10f")
            if key == "ENGINE":
                header.append(
                    r',"Grad_AeroCDrag","Grad_SolidCDrag","Grad_Radial_Distortion","Grad_Circumferential_Distortion"'
                )
                write_format.append(", %.10f, %.10f, %.10f, %.10f")
            if key == "1D_OUTPUT":
                header.append(
                    r',"Grad_Avg_TotalPress","Grad_Avg_Mach","Grad_Avg_Temperature","Grad_MassFlowRate","Grad_Avg_Pressure","Grad_Avg_Density","Grad_Avg_Velocity","Grad_Avg_Enthalpy"'
                )
                write_format.append(
                    ", %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f"
                )
            if key == "INV_DESIGN_CP":
                header.append(r',"Grad_Cp_Diff"')
                write_format.append(", %.10f")
            if key == "INV_DESIGN_HEATFLUX":
                header.append(r',"Grad_HeatFlux_Diff"')
                write_format.append(", %.10f")

    # otherwise...
    else:
        raise Exception("Unrecognized Gradient Type")

    # design variable parameters
    if kindID == "FFD_CONTROL_POINT_2D":
        header.append(r',"FFD_Box_ID","xIndex","yIndex","xAxis","yAxis"')
        write_format.append(r", %s, %s, %s, %s, %s")
    elif kindID == "FFD_CAMBER_2D":
        header.append(r',"FFD_Box_ID","xIndex"')
        write_format.append(r", %s, %s")
    elif kindID == "FFD_THICKNESS_2D":
        header.append(r',"FFD_Box_ID","xIndex"')
        write_format.append(r", %s, %s")
    elif kindID == "HICKS_HENNE":
        header.append(r',"Up/Down","Loc_Max"')
        write_format.append(r", %s, %s")
    elif kindID == "SURFACE_BUMP":
        header.append(r',"Loc_Start","Loc_End","Loc_Max"')
        write_format.append(r", %s, %s, %s")
    elif kindID == "CST":
        header.append(r',"Up/Down","Kulfan number", "Total Kulfan numbers"')
        write_format.append(r", %s, %s", "%s")
    elif kindID == "FAIRING":
        header.append(r',"ControlPoint_Index","Theta_Disp","R_Disp"')
        write_format.append(r", %s, %s, %s")
    elif kindID == "NACA_4DIGITS":
        header.append(r',"1st_digit","2nd_digit","3rd&4th_digits"')
        write_format.append(r", %s, %s, %s")
    elif kindID == "TRANSLATION":
        header.append(r',"x_Disp","y_Disp","z_Disp"')
        write_format.append(r", %s, %s, %s")
    elif kindID == "ROTATION":
        header.append(r',"x_Orig","y_Orig","z_Orig","x_End","y_End","z_End"')
        write_format.append(r", %s, %s, %s, %s, %s, %s")
    elif kindID == "FFD_CONTROL_POINT":
        header.append(
            r',"FFD_Box_ID","xIndex","yIndex","zIndex","xAxis","yAxis","zAxis"'
        )
        write_format.append(r", %s, %s, %s, %s, %s, %s, %s")
    elif kindID == "FFD_DIHEDRAL_ANGLE":
        header.append(
            r',"FFD_Box_ID","x_Orig","y_Orig","z_Orig","x_End","y_End","z_End"'
        )
        write_format.append(r", %s, %s, %s, %s, %s, %s, %s")
    elif kindID == "FFD_TWIST_ANGLE":
        header.append(
            r',"FFD_Box_ID","x_Orig","y_Orig","z_Orig","x_End","y_End","z_End"'
        )
        write_format.append(r", %s, %s, %s, %s, %s, %s, %s")
    elif kindID == "FFD_ROTATION":
        header.append(
            r',"FFD_Box_ID","x_Orig","y_Orig","z_Orig","x_End","y_End","z_End"'
        )
        write_format.append(r", %s, %s, %s, %s, %s, %s, %s")
    elif kindID == "FFD_CAMBER":
        header.append(r',"FFD_Box_ID","xIndex","yIndex"')
        write_format.append(r", %s, %s, %s")
    elif kindID == "FFD_THICKNESS":
        header.append(r',"FFD_Box_ID","xIndex","yIndex"')
        write_format.append(r", %s, %s, %s")
    elif kindID == "ANGLE_OF_ATTACK":
        pass
    elif kindID == "FFD_ANGLE_OF_ATTACK":
        pass

    # otherwise...
    else:
        raise Exception("Unrecognized Design Variable Kind")

    # finite difference step
    if grad_type == "FINITE_DIFFERENCE":
        header.append(r',"FinDiff_Step"')
        write_format.append(r", %.10f")

    # finish format
    header.append("\n")
    write_format.append("\n")

    header = "".join(header)
    write_format = "".join(write_format)

    return [header, write_format]


#: def get_gradFileFormat()


# -------------------------------------------------------------------
#  Get Optimization File Header
# -------------------------------------------------------------------


def get_optFileFormat(plot_format, special_cases=None, nZones=1):

    if special_cases is None:
        special_cases = []

    # start header, build a list of strings and join at the end
    header_list = []
    header_format = ""
    write_format = []

    # handle plot formating
    if plot_format == "TECPLOT":
        header_format = header_format + "VARIABLES="
    elif plot_format == "CSV":
        pass
    else:
        raise Exception("output plot format not recognized")

    # start header
    header_list.extend(
        [
            "Iteration",
            "CL",
            "CD",
            "CSF",
            "CMx",
            "CMy",
            "CMz",
            "CFx",
            "CFy",
            "CFz",
            "CL/CD",
            "Custom_ObjFunc",
            "HeatFlux_Total",
            "HeatFlux_Maximum",
            "Temperature_Total",
        ]
    )
    write_format.append(
        r"%4d, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f"
    )

    # special cases
    for key in special_cases:
        if key == "ROTATING_FRAME":
            header_list.extend(["CMerit", "CT", "CQ"])
            write_format.append(r", %.10f, %.10f, %.10f")
        if key == "EQUIV_AREA":
            header_list.extend(["CEquivArea", "CNearFieldOF"])
            write_format.append(r", %.10f, %.10f")
        if key == "ENGINE":
            header_list.extend(
                [
                    "AeroCDrag",
                    "SolidCDrag",
                    "Radial_Distortion",
                    "Circumferential_Distortion",
                ]
            )
            write_format.append(r", %.10f, %.10f, %.10f, %.10f")
        if key == "1D_OUTPUT":
            header_list.extend(
                [
                    "AreaAvg_TotalPress",
                    "AreaAvg_Mach",
                    "AreaAvg_Temperature",
                    "MassFlowRate",
                    "Avg_Pressure",
                    "Avg_Density",
                    "Avg_Velocity",
                    "Avg_Enthalpy",
                ]
            )
            write_format.append(
                r", %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f"
            )
        if key == "INV_DESIGN_CP":
            header_list.extend(["Cp_Diff"])
            write_format.append(r", %.10f")
        if key == "INV_DESIGN_HEATFLUX":
            header_list.extend(["HeatFlux_Diff"])
            write_format.append(r", %.10f")

    # finish formats
    header_format = (
        (header_format) + ('"') + ('","').join(header_list) + ('"') + (" \n")
    )
    write_format = "".join(write_format) + " \n"

    # build list of objective function names
    header_vars = []
    map_dict = get_headerMap(nZones)
    for variable in header_list:
        assert variable in map_dict, "unrecognized header variable"
        header_vars.append(map_dict[variable])

    # done
    return [header_format, header_vars, write_format]


#: def get_optFileFormat()


# -------------------------------------------------------------------
#  Get Extension Name
# -------------------------------------------------------------------


def get_extension(output_format):

    if output_format == "PARAVIEW":
        return ".csv"
    if output_format == "PARAVIEW_BINARY":
        return ".csv"
    if output_format == "TECPLOT":
        return ".dat"
    if output_format == "TECPLOT_BINARY":
        return ".szplt"
    if output_format == "SOLUTION":
        return ".dat"
    if output_format == "RESTART":
        return ".dat"
    if output_format == "CONFIG":
        return ".cfg"
    if output_format == "CSV":
        return ".csv"
    # otherwise
    raise Exception("Output Format Unknown")


#: def get_extension()


# -------------------------------------------------------------------
#  Check Special Case
# -------------------------------------------------------------------
def get_specialCases(config):
    """returns a list of special physical problems that were
    specified in the config file, and set to 'yes'
    """

    all_special_cases = [
        "ROTATING_FRAME",
        "EQUIV_AREA",
        "1D_OUTPUT",
        "INV_DESIGN_CP",
        "INV_DESIGN_HEATFLUX",
    ]

    special_cases = []
    for key in all_special_cases:
        if key in config and config[key] == "YES":
            special_cases.append(key)
        if "SOLVER" in config and config["SOLVER"] == key:
            special_cases.append(key)

    if config.get("TIME_MARCHING", "NO") != "NO":
        special_cases.append("TIME_MARCHING")

    # no support for more than one special case
    if len(special_cases) > 1:
        error_str = (
            "Currently cannot support " + " and ".join(special_cases) + " at once"
        )
        raise Exception(error_str)

    # Special case for harmonic balance
    if "TIME_MARCHING" in config and config["TIME_MARCHING"] == "HARMONIC_BALANCE":
        special_cases.append("HARMONIC_BALANCE")

    # Special case for rotating frame
    if (
        "GRID_MOVEMENT_KIND" in config
        and config["GRID_MOVEMENT_KIND"] == "ROTATING_FRAME"
    ):
        special_cases.append("ROTATING_FRAME")

    return special_cases


#: def get_specialCases()

# -------------------------------------------------------------------
#  Check Fluid Structure Interaction
# -------------------------------------------------------------------
def get_multizone(config):
    """returns a list of special physical problems that were
    specified in the config file, and set to 'yes'
    """

    all_multizone_problems = ["FLUID_STRUCTURE_INTERACTION"]

    multizone = []
    for key in all_multizone_problems:
        if "SOLVER" in config and config["SOLVER"] == key:
            multizone.append(key)

    return multizone


#: def get_multizone()


def next_folder(folder_format, num_format="%03d"):
    """folder = next_folder(folder_format,num_format='%03d')
    finds the next folder with given format

    Inputs:
        folder_format - folder name with wild card (*) to mark expansion
        num_format    - %d formating to expand the wild card with

    Outputs:
        folder - a folder with the next index number inserted in
        the wild card, first index is 1
    """

    assert "*" in folder_format, "wildcard (*) missing in folder_format name"

    folders = glob.glob(folder_format)
    split = folder_format.split("*")
    folder = folder_format.replace("*", num_format)

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
        folder = folder % (max_i + 1)
    else:
        # first folder, number 1
        folder = folder % 1

    return folder


def expand_part(name, config):
    names = [name]
    return names


def expand_time(name, config):
    if "TIME_MARCHING" in get_specialCases(config):
        n_time = config["UNST_ADJOINT_ITER"]
        n_start_time = 0
        if (
            config.get("TIME_DOMAIN", "NO") == "YES"
            and config.get("RESTART_SOL", "NO") == "YES"
        ):
            n_start_time = int(config["RESTART_ITER"])
        if not isinstance(name, list):
            name_pat = add_suffix(name, "%05d")
            names = [name_pat % i for i in range(n_start_time, n_time)]
        else:
            for n in range(len(name)):
                name_pat = add_suffix(name[n], "%05d")
                names = [name_pat % i for i in range(n_start_time, n_time)]
    else:
        if not isinstance(name, list):
            names = [name]
        else:
            names = name
    return names


def expand_zones(name, config):
    names = []
    if int(config.NZONES) > 1:
        if not isinstance(name, list):
            name_pat = add_suffix(name, "%d")
            names = [name_pat % i for i in range(int(config.NZONES))]
        else:
            for n in range(len(name)):
                name_pat = add_suffix(name[n], "%d")
                names.extend([name_pat % i for i in range(int(config.NZONES))])

    else:
        if not isinstance(name, list):
            names = [name]
        else:
            names = name
    return names


def expand_multipoint(name, config):
    def_objs = config["OPT_OBJECTIVE"]
    objectives = def_objs.keys()
    names = []
    n_multipoint = len(config["MULTIPOINT_WEIGHT"].split(","))

    if any(elem in optnames_multi for elem in objectives):
        if not isinstance(name, list):
            if "_point0" not in name:
                name_pat = add_suffix(name, "point%d")
                names = [name_pat % i for i in range(n_multipoint)]
            else:
                name_parts = name.split("_point0")
                name_base = name_parts[0]
                name_suff = name_parts[1]
                name_pat = name_base + "_point%d" + name_suff
                names = [name_pat % i for i in range(n_multipoint)]
        else:
            for n in range(len(name)):
                if "_point0" not in name:
                    name_pat = add_suffix(name[n], "point%d")
                    names.extend([name_pat % i for i in range(n_multipoint)])
                else:
                    name_parts = name[n].split("_point0")
                    name_base = name_parts[0]
                    name_suff = name_parts[1]
                    name_pat = name_base + "_point%d" + name_suff
                    names.extend([name_pat % i for i in range(n_multipoint)])
    else:
        if not isinstance(name, list):
            names = [name]
        else:
            names = name
    return names


def make_link(src, dst):
    """make_link(src,dst)
    makes a relative link
    Inputs:
        src - source file
        dst - destination to place link

    Windows links currently unsupported, will copy file instead
    """

    if os.path.exists(src):  # , 'source file does not exist \n%s' % src

        if os.name == "nt":
            # can't make a link in windows, need to look for other options
            if os.path.exists(dst):
                os.remove(dst)
            shutil.copy(src, dst)

        else:
            # find real file, incase source itself is a link
            src = os.path.realpath(src)

            # normalize paths
            src = os.path.normpath(src)
            dst = os.path.normpath(dst)

            # check for self referencing
            if src == dst:
                return

            # find relative folder path
            srcfolder = os.path.join(os.path.split(src)[0]) + "/"
            dstfolder = os.path.join(os.path.split(dst)[0]) + "/"
            srcfolder = os.path.relpath(srcfolder, dstfolder)
            src = os.path.join(srcfolder, os.path.split(src)[1])

            # make unix link
            if os.path.exists(dst):
                os.remove(dst)
            os.symlink(src, dst)


def restart2solution(config, state={}):
    """restart2solution(config,state={})
    moves restart file to solution file,
    optionally updates state
    direct or adjoint is read from config
    adjoint objective is read from config
    """

    # direct solution
    if config.MATH_PROBLEM == "DIRECT":
        restart = config.RESTART_FILENAME
        solution = config.SOLUTION_FILENAME
        restart = restart.split(".")[0]
        solution = solution.split(".")[0]

        if "RESTART_ASCII" in config.get("OUTPUT_FILES", ["RESTART_BINARY"]):
            restart += ".csv"
            solution += ".csv"
        else:
            restart += ".dat"
            solution += ".dat"

        # expand zones
        restarts = expand_zones(restart, config)
        solutions = expand_zones(solution, config)
        # expand unsteady time
        restarts = expand_time(restarts, config)
        solutions = expand_time(solutions, config)

        # move
        for res, sol in zip(restarts, solutions):
            if os.path.exists(res):
                shutil.move(res, sol)
        # update state
        if state:
            state.FILES.DIRECT = solution
            if os.path.exists("flow.meta"):
                state.FILES.FLOW_META = "flow.meta"

    # adjoint solution
    elif any(
        [
            config.MATH_PROBLEM == "CONTINUOUS_ADJOINT",
            config.MATH_PROBLEM == "DISCRETE_ADJOINT",
        ]
    ):
        restart = config.RESTART_ADJ_FILENAME
        solution = config.SOLUTION_ADJ_FILENAME
        restart = restart.split(".")[0]
        solution = solution.split(".")[0]

        if "RESTART_ASCII" in config.get("OUTPUT_FILES", ["RESTART_BINARY"]):
            restart += ".csv"
            solution += ".csv"
        else:
            restart += ".dat"
            solution += ".dat"
        # add suffix
        func_name = config.OBJECTIVE_FUNCTION
        suffix = get_adjointSuffix(func_name)
        restart = add_suffix(restart, suffix)
        solution = add_suffix(solution, suffix)
        # expand zones
        restarts = expand_zones(restart, config)
        solutions = expand_zones(solution, config)
        # expand unsteady time
        restarts = expand_time(restarts, config)
        solutions = expand_time(solutions, config)

        # move
        for res, sol in zip(restarts, solutions):
            shutil.move(res, sol)
        # udpate state
        if "," in func_name:
            func_name = "COMBO"
        ADJ_NAME = "ADJOINT_" + func_name
        if state:
            state.FILES[ADJ_NAME] = solution

    else:
        raise Exception("unknown math problem")
