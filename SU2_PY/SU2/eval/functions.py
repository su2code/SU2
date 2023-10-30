#!/usr/bin/env python

## \file functions.py
#  \brief python package for functions
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

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import os, sys, shutil, copy, time, subprocess
from .. import run as su2run
from .. import io as su2io
from .. import util as su2util
from ..io import redirect_folder, redirect_output


# ----------------------------------------------------------------------
#  Main Function Interface
# ----------------------------------------------------------------------


def function(func_name, config, state=None):
    """val = SU2.eval.func(func_name,config,state=None)

    Evaluates the aerodynamics and geometry functions.

    Wraps:
        SU2.eval.aerodynamics()
        SU2.eval.geometry()

    Assumptions:
        Config is already setup for deformation.
        Mesh need not be deformed.
        Updates config and state by reference.
        Redundancy if state.FUNCTIONS is not empty.

    Executes in:
        ./DIRECT or ./GEOMETRY

    Inputs:
        func_name - SU2 objective function name or 'ALL'
        config    - an SU2 config
        state     - optional, an SU2 state

    Outputs:
        If func_name is 'ALL', returns a Bunch() of
        functions with keys of objective function names
        and values of objective function floats.
        Otherwise returns a float.
    """

    # initialize
    state = su2io.State(state)

    # check for multiple objectives
    multi_objective = type(func_name) == list

    # func_name_string is only used to check whether the function has already been evaluated.
    func_name_string = func_name
    if multi_objective:
        func_name_string = func_name[0]

    # redundancy check
    if not func_name_string in state["FUNCTIONS"]:

        # Aerodynamics
        if multi_objective or func_name == "ALL":
            aerodynamics(config, state)

        elif func_name in su2io.historyOutFields:
            if (
                su2io.historyOutFields[func_name]["TYPE"] == "COEFFICIENT"
                or su2io.historyOutFields[func_name]["TYPE"] == "D_COEFFICIENT"
            ):
                aerodynamics(config, state)

        # Stability
        elif func_name in su2io.optnames_stab:
            stability(config, state)

        # Multipoint
        elif func_name in su2io.optnames_multi:
            multipoint(config, state)

        # Geometry
        elif func_name in su2io.optnames_geo:
            geometry(func_name, config, state)

        else:
            raise Exception(
                "unknown function name, %s. Please check config_template.cfg for updated list of function names"
                % func_name
            )

    #: if not redundant

    # prepare output
    if func_name == "ALL":
        func_out = state["FUNCTIONS"]
    elif multi_objective:
        # If combine_objective is true, use the 'combo' output.
        func_out = state["FUNCTIONS"]["COMBO"]
    else:
        func_out = state["FUNCTIONS"][func_name]

    if func_name_string in config["OPT_OBJECTIVE"]:
        marker = config["OPT_OBJECTIVE"][func_name_string]["MARKER"]
        if func_name_string in su2io.per_surface_map:
            name = su2io.per_surface_map[func_name_string] + "_" + marker
            if name in state["FUNCTIONS"]:
                func_out = state["FUNCTIONS"][name]

    return copy.deepcopy(func_out)


#: def function()


# ----------------------------------------------------------------------
#  Aerodynamic Functions
# ----------------------------------------------------------------------


def aerodynamics(config, state=None):
    """vals = SU2.eval.aerodynamics(config,state=None)

    Evaluates aerodynamics with the following:
              SU2.run.deform()
        SU2.run.direct()

    Assumptions:
        Config is already setup for deformation.
        Mesh may or may not be deformed.
        Updates config and state by reference.
        Redundancy if state.FUNCTIONS is not empty.

    Executes in:
        ./DIRECT

    Inputs:
        config    - an SU2 config
        state     - optional, an SU2 state

    Outputs:
        Bunch() of functions with keys of objective function names
        and values of objective function floats.
    """

    # ----------------------------------------------------
    #  Initialize
    # ----------------------------------------------------

    # initialize
    state = su2io.State(state)

    # Make sure to output aerodynamic coeff.
    if not "AERO_COEFF" in config["HISTORY_OUTPUT"]:
        config["HISTORY_OUTPUT"].append("AERO_COEFF")

    if not "MESH" in state.FILES:
        state.FILES.MESH = config["MESH_FILENAME"]
    special_cases = su2io.get_specialCases(config)

    # console output
    if config.get("CONSOLE", "VERBOSE") in ["QUIET", "CONCISE"]:
        log_direct = "log_Direct.out"
    else:
        log_direct = None

    # ----------------------------------------------------
    #  Update Mesh
    # ----------------------------------------------------

    # does decomposition and deformation
    info = update_mesh(config, state)

    # ----------------------------------------------------
    #  Adaptation (not implemented)
    # ----------------------------------------------------

    # if not state.['ADAPTED_FUNC']:
    #    config = su2run.adaptation(config)
    #    state['ADAPTED_FUNC'] = True

    # ----------------------------------------------------
    #  Direct Solution
    # ----------------------------------------------------
    opt_names = []
    for key in su2io.historyOutFields:
        if su2io.historyOutFields[key]["TYPE"] == "COEFFICIENT":
            opt_names.append(key)

    # redundancy check
    direct_done = all([key in state.FUNCTIONS for key in opt_names])
    if direct_done:
        # return aerodynamic function values
        aero = su2util.ordered_bunch()
        for key in opt_names:
            if key in state.FUNCTIONS:
                aero[key] = state.FUNCTIONS[key]
        return copy.deepcopy(aero)
    #: if redundant

    # files to pull
    files = state.FILES
    pull = []
    link = []

    # files: mesh
    name = files["MESH"]
    name = su2io.expand_part(name, config)
    link.extend(name)

    pull.extend(config.get("CONFIG_LIST", []))

    # files: restarts
    if (
        config.get("TIME_DOMAIN", "NO") == "YES"
        and config.get("RESTART_SOL", "NO") == "YES"
    ):
        if "RESTART_FILE_1" in files:  # not the case for directdiff restart
            name = files["RESTART_FILE_1"]
            name = su2io.expand_part(name, config)
            link.extend(name)
        if "RESTART_FILE_2" in files:  # not the case for 1st order time stepping
            name = files["RESTART_FILE_2"]
            name = su2io.expand_part(name, config)
            link.extend(name)

    if "FLOW_META" in files:
        pull.append(files["FLOW_META"])

    # files: direct solution
    if "DIRECT" in files:
        name = files["DIRECT"]
        name = su2io.expand_zones(name, config)
        name = su2io.expand_time(name, config)
        link.extend(name)
        ##config['RESTART_SOL'] = 'YES' # don't override config file
    else:
        if (
            config.get("TIME_DOMAIN", "NO") != "YES"
        ):  # rules out steady state optimization special cases.
            config["RESTART_SOL"] = "NO"  # for shape optimization with restart files.

    # files: target equivarea distribution
    if "EQUIV_AREA" in special_cases and "TARGET_EA" in files:
        pull.append(files["TARGET_EA"])

    # files: target pressure distribution
    if "INV_DESIGN_CP" in special_cases and "TARGET_CP" in files:
        pull.append(files["TARGET_CP"])

    # files: target heat flux distribution
    if "INV_DESIGN_HEATFLUX" in special_cases and "TARGET_HEATFLUX" in files:
        pull.append(files["TARGET_HEATFLUX"])

    # output redirection
    with redirect_folder("DIRECT", pull, link) as push:
        with redirect_output(log_direct):

            # # RUN DIRECT SOLUTION # #
            info = su2run.direct(config)

            konfig = copy.deepcopy(config)
            """
            If the time convergence criterion was activated, we have less time iterations.
            Store the changed values of TIME_ITER, ITER_AVERAGE_OBJ and UNST_ADJOINT_ITER in
            info.WND_CAUCHY_DATA"""
            if (
                konfig.get("WINDOW_CAUCHY_CRIT", "NO") == "YES"
                and konfig.TIME_MARCHING != "NO"
            ):  # Tranfer Convergence Data, if necessary
                konfig["TIME_ITER"] = info.WND_CAUCHY_DATA["TIME_ITER"]
                konfig["ITER_AVERAGE_OBJ"] = info.WND_CAUCHY_DATA["ITER_AVERAGE_OBJ"]
                konfig["UNST_ADJOINT_ITER"] = info.WND_CAUCHY_DATA["UNST_ADJOINT_ITER"]

            su2io.restart2solution(konfig, info)
            state.update(info)

            # direct files to push
            name = info.FILES["DIRECT"]
            name = su2io.expand_zones(name, konfig)
            name = su2io.expand_time(name, konfig)
            push.extend(name)

            # pressure files to push
            if "TARGET_CP" in info.FILES:
                push.append(info.FILES["TARGET_CP"])

            # heat flux files to push
            if "TARGET_HEATFLUX" in info.FILES:
                push.append(info.FILES["TARGET_HEATFLUX"])

            if "FLOW_META" in info.FILES:
                push.append(info.FILES["FLOW_META"])

    #: with output redirection
    su2io.update_persurface(konfig, state)
    # return output
    funcs = su2util.ordered_bunch()
    for key in state["FUNCTIONS"]:
        funcs[key] = state["FUNCTIONS"][key]

    return funcs


#: def aerodynamics()


# ----------------------------------------------------------------------
#  Stability Functions
# ----------------------------------------------------------------------


def stability(config, state=None, step=1e-2):

    folder = "STABILITY"  # os.path.join('STABILITY',func_name) #STABILITY/D_MOMENT_Y_D_ALPHA/

    # ----------------------------------------------------
    #  Initialize
    # ----------------------------------------------------

    # initialize
    state = su2io.State(state)
    if not "MESH" in state.FILES:
        state.FILES.MESH = config["MESH_FILENAME"]
    special_cases = su2io.get_specialCases(config)

    # console output
    if config.get("CONSOLE", "VERBOSE") in ["QUIET", "CONCISE"]:
        log_direct = "log_Direct.out"
    else:
        log_direct = None

    # ----------------------------------------------------
    #  Update Mesh
    # ----------------------------------------------------

    # does decomposition and deformation
    info = update_mesh(config, state)

    # ----------------------------------------------------
    #  CENTRAL POINT
    # ----------------------------------------------------

    # will run in DIRECT/
    func_0 = aerodynamics(config, state)

    # ----------------------------------------------------
    #  Run Forward Point
    # ----------------------------------------------------

    # files to pull
    files = state.FILES
    pull = []
    link = []

    # files: mesh
    name = files["MESH"]
    name = su2io.expand_part(name, config)
    link.extend(name)

    # files: direct solution
    if "DIRECT" in files:
        name = files["DIRECT"]
        name = su2io.expand_time(name, config)
        link.extend(name)
        ##config['RESTART_SOL'] = 'YES' # don't override config file
    else:
        config["RESTART_SOL"] = "NO"

    # files: target equivarea distribution
    if "EQUIV_AREA" in special_cases and "TARGET_EA" in files:
        pull.append(files["TARGET_EA"])

    # files: target pressure distribution
    if "INV_DESIGN_CP" in special_cases and "TARGET_CP" in files:
        pull.append(files["TARGET_CP"])

    # files: target heat flux distribution
    if "INV_DESIGN_HEATFLUX" in special_cases and "TARGET_HEATFLUX" in files:
        pull.append(files["TARGET_HEATFLUX"])

    # pull needed files, start folder
    with redirect_folder(folder, pull, link) as push:
        with redirect_output(log_direct):

            konfig = copy.deepcopy(config)
            ztate = copy.deepcopy(state)

            # TODO: GENERALIZE
            konfig.AOA = konfig.AOA + step
            ztate.FUNCTIONS.clear()

            func_1 = aerodynamics(konfig, ztate)

            ## direct files to store
            # name = ztate.FILES['DIRECT']
            # if not 'STABILITY' in state.FILES:
            # state.FILES.STABILITY = su2io.ordered_bunch()
            # state.FILES.STABILITY['DIRECT'] = name

            ## equivarea files to store
            # if 'WEIGHT_NF' in ztate.FILES:
            # state.FILES.STABILITY['WEIGHT_NF'] = ztate.FILES['WEIGHT_NF']

    # ----------------------------------------------------
    #  DIFFERENCING
    # ----------------------------------------------------

    for derv_name in su2io.optnames_stab:

        matches = [k for k in su2io.optnames_aero if k in derv_name]
        if not len(matches) == 1:
            continue
        func_name = matches[0]

        obj_func = (func_1[func_name] - func_0[func_name]) / step

        state.FUNCTIONS[derv_name] = obj_func

    # return output
    funcs = su2util.ordered_bunch()
    for key in su2io.optnames_stab:
        if key in state["FUNCTIONS"]:
            funcs[key] = state["FUNCTIONS"][key]

    return funcs


# ----------------------------------------------------------------------
#  Multipoint Functions
# ----------------------------------------------------------------------


def multipoint(config, state=None, step=1e-2):

    mach_list = (
        config["MULTIPOINT_MACH_NUMBER"].replace("(", "").replace(")", "").split(",")
    )
    reynolds_list = (
        config["MULTIPOINT_REYNOLDS_NUMBER"]
        .replace("(", "")
        .replace(")", "")
        .split(",")
    )
    freestream_temp_list = (
        config["MULTIPOINT_FREESTREAM_TEMPERATURE"]
        .replace("(", "")
        .replace(")", "")
        .split(",")
    )
    freestream_press_list = (
        config["MULTIPOINT_FREESTREAM_PRESSURE"]
        .replace("(", "")
        .replace(")", "")
        .split(",")
    )
    aoa_list = config["MULTIPOINT_AOA"].replace("(", "").replace(")", "").split(",")
    sideslip_list = (
        config["MULTIPOINT_SIDESLIP_ANGLE"].replace("(", "").replace(")", "").split(",")
    )
    target_cl_list = (
        config["MULTIPOINT_TARGET_CL"].replace("(", "").replace(")", "").split(",")
    )
    weight_list = (
        config["MULTIPOINT_WEIGHT"].replace("(", "").replace(")", "").split(",")
    )
    outlet_value_list = (
        config["MULTIPOINT_OUTLET_VALUE"].replace("(", "").replace(")", "").split(",")
    )
    solution_flow_list = su2io.expand_multipoint(config.SOLUTION_FILENAME, config)
    flow_meta_list = su2io.expand_multipoint("flow.meta", config)
    restart_sol = config["RESTART_SOL"]
    dv_value_old = config["DV_VALUE_OLD"]

    func = []
    folder = []
    for i in range(len(weight_list)):
        func.append(0)
        folder.append(0)

    for i in range(len(weight_list)):
        folder[i] = "MULTIPOINT_" + str(i)

    opt_names = []
    for key in su2io.historyOutFields:
        if su2io.historyOutFields[key]["TYPE"] == "COEFFICIENT":
            opt_names.append(key)

    # ----------------------------------------------------
    #  Initialize
    # ----------------------------------------------------

    # initialize
    state = su2io.State(state)
    if not "MESH" in state.FILES:
        state.FILES.MESH = config["MESH_FILENAME"]
    special_cases = su2io.get_specialCases(config)

    # console output
    if config.get("CONSOLE", "VERBOSE") in ["QUIET", "CONCISE"]:
        log_direct = "log_Direct.out"
    else:
        log_direct = None

    # ----------------------------------------------------
    #  Update Mesh
    # ----------------------------------------------------

    # If multiple meshes specified, use relevant mesh
    if "MULTIPOINT_MESH_FILENAME" in state.FILES:
        state.FILES.MESH = state.FILES.MULTIPOINT_MESH_FILENAME[0]
        config.MESH_FILENAME = state.FILES.MULTIPOINT_MESH_FILENAME[0]

    # does decomposition and deformation
    info = update_mesh(config, state)

    # ----------------------------------------------------
    #  FIRST POINT
    # ----------------------------------------------------

    # will run in DIRECT/

    config.AOA = aoa_list[0]
    config.SIDESLIP_ANGLE = sideslip_list[0]
    config.MACH_NUMBER = mach_list[0]
    config.REYNOLDS_NUMBER = reynolds_list[0]
    config.FREESTREAM_TEMPERATURE = freestream_temp_list[0]
    config.FREESTREAM_PRESSURE = freestream_press_list[0]
    config.TARGET_CL = target_cl_list[0]
    orig_marker_outlet = config["MARKER_OUTLET"]
    orig_marker_outlet = orig_marker_outlet.replace("(", "").replace(")", "").split(",")
    new_marker_outlet = "(" + orig_marker_outlet[0] + "," + outlet_value_list[0] + ")"
    config.MARKER_OUTLET = new_marker_outlet
    config.SOLUTION_FILENAME = solution_flow_list[0]

    # If solution file for the first point is available, use it
    if "MULTIPOINT_DIRECT" in state.FILES and state.FILES.MULTIPOINT_DIRECT[0]:
        state.FILES["DIRECT"] = state.FILES.MULTIPOINT_DIRECT[0]

    # If flow.meta file for the first point is available, rename it before using it
    if "MULTIPOINT_FLOW_META" in state.FILES and state.FILES.MULTIPOINT_FLOW_META[0]:
        os.rename(state.FILES.MULTIPOINT_FLOW_META[0], "flow.meta")
        state.FILES["FLOW_META"] = "flow.meta"

    func[0] = aerodynamics(config, state)

    # change name of flow.meta back to multipoint name
    if os.path.exists("flow.meta"):
        os.rename("flow.meta", flow_meta_list[0])
        state.FILES["FLOW_META"] = flow_meta_list[0]

    src = os.getcwd()
    src = os.path.abspath(src).rstrip("/") + "/DIRECT/"

    # files to pull
    files = state.FILES
    pull = []
    link = []

    # files: mesh
    name = files["MESH"]
    name = su2io.expand_part(name, config)
    link.extend(name)

    # files: direct solution
    if "DIRECT" in files:
        name = files["DIRECT"]
        name = su2io.expand_time(name, config)
        link.extend(name)
    else:
        config["RESTART_SOL"] = "NO"

    # files: meta data for the flow
    if "FLOW_META" in files:
        pull.append(files["FLOW_META"])

    # files: target equivarea distribution
    if "EQUIV_AREA" in special_cases and "TARGET_EA" in files:
        pull.append(files["TARGET_EA"])

    # files: target pressure distribution
    if "INV_DESIGN_CP" in special_cases and "TARGET_CP" in files:
        pull.append(files["TARGET_CP"])

    # files: target heat flux distribution
    if "INV_DESIGN_HEATFLUX" in special_cases and "TARGET_HEATFLUX" in files:
        pull.append(files["TARGET_HEATFLUX"])

    # pull needed files, start folder_0
    with redirect_folder(folder[0], pull, link) as push:
        with redirect_output(log_direct):

            konfig = copy.deepcopy(config)
            ztate = copy.deepcopy(state)
            # Reset restart to original value
            konfig["RESTART_SOL"] = restart_sol

            dst = os.getcwd()
            dst = os.path.abspath(dst).rstrip("/") + "/" + "DIRECT"

            # make unix link
            string = "ln -s " + src + " " + dst
            stringlist = string.split()
            subprocess.Popen(stringlist)

    for i in range(len(weight_list) - 1):

        konfig = copy.deepcopy(config)
        ztate = copy.deepcopy(state)

        konfig.SOLUTION_FILENAME = solution_flow_list[i + 1]

        # delete direct solution file from previous point
        if "DIRECT" in ztate.FILES:
            del ztate.FILES.DIRECT

        if "FLOW_META" in ztate.FILES:
            del ztate.FILES.FLOW_META

        # use direct solution file from relevant point
        if "MULTIPOINT_DIRECT" in state.FILES and state.FILES.MULTIPOINT_DIRECT[i + 1]:
            ztate.FILES["DIRECT"] = state.FILES.MULTIPOINT_DIRECT[i + 1]

        # use flow.meta file from relevant point
        if (
            "MULTIPOINT_FLOW_META" in state.FILES
            and state.FILES.MULTIPOINT_FLOW_META[i + 1]
        ):
            ztate.FILES["FLOW_META"] = state.FILES.MULTIPOINT_FLOW_META[i + 1]

        # use mesh file from relevant point
        if "MULTIPOINT_MESH_FILENAME" in ztate.FILES:
            ztate.FILES.MESH = ztate.FILES.MULTIPOINT_MESH_FILENAME[i + 1]
            konfig.MESH_FILENAME = ztate.FILES.MULTIPOINT_MESH_FILENAME[i + 1]
            konfig["DV_VALUE_OLD"] = dv_value_old

        files = ztate.FILES
        link = []
        pull = []

        # files: mesh
        name = files["MESH"]
        name = su2io.expand_part(name, konfig)
        link.extend(name)

        # files: direction solution
        if "DIRECT" in files:
            name = files["DIRECT"]
            name = su2io.expand_time(name, konfig)
            link.extend(name)
        else:
            konfig["RESTART_SOL"] = "NO"

        # files: meta data for the flow
        if "FLOW_META" in files:
            pull.append(files["FLOW_META"])

        # pull needed files, start folder_1
        with redirect_folder(folder[i + 1], pull, link) as push:
            with redirect_output(log_direct):

                # Perform deformation on multipoint mesh
                if "MULTIPOINT_MESH_FILENAME" in state.FILES:
                    info = update_mesh(konfig, ztate)

                # Update config values
                konfig.AOA = aoa_list[i + 1]
                konfig.SIDESLIP_ANGLE = sideslip_list[i + 1]
                konfig.MACH_NUMBER = mach_list[i + 1]
                konfig.REYNOLDS_NUMBER = reynolds_list[i + 1]
                konfig.FREESTREAM_TEMPERATURE = freestream_temp_list[i + 1]
                konfig.FREESTREAM_PRESSURE = freestream_press_list[i + 1]
                konfig.TARGET_CL = target_cl_list[i + 1]
                orig_marker_outlet = config["MARKER_OUTLET"]
                orig_marker_outlet = (
                    orig_marker_outlet.replace("(", "").replace(")", "").split(",")
                )
                new_marker_outlet = (
                    "(" + orig_marker_outlet[0] + "," + outlet_value_list[i + 1] + ")"
                )
                konfig.MARKER_OUTLET = new_marker_outlet

                ztate.FUNCTIONS.clear()

                # rename meta data to flow.meta
                if "FLOW_META" in ztate.FILES:
                    ztate.FILES["FLOW_META"] = "flow.meta"
                    os.rename(ztate.FILES.MULTIPOINT_FLOW_META[i + 1], "flow.meta")

                func[i + 1] = aerodynamics(konfig, ztate)

                dst = os.getcwd()

                # revert name of flow.meta file to multipoint name
                if os.path.exists("flow.meta"):
                    os.rename("flow.meta", flow_meta_list[i + 1])
                    ztate.FILES["FLOW_META"] = flow_meta_list[i + 1]
                    dst_flow_meta = (
                        os.path.abspath(dst).rstrip("/")
                        + "/"
                        + ztate.FILES["FLOW_META"]
                    )
                    push.append(ztate.FILES["FLOW_META"])

                # direct files to push
                dst_direct = (
                    os.path.abspath(dst).rstrip("/") + "/" + ztate.FILES["DIRECT"]
                )
                name = ztate.FILES["DIRECT"]
                name = su2io.expand_zones(name, konfig)
                name = su2io.expand_time(name, konfig)
                push.extend(name)

                if "MULTIPOINT_MESH_FILENAME" in state.FILES:
                    # Mesh files to push
                    dst_mesh = (
                        os.path.abspath(dst).rstrip("/") + "/" + ztate.FILES["MESH"]
                    )
                    name = ztate.FILES["MESH"]
                    name = su2io.expand_part(name, konfig)
                    push.extend(name)

        # Link direct solution to MULTIPOINT_# folder
        src = os.getcwd()
        src_direct = os.path.abspath(src).rstrip("/") + "/" + ztate.FILES["DIRECT"]

        # make unix link
        os.symlink(src_direct, dst_direct)

        # If the mesh doesn't already exist, link it
        if "MULTIPOINT_MESH_FILENAME" in state.FILES:
            src_mesh = os.path.abspath(src).rstrip("/") + "/" + ztate.FILES["MESH"]
            if not os.path.exists(src_mesh):
                os.symlink(src_mesh, dst_mesh)

        # link flow.meta
        if "MULTIPOINT_FLOW_META" in state.FILES:
            src_flow_meta = (
                os.path.abspath(src).rstrip("/") + "/" + ztate.FILES["FLOW_META"]
            )
            if not os.path.exists(src_flow_meta):
                os.symlink(src_flow_meta, dst_flow_meta)

    # Update MULTIPOINT_DIRECT in state.FILES
    state.FILES.MULTIPOINT_DIRECT = solution_flow_list
    if "FLOW_META" in state.FILES:
        state.FILES.MULTIPOINT_FLOW_META = flow_meta_list

    # ----------------------------------------------------
    #  WEIGHT FUNCTIONS
    # ----------------------------------------------------

    for derv_name in su2io.optnames_multi:
        matches = [k for k in opt_names if k in derv_name]
        if not len(matches) == 1:
            continue
        func_name = matches[0]
        obj_func = 0.0
        for i in range(len(weight_list)):
            obj_func = obj_func + float(weight_list[i]) * func[i][func_name]

        state.FUNCTIONS[derv_name] = obj_func

    # return output
    funcs = su2util.ordered_bunch()
    for key in su2io.optnames_multi:
        if key in state["FUNCTIONS"]:
            funcs[key] = state["FUNCTIONS"][key]

    return funcs


# ----------------------------------------------------------------------
#  Geometric Functions
# ----------------------------------------------------------------------


def geometry(func_name, config, state=None):
    """val = SU2.eval.geometry(config,state=None)

    Evaluates geometry with the following:
        SU2.run.deform()
        SU2.run.geometry()

    Assumptions:
        Config is already setup for deformation.
        Mesh may or may not be deformed.
        Updates config and state by reference.
        Redundancy if state.FUNCTIONS does not have func_name.

    Executes in:
        ./GEOMETRY

    Inputs:
        config    - an SU2 config
        state     - optional, an SU2 state

    Outputs:
        Bunch() of functions with keys of objective function names
        and values of objective function floats.
    """

    # ----------------------------------------------------
    #  Initialize
    # ----------------------------------------------------

    # initialize
    state = su2io.State(state)
    if not "MESH" in state.FILES:
        state.FILES.MESH = config["MESH_FILENAME"]
    special_cases = su2io.get_specialCases(config)

    # console output
    if config.get("CONSOLE", "VERBOSE") in ["QUIET", "CONCISE"]:
        log_geom = "log_Geometry.out"
    else:
        log_geom = None

    # ----------------------------------------------------
    #  Update Mesh (check with Trent)
    # ----------------------------------------------------

    # does decomposition and deformation
    # info = update_mesh(config,state)

    # ----------------------------------------------------
    #  Geometry Solution
    # ----------------------------------------------------

    # redundancy check
    geometry_done = func_name in state.FUNCTIONS
    # geometry_done = all([key in state.FUNCTIONS for key in su2io.optnames_geo])
    if not geometry_done:

        # files to pull
        files = state.FILES
        pull = []
        link = []

        # files: mesh
        name = files["MESH"]
        name = su2io.expand_part(name, config)
        link.extend(name)

        # update function name
        ## TODO

        # output redirection
        with redirect_folder("GEOMETRY", pull, link) as push:
            with redirect_output(log_geom):

                # setup config
                config.GEO_PARAM = func_name
                config.GEO_MODE = "FUNCTION"

                # # RUN GEOMETRY SOLUTION # #
                info = su2run.geometry(config)
                state.update(info)

                # no files to push

        #: with output redirection

    #: if not redundant

    # return output
    funcs = su2util.ordered_bunch()
    for key in su2io.optnames_geo:
        if key in state["FUNCTIONS"]:
            funcs[key] = state["FUNCTIONS"][key]
    return funcs


#: def geometry()


def update_mesh(config, state=None):
    """SU2.eval.update_mesh(config,state=None)

    updates mesh with the following:
              SU2.run.deform()

    Assumptions:
        Config is already setup for deformation.
        Mesh may or may not be deformed.
        Updates config and state by reference.

    Executes in:
        ./DECOMP and ./DEFORM

    Inputs:
        config    - an SU2 config
        state     - optional, an SU2 state

    Outputs:
        nothing

    Modifies:
        config and state by reference
    """

    # ----------------------------------------------------
    #  Initialize
    # ----------------------------------------------------

    # initialize
    state = su2io.State(state)
    if not "MESH" in state.FILES:
        state.FILES.MESH = config["MESH_FILENAME"]
    special_cases = su2io.get_specialCases(config)

    # console output
    if config.get("CONSOLE", "VERBOSE") in ["QUIET", "CONCISE"]:
        log_decomp = "log_Decomp.out"
        log_deform = "log_Deform.out"
    else:
        log_decomp = None
        log_deform = None

    # ----------------------------------------------------
    #  Deformation
    # ----------------------------------------------------

    # redundancy check
    deform_set = config["DV_KIND"] == config["DEFINITION_DV"]["KIND"]
    deform_todo = not config["DV_VALUE_NEW"] == config["DV_VALUE_OLD"]
    if deform_set and deform_todo:

        # files to pull
        pull = []
        link = config["MESH_FILENAME"]
        link = su2io.expand_part(link, config)

        pull.extend(config.get("CONFIG_LIST", []))

        # output redirection
        with redirect_folder("DEFORM", pull, link) as push:
            with redirect_output(log_deform):

                # # RUN DEFORMATION # #
                info = su2run.deform(config)
                state.update(info)

                # data to push
                meshname = info.FILES.MESH
                names = su2io.expand_part(meshname, config)
                push.extend(names)

        #: with redirect output

    elif deform_set and not deform_todo:
        state.VARIABLES.DV_VALUE_NEW = config.DV_VALUE_NEW

    #: if not redundant

    return
