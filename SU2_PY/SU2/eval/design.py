#!/usr/bin/env python

## \file design.py
#  \brief python package for designs
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

import os, copy
from .. import io as su2io
from . import func as su2func
from . import grad as su2grad
from ..io import redirect_folder, save_data

# todo:
# shouldnt be needed, but self.append_state() (ie after initialization)


# ----------------------------------------------------------------------
#  Design Class
# ----------------------------------------------------------------------


class Design(object):
    """SU2.eval.Design(config,state=None,folder='DESIGNS/DSN_*')

    Starts a design class, which manages a config and state.
    Will run design in folder, and with self indexing name if '*' is
    included in the folder name.
    Methods are wrappers for SU2.eval.func() and SU2.eval.grad()

    Attributes:
        state  - design state
        config - design config
        files  - design files
        folder - design folder
        funcs  - design function value bunch
        grads  - design gradient values bunch

    Methods:
        Optimizer Interface
        The following methods take a design vector for input
        as a list (shape n) or numpy array (shape n or nx1 or 1xn).
        Values are returned as floats or lists or lists of lists.
        See SU2.eval.obj_f, etc for more detail.

        obj_f(dvs)     - objective function              : float
        obj_df(dvs)    - objective function derivatives  : list
        con_ceq(dvs)   - equality constraints            : list
        con_dceq(dvs)  - equality constraint derivatives : list[list]
        con_cieq(dvs)  - inequality constraints          : list
        con_dcieq(dvs) - inequality constraint gradients : list[list]

        Functional Interface
        The following methods take an objective function name for input.
        func(func_name)                  - function of specified name
        grad(func_name,method='CONTINUOUS_ADJOINT') - gradient of specified name
    """

    def __init__(self, config, state=None, folder="DESIGNS/DSN_*"):
        """Initializes an SU2 Design"""

        ## ???: Move to Project, no next folder here

        if "*" in folder:
            folder = su2io.next_folder(folder)

        config = copy.deepcopy(config)
        state = copy.deepcopy(state)
        state = su2io.State(state)
        state.find_files(config)

        self.config = config
        self.state = state
        self.files = state.FILES
        self.funcs = state.FUNCTIONS
        self.grads = state.GRADIENTS
        self.folder = folder

        self.filename = "design.pkl"

        # initialize folder with files
        pull, link = state.pullnlink(config)
        with redirect_folder(folder, pull, link, force=True):
            # save design, config
            save_data(self.filename, self)
            config.dump("config_DSN.cfg")

    def _eval(self, eval_func, *args):
        """Evaluates an SU2 Design
        always adds config and state to the inputs list
        """

        config = self.config
        state = self.state
        files = self.files
        folder = self.folder

        filename = self.filename

        # check folder
        assert os.path.exists(folder), "cannot find design folder %s" % folder

        konfig = copy.deepcopy(config)

        """
        If the time convergence criterion was activated, we have less time iterations.
        Store the changed values of TIME_ITER, ITER_AVERAGE_OBJ and UNST_ADJOINT_ITER in
        state.WND_CAUCHY_DATA"""
        if (
            "TIME_ITER" in state.WND_CAUCHY_DATA
        ):  # Use Convergence data, if we have already a direct run
            konfig["TIME_ITER"] = state.WND_CAUCHY_DATA["TIME_ITER"]
            konfig["ITER_AVERAGE_OBJ"] = state.WND_CAUCHY_DATA["ITER_AVERAGE_OBJ"]
            konfig["UNST_ADJOINT_ITER"] = state.WND_CAUCHY_DATA["UNST_ADJOINT_ITER"]

        # list files to pull and link
        pull, link = state.pullnlink(konfig)

        # output redirection, don't re-pull files
        with redirect_folder(folder, pull, link, force=False) as push:

            # get timestamp
            timestamp = state.tic()

            # run
            inputs = args + (config, state)
            vals = eval_func(*inputs)

            # save design
            if state.toc(timestamp):
                save_data(filename, self)

        #: with redirect folder

        # update files
        files.update(state["FILES"])

        return vals

    def obj_f(self, dvs):
        """Evaluates SU2 Design Objectives"""
        return self._eval(obj_f, dvs)

    def obj_df(self, dvs):
        """Evaluates SU2 Design Objective Gradients"""
        return self._eval(obj_df, dvs)

    def con_ceq(self, dvs):
        """Evaluates SU2 Design Equality Constraints"""
        return self._eval(con_ceq, dvs)

    def con_dceq(self, dvs):
        """Evaluates SU2 Design Equality Constraint Gradients"""
        return self._eval(con_dceq, dvs)

    def con_cieq(self, dvs):
        """Evaluates SU2 Design Inequality Constraints"""
        return self._eval(con_cieq, dvs)

    def con_dcieq(self, dvs):
        """Evaluates SU2 Design Inequality Constraint Gradients"""
        return self._eval(con_dcieq, dvs)

    def func(self, func_name):
        """Evaluates SU2 Design Functions by Name"""
        return self._eval(su2func, func_name)

    def grad(self, func_name, method="CONTINUOUS_ADJOINT"):
        """Evaluates SU2 Design Gradients by Name"""
        return self._eval(su2grad, func_name, method)

    def touch(self):
        return self._eval(touch)

    def skip(self, *args, **kwarg):
        return self._eval(skip)

    def __repr__(self):
        return "<Design> %s" % self.folder

    def __str__(self):
        output = self.__repr__()
        output += "\n%s" % self.state
        return output


#: class Design()


# ----------------------------------------------------------------------
#  Optimization Interface Functions
# ----------------------------------------------------------------------


def obj_f(dvs, config, state=None):
    """val = SU2.eval.obj_f(dvs,config,state=None)

    Evaluates SU2 Objectives
    Wraps SU2.eval.func()

    Takes a design vector for input as a list (shape n)
    or numpy array (shape n or nx1 or 1xn), a config
    and optionally a state.

    Outputs a float.
    """

    # unpack config and state
    config.unpack_dvs(dvs)
    state = su2io.State(state)

    def_objs = config["OPT_OBJECTIVE"]
    objectives = def_objs.keys()

    # evaluate each objective
    vals_out = []
    func = 0.0
    for i_obj, this_obj in enumerate(objectives):
        scale = def_objs[this_obj]["SCALE"]
        global_factor = float(config["OPT_GRADIENT_FACTOR"])
        sign = su2io.get_objectiveSign(this_obj)

        # Evaluate Objective Function scaling and sign
        # If default evaluate as normal,
        if def_objs[this_obj]["OBJTYPE"] == "DEFAULT":
            func += su2func(this_obj, config, state) * sign * scale * global_factor
        # otherwise evaluate the penalty function (OBJTYPE = '>','<', or '=')
        else:
            func += obj_p(config, state, this_obj, def_objs) * scale
    vals_out.append(func)

    #: for each objective
    # If evaluating the combined function is desired, update it here.
    # This is only used when OPT_COMBINE_OBJECTIVE = YES
    if "COMBO" in state.FUNCTIONS:
        state["FUNCTIONS"]["COMBO"] = func

    return vals_out


#: def obj_f()


def obj_p(config, state, this_obj, def_objs):
    # Penalty function: square of the difference between value and limit
    # This function is used when a constraint-type term is added to OPT_OBJECTIVE
    # This code, and obj_dp, must be changed to use a non-quadratic penalty function
    funcval = su2func(this_obj, config, state)
    constraint = float(def_objs[this_obj]["VALUE"])
    penalty = 0.0
    if (
        def_objs[this_obj]["OBJTYPE"] == "="
        or (def_objs[this_obj]["OBJTYPE"] == ">" and funcval < constraint)
        or (def_objs[this_obj]["OBJTYPE"] == "<" and funcval > constraint)
    ):
        penalty = (constraint - funcval) ** 2.0
    # If 'DEFAULT' objtype this returns the function value.
    elif def_objs[this_obj]["OBJTYPE"] == "DEFAULT":
        penalty = funcval
    return penalty


#: def obj_p()


def obj_dp(config, state, this_obj, def_objs):
    # Partial Derivative of Penalty function: square of the difference between value and limit
    # This function is used when a constraint-type term is added to OPT_OBJECTIVE
    # This code, and obj_p, must be changed to use a non-quadratic penalty function
    funcval = su2func(this_obj, config, state)
    constraint = float(def_objs[this_obj]["VALUE"])
    dpenalty = 0.0

    # Inequalities will be 0 or a positive value
    if (def_objs[this_obj]["OBJTYPE"] == ">" and funcval < constraint) or (
        def_objs[this_obj]["OBJTYPE"] == "<" and funcval > constraint
    ):
        dpenalty = 2.0 * abs(constraint - funcval)
    # Equalities dp will be positive if value>constraint, negative if value<constraint
    elif def_objs[this_obj]["OBJTYPE"] == "=":
        dpenalty = 2.0 * (funcval - constraint)
    # If 'DEFAULT' objtype, this will return 1.0
    elif def_objs[this_obj]["OBJTYPE"] == "DEFAULT":
        dpenalty = 1.0

    return dpenalty


#: def obj_dp()


def obj_df(dvs, config, state=None):
    """vals = SU2.eval.obj_df(dvs,config,state=None)

    Evaluates SU2 Objective Gradients
    Wraps SU2.eval.grad()

    Takes a design vector for input as a list (shape n)
    or numpy array (shape n or nx1 or 1xn), a config
    and optionally a state.

    Outputs a list of gradients.
    """

    # unpack config and state
    config.unpack_dvs(dvs)
    state = su2io.State(state)
    grad_method = config.get("GRADIENT_METHOD", "CONTINUOUS_ADJOINT")

    def_objs = config["OPT_OBJECTIVE"]
    objectives = def_objs.keys()

    # Number of objective functionals
    n_obj = len(objectives)
    # Whether to calculate gradients one-by-one or all-at-once
    combine_obj = config["OPT_COMBINE_OBJECTIVE"] == "YES"

    dv_scales = config["DEFINITION_DV"]["SCALE"]
    dv_size = config["DEFINITION_DV"]["SIZE"]

    # evaluate each objective
    vals_out = []
    if combine_obj and n_obj > 1:
        # Evaluate objectives all-at-once; for adjoint methods this results in a
        # single, combined objective.
        scale = [1.0] * n_obj
        obj_list = ["DRAG"] * n_obj
        for i_obj, this_obj in enumerate(objectives):
            obj_list[i_obj] = this_obj
            scale[i_obj] = def_objs[this_obj]["SCALE"]
            if def_objs[this_obj]["OBJTYPE"] == "DEFAULT":
                # Standard case
                sign = su2io.get_objectiveSign(this_obj)
                scale[i_obj] *= sign
            else:
                # For a penalty function, the term is scaled by the partial derivative
                # d p(j) / dx = (dj / dx) * ( dp / dj)
                scale[i_obj] *= obj_dp(config, state, this_obj, def_objs)

        config["OBJECTIVE_WEIGHT"] = ",".join(map(str, scale))
        grad = su2grad(obj_list, grad_method, config, state)
        # scaling : obj scale  and sign are accounted for in combo gradient, dv scale now applied
        global_factor = float(config["OPT_GRADIENT_FACTOR"])
        k = 0
        for i_dv, dv_scl in enumerate(dv_scales):
            for i_grd in range(dv_size[i_dv]):
                grad[k] = grad[k] * global_factor / dv_scl
                k = k + 1

        vals_out.append(grad)
    else:
        # Evaluate objectives one-by-one
        marker_monitored = config["MARKER_MONITORING"]
        for i_obj, this_obj in enumerate(objectives):
            # For multiple objectives are evaluated one-by-one rather than combined
            # MARKER_MONITORING should be updated to only include the marker for i_obj
            # For single objectives, multiple markers can be used
            if n_obj > 1:
                config["MARKER_MONITORING"] = marker_monitored[i_obj]
            scale = def_objs[this_obj]["SCALE"]
            global_factor = float(config["OPT_GRADIENT_FACTOR"])
            sign = su2io.get_objectiveSign(this_obj)
            if def_objs[this_obj]["OBJTYPE"] != "DEFAULT":
                # For a penalty function, the term is scaled by the partial derivative
                # and the sign is always positive
                # d p(j) / dx = (dj / dx) * ( dp / dj)
                scale *= obj_dp(config, state, this_obj, def_objs)
                sign = 1.0

            # Evaluate Objective Gradient
            grad = su2grad(this_obj, grad_method, config, state)

            # scaling and sign
            k = 0
            for i_dv, dv_scl in enumerate(dv_scales):
                for i_grd in range(dv_size[i_dv]):
                    grad[k] = grad[k] * sign * scale * global_factor / dv_scl
                    k = k + 1

            vals_out.append(grad)

    #: for each objective

    return vals_out


#: def obj_df()


def con_ceq(dvs, config, state=None):
    """vals = SU2.eval.con_ceq(dvs,config,state=None)

    Evaluates SU2 Equality Constraints
    Wraps SU2.eval.func()

    Takes a design vector for input as a list (shape n)
    or numpy array (shape n or nx1 or 1xn), a config
    and optionally a state.

    Returns: a list of constraint values, ordered
    by the OPT_CONSTRAINT config parameter.
    """

    # unpack state and config
    config.unpack_dvs(dvs)
    state = su2io.State(state)

    def_cons = config["OPT_CONSTRAINT"]["EQUALITY"]
    constraints = def_cons.keys()

    # evaluate each constraint
    vals_out = []
    for i_obj, this_con in enumerate(constraints):
        global_factor = float(config["OPT_GRADIENT_FACTOR"])
        push = def_cons[this_con]["SCALE"]
        value = def_cons[this_con]["VALUE"]

        # Evaluate Constraint Function
        func = su2func(this_con, config, state)

        # scaling and centering
        func = (func - value) * global_factor * push

        vals_out.append(func)

    #: for each constraint

    return vals_out


#: def obj_ceq()


def con_dceq(dvs, config, state=None):
    """vals = SU2.eval.con_dceq(dvs,config,state=None)

    Evaluates SU2 Equality Constraint Gradients
    Wraps SU2.eval.grad()

    Takes a design vector for input as a list (shape n)
    or numpy array (shape n or nx1 or 1xn), a config
    and optionally a state.

    Returns a list of lists of constraint gradients,
    ordered by the OPT_CONSTRAINT config parameter.
    """

    # unpack state and config
    config.unpack_dvs(dvs)
    state = su2io.State(state)
    grad_method = config.get("GRADIENT_METHOD", "CONTINUOUS_ADJOINT")

    def_cons = config["OPT_CONSTRAINT"]["EQUALITY"]
    constraints = def_cons.keys()

    dv_scales = config["DEFINITION_DV"]["SCALE"]
    dv_size = config["DEFINITION_DV"]["SIZE"]

    # evaluate each constraint
    vals_out = []
    for i_obj, this_con in enumerate(constraints):
        global_factor = float(config["OPT_GRADIENT_FACTOR"])
        value = def_cons[this_con]["VALUE"]

        # Evaluate Constraint Gradient
        grad = su2grad(this_con, grad_method, config, state)

        # scaling
        k = 0
        for i_dv, dv_scl in enumerate(dv_scales):
            for i_grd in range(dv_size[i_dv]):
                grad[k] = grad[k] * global_factor / dv_scl
                k = k + 1

        vals_out.append(grad)

    #: for each constraint

    return vals_out


#: def obj_dceq()


def con_cieq(dvs, config, state=None):
    """vals = SU2.eval.con_cieq(dvs,config,state=None)

    Evaluates SU2 Inequality Constraints
    Wraps SU2.eval.func()
    Convention is con(x)<=0

    Takes a design vector for input as a list (shape n)
    or numpy array (shape n or nx1 or 1xn), a config
    and optionally a state.

    Returns a list of constraint gradients, ordered
    by the OPT_CONSTRAINT config parameter.
    """

    # unpack state and config
    config.unpack_dvs(dvs)
    state = su2io.State(state)

    def_cons = config["OPT_CONSTRAINT"]["INEQUALITY"]
    constraints = def_cons.keys()

    # evaluate each constraint
    vals_out = []
    for i_obj, this_con in enumerate(constraints):
        global_factor = float(config["OPT_GRADIENT_FACTOR"])
        push = def_cons[this_con]["SCALE"]
        value = def_cons[this_con]["VALUE"]
        sign = def_cons[this_con]["SIGN"]
        sign = su2io.get_constraintSign(sign)

        # Evaluate Constraint Function
        func = su2func(this_con, config, state)

        # scaling and centering
        func = (func - value) * sign * global_factor * push

        vals_out.append(func)

    #: for each constraint

    return vals_out


#: def obj_cieq()


def con_dcieq(dvs, config, state=None):
    """vals = SU2.eval.con_dceq(dvs,config,state=None)

    Evaluates SU2 Inequality Constraint Gradients
    Wraps SU2.eval.grad()
    Convention is con(x)<=0

    Takes a design vector for input as a list (shape n)
    or numpy array (shape n or nx1 or 1xn), a config
    and optionally a state.

    Returns a list of lists of constraint gradients,
    ordered by the OPT_CONSTRAINT config parameter.
    """

    # unpack state and config
    config.unpack_dvs(dvs)
    state = su2io.State(state)
    grad_method = config.get("GRADIENT_METHOD", "CONTINUOUS_ADJOINT")

    def_cons = config["OPT_CONSTRAINT"]["INEQUALITY"]
    constraints = def_cons.keys()

    dv_scales = config["DEFINITION_DV"]["SCALE"]
    dv_size = config["DEFINITION_DV"]["SIZE"]

    # evaluate each constraint
    vals_out = []
    for i_obj, this_con in enumerate(constraints):
        global_factor = float(config["OPT_GRADIENT_FACTOR"])
        value = def_cons[this_con]["VALUE"]
        sign = def_cons[this_con]["SIGN"]
        sign = su2io.get_constraintSign(sign)

        # Evaluate Constraint Gradient
        grad = su2grad(this_con, grad_method, config, state)

        # scaling and sign
        k = 0
        for i_dv, dv_scl in enumerate(dv_scales):
            for i_grd in range(dv_size[i_dv]):
                grad[k] = grad[k] * sign * global_factor / dv_scl
                k = k + 1

        vals_out.append(grad)

    #: for each constraint

    return vals_out


#: def obj_dcieq()


def touch(config, state):
    """SU2.eval.touch(config,state)
    resets state timestamp
    """
    state.set_timestamp()


def skip(config, state):
    """SU2.eval.skip(config,state)
    does nothing
    """
    pass
