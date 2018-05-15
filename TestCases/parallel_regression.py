#!/usr/bin/env python

## \file parallel_regression.py
#  \brief Python script for automated regression testing of SU2 examples
#  \author A. Aranake, A. Campos, T. Economon, T. Lukaczyk, S. Padron
#  \version 6.0.1 "Falcon"
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
# Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
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

# make print(*args) function available in PY2.6+, does'nt work on PY < 2.6
from __future__ import print_function

import sys
from TestCase import TestCase    

def main():
    '''This program runs SU2 and ensures that the output matches specified values. 
       This will be used to do checks when code is pushed to github 
       to make sure nothing is broken. '''

    test_list = []

    ##########################
    ### Compressible Euler ###
    ##########################

    # Channel
    channel           = TestCase('channel')
    channel.cfg_dir   = "euler/channel"
    channel.cfg_file  = "inv_channel_RK.cfg"
    channel.test_iter = 20
    channel.test_vals = [-2.637944, 2.830321, 0.040796, 0.001679] #last 4 columns
    channel.su2_exec  = "parallel_computation.py -f"
    channel.timeout   = 1600
    channel.tol       = 0.00001
    test_list.append(channel)

    # NACA0012 
    naca0012           = TestCase('naca0012')
    naca0012.cfg_dir   = "euler/naca0012"
    naca0012.cfg_file  = "inv_NACA0012_Roe.cfg"
    naca0012.test_iter = 20
    naca0012.test_vals = [-4.098259, -3.600407, 0.335045, 0.022340] #last 4 columns
    naca0012.su2_exec  = "parallel_computation.py -f"
    naca0012.timeout   = 1600
    naca0012.tol       = 0.00001
    test_list.append(naca0012)

    # Supersonic wedge 
    wedge           = TestCase('wedge')
    wedge.cfg_dir   = "euler/wedge"
    wedge.cfg_file  = "inv_wedge_HLLC.cfg"
    wedge.test_iter = 20
    wedge.test_vals = [-0.804324, 4.937633, -0.250431, 0.044119] #last 4 columns
    wedge.su2_exec  = "parallel_computation.py -f"
    wedge.timeout   = 1600
    wedge.tol       = 0.00001
    test_list.append(wedge)

    # ONERA M6 Wing
    oneram6           = TestCase('oneram6')
    oneram6.cfg_dir   = "euler/oneram6"
    oneram6.cfg_file  = "inv_ONERAM6.cfg"
    oneram6.test_iter = 10
    oneram6.test_vals = [-13.398436, -12.930311, 0.282557, 0.012706] #last 4 columns
    oneram6.su2_exec  = "parallel_computation.py -f"
    oneram6.timeout   = 3200
    oneram6.tol       = 0.00001
    test_list.append(oneram6)

    # Fixed CL NACA0012
    fixedCL_naca0012           = TestCase('fixedcl_naca0012')
    fixedCL_naca0012.cfg_dir   = "fixed_cl/naca0012"
    fixedCL_naca0012.cfg_file  = "inv_NACA0012.cfg"
    fixedCL_naca0012.test_iter = 100
    fixedCL_naca0012.test_vals = [-2.510935, 2.891414, 0.292482, 0.019181] #last 4 columns
    fixedCL_naca0012.su2_exec  = "parallel_computation.py -f"
    fixedCL_naca0012.timeout   = 1600
    fixedCL_naca0012.tol       = 0.00001
    test_list.append(fixedCL_naca0012)
    
    # Polar sweep of the inviscid NACA0012
    polar_naca0012           = TestCase('polar_naca0012')
    polar_naca0012.cfg_dir   = "polar/naca0012"
    polar_naca0012.cfg_file  = "inv_NACA0012.cfg"
    polar_naca0012.polar     = True
    polar_naca0012.test_iter = 10
    polar_naca0012.test_vals = [-1.298934, 4.139523, 0.013808, 0.009041] #last 4 columns
    polar_naca0012.su2_exec  = "compute_polar.py -i 11"
    polar_naca0012.timeout   = 1600
    polar_naca0012.tol       = 0.00001
    test_list.append(polar_naca0012)

    ##########################
    ###  Compressible N-S  ###
    ##########################

    # Laminar flat plate
    flatplate           = TestCase('flatplate')
    flatplate.cfg_dir   = "navierstokes/flatplate"
    flatplate.cfg_file  = "lam_flatplate.cfg"
    flatplate.test_iter = 20
    flatplate.test_vals = [-4.649439, 0.812049, -0.130630, 0.024358] #last 4 columns
    flatplate.su2_exec  = "parallel_computation.py -f"
    flatplate.timeout   = 1600
    flatplate.tol       = 0.00001
    test_list.append(flatplate)

    # Laminar cylinder (steady)
    cylinder           = TestCase('cylinder')
    cylinder.cfg_dir   = "navierstokes/cylinder"
    cylinder.cfg_file  = "lam_cylinder.cfg"
    cylinder.test_iter = 25
    cylinder.test_vals = [-6.756603, -1.288658, -0.131874, 0.509674] #last 4 columns
    cylinder.su2_exec  = "parallel_computation.py -f"
    cylinder.timeout   = 1600
    cylinder.tol       = 0.00001
    test_list.append(cylinder)

    # Laminar cylinder (low Mach correction)
    cylinder_lowmach           = TestCase('cylinder_lowmach')
    cylinder_lowmach.cfg_dir   = "navierstokes/cylinder"
    cylinder_lowmach.cfg_file  = "cylinder_lowmach.cfg"
    cylinder_lowmach.test_iter = 25
    cylinder_lowmach.test_vals = [-6.858759, -1.396741, -1.433192, 110.131583] #last 4 columns
    cylinder_lowmach.su2_exec  = "parallel_computation.py -f"
    cylinder_lowmach.timeout   = 1600
    cylinder_lowmach.tol       = 0.00001
    test_list.append(cylinder_lowmach)

    # 2D Poiseuille flow (body force driven with periodic inlet / outlet)
    poiseuille           = TestCase('poiseuille')
    poiseuille.cfg_dir   = "navierstokes/poiseuille"
    poiseuille.cfg_file  = "lam_poiseuille.cfg"
    poiseuille.test_iter = 10
    poiseuille.test_vals = [-12.027716, -3.326792, -0.000000, 2.351005] #last 4 columns
    poiseuille.su2_exec  = "parallel_computation.py -f"
    poiseuille.timeout   = 1600
    poiseuille.tol       = 0.001
    test_list.append(poiseuille)

    ##########################
    ### Compressible RANS  ###
    ##########################

    # RAE2822 SA
    rae2822_sa           = TestCase('rae2822_sa')
    rae2822_sa.cfg_dir   = "rans/rae2822"
    rae2822_sa.cfg_file  = "turb_SA_RAE2822.cfg"
    rae2822_sa.test_iter = 20
    rae2822_sa.test_vals = [-2.003666, -5.226664, 0.832084, 0.053819] #last 4 columns
    rae2822_sa.su2_exec  = "parallel_computation.py -f"
    rae2822_sa.timeout   = 1600
    rae2822_sa.tol       = 0.00001
    test_list.append(rae2822_sa)
    
    # RAE2822 SST
    rae2822_sst           = TestCase('rae2822_sst')
    rae2822_sst.cfg_dir   = "rans/rae2822"
    rae2822_sst.cfg_file  = "turb_SST_RAE2822.cfg"
    rae2822_sst.test_iter = 20
    rae2822_sst.test_vals = [-0.510807, 4.909622, 0.835523, 0.054264] #last 4 columns
    rae2822_sst.su2_exec  = "parallel_computation.py -f"
    rae2822_sst.timeout   = 1600
    rae2822_sst.tol       = 0.00001
    test_list.append(rae2822_sst)

    # Flat plate
    turb_flatplate           = TestCase('turb_flatplate')
    turb_flatplate.cfg_dir   = "rans/flatplate"
    turb_flatplate.cfg_file  = "turb_SA_flatplate.cfg"
    turb_flatplate.test_iter = 20
    turb_flatplate.test_vals = [-4.146902, -6.736668, -0.176251, 0.057569] #last 4 columns
    turb_flatplate.su2_exec  = "parallel_computation.py -f"
    turb_flatplate.timeout   = 1600
    turb_flatplate.tol       = 0.00001
    test_list.append(turb_flatplate)

    # ONERA M6 Wing
    turb_oneram6           = TestCase('turb_oneram6')
    turb_oneram6.cfg_dir   = "rans/oneram6"
    turb_oneram6.cfg_file  = "turb_ONERAM6.cfg"
    turb_oneram6.test_iter = 10
    turb_oneram6.test_vals = [-2.327509, -6.563372, 0.230438, 0.155815] #last 4 columns
    turb_oneram6.su2_exec  = "parallel_computation.py -f"
    turb_oneram6.timeout   = 3200
    turb_oneram6.tol       = 0.00001
    test_list.append(turb_oneram6)

    # NACA0012 (SA, FUN3D finest grid results: CL=1.0983, CD=0.01242)
    turb_naca0012_sa           = TestCase('turb_naca0012_sa')
    turb_naca0012_sa.cfg_dir   = "rans/naca0012"
    turb_naca0012_sa.cfg_file  = "turb_NACA0012_sa.cfg"
    turb_naca0012_sa.test_iter = 10
    turb_naca0012_sa.test_vals = [-12.000763, -9.145363, 1.070528, 0.019417] #last 4 columns
    turb_naca0012_sa.su2_exec  = "parallel_computation.py -f"
    turb_naca0012_sa.timeout   = 3200
    turb_naca0012_sa.tol       = 0.00001
    test_list.append(turb_naca0012_sa)

    # NACA0012 (SA, FUN3D results for finest grid: CL=1.0983, CD=0.01242) with binary restart
    turb_naca0012_sa_bin           = TestCase('turb_naca0012_sa_bin')
    turb_naca0012_sa_bin.cfg_dir   = "rans/naca0012"
    turb_naca0012_sa_bin.cfg_file  = "turb_NACA0012_sa_binary.cfg"
    turb_naca0012_sa_bin.test_iter = 10
    turb_naca0012_sa_bin.test_vals = [-12.000899, -9.145363, 1.070528, 0.019417] #last 4 columns
    turb_naca0012_sa_bin.su2_exec  = "parallel_computation.py -f"
    turb_naca0012_sa_bin.timeout   = 3200
    turb_naca0012_sa_bin.tol       = 0.00001
    test_list.append(turb_naca0012_sa_bin)
    
    # NACA0012 (SST, FUN3D finest grid results: CL=1.0840, CD=0.01253)
    turb_naca0012_sst           = TestCase('turb_naca0012_sst')
    turb_naca0012_sst.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst.cfg_file  = "turb_NACA0012_sst.cfg"
    turb_naca0012_sst.test_iter = 10
    turb_naca0012_sst.test_vals = [-15.039645, -7.220177, 1.059622, 0.019138] #last 4 columns
    turb_naca0012_sst.su2_exec  = "parallel_computation.py -f"
    turb_naca0012_sst.timeout   = 3200
    turb_naca0012_sst.tol       = 0.00001
    test_list.append(turb_naca0012_sst)

    # PROPELLER
    propeller           = TestCase('propeller')
    propeller.cfg_dir   = "rans/propeller"
    propeller.cfg_file  = "propeller.cfg"
    propeller.test_iter = 10
    propeller.test_vals = [-3.378804, -8.130406, 0.000051, 0.055764] #last 4 columns
    propeller.su2_exec  = "parallel_computation.py -f"
    propeller.timeout   = 3200
    propeller.tol       = 0.00001
    test_list.append(propeller)
    
    #################################
    ## Compressible RANS Restart  ###
    #################################
    
    # NACA0012 SST Multigrid restart
    turb_naca0012_sst_restart_mg           = TestCase('turb_naca0012_sst_restart_mg')
    turb_naca0012_sst_restart_mg.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst_restart_mg.cfg_file  = "turb_NACA0012_sst_multigrid_restart.cfg"
    turb_naca0012_sst_restart_mg.test_iter = 20
    turb_naca0012_sst_restart_mg.ntest_vals = 5
    turb_naca0012_sst_restart_mg.test_vals = [-6.436237, -4.558626, 1.231779, -0.007019, 0.081472] #last 5 columns
    turb_naca0012_sst_restart_mg.su2_exec  = "parallel_computation.py -f"
    turb_naca0012_sst_restart_mg.timeout   = 3200
    turb_naca0012_sst_restart_mg.tol       = 0.000001
    test_list.append(turb_naca0012_sst_restart_mg)

    #############################
    ### Incompressible Euler  ###
    #############################

    # NACA0012 Hydrofoil
    inc_euler_naca0012           = TestCase('inc_euler_naca0012')
    inc_euler_naca0012.cfg_dir   = "incomp_euler/naca0012"
    inc_euler_naca0012.cfg_file  = "incomp_NACA0012.cfg"
    inc_euler_naca0012.test_iter = 20
    inc_euler_naca0012.test_vals = [-4.843880, -3.769636, 0.505456, 0.007283] #last 4 columns
    inc_euler_naca0012.su2_exec  = "parallel_computation.py -f"
    inc_euler_naca0012.timeout   = 1600
    inc_euler_naca0012.tol       = 0.00001
    test_list.append(inc_euler_naca0012)

    #############################
    ### Incompressible N-S    ###
    #############################

    # Laminar cylinder
    inc_lam_cylinder          = TestCase('inc_lam_cylinder')
    inc_lam_cylinder.cfg_dir   = "incomp_navierstokes/cylinder"
    inc_lam_cylinder.cfg_file  = "incomp_cylinder.cfg"
    inc_lam_cylinder.test_iter = 10
    inc_lam_cylinder.test_vals = [-3.989654, -3.177372, 0.002595, 7.748101] #last 4 columns
    inc_lam_cylinder.su2_exec  = "parallel_computation.py -f"
    inc_lam_cylinder.timeout   = 1600
    inc_lam_cylinder.tol       = 0.00001
    test_list.append(inc_lam_cylinder)

    ############################
    ### Incompressible RANS  ###
    ############################

    # NACA0012
    inc_turb_naca0012           = TestCase('inc_turb_naca0012')
    inc_turb_naca0012.cfg_dir   = "incomp_rans/naca0012"
    inc_turb_naca0012.cfg_file  = "naca0012.cfg"
    inc_turb_naca0012.test_iter = 20
    inc_turb_naca0012.test_vals = [-4.788498, -11.040505, 0.000003, 0.309508] #last 4 columns
    inc_turb_naca0012.su2_exec  = "parallel_computation.py -f"
    inc_turb_naca0012.timeout   = 1600
    inc_turb_naca0012.tol       = 0.00001
    test_list.append(inc_turb_naca0012)

    ############################
    ###      Transition      ###
    ############################

    # Schubauer-Klebanoff Natural Transition Case
    schubauer_klebanoff_transition              = TestCase('Schubauer_Klebanoff')
    schubauer_klebanoff_transition.cfg_dir      = "transition/Schubauer_Klebanoff"
    schubauer_klebanoff_transition.cfg_file     = "transitional_BC_model_ConfigFile.cfg"
    schubauer_klebanoff_transition.test_iter    = 10
    schubauer_klebanoff_transition.test_vals    = [-8.225238, -14.278195, 0.000045, 0.007987] #last 4 columns
    schubauer_klebanoff_transition.su2_exec     = "parallel_computation.py -f"
    schubauer_klebanoff_transition.timeout      = 1600
    schubauer_klebanoff_transition.tol          = 0.00001
    test_list.append(schubauer_klebanoff_transition)

    #####################################
    ### Cont. adj. compressible Euler ###
    #####################################

    # Inviscid NACA0012
    contadj_naca0012           = TestCase('contadj_naca0012')
    contadj_naca0012.cfg_dir   = "cont_adj_euler/naca0012"
    contadj_naca0012.cfg_file  = "inv_NACA0012.cfg"
    contadj_naca0012.test_iter = 5
    contadj_naca0012.test_vals = [-9.783705, -15.193027, 0.300920, 0.019552] #last 4 columns
    contadj_naca0012.su2_exec  = "parallel_computation.py -f"
    contadj_naca0012.timeout   = 1600
    contadj_naca0012.tol       = 0.00001
    test_list.append(contadj_naca0012)

    # Inviscid ONERA M6
    contadj_oneram6           = TestCase('contadj_oneram6')
    contadj_oneram6.cfg_dir   = "cont_adj_euler/oneram6"
    contadj_oneram6.cfg_file  = "inv_ONERAM6.cfg"
    contadj_oneram6.test_iter = 10
    contadj_oneram6.test_vals = [-12.130871, -12.702611, 0.685900, 0.007594] #last 4 columns
    contadj_oneram6.su2_exec  = "parallel_computation.py -f"
    contadj_oneram6.timeout   = 1600
    contadj_oneram6.tol       = 0.00001
    test_list.append(contadj_oneram6)

    # Inviscid WEDGE: tests averaged outflow total pressure adjoint
    contadj_wedge             = TestCase('contadj_wedge')
    contadj_wedge.cfg_dir   = "cont_adj_euler/wedge"
    contadj_wedge.cfg_file  = "inv_wedge_ROE.cfg"
    contadj_wedge.test_iter = 10  
    contadj_wedge.test_vals = [2.856008, -2.767216, 1.0029e+06, 7.0328e-14] #last 4 columns
    contadj_wedge.su2_exec  = "parallel_computation.py -f"
    contadj_wedge.timeout   = 1600
    contadj_wedge.tol       = 0.00001
    test_list.append(contadj_wedge)

    # Inviscid fixed CL NACA0012
    contadj_fixed_CL_naca0012           = TestCase('contadj_fixedcl_naca0012')
    contadj_fixed_CL_naca0012.cfg_dir   = "fixed_cl/naca0012"
    contadj_fixed_CL_naca0012.cfg_file  = "inv_NACA0012_ContAdj.cfg"
    contadj_fixed_CL_naca0012.test_iter = 100
    contadj_fixed_CL_naca0012.test_vals = [0.326574, -5.209644, 0.252580, -0.000078] #last 4 columns
    contadj_fixed_CL_naca0012.su2_exec  = "parallel_computation.py -f"
    contadj_fixed_CL_naca0012.timeout   = 1600
    contadj_fixed_CL_naca0012.tol       = 0.00001
    test_list.append(contadj_fixed_CL_naca0012)

    ###################################
    ### Cont. adj. compressible N-S ###
    ###################################

    # Adjoint laminar cylinder
    contadj_ns_cylinder           = TestCase('contadj_ns_cylinder')
    contadj_ns_cylinder.cfg_dir   = "cont_adj_navierstokes/cylinder"
    contadj_ns_cylinder.cfg_file  = "lam_cylinder.cfg"
    contadj_ns_cylinder.test_iter = 20
    contadj_ns_cylinder.test_vals = [-3.648105, -9.104817, 2.056700, -0.000000] #last 4 columns
    contadj_ns_cylinder.su2_exec  = "parallel_computation.py -f"
    contadj_ns_cylinder.timeout   = 1600
    contadj_ns_cylinder.tol       = 0.00001
    test_list.append(contadj_ns_cylinder)

    # Adjoint laminar naca0012 subsonic
    contadj_ns_naca0012_sub           = TestCase('contadj_ns_naca0012_sub')
    contadj_ns_naca0012_sub.cfg_dir   = "cont_adj_navierstokes/naca0012_sub"
    contadj_ns_naca0012_sub.cfg_file  = "lam_NACA0012.cfg"
    contadj_ns_naca0012_sub.test_iter = 20
    contadj_ns_naca0012_sub.test_vals = [-2.743268, -8.215193, 0.518810, 0.001210] #last 4 columns
    contadj_ns_naca0012_sub.su2_exec  = "parallel_computation.py -f"
    contadj_ns_naca0012_sub.timeout   = 1600
    contadj_ns_naca0012_sub.tol       = 0.00001
    test_list.append(contadj_ns_naca0012_sub)
    
    # Adjoint laminar naca0012 transonic
    contadj_ns_naca0012_trans           = TestCase('contadj_ns_naca0012_trans')
    contadj_ns_naca0012_trans.cfg_dir   = "cont_adj_navierstokes/naca0012_trans"
    contadj_ns_naca0012_trans.cfg_file  = "lam_NACA0012.cfg"
    contadj_ns_naca0012_trans.test_iter = 20
    contadj_ns_naca0012_trans.test_vals = [ -1.039664, -6.575019, 1.772300, 0.012495] #last 4 columns
    contadj_ns_naca0012_trans.su2_exec  = "parallel_computation.py -f"
    contadj_ns_naca0012_trans.timeout   = 1600
    contadj_ns_naca0012_trans.tol       = 0.00001
    test_list.append(contadj_ns_naca0012_trans)
    
    #######################################################
    ### Cont. adj. compressible RANS (frozen viscosity) ###
    #######################################################

    # Adjoint turbulent NACA0012
    contadj_rans_naca0012           = TestCase('contadj_rans_naca0012')
    contadj_rans_naca0012.cfg_dir   = "cont_adj_rans/naca0012"
    contadj_rans_naca0012.cfg_file  = "turb_nasa.cfg"
    contadj_rans_naca0012.test_iter = 20
    contadj_rans_naca0012.test_vals = [ -0.794162, -5.761722, 19.214000, -0.000000] #last 4 columns
    contadj_rans_naca0012.su2_exec  = "parallel_computation.py -f"
    contadj_rans_naca0012.timeout   = 1600
    contadj_rans_naca0012.tol       = 0.00001
    test_list.append(contadj_rans_naca0012)
   
    # Adjoint turbulent NACA0012 with binary restarts
    contadj_rans_naca0012_bin           = TestCase('contadj_rans_naca0012_bin')
    contadj_rans_naca0012_bin.cfg_dir   = "cont_adj_rans/naca0012"
    contadj_rans_naca0012_bin.cfg_file  = "turb_nasa_binary.cfg"
    contadj_rans_naca0012_bin.test_iter = 20
    contadj_rans_naca0012_bin.test_vals = [-0.794169, -5.761671, 19.214000, -0.000000] #last 4 columns
    contadj_rans_naca0012_bin.su2_exec  = "parallel_computation.py -f"
    contadj_rans_naca0012_bin.timeout   = 1600
    contadj_rans_naca0012_bin.tol       = 0.00001
    test_list.append(contadj_rans_naca0012_bin)
 
    # Adjoint turbulent RAE2822
    contadj_rans_rae2822           = TestCase('contadj_rans_rae822')
    contadj_rans_rae2822.cfg_dir   = "cont_adj_rans/rae2822"
    contadj_rans_rae2822.cfg_file  = "turb_SA_RAE2822.cfg"
    contadj_rans_rae2822.test_iter = 20
    contadj_rans_rae2822.test_vals = [-5.369195, -10.871609, -0.212470, 0.005448] #last 4 columns
    contadj_rans_rae2822.su2_exec  = "parallel_computation.py -f"
    contadj_rans_rae2822.timeout   = 1600
    contadj_rans_rae2822.tol       = 0.00001
    test_list.append(contadj_rans_rae2822)

    ######################################                                                                                  
    ### Harmonic Balance               ###                                                                                  
    ######################################                                                                                    

    # Description of the regression test 
    harmonic_balance           = TestCase('harmonic_balance')
    harmonic_balance.cfg_dir   = "harmonic_balance"
    harmonic_balance.cfg_file  = "HB.cfg"
    harmonic_balance.test_iter = 25
    harmonic_balance.test_vals = [-1.569573, 3.941896, 0.008780, 0.079775] #last 4 columns
    harmonic_balance.su2_exec  = "parallel_computation.py -f"
    harmonic_balance.timeout   = 1600
    harmonic_balance.tol       = 0.00001
    test_list.append(harmonic_balance)

    # Turbulent pitching NACA 64a010 airfoil
    hb_rans_preconditioning           = TestCase('hb_rans_preconditioning')
    hb_rans_preconditioning.cfg_dir   = "harmonic_balance/hb_rans_preconditioning"
    hb_rans_preconditioning.cfg_file  = "davis.cfg"
    hb_rans_preconditioning.test_iter = 25
    hb_rans_preconditioning.test_vals = [-1.900976,-5.876722, 0.007761, 0.125934] #last 4 columns
    hb_rans_preconditioning.su2_exec  = "parallel_computation.py -f"
    hb_rans_preconditioning.timeout   = 1600
    hb_rans_preconditioning.tol       = 0.00001
    test_list.append(hb_rans_preconditioning)

    ######################################
    ### Moving Wall                    ###
    ######################################

    # Lid-driven cavity
    cavity           = TestCase('cavity')
    cavity.cfg_dir   = "moving_wall/cavity"
    cavity.cfg_file  = "lam_cavity.cfg"
    cavity.test_iter = 25
    cavity.test_vals = [-5.596392, -0.132416, 0.140287, 0.832445] #last 4 columns
    cavity.su2_exec  = "parallel_computation.py -f"
    cavity.timeout   = 1600
    cavity.tol       = 0.00001
    test_list.append(cavity)

    # Spinning cylinder
    spinning_cylinder           = TestCase('spinning_cylinder')
    spinning_cylinder.cfg_dir   = "moving_wall/spinning_cylinder"
    spinning_cylinder.cfg_file  = "spinning_cylinder.cfg"
    spinning_cylinder.test_iter = 25
    spinning_cylinder.test_vals = [-7.710702, -2.264970, 1.596921, 1.584878] #last 4 columns
    spinning_cylinder.su2_exec  = "parallel_computation.py -f"
    spinning_cylinder.timeout   = 1600
    spinning_cylinder.tol       = 0.00001
    test_list.append(spinning_cylinder)

    ######################################
    ### Unsteady                       ###
    ######################################

    # Square cylinder
    square_cylinder           = TestCase('square_cylinder')
    square_cylinder.cfg_dir   = "unsteady/square_cylinder"
    square_cylinder.cfg_file  = "turb_square.cfg"
    square_cylinder.test_iter = 3
    square_cylinder.test_vals = [-1.166475, 0.076766, 1.398549, 2.197049] #last 4 columns
    square_cylinder.su2_exec  = "parallel_computation.py -f"
    square_cylinder.timeout   = 1600
    square_cylinder.tol       = 0.00001
    square_cylinder.unsteady  = True
    test_list.append(square_cylinder)

    # Gust
    sine_gust           = TestCase('sine_gust')
    sine_gust.cfg_dir   = "gust"
    sine_gust.cfg_file  = "inv_gust_NACA0012.cfg"
    sine_gust.test_iter = 5
    sine_gust.test_vals = [-1.977531, 3.481790, -0.008934, -0.004244] #last 4 columns
    sine_gust.su2_exec  = "parallel_computation.py -f"
    sine_gust.timeout   = 1600
    sine_gust.tol       = 0.00001
    sine_gust.unsteady  = True
    test_list.append(sine_gust)

    # Aeroelastic
    aeroelastic           = TestCase('aeroelastic')
    aeroelastic.cfg_dir   = "aeroelastic"
    aeroelastic.cfg_file  = "aeroelastic_NACA64A010.cfg"
    aeroelastic.test_iter = 2
    aeroelastic.test_vals = [0.077169, 0.036452, -0.001685, -0.000113] #last 4 columns
    aeroelastic.su2_exec  = "parallel_computation.py -f"
    aeroelastic.timeout   = 1600
    aeroelastic.tol       = 0.000001
    aeroelastic.unsteady  = True
    test_list.append(aeroelastic)

    # Delayed Detached Eddy Simulation
    ddes_flatplate        = TestCase('ddes_flatplate')
    ddes_flatplate.cfg_dir   = "ddes/flatplate"
    ddes_flatplate.cfg_file  = "ddes_flatplate.cfg"
    ddes_flatplate.test_iter = 10
    ddes_flatplate.test_vals = [-2.714721, -5.883008, -0.214968, 0.023783] #last 4 columns
    ddes_flatplate.su2_exec  = "parallel_computation.py -f"
    ddes_flatplate.timeout   = 1600
    ddes_flatplate.tol       = 0.00001
    ddes_flatplate.unsteady  = True
    test_list.append(ddes_flatplate)    

    ######################################
    ### NICFD                          ###
    ######################################	

    # Rarefaction shock wave edge_VW
    edge_VW           = TestCase('edge_VW')
    edge_VW.cfg_dir   = "nicf/edge"
    edge_VW.cfg_file  = "edge_VW.cfg"
    edge_VW.test_iter = 100
    edge_VW.test_vals = [-5.063606, 1.110059, -0.000009, 0.000000] #last 4 columns
    edge_VW.su2_exec  = "parallel_computation.py -f"
    edge_VW.timeout   = 1600
    edge_VW.tol       = 0.00001
    test_list.append(edge_VW)

    # Rarefaction shock wave edge_PPR
    edge_PPR           = TestCase('edge_PPR')
    edge_PPR.cfg_dir   = "nicf/edge"
    edge_PPR.cfg_file  = "edge_PPR.cfg"
    edge_PPR.test_iter = 100
    edge_PPR.test_vals = [-5.485000, 0.656415, -0.000037, 0.000000] #last 4 columns
    edge_PPR.su2_exec  = "parallel_computation.py -f"
    edge_PPR.timeout   = 1600
    edge_PPR.tol       = 0.00001
    test_list.append(edge_PPR)
    
    ######################################
    ### Turbomachinery                 ###
    ######################################	

    # Jones APU Turbocharger
    Jones_tc           = TestCase('jones_turbocharger')
    Jones_tc.cfg_dir   = "turbomachinery/APU_turbocharger"
    Jones_tc.cfg_file  = "Jones.cfg"
    Jones_tc.test_iter = 5
    Jones_tc.test_vals = [-5.294493, 0.429399, 80.085030, 1.034808] #last 4 columns
    Jones_tc.su2_exec  = "parallel_computation.py -f"
    Jones_tc.timeout   = 1600
    Jones_tc.tol       = 0.00001
    test_list.append(Jones_tc)

	# Jones APU Turbocharger restart
    Jones_tc_rst           = TestCase('jones_turbocharger_restart')
    Jones_tc_rst.cfg_dir   = "turbomachinery/APU_turbocharger"
    Jones_tc_rst.cfg_file  = "Jones_rst.cfg"
    Jones_tc_rst.test_iter = 5
    Jones_tc_rst.test_vals = [-4.409575, -1.605695, 82.250580, 2.791302] #last 4 columns
    Jones_tc_rst.su2_exec  = "parallel_computation.py -f"
    Jones_tc_rst.timeout   = 1600
    Jones_tc_rst.tol       = 0.00001
    test_list.append(Jones_tc_rst)

    # 2D axial stage
    axial_stage2D           = TestCase('axial_stage2D')
    axial_stage2D.cfg_dir   = "turbomachinery/axial_stage_2D"
    axial_stage2D.cfg_file  = "Axial_stage2D.cfg"
    axial_stage2D.test_iter = 20
    axial_stage2D.test_vals = [-1.835157, 5.788694, 73.679930, 0.888906] #last 4 columns
    axial_stage2D.su2_exec  = "parallel_computation.py -f"
    axial_stage2D.timeout   = 1600
    axial_stage2D.tol       = 0.00001
    test_list.append(axial_stage2D)
    
    # 2D transonic stator
    transonic_stator           = TestCase('transonic_stator')
    transonic_stator.cfg_dir   = "turbomachinery/transonic_stator_2D"
    transonic_stator.cfg_file  = "transonic_stator.cfg"
    transonic_stator.test_iter = 20
    transonic_stator.test_vals = [-1.102837, 6.115650, 67.308710, 0.071862] #last 4 columns
    transonic_stator.su2_exec  = "parallel_computation.py -f"
    transonic_stator.timeout   = 1600
    transonic_stator.tol       = 0.00001
    test_list.append(transonic_stator)
    
    # 2D transonic stator restart
    transonic_stator_rst           = TestCase('transonic_stator_restart')
    transonic_stator_rst.cfg_dir   = "turbomachinery/transonic_stator_2D"
    transonic_stator_rst.cfg_file  = "transonic_stator_rst.cfg"
    transonic_stator_rst.test_iter = 20
    transonic_stator_rst.test_vals = [-0.484778, 4.462654, 6.460856, 0.004009] #last 4 columns
    transonic_stator_rst.su2_exec  = "parallel_computation.py -f"
    transonic_stator_rst.timeout   = 1600
    transonic_stator_rst.tol       = 0.00001
    test_list.append(transonic_stator_rst)

    ######################################
    ### Sliding Mesh                   ###
    ######################################

    # Uniform flow
    uniform_flow         = TestCase('uniform_flow')
    uniform_flow.cfg_dir   = "sliding_interface/uniform_flow"
    uniform_flow.cfg_file  = "uniform_NN.cfg"
    uniform_flow.test_iter = 50
    uniform_flow.test_vals = [-0.368836, 5.156090, 0.000000, 0.000000] #last 4 columns
    uniform_flow.su2_exec  = "parallel_computation.py -f"
    uniform_flow.timeout   = 1600
    uniform_flow.tol       = 0.000001
    uniform_flow.unsteady  = True
    test_list.append(uniform_flow) 

    # Channel_2D
    channel_2D           = TestCase('channel_2D')
    channel_2D.cfg_dir   = "sliding_interface/channel_2D"
    channel_2D.cfg_file  = "channel_2D_WA.cfg"
    channel_2D.test_iter = 4
    channel_2D.test_vals = [-1.656732, 4.263236, 0.000000, 0.000000] #last 4 columns
    channel_2D.su2_exec  = "parallel_computation.py -f"
    channel_2D.timeout   = 100
    channel_2D.tol       = 0.00001
    channel_2D.unsteady  = True
    test_list.append(channel_2D)

    # Channel_3D
    channel_3D           = TestCase('channel_3D')
    channel_3D.cfg_dir   = "sliding_interface/channel_3D"
    channel_3D.cfg_file  = "channel_3D_WA.cfg"
    channel_3D.test_iter = 2
    channel_3D.test_vals = [-1.999171, 3.956649, 0.000000, 0.000000] #last 4 columns
    channel_3D.su2_exec  = "parallel_computation.py -f"
    channel_3D.timeout   = 1600
    channel_3D.tol       = 0.00001
    channel_3D.unsteady  = True
    test_list.append(channel_3D)

    # Pipe
    pipe           = TestCase('pipe')
    pipe.cfg_dir   = "sliding_interface/pipe"
    pipe.cfg_file  = "pipe_NN.cfg"
    pipe.test_iter = 2
    pipe.test_vals = [-3.503708, 3.194241, 0.000000, 0.000000] #last 4 columns
    pipe.su2_exec  = "parallel_computation.py -f"
    pipe.timeout   = 1600
    pipe.tol       = 0.00001
    pipe.unsteady  = True
    test_list.append(pipe)

    # Rotating cylinders
    rotating_cylinders           = TestCase('rotating_cylinders')
    rotating_cylinders.cfg_dir   = "sliding_interface/rotating_cylinders"
    rotating_cylinders.cfg_file  = "rot_cylinders_WA.cfg"
    rotating_cylinders.test_iter = 3
    rotating_cylinders.test_vals = [-0.888640, 4.881651, 0.000000, 0.000000] #last 4 columns
    rotating_cylinders.su2_exec  = "parallel_computation.py -f"
    rotating_cylinders.timeout   = 1600
    rotating_cylinders.tol       = 0.00001
    rotating_cylinders.unsteady  = True
    test_list.append(rotating_cylinders)

    # Supersonic vortex shedding
    supersonic_vortex_shedding           = TestCase('supersonic_vortex_shedding')
    supersonic_vortex_shedding.cfg_dir   = "sliding_interface/supersonic_vortex_shedding"
    supersonic_vortex_shedding.cfg_file  = "sup_vor_shed_WA.cfg"
    supersonic_vortex_shedding.test_iter = 5
    supersonic_vortex_shedding.test_vals = [-1.142212, 4.585974, 0.000000, 0.000000] #last 4 columns
    supersonic_vortex_shedding.su2_exec  = "parallel_computation.py -f"
    supersonic_vortex_shedding.timeout   = 1600
    supersonic_vortex_shedding.tol       = 0.00001
    supersonic_vortex_shedding.unsteady  = True
    test_list.append(supersonic_vortex_shedding)

    # Bars_SST_2D
    bars_SST_2D           = TestCase('bars_SST_2D')
    bars_SST_2D.cfg_dir   = "sliding_interface/bars_SST_2D"
    bars_SST_2D.cfg_file  = "bars.cfg"
    bars_SST_2D.test_iter = 13
    bars_SST_2D.test_vals = [-2.133553, 1.644733, -0.000831, 0.117497] #last 4 columns
    bars_SST_2D.su2_exec  = "SU2_CFD"
    bars_SST_2D.timeout   = 1600
    bars_SST_2D.tol       = 0.00001
    test_list.append(bars_SST_2D)

    ##########################
    ### FEA - FSI          ###
    ##########################   

    # Static beam, 3d
    statbeam3d           = TestCase('statbeam3d')
    statbeam3d.cfg_dir   = "fea_fsi/StatBeam_3d"
    statbeam3d.cfg_file  = "configBeam_3d.cfg"
    statbeam3d.test_iter = 0
    statbeam3d.test_vals = [-8.326984, -8.110072, -8.051135, 64095.000000] #last 4 columns
    statbeam3d.su2_exec  = "parallel_computation_fsi.py -f"
    statbeam3d.timeout   = 1600
    statbeam3d.tol       = 0.00001
    test_list.append(statbeam3d)

    # Dynamic beam, 2d
    dynbeam2d           = TestCase('dynbeam2d')
    dynbeam2d.cfg_dir   = "fea_fsi/DynBeam_2d"
    dynbeam2d.cfg_file  = "configBeam_2d.cfg"
    dynbeam2d.test_iter = 6
    dynbeam2d.test_vals = [-9.420641, -5.365871, -12.430382, 6.5210e+04] #last 4 columns
    dynbeam2d.su2_exec  = "parallel_computation_fsi.py -f"
    dynbeam2d.timeout   = 1600
    dynbeam2d.tol       = 0.00001
    test_list.append(dynbeam2d)

    # FSI, 2d
    fsi2d           = TestCase('fsi2d')
    fsi2d.cfg_dir   = "fea_fsi/WallChannel_2d"
    fsi2d.cfg_file  = "configFSI_2D.cfg"
    fsi2d.test_iter = 4
    fsi2d.test_vals = [2.000000, 0.500000, -7.780230, -1.142095] #last 4 columns
    fsi2d.su2_exec  = "parallel_computation_fsi.py -f"
    fsi2d.timeout   = 1600
    fsi2d.tol       = 0.00001
    test_list.append(fsi2d)

    ##########################
    ### Zonal multiphysics ###
    ##########################

    # CHT incompressible
    cht_incompressible           = TestCase('cht_incompressible')
    cht_incompressible.cfg_dir   = "coupled_cht/incompressible"
    cht_incompressible.cfg_file  = "coupled_cht_incompressible.cfg"
    cht_incompressible.test_iter = 10
    cht_incompressible.test_vals = [0.000000, 0.000000, -7.962470, -4634.558027] #last 4 columns
    cht_incompressible.su2_exec  = "parallel_computation.py -f"
    cht_incompressible.timeout   = 1600
    cht_incompressible.tol       = 0.00001
    test_list.append(cht_incompressible)

    ##########################
    ###   Python wrapper   ###
    ##########################
    
    # NACA0012 
    pywrapper_naca0012           = TestCase('pywrapper_naca0012')
    pywrapper_naca0012.cfg_dir   = "euler/naca0012"
    pywrapper_naca0012.cfg_file  = "inv_NACA0012_Roe.cfg"
    pywrapper_naca0012.test_iter = 100
    pywrapper_naca0012.test_vals = [-6.192069, -5.593575, 0.334797, 0.022183] #last 4 columns
    pywrapper_naca0012.su2_exec  = "mpirun -np 2 SU2_CFD.py --parallel -f"
    pywrapper_naca0012.timeout   = 1600
    pywrapper_naca0012.tol       = 0.00001
    test_list.append(pywrapper_naca0012)

    # NACA0012 (SST, FUN3D results for finest grid: CL=1.0840, CD=0.01253)
    pywrapper_turb_naca0012_sst           = TestCase('pywrapper_turb_naca0012_sst')
    pywrapper_turb_naca0012_sst.cfg_dir   = "rans/naca0012"
    pywrapper_turb_naca0012_sst.cfg_file  = "turb_NACA0012_sst.cfg"
    pywrapper_turb_naca0012_sst.test_iter = 10
    pywrapper_turb_naca0012_sst.test_vals = [-15.039645, -7.220177, 1.059622, 0.019138] #last 4 columns
    pywrapper_turb_naca0012_sst.su2_exec  = "mpirun -np 2 SU2_CFD.py --parallel -f"
    pywrapper_turb_naca0012_sst.timeout   = 3200
    pywrapper_turb_naca0012_sst.tol       = 0.00001
    test_list.append(pywrapper_turb_naca0012_sst)

    # Square cylinder
    pywrapper_square_cylinder           = TestCase('pywrapper_square_cylinder')
    pywrapper_square_cylinder.cfg_dir   = "unsteady/square_cylinder"
    pywrapper_square_cylinder.cfg_file  = "turb_square.cfg"
    pywrapper_square_cylinder.test_iter = 3
    pywrapper_square_cylinder.test_vals = [-1.166475, 0.076766, 1.398549, 2.197049] #last 4 columns
    pywrapper_square_cylinder.su2_exec  = "mpirun -np 2 SU2_CFD.py --parallel -f"
    pywrapper_square_cylinder.timeout   = 1600
    pywrapper_square_cylinder.tol       = 0.00001
    pywrapper_square_cylinder.unsteady  = True
    test_list.append(pywrapper_square_cylinder)

    # Aeroelastic
    pywrapper_aeroelastic         = TestCase('pywrapper_aeroelastic')
    pywrapper_aeroelastic.cfg_dir   = "aeroelastic"
    pywrapper_aeroelastic.cfg_file  = "aeroelastic_NACA64A010.cfg"
    pywrapper_aeroelastic.test_iter = 2
    pywrapper_aeroelastic.test_vals = [0.077169, 0.036452, -0.001685, -0.000113] #last 4 columns
    pywrapper_aeroelastic.su2_exec  = "mpirun -np 2 SU2_CFD.py --parallel -f"
    pywrapper_aeroelastic.timeout   = 1600
    pywrapper_aeroelastic.tol       = 0.000001
    pywrapper_aeroelastic.unsteady  = True
    test_list.append(pywrapper_aeroelastic)

    # FSI, 2d
    pywrapper_fsi2d           = TestCase('pywrapper_fsi2d')
    pywrapper_fsi2d.cfg_dir   = "fea_fsi/WallChannel_2d"
    pywrapper_fsi2d.cfg_file  = "configFSI_2D.cfg"
    pywrapper_fsi2d.test_iter = 4
    pywrapper_fsi2d.test_vals = [2.000000, 0.500000, -7.780230, -1.142095] #last 4 columns
    pywrapper_fsi2d.su2_exec  = "mpirun -np 2 SU2_CFD.py --nZone 2 --fsi True --parallel -f"
    pywrapper_fsi2d.timeout   = 1600
    pywrapper_fsi2d.tol       = 0.00001
    test_list.append(pywrapper_fsi2d)

    # Unsteady CHT
    pywrapper_unsteadyCHT               = TestCase('pywrapper_unsteadyCHT')
    pywrapper_unsteadyCHT.cfg_dir       = "py_wrapper/flatPlate_unsteady_CHT"
    pywrapper_unsteadyCHT.cfg_file      = "unsteady_CHT_FlatPlate_Conf.cfg"
    pywrapper_unsteadyCHT.test_iter     = 5
    pywrapper_unsteadyCHT.test_vals     = [-1.598116, 2.263342, 0.001088, 0.145901] #last 4 columns
    pywrapper_unsteadyCHT.su2_exec      = "mpirun -np 2 python launch_unsteady_CHT_FlatPlate.py --parallel -f"
    pywrapper_unsteadyCHT.timeout       = 1600
    pywrapper_unsteadyCHT.tol           = 0.00001
    pywrapper_unsteadyCHT.unsteady      = True
    test_list.append(pywrapper_unsteadyCHT)

    # Rigid motion
    pywrapper_rigidMotion               = TestCase('pywrapper_rigidMotion')
    pywrapper_rigidMotion.cfg_dir       = "py_wrapper/flatPlate_rigidMotion"
    pywrapper_rigidMotion.cfg_file      = "flatPlate_rigidMotion_Conf.cfg"
    pywrapper_rigidMotion.test_iter     = 5
    pywrapper_rigidMotion.test_vals     = [-1.598116, 2.259705, -0.040360, 0.144232] #last 4 columns
    pywrapper_rigidMotion.su2_exec      = "mpirun -np 2 python launch_flatPlate_rigidMotion.py --parallel -f"
    pywrapper_rigidMotion.timeout       = 1600
    pywrapper_rigidMotion.tol           = 0.00001
    pywrapper_rigidMotion.unsteady      = True
    test_list.append(pywrapper_rigidMotion)
    
    ######################################
    ### RUN TUTORIAL CASES             ###
    ######################################
    
    # Inviscid Bump
    tutorial_inv_bump            = TestCase('inviscid_bump_tutorial')
    tutorial_inv_bump.cfg_dir    = "../Tutorials/Inviscid_Bump"
    tutorial_inv_bump.cfg_file   = "inv_channel.cfg"
    tutorial_inv_bump.test_iter  = 0
    tutorial_inv_bump.test_vals  = [-1.437425, 4.075857, -0.261018, 0.060652] #last 4 columns
    tutorial_inv_bump.su2_exec   = "mpirun -np 2 SU2_CFD"
    tutorial_inv_bump.timeout    = 1600
    tutorial_inv_bump.tol        = 0.00001
    tutorial_inv_bump.no_restart = True
    test_list.append(tutorial_inv_bump)
    
    # Inviscid Wedge
    tutorial_inv_wedge            = TestCase('inviscid_wedge_tutorial')
    tutorial_inv_wedge.cfg_dir    = "../Tutorials/Inviscid_Wedge"
    tutorial_inv_wedge.cfg_file   = "inv_wedge_HLLC.cfg"
    tutorial_inv_wedge.test_iter  = 0
    tutorial_inv_wedge.test_vals  = [-0.481460, 5.253008, -0.241559, 0.042452] #last 4 columns
    tutorial_inv_wedge.su2_exec   = "mpirun -np 2 SU2_CFD"
    tutorial_inv_wedge.timeout    = 1600
    tutorial_inv_wedge.tol        = 0.00001
    tutorial_inv_wedge.no_restart = True
    test_list.append(tutorial_inv_wedge)
    
    # Inviscid ONERA M6
    tutorial_inv_onera            = TestCase('inviscid_onera_tutorial')
    tutorial_inv_onera.cfg_dir    = "../Tutorials/Inviscid_ONERAM6"
    tutorial_inv_onera.cfg_file   = "inv_ONERAM6.cfg"
    tutorial_inv_onera.test_iter  = 0
    tutorial_inv_onera.test_vals  = [-5.204928, -4.597762, 0.167072, 0.053113] #last 4 columns
    tutorial_inv_onera.su2_exec   = "mpirun -np 2 SU2_CFD"
    tutorial_inv_onera.timeout    = 1600
    tutorial_inv_onera.tol        = 0.00001
    tutorial_inv_onera.no_restart = True
    test_list.append(tutorial_inv_onera)
    
    # Laminar Cylinder
    tutorial_lam_cylinder            = TestCase('laminar_cylinder_tutorial')
    tutorial_lam_cylinder.cfg_dir    = "../Tutorials/Laminar_Cylinder"
    tutorial_lam_cylinder.cfg_file   = "lam_cylinder.cfg"
    tutorial_lam_cylinder.test_iter  = 0
    tutorial_lam_cylinder.test_vals  = [-6.162141, -0.699617, 0.130622, 60.172836] #last 4 columns
    tutorial_lam_cylinder.su2_exec   = "mpirun -np 2 SU2_CFD"
    tutorial_lam_cylinder.timeout    = 1600
    tutorial_lam_cylinder.tol        = 0.00001
    tutorial_lam_cylinder.no_restart = True
    test_list.append(tutorial_lam_cylinder)

    # Laminar Flat Plate
    tutorial_lam_flatplate            = TestCase('laminar_flatplate_tutorial')
    tutorial_lam_flatplate.cfg_dir    = "../Tutorials/Laminar_Flat_Plate"
    tutorial_lam_flatplate.cfg_file   = "lam_flatplate.cfg"
    tutorial_lam_flatplate.test_iter  = 0
    tutorial_lam_flatplate.test_vals  = [-2.821818, 2.657591, -0.683921, 0.028640] #last 4 columns
    tutorial_lam_flatplate.su2_exec   = "mpirun -np 2 SU2_CFD"
    tutorial_lam_flatplate.timeout    = 1600
    tutorial_lam_flatplate.tol        = 0.00001
    tutorial_lam_flatplate.no_restart = True
    test_list.append(tutorial_lam_flatplate)
    
    # Turbulent Flat Plate
    tutorial_turb_flatplate            = TestCase('turbulent_flatplate_tutorial')
    tutorial_turb_flatplate.cfg_dir    = "../Tutorials/Turbulent_Flat_Plate"
    tutorial_turb_flatplate.cfg_file   = "turb_SA_flatplate.cfg"
    tutorial_turb_flatplate.test_iter  = 0
    tutorial_turb_flatplate.test_vals  = [-2.258584, -4.899476, -0.786721, 0.200332] #last 4 columns
    tutorial_turb_flatplate.su2_exec   = "mpirun -np 2 SU2_CFD"
    tutorial_turb_flatplate.timeout    = 1600
    tutorial_turb_flatplate.tol        = 0.00001
    tutorial_turb_flatplate.no_restart = True
    test_list.append(tutorial_turb_flatplate)
    
    # Transitional FlatPlate
    tutorial_trans_flatplate            = TestCase('transitional_flatplate_tutorial')
    tutorial_trans_flatplate.cfg_dir    = "../Tutorials/Transitional_Flat_Plate"
    tutorial_trans_flatplate.cfg_file   = "transitional_BC_model_ConfigFile.cfg"
    tutorial_trans_flatplate.test_iter  = 0
    tutorial_trans_flatplate.test_vals  = [-22.021786, -15.330906, 0.000000, 0.023952] #last 4 columns
    tutorial_trans_flatplate.su2_exec   = "mpirun -np 2 SU2_CFD"
    tutorial_trans_flatplate.timeout    = 1600
    tutorial_trans_flatplate.tol        = 0.00001
    tutorial_trans_flatplate.no_restart = True
    test_list.append(tutorial_trans_flatplate)

    # Turbulent ONERA M6
    tutorial_turb_oneram6            = TestCase('turbulent_oneram6_tutorial')
    tutorial_turb_oneram6.cfg_dir    = "../Tutorials/Turbulent_ONERAM6"
    tutorial_turb_oneram6.cfg_file   = "turb_ONERAM6.cfg"
    tutorial_turb_oneram6.test_iter  = 0
    tutorial_turb_oneram6.test_vals  = [-4.499497, -11.517143, 0.391850, 0.343722] #last 4 columns
    tutorial_turb_oneram6.su2_exec   = "mpirun -np 2 SU2_CFD"
    tutorial_turb_oneram6.timeout    = 1600
    tutorial_turb_oneram6.tol        = 0.00001
    tutorial_turb_oneram6.no_restart = True
    test_list.append(tutorial_turb_oneram6)

    # Inviscid NACA 0012 Design
    tutorial_design_inv_naca0012            = TestCase('design_inv_naca0012')
    tutorial_design_inv_naca0012.cfg_dir    = "../Tutorials/Inviscid_2D_Unconstrained_NACA0012"
    tutorial_design_inv_naca0012.cfg_file   = "inv_NACA0012_basic.cfg"
    tutorial_design_inv_naca0012.test_iter  = 0
    tutorial_design_inv_naca0012.test_vals  = [-3.585391, -2.989014, 0.095708, 0.174500] #last 4 columns
    tutorial_design_inv_naca0012.su2_exec   = "mpirun -np 2 SU2_CFD"
    tutorial_design_inv_naca0012.timeout    = 1600
    tutorial_design_inv_naca0012.tol        = 0.00001
    tutorial_design_inv_naca0012.no_restart = True
    test_list.append(tutorial_design_inv_naca0012)

    # Turbulent RAE 2822 Design
    tutorial_design_turb_rae2822            = TestCase('design_turb_rae2822')
    tutorial_design_turb_rae2822.cfg_dir    = "../Tutorials/Turbulent_2D_Constrained_RAE2822"
    tutorial_design_turb_rae2822.cfg_file   = "turb_SA_RAE2822.cfg"
    tutorial_design_turb_rae2822.test_iter  = 0
    tutorial_design_turb_rae2822.test_vals  = [-1.700114, -4.931315, 0.293884, 0.331019] #last 4 columns
    tutorial_design_turb_rae2822.su2_exec   = "mpirun -np 2 SU2_CFD"
    tutorial_design_turb_rae2822.timeout    = 1600
    tutorial_design_turb_rae2822.tol        = 0.00001
    tutorial_design_turb_rae2822.no_restart = True
    test_list.append(tutorial_design_turb_rae2822)

    # Multi Objective Design
    tutorial_design_multiobj            = TestCase('design_multiobj')
    tutorial_design_multiobj.cfg_dir    = "../Tutorials/Multi_Objective_Shape_Design"
    tutorial_design_multiobj.cfg_file   = "inv_wedge_ROE_multiobj_combo.cfg"
    tutorial_design_multiobj.test_iter  = 0
    tutorial_design_multiobj.test_vals  = [2.657333, -3.020635, 324840.000000, 0.000000] #last 4 columns
    tutorial_design_multiobj.su2_exec   = "mpirun -np 2 SU2_CFD"
    tutorial_design_multiobj.timeout    = 1600
    tutorial_design_multiobj.tol        = 0.00001
    tutorial_design_multiobj.no_restart = True
    test_list.append(tutorial_design_multiobj)

    ######################################
    ### RUN TESTS                      ###
    ######################################
    
    pass_list = [ test.run_test() for test in test_list ]


    ######################################
    ### RUN SU2_DEF TESTS              ###
    ######################################
    
    # Inviscid NACA0012 (triangles)
    naca0012_def            = TestCase('naca0012_def')
    naca0012_def.cfg_dir   = "deformation/naca0012"
    naca0012_def.cfg_file  = "def_NACA0012.cfg"
    naca0012_def.test_iter = 10
    naca0012_def.test_vals = [0.0035688] #residual
    naca0012_def.su2_exec  = "mpirun -n 2 SU2_DEF"
    naca0012_def.timeout   = 1600
    naca0012_def.tol       = 1e-8
    
    pass_list.append(naca0012_def.run_def())
    test_list.append(naca0012_def)

    # RAE2822 (mixed tris + quads)
    rae2822_def            = TestCase('rae2822_def')
    rae2822_def.cfg_dir   = "deformation/rae2822"
    rae2822_def.cfg_file  = "def_RAE2822.cfg"
    rae2822_def.test_iter = 10
    rae2822_def.test_vals = [8.24031e-09] #residual
    rae2822_def.su2_exec  = "mpirun -n 2 SU2_DEF"
    rae2822_def.timeout   = 1600
    rae2822_def.tol       = 1e-13
    
    pass_list.append(rae2822_def.run_def())
    test_list.append(rae2822_def)
    
    # Turb NACA4412 (quads, wall distance)
    naca4412_def            = TestCase('naca4412_def')
    naca4412_def.cfg_dir   = "deformation/naca4412"
    naca4412_def.cfg_file  = "def_NACA4412.cfg"
    naca4412_def.test_iter = 10
    naca4412_def.test_vals = [2.210380e-12] #residual
    naca4412_def.su2_exec  = "mpirun -n 2 SU2_DEF"
    naca4412_def.timeout   = 1600
    naca4412_def.tol       = 1e-12
    
    pass_list.append(naca4412_def.run_def())
    test_list.append(naca4412_def)
    
    # Brick of tets (inverse volume)
    brick_tets_def            = TestCase('brick_tets_def')
    brick_tets_def.cfg_dir   = "deformation/brick_tets"
    brick_tets_def.cfg_file  = "def_brick_tets.cfg"
    brick_tets_def.test_iter = 10
    brick_tets_def.test_vals = [9.584000e-04] #residual
    brick_tets_def.su2_exec  = "mpirun -n 2 SU2_DEF"
    brick_tets_def.timeout   = 1600
    brick_tets_def.tol       = 1e-9
    
    pass_list.append(brick_tets_def.run_def())
    test_list.append(brick_tets_def)
    
    # Brick of isotropic hexas (inverse volume)
    brick_hex_def           = TestCase('brick_hex_def')
    brick_hex_def.cfg_dir   = "deformation/brick_hex"
    brick_hex_def.cfg_file  = "def_brick_hex.cfg"
    brick_hex_def.test_iter = 10
    brick_hex_def.test_vals = [2.024230e-04] #residual
    brick_hex_def.su2_exec  = "mpirun -n 2 SU2_DEF"
    brick_hex_def.timeout   = 1600
    brick_hex_def.tol       = 1e-9
    
    pass_list.append(brick_hex_def.run_def())
    test_list.append(brick_hex_def)
    
    # Brick with a pyramid layer (inverse volume)
    brick_pyra_def           = TestCase('brick_pyra_def')
    brick_pyra_def.cfg_dir   = "deformation/brick_pyra"
    brick_pyra_def.cfg_file  = "def_brick_pyra.cfg"
    brick_pyra_def.test_iter = 10
    brick_pyra_def.test_vals = [0.00194185] #residual
    brick_pyra_def.su2_exec  = "mpirun -n 2 SU2_DEF"
    brick_pyra_def.timeout   = 1600
    brick_pyra_def.tol       = 1e-8
    
    pass_list.append(brick_pyra_def.run_def())
    test_list.append(brick_pyra_def)
    
    # Brick of isotropic prisms (inverse volume)
    brick_prism_def           = TestCase('brick_prism_def')
    brick_prism_def.cfg_dir   = "deformation/brick_prism"
    brick_prism_def.cfg_file  = "def_brick_prism.cfg"
    brick_prism_def.test_iter = 10
    brick_prism_def.test_vals = [6.210870e-03] #residual
    brick_prism_def.su2_exec  = "mpirun -n 2 SU2_DEF"
    brick_prism_def.timeout   = 1600
    brick_prism_def.tol       = 1e-8
    
    pass_list.append(brick_prism_def.run_def())
    test_list.append(brick_prism_def)
    
    # Brick of prisms with high aspect ratio cells near the wall (wall distance)
    brick_prism_rans_def           = TestCase('brick_prism_rans_def')
    brick_prism_rans_def.cfg_dir   = "deformation/brick_prism_rans"
    brick_prism_rans_def.cfg_file  = "def_brick_prism_rans.cfg"
    brick_prism_rans_def.test_iter = 10
    brick_prism_rans_def.test_vals = [2.788570e-07] #residual
    brick_prism_rans_def.su2_exec  = "mpirun -n 2 SU2_DEF"
    brick_prism_rans_def.timeout   = 1600
    brick_prism_rans_def.tol       = 1e-12
    
    pass_list.append(brick_prism_rans_def.run_def())
    test_list.append(brick_prism_rans_def)
    
    # Brick of hexas with high aspect ratio cells near the wall (inverse volume)
    brick_hex_rans_def           = TestCase('brick_hex_rans_def')
    brick_hex_rans_def.cfg_dir   = "deformation/brick_hex_rans"
    brick_hex_rans_def.cfg_file  = "def_brick_hex_rans.cfg"
    brick_hex_rans_def.test_iter = 10
    brick_hex_rans_def.test_vals = [3.566260e-06] #residual
    brick_hex_rans_def.su2_exec  = "mpirun -n 2 SU2_DEF"
    brick_hex_rans_def.timeout   = 1600
    brick_hex_rans_def.tol       = 1e-11
    
    pass_list.append(brick_hex_rans_def.run_def())
    test_list.append(brick_hex_rans_def)

    # Cylindrical FFD test
    cylinder_ffd_def           = TestCase('cylinder_ffd_def')
    cylinder_ffd_def.cfg_dir   = "deformation/cylindrical_ffd"
    cylinder_ffd_def.cfg_file  = "def_cylindrical.cfg"
    cylinder_ffd_def.test_iter = 10
    cylinder_ffd_def.test_vals = [0.000901831] #residual
    cylinder_ffd_def.su2_exec  = "mpirun -n 2 SU2_DEF"
    cylinder_ffd_def.timeout   = 1600
    cylinder_ffd_def.tol       = 1e-9

    pass_list.append(cylinder_ffd_def.run_def())
    test_list.append(cylinder_ffd_def)

    # Spherical FFD test
    sphere_ffd_def           = TestCase('sphere_ffd_def')
    sphere_ffd_def.cfg_dir   = "deformation/spherical_ffd"
    sphere_ffd_def.cfg_file  = "def_spherical.cfg"
    sphere_ffd_def.test_iter = 10
    sphere_ffd_def.test_vals = [0.00360079] #residual
    sphere_ffd_def.su2_exec  = "mpirun -n 2 SU2_DEF"
    sphere_ffd_def.timeout   = 1600
    sphere_ffd_def.tol       = 1e-8

    pass_list.append(sphere_ffd_def.run_def())
    test_list.append(sphere_ffd_def)

    # Spherical FFD test using BSplines
    sphere_ffd_def_bspline           = TestCase('sphere_ffd_def_bspline')
    sphere_ffd_def_bspline.cfg_dir   = "deformation/spherical_ffd"
    sphere_ffd_def_bspline.cfg_file  = "def_spherical_bspline.cfg"
    sphere_ffd_def_bspline.test_iter = 10
    sphere_ffd_def_bspline.test_vals = [0.00208253] #residual
    sphere_ffd_def_bspline.su2_exec  = "mpirun -n 2 SU2_DEF"
    sphere_ffd_def_bspline.timeout   = 1600
    sphere_ffd_def_bspline.tol       = 1e-8

    pass_list.append(sphere_ffd_def_bspline.run_def())
    test_list.append(sphere_ffd_def_bspline)


    # Tests summary
    print('==================================================================')
    print('Summary of the parallel tests')
    print('python version:', sys.version)
    for i, test in enumerate(test_list):
        if (pass_list[i]):
            print('  passed - %s'%test.tag)
        else:
            print('* FAILED - %s'%test.tag)

    if all(pass_list):
        sys.exit(0)
    else:
        sys.exit(1)
    # done

if __name__ == '__main__':
    main()
