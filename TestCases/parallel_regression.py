#!/usr/bin/env python

## \file parallel_regression.py
#  \brief Python script for automated regression testing of SU2 examples
#  \author A. Aranake, A. Campos, T. Economon, T. Lukaczyk, S. Padron
#  \version 5.0.0 "Raven"
#
# SU2 Original Developers: Dr. Francisco D. Palacios.
#                          Dr. Thomas D. Economon.
#
# SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
#                 Prof. Piero Colonna's group at Delft University of Technology.
#                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
#                 Prof. Rafael Palacios' group at Imperial College London.
#                 Prof. Edwin van der Weide's group at the University of Twente.
#                 Prof. Vincent Terrapon's group at the University of Liege.
#
# Copyright (C) 2012-2017 SU2, the open-source CFD code.
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
    channel.test_iter = 100
    channel.test_vals = [-3.071007, 2.301191, 0.008562, 0.028922] #last 4 columns
    channel.su2_exec  = "parallel_computation.py -f"
    channel.timeout   = 1600
    channel.tol       = 0.00001
    test_list.append(channel)

    # NACA0012 
    naca0012           = TestCase('naca0012')
    naca0012.cfg_dir   = "euler/naca0012"
    naca0012.cfg_file  = "inv_NACA0012_Roe.cfg"
    naca0012.test_iter = 100
    naca0012.test_vals = [-6.237188, -5.641250, 0.334843, 0.022206] #last 4 columns
    naca0012.su2_exec  = "parallel_computation.py -f"
    naca0012.timeout   = 1600
    naca0012.tol       = 0.00001
    test_list.append(naca0012)

    # Supersonic wedge 
    wedge           = TestCase('wedge')
    wedge.cfg_dir   = "euler/wedge"
    wedge.cfg_file  = "inv_wedge_HLLC.cfg"
    wedge.test_iter = 100
    wedge.test_vals = [-1.690232, 3.924432, -0.252221, 0.044419] #last 4 columns
    wedge.su2_exec  = "parallel_computation.py -f"
    wedge.timeout   = 1600
    wedge.tol       = 0.00001
    test_list.append(wedge)

    # ONERA M6 Wing
    oneram6           = TestCase('oneram6')
    oneram6.cfg_dir   = "euler/oneram6"
    oneram6.cfg_file  = "inv_ONERAM6.cfg"
    oneram6.test_iter = 10
    oneram6.test_vals = [-13.400678, -12.932056, 0.282557, 0.012706] #last 4 columns
    oneram6.su2_exec  = "parallel_computation.py -f"
    oneram6.timeout   = 3200
    oneram6.tol       = 0.00001
    test_list.append(oneram6)

    # Fixed CL NACA0012
    fixedCL_naca0012           = TestCase('fixedcl_naca0012')
    fixedCL_naca0012.cfg_dir   = "fixed_cl/naca0012"
    fixedCL_naca0012.cfg_file  = "inv_NACA0012.cfg"
    fixedCL_naca0012.test_iter = 500
    fixedCL_naca0012.test_vals = [-5.930145,-0.546310,0.300089,0.019393] #last 4 columns
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
    polar_naca0012.test_vals = [-1.298373, 4.139860, 0.010703, 0.008932] #last 4 columns
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
    flatplate.test_iter = 100
    flatplate.test_vals = [-5.231727, 0.261637, -0.166869, 0.012707] #last 4 columns
    flatplate.su2_exec  = "parallel_computation.py -f"
    flatplate.timeout   = 1600
    flatplate.tol       = 0.00001
    test_list.append(flatplate)

    # Laminar cylinder (steady)
    cylinder           = TestCase('cylinder')
    cylinder.cfg_dir   = "navierstokes/cylinder"
    cylinder.cfg_file  = "lam_cylinder.cfg"
    cylinder.test_iter = 25
    cylinder.test_vals = [-6.757291, -1.289309, -0.125948, 0.625438] #last 4 columns
    cylinder.su2_exec  = "parallel_computation.py -f"
    cylinder.timeout   = 1600
    cylinder.tol       = 0.00001
    test_list.append(cylinder)

    # Laminar cylinder (low Mach correction)
    cylinder_lowmach           = TestCase('cylinder_lowmach')
    cylinder_lowmach.cfg_dir   = "navierstokes/cylinder"
    cylinder_lowmach.cfg_file  = "cylinder_lowmach.cfg"
    cylinder_lowmach.test_iter = 25
    cylinder_lowmach.test_vals = [-6.861860, -1.399846, -1.557250, 110.230719] #last 4 columns
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
    poiseuille.tol       = 0.00001
    test_list.append(poiseuille)

    ##########################
    ### Compressible RANS  ###
    ##########################

    # RAE2822 SA
    rae2822_sa           = TestCase('rae2822_sa')
    rae2822_sa.cfg_dir   = "rans/rae2822"
    rae2822_sa.cfg_file  = "turb_SA_RAE2822.cfg"
    rae2822_sa.test_iter = 100
    rae2822_sa.test_vals = [-3.647743, -5.442583, 0.903159, 0.025800] #last 4 columns
    rae2822_sa.su2_exec  = "parallel_computation.py -f"
    rae2822_sa.timeout   = 1600
    rae2822_sa.tol       = 0.00001
    test_list.append(rae2822_sa)
    
    # RAE2822 SST
    rae2822_sst           = TestCase('rae2822_sst')
    rae2822_sst.cfg_dir   = "rans/rae2822"
    rae2822_sst.cfg_file  = "turb_SST_RAE2822.cfg"
    rae2822_sst.test_iter = 100
    rae2822_sst.test_vals = [-1.198777, 4.013989, 0.904265, 0.025941] #last 4 columns
    rae2822_sst.su2_exec  = "parallel_computation.py -f"
    rae2822_sst.timeout   = 1600
    rae2822_sst.tol       = 0.00001
    test_list.append(rae2822_sst)

    # Flat plate
    turb_flatplate           = TestCase('turb_flatplate')
    turb_flatplate.cfg_dir   = "rans/flatplate"
    turb_flatplate.cfg_file  = "turb_SA_flatplate.cfg"
    turb_flatplate.test_iter = 100
    turb_flatplate.test_vals = [-5.071289, -7.366446, -0.187145, 0.010957] #last 4 columns
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
    turb_naca0012_sa_bin.test_vals = [-12.000763, -9.145363, 1.070528, 0.019417] #last 4 columns
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
    propeller.test_vals = [-3.380312, -8.129519, 0.000043, 0.056159] #last 4 columns
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
    turb_naca0012_sst_restart_mg.test_iter = 600
    turb_naca0012_sst_restart_mg.ntest_vals = 5
    turb_naca0012_sst_restart_mg.test_vals = [-6.482396, -4.653908, 1.155000, -0.006296, 0.078747] #last 5 columns
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
    inc_euler_naca0012.test_vals = [-3.544713,-3.135163,0.968231,0.010161] #last 4 columns
    inc_euler_naca0012.su2_exec  = "parallel_computation.py -f"
    inc_euler_naca0012.timeout   = 1600
    inc_euler_naca0012.tol       = 0.00001
    test_list.append(inc_euler_naca0012)

    ############################
    ### Incompressible RANS  ###
    ############################

    # NACA0012
    inc_turb_naca0012           = TestCase('inc_turb_naca0012')
    inc_turb_naca0012.cfg_dir   = "incomp_rans/naca0012"
    inc_turb_naca0012.cfg_file  = "naca0012.cfg"
    inc_turb_naca0012.test_iter = 20
    inc_turb_naca0012.test_vals = [-4.710048, -11.007498, 0.000002, 0.210441] #last 4 columns
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
    schubauer_klebanoff_transition.test_iter    = 250
    schubauer_klebanoff_transition.test_vals    = [-9.467876, -15.630835, 0.000520, 0.001062] #last 4 columns
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
    contadj_naca0012.test_vals = [-9.783712, -15.192812, 0.300920, 0.019552] #last 4 columns
    contadj_naca0012.su2_exec  = "parallel_computation.py -f"
    contadj_naca0012.timeout   = 1600
    contadj_naca0012.tol       = 0.00001
    test_list.append(contadj_naca0012)

    # Inviscid ONERA M6
    contadj_oneram6           = TestCase('contadj_oneram6')
    contadj_oneram6.cfg_dir   = "cont_adj_euler/oneram6"
    contadj_oneram6.cfg_file  = "inv_ONERAM6.cfg"
    contadj_oneram6.test_iter = 10
    contadj_oneram6.test_vals = [-12.131587, -12.703243, 0.685900, 0.007594] #last 4 columns
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
    contadj_fixed_CL_naca0012.test_iter = 500
    contadj_fixed_CL_naca0012.test_vals = [-2.921695, -8.474107, 0.250530, 0.000003] #last 4 columns
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
    contadj_ns_cylinder.test_iter = 100
    contadj_ns_cylinder.test_vals = [-3.659040, -9.117773, 2.056700, -0.000000] #last 4 columns
    contadj_ns_cylinder.su2_exec  = "parallel_computation.py -f"
    contadj_ns_cylinder.timeout   = 1600
    contadj_ns_cylinder.tol       = 0.00001
    test_list.append(contadj_ns_cylinder)

    # Adjoint laminar naca0012 subsonic
    contadj_ns_naca0012_sub           = TestCase('contadj_ns_naca0012_sub')
    contadj_ns_naca0012_sub.cfg_dir   = "cont_adj_navierstokes/naca0012_sub"
    contadj_ns_naca0012_sub.cfg_file  = "lam_NACA0012.cfg"
    contadj_ns_naca0012_sub.test_iter = 100
    contadj_ns_naca0012_sub.test_vals = [-2.744551, -8.216469, 0.518810, 0.001210] #last 4 columns
    contadj_ns_naca0012_sub.su2_exec  = "parallel_computation.py -f"
    contadj_ns_naca0012_sub.timeout   = 1600
    contadj_ns_naca0012_sub.tol       = 0.00001
    test_list.append(contadj_ns_naca0012_sub)
    
    # Adjoint laminar naca0012 transonic
    contadj_ns_naca0012_trans           = TestCase('contadj_ns_naca0012_trans')
    contadj_ns_naca0012_trans.cfg_dir   = "cont_adj_navierstokes/naca0012_trans"
    contadj_ns_naca0012_trans.cfg_file  = "lam_NACA0012.cfg"
    contadj_ns_naca0012_trans.test_iter = 100
    contadj_ns_naca0012_trans.test_vals = [-1.041539, -6.578524, 1.772300, 0.012494] #last 4 columns
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
    contadj_rans_naca0012.test_iter = 100
    contadj_rans_naca0012.test_vals = [-0.814757, -5.726517, 19.169000, -0.000000] #last 4 columns
    contadj_rans_naca0012.su2_exec  = "parallel_computation.py -f"
    contadj_rans_naca0012.timeout   = 1600
    contadj_rans_naca0012.tol       = 0.00001
    test_list.append(contadj_rans_naca0012)
   
    # Adjoint turbulent NACA0012 with binary restarts
    contadj_rans_naca0012_bin           = TestCase('contadj_rans_naca0012_bin')
    contadj_rans_naca0012_bin.cfg_dir   = "cont_adj_rans/naca0012"
    contadj_rans_naca0012_bin.cfg_file  = "turb_nasa_binary.cfg"
    contadj_rans_naca0012_bin.test_iter = 100
    contadj_rans_naca0012_bin.test_vals = [-0.814736,-5.726472,1.9169e+01,-4.6176e-07] #last 4 columns
    contadj_rans_naca0012_bin.su2_exec  = "parallel_computation.py -f"
    contadj_rans_naca0012_bin.timeout   = 1600
    contadj_rans_naca0012_bin.tol       = 0.00001
    test_list.append(contadj_rans_naca0012_bin)
 
    # Adjoint turbulent RAE2822
    contadj_rans_rae2822           = TestCase('contadj_rans_rae822')
    contadj_rans_rae2822.cfg_dir   = "cont_adj_rans/rae2822"
    contadj_rans_rae2822.cfg_file  = "turb_SA_RAE2822.cfg"
    contadj_rans_rae2822.test_iter = 100
    contadj_rans_rae2822.test_vals = [-5.380910, -10.886057, -0.212470, 0.005448] #last 4 columns
    contadj_rans_rae2822.su2_exec  = "parallel_computation.py -f"
    contadj_rans_rae2822.timeout   = 1600
    contadj_rans_rae2822.tol       = 0.00001
    test_list.append(contadj_rans_rae2822)

    #######################################
    ### Cont. adj. incompressible Euler ###
    #######################################

    # Adjoint Incompressible Inviscid NACA0012
    contadj_incomp_NACA0012           = TestCase('contadj_incomp_NACA0012')
    contadj_incomp_NACA0012.cfg_dir   = "cont_adj_incomp_euler/naca0012"
    contadj_incomp_NACA0012.cfg_file  = "incomp_NACA0012.cfg"
    contadj_incomp_NACA0012.test_iter = 5
    contadj_incomp_NACA0012.test_vals = [-11.968536, -12.133235, 1.939900, 0.000000] #last 4 columns
    contadj_incomp_NACA0012.su2_exec  = "parallel_computation.py -f"
    contadj_incomp_NACA0012.timeout   = 1600
    contadj_incomp_NACA0012.tol       = 0.00001
    test_list.append(contadj_incomp_NACA0012)

    #####################################
    ### Cont. adj. incompressible N-S ###
    #####################################

    # Adjoint Incompressible Viscous Cylinder
    contadj_incomp_cylinder           = TestCase('contadj_incomp_cylinder')
    contadj_incomp_cylinder.cfg_dir   = "cont_adj_incomp_navierstokes/cylinder"
    contadj_incomp_cylinder.cfg_file  = "lam_incomp_cylinder.cfg"
    contadj_incomp_cylinder.test_iter = 25
    contadj_incomp_cylinder.test_vals = [-5.718840, -7.012324, 2.932100, 0.000000] #last 4 columns
    contadj_incomp_cylinder.su2_exec  = "parallel_computation.py -f"
    contadj_incomp_cylinder.timeout   = 1600
    contadj_incomp_cylinder.tol       = 0.00001
    test_list.append(contadj_incomp_cylinder)

    ######################################
    ### Moving Wall                    ###
    ######################################

    # Lid-driven cavity
    cavity           = TestCase('cavity')
    cavity.cfg_dir   = "moving_wall/cavity"
    cavity.cfg_file  = "lam_cavity.cfg"
    cavity.test_iter = 25
    cavity.test_vals = [-5.597018, -0.133158, 0.169034, 0.821371] #last 4 columns
    cavity.su2_exec  = "parallel_computation.py -f"
    cavity.timeout   = 1600
    cavity.tol       = 0.00001
    test_list.append(cavity)

    # Spinning cylinder
    spinning_cylinder           = TestCase('spinning_cylinder')
    spinning_cylinder.cfg_dir   = "moving_wall/spinning_cylinder"
    spinning_cylinder.cfg_file  = "spinning_cylinder.cfg"
    spinning_cylinder.test_iter = 25
    spinning_cylinder.test_vals = [-7.648320, -2.202753, 1.236836, 1.609040] #last 4 columns
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
    square_cylinder.test_vals = [-1.166422,0.076751,1.398549,2.197047] #last 4 columns
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
    sine_gust.test_vals = [-1.977531, 3.481790, -0.014552, -0.004969] #last 4 columns
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

    ######################################
    ### NICFD                          ###
    ######################################	

    # Rarefaction shock wave edge_VW
    edge_VW           = TestCase('edge_VW')
    edge_VW.cfg_dir   = "nicf/edge"
    edge_VW.cfg_file  = "edge_VW.cfg"
    edge_VW.test_iter = 100
    edge_VW.test_vals = [-5.059867, 1.114232, -0.000009, 0.000000] #last 4 columns
    edge_VW.su2_exec  = "parallel_computation.py -f"
    edge_VW.timeout   = 1600
    edge_VW.tol       = 0.00001
    test_list.append(edge_VW)

    # Rarefaction shock wave edge_PPR
    edge_PPR           = TestCase('edge_PPR')
    edge_PPR.cfg_dir   = "nicf/edge"
    edge_PPR.cfg_file  = "edge_PPR.cfg"
    edge_PPR.test_iter = 100
    edge_PPR.test_vals = [-5.487473, 0.653442, -0.000037, 0.000000] #last 4 columns
    edge_PPR.su2_exec  = "parallel_computation.py -f"
    edge_PPR.timeout   = 1600
    edge_PPR.tol       = 0.00001
    test_list.append(edge_PPR)
    
    ##########################
    ###   Python wrapper   ###
    ##########################
    
    # NACA0012 
    pywrapper_naca0012           = TestCase('pywrapper_naca0012')
    pywrapper_naca0012.cfg_dir   = "euler/naca0012"
    pywrapper_naca0012.cfg_file  = "inv_NACA0012_Roe.cfg"
    pywrapper_naca0012.test_iter = 100
    pywrapper_naca0012.test_vals = [-6.237188, -5.641250, 0.334843, 0.022206] #last 4 columns
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
    pywrapper_square_cylinder.test_vals = [-1.166422,0.076751,1.398549,2.197047] #last 4 columns
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

    
    ######################################
    ### RUN TESTS                      ###
    ######################################
    
    pass_list = [ test.run_test() for test in test_list ]


    # Tests summary
    print '=================================================================='
    print 'Summary of the parallel tests'
    for i, test in enumerate(test_list):
        if (pass_list[i]):
            print '  passed - %s'%test.tag
        else:
            print '* FAILED - %s'%test.tag

    if all(pass_list):
        sys.exit(0)
    else:
        sys.exit(1)
    # done

if __name__ == '__main__':
    main()
