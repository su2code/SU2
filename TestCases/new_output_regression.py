#!/usr/bin/env python

## \file serial_regression.py
#  \brief Python script for automated regression testing of SU2 examples
#  \author A. Aranake, A. Campos, T. Economon, T. Lukaczyk, S. Padron
#  \version 6.1.0 "Falcon"
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

from __future__ import print_function, division, absolute_import
import sys
from TestCase import TestCase

def main():
    '''This program runs SU2 and ensures that the output matches specified values. 
       This will be used to do checks when code is pushed to github 
       to make sure nothing is broken. '''

    test_list = []

    #########################
    ## Compressible Euler ###
    #########################

    # Channel
    channel           = TestCase('channel')
    channel.cfg_dir   = "euler/channel"
    channel.cfg_file  = "inv_channel_RK.cfg"
    channel.test_iter = 10
    channel.test_vals = [-2.454049, 3.065639, -0.200679, 0.036298] #last 4 columns
    channel.su2_exec  = "SU2_CFD"
    channel.timeout   = 1600
    channel.new_output = True
    channel.tol       = 0.00001
    test_list.append(channel)

    # NACA0012 
    naca0012           = TestCase('naca0012')
    naca0012.cfg_dir   = "euler/naca0012"
    naca0012.cfg_file  = "inv_NACA0012_Roe.cfg"
    naca0012.test_iter = 20
    naca0012.test_vals = [-4.047448, -3.538057, 0.338691, 0.023131] #last 4 columns
    naca0012.su2_exec  = "SU2_CFD"
    naca0012.timeout   = 1600
    naca0012.new_output= True
    naca0012.tol       = 0.00001
    test_list.append(naca0012)

    # Supersonic wedge 
    wedge           = TestCase('wedge')
    wedge.cfg_dir   = "euler/wedge"
    wedge.cfg_file  = "inv_wedge_HLLC.cfg"
    wedge.test_iter = 20
    wedge.test_vals = [-0.804690, 4.936631, -0.251188, 0.044255] #last 4 columns
    wedge.su2_exec  = "SU2_CFD"
    wedge.timeout   = 1600
    wedge.new_output= True
    wedge.tol       = 0.00001
    test_list.append(wedge)

    # ONERA M6 Wing
    oneram6           = TestCase('oneram6')
    oneram6.cfg_dir   = "euler/oneram6"
    oneram6.cfg_file  = "inv_ONERAM6.cfg"
    oneram6.test_iter = 10
    oneram6.test_vals = [-10.384532, -9.835738, 0.282580, 0.012694] #last 4 columns
    oneram6.su2_exec  = "SU2_CFD"
    oneram6.timeout   = 9600
    oneram6.new_output = True
    oneram6.tol       = 0.00001
    test_list.append(oneram6)
    
    ##########################
    ###  Compressible N-S  ###
    ##########################

    # Laminar flat plate
    flatplate           = TestCase('flatplate')
    flatplate.cfg_dir   = "navierstokes/flatplate"
    flatplate.cfg_file  = "lam_flatplate.cfg"
    flatplate.test_iter = 20
    flatplate.test_vals = [-4.680777, 0.781234, -0.135957, 0.022977] #last 4 columns
    flatplate.su2_exec  = "SU2_CFD"
    flatplate.new_output = True
    flatplate.timeout   = 1600
    flatplate.tol       = 0.00001
    test_list.append(flatplate)

    # Laminar cylinder (steady)
    cylinder           = TestCase('cylinder')
    cylinder.cfg_dir   = "navierstokes/cylinder"
    cylinder.cfg_file  = "lam_cylinder.cfg"
    cylinder.test_iter = 25
    cylinder.test_vals = [-6.765432, -1.297428, 0.019508, 0.310040] #last 4 columns
    cylinder.su2_exec  = "SU2_CFD"
    cylinder.new_output = True
    cylinder.timeout   = 1600
    cylinder.tol       = 0.00001
    test_list.append(cylinder)

    # Laminar cylinder (low Mach correction)
    cylinder_lowmach           = TestCase('cylinder_lowmach')
    cylinder_lowmach.cfg_dir   = "navierstokes/cylinder"
    cylinder_lowmach.cfg_file  = "cylinder_lowmach.cfg"
    cylinder_lowmach.test_iter = 25
    cylinder_lowmach.test_vals = [-6.850123, -1.388088, -0.056090, 108.140177] #last 4 columns
    cylinder_lowmach.su2_exec  = "SU2_CFD"
    cylinder_lowmach.new_output = True
    cylinder_lowmach.timeout   = 1600
    cylinder_lowmach.tol       = 0.00001
    test_list.append(cylinder_lowmach)

    # 2D Poiseuille flow (body force driven with periodic inlet / outlet)
   # poiseuille           = TestCase('poiseuille')
   # poiseuille.cfg_dir   = "navierstokes/poiseuille"
   # poiseuille.cfg_file  = "lam_poiseuille.cfg"
   # poiseuille.test_iter = 10
   # poiseuille.test_vals = [-12.272146, -3.335311, 0.000001, 2.351005] #last 4 columns
   # poiseuille.su2_exec  = "SU2_CFD"
   # poiseuille.new_output = True
   # poiseuille.timeout   = 1600
   # poiseuille.tol       = 0.00001
   # test_list.append(poiseuille)

    # 2D Poiseuille flow (inlet profile file)
    poiseuille_profile           = TestCase('poiseuille_profile')
    poiseuille_profile.cfg_dir   = "navierstokes/poiseuille"
    poiseuille_profile.cfg_file  = "profile_poiseuille.cfg"
    poiseuille_profile.test_iter = 10
    poiseuille_profile.test_vals = [-12.494724, -7.712336, -0.000000, 2.085796] #last 4 columns
    poiseuille_profile.su2_exec  = "SU2_CFD"
    poiseuille_profile.new_output = True
    poiseuille_profile.timeout   = 1600
    poiseuille_profile.tol       = 0.00001
    test_list.append(poiseuille_profile)
    
    ##########################
    ### Compressible RANS  ###
    ##########################

    # RAE2822 SA
    rae2822_sa           = TestCase('rae2822_sa')
    rae2822_sa.cfg_dir   = "rans/rae2822"
    rae2822_sa.cfg_file  = "turb_SA_RAE2822.cfg"
    rae2822_sa.test_iter = 20
    rae2822_sa.test_vals = [-2.000469, -5.228296, 0.820188, 0.052004] #last 4 columns
    rae2822_sa.su2_exec  = "SU2_CFD"
    rae2822_sa.timeout   = 1600
    rae2822_sa.new_output = True
    rae2822_sa.tol       = 0.00001
    test_list.append(rae2822_sa)
    
    # RAE2822 SST
    rae2822_sst           = TestCase('rae2822_sst')
    rae2822_sst.cfg_dir   = "rans/rae2822"
    rae2822_sst.cfg_file  = "turb_SST_RAE2822.cfg"
    rae2822_sst.test_iter = 20
    rae2822_sst.test_vals = [-0.510826, 4.909714, 0.825023, 0.052675] #last 4 columns
    rae2822_sst.su2_exec  = "SU2_CFD"
    rae2822_sst.new_output = True
    rae2822_sst.timeout   = 1600
    rae2822_sst.tol       = 0.00001
    test_list.append(rae2822_sst)

    # Flat plate
    turb_flatplate           = TestCase('turb_flatplate')
    turb_flatplate.cfg_dir   = "rans/flatplate"
    turb_flatplate.cfg_file  = "turb_SA_flatplate.cfg"
    turb_flatplate.test_iter = 20
    turb_flatplate.test_vals = [-4.157169, -6.737133, -0.176253, 0.057446] #last 4 columns
    turb_flatplate.su2_exec  = "SU2_CFD"
    turb_flatplate.new_output  = True
    turb_flatplate.timeout   = 1600
    turb_flatplate.tol       = 0.00001
    test_list.append(turb_flatplate)

    # ONERA M6 Wing
    turb_oneram6           = TestCase('turb_oneram6')
    turb_oneram6.cfg_dir   = "rans/oneram6"
    turb_oneram6.cfg_file  = "turb_ONERAM6.cfg"
    turb_oneram6.test_iter = 10
    turb_oneram6.test_vals = [-2.327431, -6.564331, 0.230257, 0.155839]#last 4 columns
    turb_oneram6.su2_exec  = "SU2_CFD"
    turb_oneram6.new_output = True
    turb_oneram6.timeout   = 3200
    turb_oneram6.tol       = 0.00001
    test_list.append(turb_oneram6)

    # NACA0012 (SA, FUN3D results for finest grid: CL=1.0983, CD=0.01242)
    turb_naca0012_sa           = TestCase('turb_naca0012_sa')
    turb_naca0012_sa.cfg_dir   = "rans/naca0012"
    turb_naca0012_sa.cfg_file  = "turb_NACA0012_sa.cfg"
    turb_naca0012_sa.test_iter = 10
    turb_naca0012_sa.test_vals = [-11.981166, -9.145363, 1.070528, 0.019417] #last 4 columns
    turb_naca0012_sa.su2_exec  = "SU2_CFD"
    turb_naca0012_sa.new_output = True
    turb_naca0012_sa.timeout   = 3200
    turb_naca0012_sa.tol       = 0.00001
    test_list.append(turb_naca0012_sa)
   
    # NACA0012 (SA, FUN3D results for finest grid: CL=1.0983, CD=0.01242) with binary restart
    turb_naca0012_sa_bin           = TestCase('turb_naca0012_sa_bin')
    turb_naca0012_sa_bin.cfg_dir   = "rans/naca0012"
    turb_naca0012_sa_bin.cfg_file  = "turb_NACA0012_sa_binary.cfg"
    turb_naca0012_sa_bin.test_iter = 10
    turb_naca0012_sa_bin.test_vals = [-11.981289, -9.145363, 1.070528, 0.019417] #last 4 columns
    turb_naca0012_sa_bin.su2_exec  = "SU2_CFD"
    turb_naca0012_sa_bin.new_output  = True
    turb_naca0012_sa_bin.timeout   = 3200
    turb_naca0012_sa_bin.tol       = 0.00001
    test_list.append(turb_naca0012_sa_bin)
 
    # NACA0012 (SST, FUN3D results for finest grid: CL=1.0840, CD=0.01253)
    turb_naca0012_sst           = TestCase('turb_naca0012_sst')
    turb_naca0012_sst.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst.cfg_file  = "turb_NACA0012_sst.cfg"
    turb_naca0012_sst.test_iter = 10
    turb_naca0012_sst.test_vals = [-12.445710, -6.918658, 1.059622, 0.019138] #last 4 columns
    turb_naca0012_sst.su2_exec  = "SU2_CFD"
    turb_naca0012_sst.new_output  = True
    turb_naca0012_sst.timeout   = 3200
    turb_naca0012_sst.tol       = 0.00001
    test_list.append(turb_naca0012_sst)

    # PROPELLER 
    propeller           = TestCase('propeller')
    propeller.cfg_dir   = "rans/propeller"
    propeller.cfg_file  = "propeller.cfg"
    propeller.test_iter = 10
    propeller.test_vals = [-3.378876, -8.396837, 0.000047, 0.055591] #last 4 columns
    propeller.su2_exec  = "SU2_CFD"
    propeller.new_output = True
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
    turb_naca0012_sst_restart_mg.test_iter = 50
    turb_naca0012_sst_restart_mg.ntest_vals = 5
    turb_naca0012_sst_restart_mg.test_vals = [-6.459444, -4.595710, 1.201844, -0.007146, 0.080517] #last 5 columns
    turb_naca0012_sst_restart_mg.su2_exec  = "SU2_CFD"
    turb_naca0012_sst_restart_mg.new_output  = True
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
    inc_euler_naca0012.test_vals = [-4.880274, -3.797906, 0.501143, 0.007051] #last 4 columns
    inc_euler_naca0012.su2_exec  = "SU2_CFD"
    inc_euler_naca0012.new_output = True
    inc_euler_naca0012.timeout   = 1600
    inc_euler_naca0012.tol       = 0.00001
    test_list.append(inc_euler_naca0012)

    # C-D nozzle with pressure inlet and mass flow outlet
    inc_nozzle           = TestCase('inc_nozzle')
    inc_nozzle.cfg_dir   = "incomp_euler/nozzle"
    inc_nozzle.cfg_file  = "inv_nozzle.cfg"
    inc_nozzle.test_iter = 20
    inc_nozzle.test_vals = [-5.799445, -4.785945, -0.000443, 0.124533] #last 4 columns
    inc_nozzle.su2_exec  = "SU2_CFD"
    inc_nozzle.new_output = True
    inc_nozzle.timeout   = 1600
    inc_nozzle.tol       = 0.00001
    test_list.append(inc_nozzle)

    #############################
    ### Incompressible N-S    ###
    #############################

    # Laminar cylinder
    inc_lam_cylinder          = TestCase('inc_lam_cylinder')
    inc_lam_cylinder.cfg_dir   = "incomp_navierstokes/cylinder"
    inc_lam_cylinder.cfg_file  = "incomp_cylinder.cfg"
    inc_lam_cylinder.test_iter = 10
    inc_lam_cylinder.test_vals = [-4.004277, -3.227956, 0.003852, 7.626578] #last 4 columns
    inc_lam_cylinder.new_output  = True
    inc_lam_cylinder.su2_exec  = "SU2_CFD"
    inc_lam_cylinder.timeout   = 1600
    inc_lam_cylinder.tol       = 0.00001
    test_list.append(inc_lam_cylinder)

    # Buoyancy-driven cavity
    inc_buoyancy          = TestCase('inc_buoyancy')
    inc_buoyancy.cfg_dir   = "incomp_navierstokes/buoyancy_cavity"
    inc_buoyancy.cfg_file  = "lam_buoyancy_cavity.cfg"
    inc_buoyancy.test_iter = 20
    inc_buoyancy.test_vals = [-4.436657, 0.507847, 0.000000, 0.000000] #last 4 columns
    inc_buoyancy.new_output  = True
    inc_buoyancy.su2_exec  = "SU2_CFD"
    inc_buoyancy.timeout   = 1600
    inc_buoyancy.tol       = 0.00001
    test_list.append(inc_buoyancy)

    # Laminar heated cylinder with polynomial fluid model
    inc_poly_cylinder          = TestCase('inc_poly_cylinder')
    inc_poly_cylinder.cfg_dir   = "incomp_navierstokes/cylinder"
    inc_poly_cylinder.cfg_file  = "poly_cylinder.cfg"
    inc_poly_cylinder.test_iter = 20
    inc_poly_cylinder.test_vals = [-8.108218, -2.158606, 0.019142, 1.902461] #last 4 columns
    inc_poly_cylinder.new_output  = True
    inc_poly_cylinder.su2_exec  = "SU2_CFD"
    inc_poly_cylinder.timeout   = 1600
    inc_poly_cylinder.tol       = 0.00001
    test_list.append(inc_poly_cylinder)

    ############################
    ### Incompressible RANS  ###
    ############################

    # NACA0012
    inc_turb_naca0012           = TestCase('inc_turb_naca0012')
    inc_turb_naca0012.cfg_dir   = "incomp_rans/naca0012"
    inc_turb_naca0012.cfg_file  = "naca0012.cfg"
    inc_turb_naca0012.test_iter = 20
    inc_turb_naca0012.test_vals = [-4.788495, -11.040511, 0.000023, 0.309503] #last 4 columns
    inc_turb_naca0012.new_output  = True
    inc_turb_naca0012.su2_exec  = "SU2_CFD"
    inc_turb_naca0012.timeout   = 1600
    inc_turb_naca0012.tol       = 0.00001
    test_list.append(inc_turb_naca0012)
    
    ####################
    ### DG-FEM Euler ###
    ####################
    
    # NACA0012
    fem_euler_naca0012           = TestCase('fem_euler_naca0012')
    fem_euler_naca0012.cfg_dir   = "hom_euler/NACA0012_5thOrder"
    fem_euler_naca0012.cfg_file  = "fem_NACA0012_reg.cfg"
    fem_euler_naca0012.test_iter = 10
    fem_euler_naca0012.test_vals = [-6.519946,-5.976944,0.255551,0.000028] #last 4 columns
    fem_euler_naca0012.su2_exec  = "SU2_CFD"
    fem_euler_naca0012.new_output = True
    fem_euler_naca0012.timeout   = 1600
    fem_euler_naca0012.tol       = 0.00001
    test_list.append(fem_euler_naca0012)
    
    ############################
    ### DG-FEM Navier-Stokes ###
    ############################
    
    # Flat plate
    fem_ns_flatplate           = TestCase('fem_ns_flatplate')
    fem_ns_flatplate.cfg_dir   = "hom_navierstokes/FlatPlate/nPoly4"
    fem_ns_flatplate.cfg_file  = "lam_flatplate_reg.cfg"
    fem_ns_flatplate.test_iter = 25
    fem_ns_flatplate.test_vals = [1.383727,3.175247,0.058387,0.257951] #last 4 columns
    fem_ns_flatplate.su2_exec  = "SU2_CFD"
    fem_ns_flatplate.new_output = True
    fem_ns_flatplate.timeout   = 1600
    fem_ns_flatplate.tol       = 0.00001
    test_list.append(fem_ns_flatplate)
    
    # Steady cylinder
    fem_ns_cylinder           = TestCase('fem_ns_cylinder')
    fem_ns_cylinder.cfg_dir   = "hom_navierstokes/CylinderViscous/nPoly3"
    fem_ns_cylinder.cfg_file  = "fem_Cylinder_reg.cfg"
    fem_ns_cylinder.test_iter = 10
    fem_ns_cylinder.test_vals = [0.454960,0.979123,-0.000028,79.984799] #last 4 columns
    fem_ns_cylinder.su2_exec  = "SU2_CFD"
    fem_ns_cylinder.new_output  = True
    fem_ns_cylinder.timeout   = 1600
    fem_ns_cylinder.tol       = 0.00001
    test_list.append(fem_ns_cylinder)

    # Steady sphere
    fem_ns_sphere           = TestCase('fem_ns_sphere')
    fem_ns_sphere.cfg_dir   = "hom_navierstokes/SphereViscous/nPoly3_QuadDominant"
    fem_ns_sphere.cfg_file  = "fem_Sphere_reg.cfg"
    fem_ns_sphere.test_iter = 10
    fem_ns_sphere.test_vals = [-0.288121,0.240324,0.000258,21.797363] #last 4 columns
    fem_ns_sphere.su2_exec  = "SU2_CFD"
    fem_ns_sphere.new_output = True
    fem_ns_sphere.timeout   = 1600
    fem_ns_sphere.tol       = 0.00001
    test_list.append(fem_ns_sphere)

    # Unsteady sphere ADER
    fem_ns_sphere_ader           = TestCase('fem_ns_sphere_ader')
    fem_ns_sphere_ader.cfg_dir   = "hom_navierstokes/SphereViscous/nPoly3_QuadDominant"
    fem_ns_sphere_ader.cfg_file  = "fem_Sphere_reg_ADER.cfg"
    fem_ns_sphere_ader.test_iter = 10
    fem_ns_sphere_ader.test_vals = [-35.000000,-35.000000,0.000047,31.110911] #last 4 columns
    fem_ns_sphere_ader.su2_exec  = "SU2_CFD"
    fem_ns_sphere_ader.new_output = True
    fem_ns_sphere_ader.unsteady   = True
    fem_ns_sphere_ader.timeout   = 1600
    fem_ns_sphere_ader.tol       = 0.00001
    test_list.append(fem_ns_sphere_ader)

    # Unsteady cylinder
    fem_ns_unsteady_cylinder           = TestCase('fem_ns_unsteady_cylinder')
    fem_ns_unsteady_cylinder.cfg_dir   = "hom_navierstokes/UnsteadyCylinder/nPoly4"
    fem_ns_unsteady_cylinder.cfg_file  = "fem_unst_cylinder.cfg"
    fem_ns_unsteady_cylinder.test_iter = 11
    fem_ns_unsteady_cylinder.test_vals = [-3.558582,-3.014464,-0.038927,1.383983] #last 4 columns
    fem_ns_unsteady_cylinder.su2_exec  = "SU2_CFD"
    fem_ns_unsteady_cylinder.new_output = True
    fem_ns_unsteady_cylinder.unsteady   = True
    fem_ns_unsteady_cylinder.timeout   = 1600
    fem_ns_unsteady_cylinder.tol       = 0.00001
    test_list.append(fem_ns_unsteady_cylinder)

    # Unsteady cylinder ADER
    fem_ns_unsteady_cylinder_ader           = TestCase('fem_ns_unsteady_cylinder_ader')
    fem_ns_unsteady_cylinder_ader.cfg_dir   = "hom_navierstokes/UnsteadyCylinder/nPoly4"
    fem_ns_unsteady_cylinder_ader.cfg_file  = "fem_unst_cylinder_ADER.cfg"
    fem_ns_unsteady_cylinder_ader.test_iter = 11
    fem_ns_unsteady_cylinder_ader.test_vals = [-35.000000,-35.000000,-0.041003,1.391339] #last 4 columns
    fem_ns_unsteady_cylinder_ader.su2_exec  = "SU2_CFD"
    fem_ns_unsteady_cylinder_ader.new_output = True
    fem_ns_unsteady_cylinder_ader.unsteady   = True
    fem_ns_unsteady_cylinder_ader.timeout   = 1600
    fem_ns_unsteady_cylinder_ader.tol       = 0.00001
    test_list.append(fem_ns_unsteady_cylinder_ader)

    ######################################
    ### Sliding Mesh                   ###
    ######################################

    # Uniform flow
    uniform_flow         = TestCase('uniform_flow')
    uniform_flow.cfg_dir   = "sliding_interface/uniform_flow"
    uniform_flow.cfg_file  = "uniform_NN.cfg"
    uniform_flow.test_iter = 10
    uniform_flow.test_vals = [10.000000, 0.000000, 0.432864, -5.596207] #last 4 columns
    uniform_flow.su2_exec  = "SU2_CFD"
    uniform_flow.timeout   = 1600
    uniform_flow.tol       = 0.000001
    uniform_flow.unsteady  = True
    uniform_flow.multizone = True
    test_list.append(uniform_flow) 

   # Channel_2D
    channel_2D           = TestCase('channel_2D')
    channel_2D.cfg_dir   = "sliding_interface/channel_2D"
    channel_2D.cfg_file  = "channel_2D_WA.cfg"
    channel_2D.test_iter = 4
    channel_2D.test_vals = [0.000000, 1.497633, 1.199455, 0.904683] #last 4 columns
    channel_2D.su2_exec  = "SU2_CFD"
    channel_2D.timeout   = 100
    channel_2D.tol       = 0.00001
    channel_2D.unsteady  = True
    channel_2D.multizone = True
    test_list.append(channel_2D)

    # Channel_3D
    channel_3D           = TestCase('channel_3D')
    channel_3D.cfg_dir   = "sliding_interface/channel_3D"
    channel_3D.cfg_file  = "channel_3D_WA.cfg"
    channel_3D.test_iter = 1
    channel_3D.test_vals = [1.000000, 0.000000, 1.834878, 1.892644, 1.606895] #last 4 columns
    channel_3D.su2_exec  = "SU2_CFD"
    channel_3D.timeout   = 1600
    channel_3D.tol       = 0.00001
    channel_3D.unsteady  = True
    channel_3D.multizone = True
    test_list.append(channel_3D)

    # Pipe
    pipe           = TestCase('pipe')
    pipe.cfg_dir   = "sliding_interface/pipe"
    pipe.cfg_file  = "pipe_NN.cfg"
    pipe.test_iter = 2
    pipe.test_vals = [ -0.341198, -0.045916, 0.800374, 0.699486] #last 4 columns
    pipe.su2_exec  = "SU2_CFD"
    pipe.timeout   = 1600
    pipe.tol       = 0.00001
    pipe.unsteady  = True
    pipe.multizone = True
    test_list.append(pipe)
    
    # Rotating cylinders
    rotating_cylinders           = TestCase('rotating_cylinders')
    rotating_cylinders.cfg_dir   = "sliding_interface/rotating_cylinders"
    rotating_cylinders.cfg_file  = "rot_cylinders_WA.cfg"
    rotating_cylinders.test_iter = 3
    rotating_cylinders.test_vals = [0.000000,2.018268, 1.942868, 1.696582] #last 4 columns
    rotating_cylinders.su2_exec  = "SU2_CFD"
    rotating_cylinders.timeout   = 1600
    rotating_cylinders.tol       = 0.00001
    rotating_cylinders.unsteady  = True
    rotating_cylinders.multizone = True
    test_list.append(rotating_cylinders)
    
    # Supersonic vortex shedding
    supersonic_vortex_shedding           = TestCase('supersonic_vortex_shedding')
    supersonic_vortex_shedding.cfg_dir   = "sliding_interface/supersonic_vortex_shedding"
    supersonic_vortex_shedding.cfg_file  = "sup_vor_shed_WA.cfg"
    supersonic_vortex_shedding.test_iter = 5
    supersonic_vortex_shedding.test_vals = [5.000000, 0.000000, 2.162202, 2.304635] #last 4 columns
    supersonic_vortex_shedding.su2_exec  = "SU2_CFD"
    supersonic_vortex_shedding.timeout   = 1600
    supersonic_vortex_shedding.tol       = 0.00001
    supersonic_vortex_shedding.unsteady  = True
    supersonic_vortex_shedding.multizone = True
    test_list.append(supersonic_vortex_shedding)
    
    # Bars_SST_2D
    bars_SST_2D           = TestCase('bars_SST_2D')
    bars_SST_2D.cfg_dir   = "sliding_interface/bars_SST_2D"
    bars_SST_2D.cfg_file  = "bars.cfg"
    bars_SST_2D.test_iter = 13
    bars_SST_2D.test_vals = [-0.941341, -0.270420, -1.840081,-1.159466] #last 4 columns
    bars_SST_2D.su2_exec  = "SU2_CFD"
    bars_SST_2D.timeout   = 1600
    bars_SST_2D.tol       = 0.00001
    bars_SST_2D.multizone = True
    test_list.append(bars_SST_2D)
    
    # Sliding mesh with incompressible flows (steady)
    slinc_steady           = TestCase('slinc_steady')
    slinc_steady.cfg_dir   = "sliding_interface/incompressible_steady"
    slinc_steady.cfg_file  = "config.cfg"
    slinc_steady.test_iter = 19
    slinc_steady.test_vals = [19.000000, -1.766116,  -2.206522] #last 3 columns
    slinc_steady.su2_exec  = "SU2_CFD"
    slinc_steady.timeout   = 100
    slinc_steady.tol       = 0.00001
    slinc_steady.multizone = True
    test_list.append(slinc_steady)
    
    # Sliding mesh with incompressible flows (unsteady)
    # slinc_unsteady           = TestCase('slinc_unsteady')
    # slinc_unsteady.cfg_dir   = "sliding_interface/incompressible_unsteady"
    # slinc_unsteady.cfg_file  = "config.cfg"
    # slinc_unsteady.test_iter = 19
    # slinc_unsteady.test_vals = [-3.515218,1.930028,0.000000,0.000000] #last 4 columns
    # slinc_unsteady.su2_exec  = "SU2_CFD"
    # slinc_unsteady.timeout   = 100
    # slinc_unsteady.tol       = 0.00001
    # slinc_unsteady.unsteady  = True
    # test_list.append(slinc_unsteady)

    ##########################
    ### FEA - FSI          ###
    ##########################

#    # Static beam, 3d
#    statbeam3d           = TestCase('statbeam3d')
#    statbeam3d.cfg_dir   = "fea_fsi/StatBeam_3d"
#    statbeam3d.cfg_file  = "configBeam_3d.cfg"
#    statbeam3d.new_output= True
#    statbeam3d.test_iter = 0
#    statbeam3d.test_vals = [-8.498274, -8.230638, -8.123824, 6.4095e+04] #last 4 columns
#    statbeam3d.su2_exec  = "SU2_CFD"
#    statbeam3d.timeout   = 1600
#    statbeam3d.tol       = 0.00001
#    test_list.append(statbeam3d)
#
#    # # Dynamic beam, 2d
#    dynbeam2d           = TestCase('dynbeam2d')
#    dynbeam2d.cfg_dir   = "fea_fsi/DynBeam_2d"
#    dynbeam2d.cfg_file  = "configBeam_2d.cfg"
#    dynbeam2d.unsteady  = True
#    dynbeam2d.new_output= True  
#    dynbeam2d.test_iter = 6 
#    dynbeam2d.test_vals = [-9.420640, -5.365872, -12.430382, 6.5210e+04] #last 4 columns
#    dynbeam2d.su2_exec  = "SU2_CFD"
#    dynbeam2d.timeout   = 1600
#    dynbeam2d.tol       = 0.00001
#    test_list.append(dynbeam2d)

    # # FSI, 2d
    fsi2d           = TestCase('fsi2d')
    fsi2d.cfg_dir   = "fea_fsi/WallChannel_2d"
    fsi2d.cfg_file  = "configFSI.cfg"
    fsi2d.test_iter = 4
    fsi2d.test_vals = [4.000000, 0.000000, -3.578187, -5.758338] #last 4 columns
    fsi2d.su2_exec  = "SU2_CFD"
    fsi2d.timeout   = 1600
    fsi2d.multizone = True
    fsi2d.unsteady  = True
    fsi2d.tol       = 0.00001
    test_list.append(fsi2d)    

    # # FSI, 2D airfoil with RBF interpolation
    # airfoilRBF           = TestCase('airfoil_fsi_rbf')
    # airfoilRBF.cfg_dir   = "fea_fsi/Airfoil_RBF"
    # airfoilRBF.cfg_file  = "config.cfg"
    # airfoilRBF.test_iter = 50
    # airfoilRBF.test_vals = [-8.000964, -2.600088, 0.276433, 0.000824] #last 4 columns
    # airfoilRBF.su2_exec  = "SU2_CFD"
    # airfoilRBF.timeout   = 1600
    # airfoilRBF.tol       = 0.00001
    # test_list.append(airfoilRBF)
   
    # ##########################
    # ### Zonal multiphysics ###
    # ##########################

    # CHT incompressible
    cht_incompressible           = TestCase('cht_incompressible')
    cht_incompressible.cfg_dir   = "coupled_cht/incompressible"
    cht_incompressible.cfg_file  = "config.cfg"
    cht_incompressible.test_iter = 10
    cht_incompressible.test_vals = [ 10, -1.552349, -6.795633, -3.342894] #last 4 columns
    cht_incompressible.su2_exec  = "SU2_CFD"
    cht_incompressible.timeout   = 1600
    cht_incompressible.multizone = True
    cht_incompressible.tol       = 0.00001
    test_list.append(cht_incompressible)

    ######################################
    ### RUN TESTS                      ###
    ######################################  

    pass_list = [ test.run_test() for test in test_list ]

    # Tests summary
    print('==================================================================')
    print('Summary of the serial tests')
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
