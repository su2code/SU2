#!/usr/bin/env python

## \file parallel_regression.py
#  \brief Python script for automated regression testing of SU2 examples
#  \author A. Aranake, A. Campos, T. Economon, T. Lukaczyk, S. Padron
#  \version 7.0.0 "Blackbird"
#
# SU2 Project Website: https://su2code.github.io
# 
# The SU2 Project is maintained by the SU2 Foundation 
# (http://su2foundation.org)
#
# Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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
    channel.test_vals = [-2.651955, 2.814536, 0.031770, 0.002870] #last 4 columns
    channel.su2_exec  = "parallel_computation.py -f"
    channel.timeout   = 1600
    channel.tol       = 0.00001
    test_list.append(channel)

    # NACA0012 
    naca0012           = TestCase('naca0012')
    naca0012.cfg_dir   = "euler/naca0012"
    naca0012.cfg_file  = "inv_NACA0012_Roe.cfg"
    naca0012.test_iter = 20
    naca0012.test_vals = [-4.055696, -3.564675, 0.336752, 0.021541] #last 4 columns
    naca0012.su2_exec  = "parallel_computation.py -f"
    naca0012.timeout   = 1600
    naca0012.tol       = 0.00001
    test_list.append(naca0012)

    # Supersonic wedge 
    wedge           = TestCase('wedge')
    wedge.cfg_dir   = "euler/wedge"
    wedge.cfg_file  = "inv_wedge_HLLC.cfg"
    wedge.test_iter = 20
    wedge.test_vals = [-0.941371, 4.787744, -0.208777, 0.036781] #last 4 columns
    wedge.su2_exec  = "parallel_computation.py -f"
    wedge.timeout   = 1600
    wedge.tol       = 0.00001
    test_list.append(wedge)

    # ONERA M6 Wing
    oneram6           = TestCase('oneram6')
    oneram6.cfg_dir   = "euler/oneram6"
    oneram6.cfg_file  = "inv_ONERAM6.cfg"
    oneram6.test_iter = 10
    oneram6.test_vals = [-10.208444, -9.625586, 0.281704, 0.011821] #last 4 columns
    oneram6.su2_exec  = "parallel_computation.py -f"
    oneram6.timeout   = 3200
    oneram6.tol       = 0.00001
    test_list.append(oneram6)

    # Fixed CL NACA0012
    fixedCL_naca0012           = TestCase('fixedcl_naca0012')
    fixedCL_naca0012.cfg_dir   = "fixed_cl/naca0012"
    fixedCL_naca0012.cfg_file  = "inv_NACA0012.cfg"
    fixedCL_naca0012.test_iter = 10
    fixedCL_naca0012.test_vals = [-12.137680, -6.705671, 0.300000, 0.019470] #last 4 columns
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
    polar_naca0012.test_vals = [-1.216242, 4.253547, 0.001423, 0.015977] #last 4 columns
    polar_naca0012.su2_exec  = "compute_polar.py -i 11"
    polar_naca0012.timeout   = 1600
    polar_naca0012.tol       = 0.00001
    test_list.append(polar_naca0012)

    # HYPERSONIC FLOW PAST BLUNT BODY          
    bluntbody           = TestCase('bluntbody')
    bluntbody.cfg_dir   = "euler/bluntbody"
    bluntbody.cfg_file  = "blunt.cfg"
    bluntbody.test_iter = 20
    bluntbody.test_vals = [0.540009, 6.916653, 0.000000, 1.868976] #last 4 columns
    bluntbody.su2_exec  = "parallel_computation.py -f"
    bluntbody.timeout   = 1600
    bluntbody.tol       = 0.00001
    test_list.append(bluntbody)

    ##########################
    ###  Compressible N-S  ###
    ##########################

    # Laminar flat plate
    flatplate           = TestCase('flatplate')
    flatplate.cfg_dir   = "navierstokes/flatplate"
    flatplate.cfg_file  = "lam_flatplate.cfg"
    flatplate.test_iter = 20
    flatplate.test_vals = [-4.648252, 0.813253, -0.130643, 0.024357] #last 4 columns
    flatplate.su2_exec  = "parallel_computation.py -f"
    flatplate.timeout   = 1600
    flatplate.tol       = 0.00001
    test_list.append(flatplate)

    # Laminar cylinder (steady)
    cylinder           = TestCase('cylinder')
    cylinder.cfg_dir   = "navierstokes/cylinder"
    cylinder.cfg_file  = "lam_cylinder.cfg"
    cylinder.test_iter = 25
    cylinder.test_vals = [-6.759136, -1.291221, 0.107150, 0.853257] #last 4 columns
    cylinder.su2_exec  = "parallel_computation.py -f"
    cylinder.timeout   = 1600
    cylinder.tol       = 0.00001
    test_list.append(cylinder)

    # Laminar cylinder (low Mach correction)
    cylinder_lowmach           = TestCase('cylinder_lowmach')
    cylinder_lowmach.cfg_dir   = "navierstokes/cylinder"
    cylinder_lowmach.cfg_file  = "cylinder_lowmach.cfg"
    cylinder_lowmach.test_iter = 25
    cylinder_lowmach.test_vals = [-6.870761, -1.408778, -0.228736, 112.418622] #last 4 columns
    cylinder_lowmach.su2_exec  = "parallel_computation.py -f"
    cylinder_lowmach.timeout   = 1600
    cylinder_lowmach.tol       = 0.00001
    test_list.append(cylinder_lowmach)

    # 2D Poiseuille flow (body force driven with periodic inlet / outlet)
    poiseuille           = TestCase('poiseuille')
    poiseuille.cfg_dir   = "navierstokes/poiseuille"
    poiseuille.cfg_file  = "lam_poiseuille.cfg"
    poiseuille.test_iter = 10
    poiseuille.test_vals = [-5.050864, 0.648220, 0.000349, 13.639525] #last 4 columns
    poiseuille.su2_exec  = "parallel_computation.py -f"
    poiseuille.timeout   = 1600
    poiseuille.tol       = 0.001
    test_list.append(poiseuille)

    # 2D Poiseuille flow (inlet profile file)
    poiseuille_profile           = TestCase('poiseuille_profile')
    poiseuille_profile.cfg_dir   = "navierstokes/poiseuille"
    poiseuille_profile.cfg_file  = "profile_poiseuille.cfg"
    poiseuille_profile.test_iter = 10
    poiseuille_profile.test_vals = [-12.493462, -7.671815, -0.000000, 2.085796] #last 4 columns
    poiseuille_profile.su2_exec  = "parallel_computation.py -f"
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
    rae2822_sa.test_vals = [-2.013881, -5.271311, 0.814981, 0.061858] #last 4 columns
    rae2822_sa.su2_exec  = "parallel_computation.py -f"
    rae2822_sa.timeout   = 1600
    rae2822_sa.tol       = 0.00001
    test_list.append(rae2822_sa)
    
    # RAE2822 SST
    rae2822_sst           = TestCase('rae2822_sst')
    rae2822_sst.cfg_dir   = "rans/rae2822"
    rae2822_sst.cfg_file  = "turb_SST_RAE2822.cfg"
    rae2822_sst.test_iter = 20
    rae2822_sst.test_vals = [-0.510623, 4.874845, 0.817046, 0.062059] #last 4 columns
    rae2822_sst.su2_exec  = "parallel_computation.py -f"
    rae2822_sst.timeout   = 1600
    rae2822_sst.tol       = 0.00001
    test_list.append(rae2822_sst)

    # RAE2822 SST_SUST
    rae2822_sst_sust           = TestCase('rae2822_sst_sust')
    rae2822_sst_sust.cfg_dir   = "rans/rae2822"
    rae2822_sst_sust.cfg_file  = "turb_SST_SUST_RAE2822.cfg"
    rae2822_sst_sust.test_iter = 20
    rae2822_sst_sust.test_vals = [-2.427686, 4.874845, 0.817046, 0.062059] #last 4 columns
    rae2822_sst_sust.su2_exec  = "parallel_computation.py -f"
    rae2822_sst_sust.timeout   = 1600
    rae2822_sst_sust.tol       = 0.00001
    test_list.append(rae2822_sst_sust)

    # Flat plate
    turb_flatplate           = TestCase('turb_flatplate')
    turb_flatplate.cfg_dir   = "rans/flatplate"
    turb_flatplate.cfg_file  = "turb_SA_flatplate.cfg"
    turb_flatplate.test_iter = 20
    turb_flatplate.test_vals = [-4.145487, -6.734014, -0.176490, 0.057451] #last 4 columns
    turb_flatplate.su2_exec  = "parallel_computation.py -f"
    turb_flatplate.timeout   = 1600
    turb_flatplate.tol       = 0.00001
    test_list.append(turb_flatplate)

    # ONERA M6 Wing
    turb_oneram6           = TestCase('turb_oneram6')
    turb_oneram6.cfg_dir   = "rans/oneram6"
    turb_oneram6.cfg_file  = "turb_ONERAM6.cfg"
    turb_oneram6.test_iter = 10
    turb_oneram6.test_vals = [-2.372346, -6.579371, 0.229867, 0.147637] #last 4 columns
    turb_oneram6.su2_exec  = "parallel_computation.py -f"
    turb_oneram6.timeout   = 3200
    turb_oneram6.tol       = 0.00001
    test_list.append(turb_oneram6)

    # NACA0012 (SA, FUN3D finest grid results: CL=1.0983, CD=0.01242)
    turb_naca0012_sa           = TestCase('turb_naca0012_sa')
    turb_naca0012_sa.cfg_dir   = "rans/naca0012"
    turb_naca0012_sa.cfg_file  = "turb_NACA0012_sa.cfg"
    turb_naca0012_sa.test_iter = 10
    turb_naca0012_sa.test_vals = [-12.078482, -16.147828, 1.064326, 0.019770] #last 4 columns
    turb_naca0012_sa.su2_exec  = "parallel_computation.py -f"
    turb_naca0012_sa.timeout   = 3200
    turb_naca0012_sa.tol       = 0.00001
    test_list.append(turb_naca0012_sa)
    
    # NACA0012 (SST, FUN3D finest grid results: CL=1.0840, CD=0.01253)
    turb_naca0012_sst           = TestCase('turb_naca0012_sst')
    turb_naca0012_sst.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst.cfg_file  = "turb_NACA0012_sst.cfg"
    turb_naca0012_sst.test_iter = 10
    turb_naca0012_sst.test_vals = [-15.274205, -6.250384, 1.049989, 0.019166] #last 4 columns
    turb_naca0012_sst.su2_exec  = "parallel_computation.py -f"
    turb_naca0012_sst.timeout   = 3200
    turb_naca0012_sst.tol       = 0.00001
    test_list.append(turb_naca0012_sst)

    # NACA0012 (SST_SUST, FUN3D finest grid results: CL=1.0840, CD=0.01253)
    turb_naca0012_sst_sust           = TestCase('turb_naca0012_sst_sust')
    turb_naca0012_sst_sust.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst_sust.cfg_file  = "turb_NACA0012_sst_sust.cfg"
    turb_naca0012_sst_sust.test_iter = 10
    turb_naca0012_sst_sust.test_vals = [-14.851291, -5.980663, 1.005228, 0.019014] #last 4 columns
    turb_naca0012_sst_sust.su2_exec  = "parallel_computation.py -f"
    turb_naca0012_sst_sust.timeout   = 3200
    turb_naca0012_sst_sust.tol       = 0.00001
    test_list.append(turb_naca0012_sst_sust)

    # PROPELLER
    propeller           = TestCase('propeller')
    propeller.cfg_dir   = "rans/propeller"
    propeller.cfg_file  = "propeller.cfg"
    propeller.test_iter = 10
    propeller.test_vals = [-3.389575, -8.409529, 0.000048, 0.056329] #last 4 columns
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
    turb_naca0012_sst_restart_mg.test_vals = [-7.826679, -8.082385, -3.174847, -0.000015, 0.081277] #last 5 columns
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
    inc_euler_naca0012.test_vals = [-4.802192, -3.799720, 0.497093, 0.007515] #last 4 columns
    inc_euler_naca0012.su2_exec  = "parallel_computation.py -f"
    inc_euler_naca0012.timeout   = 1600
    inc_euler_naca0012.tol       = 0.00001
    test_list.append(inc_euler_naca0012)

    # C-D nozzle with pressure inlet and mass flow outlet
    inc_nozzle           = TestCase('inc_nozzle')
    inc_nozzle.cfg_dir   = "incomp_euler/nozzle"
    inc_nozzle.cfg_file  = "inv_nozzle.cfg"
    inc_nozzle.test_iter = 20
    inc_nozzle.test_vals = [-5.959953, -4.932037, -0.000109, 0.121192] #last 4 columns
    inc_nozzle.su2_exec  = "parallel_computation.py -f"
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
    inc_lam_cylinder.test_vals = [-3.987107, -3.176477, -0.057813, 7.691336] #last 4 columns
    inc_lam_cylinder.su2_exec  = "parallel_computation.py -f"
    inc_lam_cylinder.timeout   = 1600
    inc_lam_cylinder.tol       = 0.00001
    test_list.append(inc_lam_cylinder)

    # Buoyancy-driven cavity
    inc_buoyancy          = TestCase('inc_buoyancy')
    inc_buoyancy.cfg_dir   = "incomp_navierstokes/buoyancy_cavity"
    inc_buoyancy.cfg_file  = "lam_buoyancy_cavity.cfg"
    inc_buoyancy.test_iter = 20
    inc_buoyancy.test_vals = [-4.436564, 0.508012, 0.000000, 0.000000] #last 4 columns
    inc_buoyancy.su2_exec  = "parallel_computation.py -f"
    inc_buoyancy.timeout   = 1600
    inc_buoyancy.tol       = 0.00001
    test_list.append(inc_buoyancy)

    # Laminar heated cylinder with polynomial fluid model
    inc_poly_cylinder          = TestCase('inc_poly_cylinder')
    inc_poly_cylinder.cfg_dir   = "incomp_navierstokes/cylinder"
    inc_poly_cylinder.cfg_file  = "poly_cylinder.cfg"
    inc_poly_cylinder.test_iter = 20
    inc_poly_cylinder.test_vals = [-7.812686, -2.080649, 0.016033, 1.912088] #last 4 columns
    inc_poly_cylinder.su2_exec  = "parallel_computation.py -f"
    inc_poly_cylinder.timeout   = 1600
    inc_poly_cylinder.tol       = 0.00001
    test_list.append(inc_poly_cylinder)
    
    # X-coarse laminar bend as a mixed element CGNS test
    inc_lam_bend          = TestCase('inc_lam_bend')
    inc_lam_bend.cfg_dir   = "incomp_navierstokes/bend"
    inc_lam_bend.cfg_file  = "lam_bend.cfg"
    inc_lam_bend.test_iter = 10
    inc_lam_bend.test_vals = [-3.437567, -3.088005, -0.022291, -0.172738] #last 4 columns
    inc_lam_bend.su2_exec  = "mpirun -n 2 SU2_CFD"
    inc_lam_bend.timeout   = 1600
    inc_lam_bend.tol       = 0.00001
    test_list.append(inc_lam_bend)

    ############################
    ### Incompressible RANS  ###
    ############################

    # NACA0012
    inc_turb_naca0012           = TestCase('inc_turb_naca0012')
    inc_turb_naca0012.cfg_dir   = "incomp_rans/naca0012"
    inc_turb_naca0012.cfg_file  = "naca0012.cfg"
    inc_turb_naca0012.test_iter = 20
    inc_turb_naca0012.test_vals = [-4.788584, -11.040550, -0.000000, 0.309515] #last 4 columns
    inc_turb_naca0012.su2_exec  = "parallel_computation.py -f"
    inc_turb_naca0012.timeout   = 1600
    inc_turb_naca0012.tol       = 0.00001
    test_list.append(inc_turb_naca0012)

    # NACA0012, SST_SUST
    inc_turb_naca0012_sst_sust           = TestCase('inc_turb_naca0012_sst_sust')
    inc_turb_naca0012_sst_sust.cfg_dir   = "incomp_rans/naca0012"
    inc_turb_naca0012_sst_sust.cfg_file  = "naca0012_SST_SUST.cfg"
    inc_turb_naca0012_sst_sust.test_iter = 20
    inc_turb_naca0012_sst_sust.test_vals = [-7.276357, 0.146097, -0.000000, 0.312019] #last 4 columns
    inc_turb_naca0012_sst_sust.su2_exec  = "parallel_computation.py -f"
    inc_turb_naca0012_sst_sust.timeout   = 1600
    inc_turb_naca0012_sst_sust.tol       = 0.00001
    test_list.append(inc_turb_naca0012_sst_sust)
    
    ####################
    ### DG-FEM Euler ###
    ####################
    
    # NACA0012
    fem_euler_naca0012           = TestCase('fem_euler_naca0012')
    fem_euler_naca0012.cfg_dir   = "hom_euler/NACA0012_5thOrder"
    fem_euler_naca0012.cfg_file  = "fem_NACA0012_reg.cfg"
    fem_euler_naca0012.test_iter = 10
    fem_euler_naca0012.test_vals = [-6.519946,-5.976944,0.255551,0.000028] #last 4 columns
    fem_euler_naca0012.su2_exec  = "mpirun -n 2 SU2_CFD"
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
    fem_ns_flatplate.su2_exec  = "mpirun -n 2 SU2_CFD"
    fem_ns_flatplate.timeout   = 1600
    fem_ns_flatplate.tol       = 0.00001
    test_list.append(fem_ns_flatplate)
    
    # Steady cylinder
    fem_ns_cylinder           = TestCase('fem_ns_cylinder')
    fem_ns_cylinder.cfg_dir   = "hom_navierstokes/CylinderViscous/nPoly3"
    fem_ns_cylinder.cfg_file  = "fem_Cylinder_reg.cfg"
    fem_ns_cylinder.test_iter = 10
    fem_ns_cylinder.test_vals = [0.454960,0.979123,-0.000028,79.984799] #last 4 columns
    fem_ns_cylinder.su2_exec  = "mpirun -n 2 SU2_CFD"
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
    fem_ns_unsteady_cylinder.timeout   = 1600
    fem_ns_unsteady_cylinder.unsteady  = True
    fem_ns_unsteady_cylinder.tol       = 0.00001
    test_list.append(fem_ns_unsteady_cylinder)

    # Unsteady cylinder ADER
    fem_ns_unsteady_cylinder_ader           = TestCase('fem_ns_unsteady_cylinder_ader')
    fem_ns_unsteady_cylinder_ader.cfg_dir   = "hom_navierstokes/UnsteadyCylinder/nPoly4"
    fem_ns_unsteady_cylinder_ader.cfg_file  = "fem_unst_cylinder_ADER.cfg"
    fem_ns_unsteady_cylinder_ader.test_iter = 11
    fem_ns_unsteady_cylinder_ader.test_vals = [-35.000000,-35.000000,-0.041003,1.391339] #last 4 columns
    fem_ns_unsteady_cylinder_ader.su2_exec  = "SU2_CFD"
    fem_ns_unsteady_cylinder_ader.timeout   = 1600
    fem_ns_unsteady_cylinder_ader.unsteady  = True
    fem_ns_unsteady_cylinder_ader.tol       = 0.00001
    test_list.append(fem_ns_unsteady_cylinder_ader)

    ############################
    ###      Transition      ###
    ############################

    # Schubauer-Klebanoff Natural Transition Case
    schubauer_klebanoff_transition              = TestCase('Schubauer_Klebanoff')
    schubauer_klebanoff_transition.cfg_dir      = "transition/Schubauer_Klebanoff"
    schubauer_klebanoff_transition.cfg_file     = "transitional_BC_model_ConfigFile.cfg"
    schubauer_klebanoff_transition.test_iter    = 10
    schubauer_klebanoff_transition.test_vals    = [-7.994738, -14.268367, 0.000046, 0.007987] #last 4 columns
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
    contadj_naca0012.test_vals = [-9.300829, -14.583240, 0.300920, 0.019552] #last 4 columns
    contadj_naca0012.su2_exec  = "parallel_computation.py -f"
    contadj_naca0012.timeout   = 1600
    contadj_naca0012.tol       = 0.00001
    test_list.append(contadj_naca0012)

    # Inviscid ONERA M6
    contadj_oneram6           = TestCase('contadj_oneram6')
    contadj_oneram6.cfg_dir   = "cont_adj_euler/oneram6"
    contadj_oneram6.cfg_file  = "inv_ONERAM6.cfg"
    contadj_oneram6.test_iter = 10
    contadj_oneram6.test_vals = [-12.132202, -12.702416, 0.685900, 0.007594] #last 4 columns
    contadj_oneram6.su2_exec  = "parallel_computation.py -f"
    contadj_oneram6.timeout   = 1600
    contadj_oneram6.tol       = 0.00001
    test_list.append(contadj_oneram6)

    # Inviscid WEDGE: tests averaged outflow total pressure adjoint
    contadj_wedge             = TestCase('contadj_wedge')
    contadj_wedge.cfg_dir   = "cont_adj_euler/wedge"
    contadj_wedge.cfg_file  = "inv_wedge_ROE.cfg"
    contadj_wedge.test_iter = 10  
    contadj_wedge.test_vals = [2.872691, -2.755572, 853010.000000, 0.000000] #last 4 columns
    contadj_wedge.su2_exec  = "parallel_computation.py -f"
    contadj_wedge.timeout   = 1600
    contadj_wedge.tol       = 0.00001
    test_list.append(contadj_wedge)

    # Inviscid fixed CL NACA0012
    contadj_fixed_CL_naca0012           = TestCase('contadj_fixedcl_naca0012')
    contadj_fixed_CL_naca0012.cfg_dir   = "fixed_cl/naca0012"
    contadj_fixed_CL_naca0012.cfg_file  = "inv_NACA0012_ContAdj.cfg"
    contadj_fixed_CL_naca0012.test_iter = 100
    contadj_fixed_CL_naca0012.test_vals = [0.317289, -5.193066, 0.360470, 0.000069] #last 4 columns
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
    contadj_ns_cylinder.test_vals = [-3.644966, -9.102024, 2.056700, -0.000000] #last 4 columns
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
    contadj_rans_naca0012.test_vals = [-0.794162, -5.761722, 19.214000, -0.000000] #last 4 columns
    contadj_rans_naca0012.su2_exec  = "parallel_computation.py -f"
    contadj_rans_naca0012.timeout   = 1600
    contadj_rans_naca0012.tol       = 0.00001
    test_list.append(contadj_rans_naca0012)
   
    # Adjoint turbulent NACA0012 with binary restarts
    contadj_rans_naca0012_bin           = TestCase('contadj_rans_naca0012_bin')
    contadj_rans_naca0012_bin.cfg_dir   = "cont_adj_rans/naca0012"
    contadj_rans_naca0012_bin.cfg_file  = "turb_nasa_binary.cfg"
    contadj_rans_naca0012_bin.test_iter = 18
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
    contadj_rans_rae2822.test_vals = [-5.381383, -10.883812, -0.212470, 0.005448] #last 4 columns
    contadj_rans_rae2822.su2_exec  = "parallel_computation.py -f"
    contadj_rans_rae2822.timeout   = 1600
    contadj_rans_rae2822.tol       = 0.00001
    test_list.append(contadj_rans_rae2822)

    #############################
    ### Compressibele RANS UQ ###
    #############################

    # NACA0012 1c
    turb_naca0012_1c           = TestCase('turb_naca0012_1c')
    turb_naca0012_1c.cfg_dir   = "rans_uq/naca0012"
    turb_naca0012_1c.cfg_file  = "turb_NACA0012_uq_1c.cfg"
    turb_naca0012_1c.test_iter = 10
    turb_naca0012_1c.test_vals = [-4.906243, 1.337959, 6.052221, 2.395629] #last 4 columns
    turb_naca0012_1c.su2_exec  = "parallel_computation.py -f"
    turb_naca0012_1c.timeout   = 1600
    turb_naca0012_1c.tol       = 0.00001
    test_list.append(turb_naca0012_1c)

    # NACA0012 2c
    turb_naca0012_2c           = TestCase('turb_naca0012_2c')
    turb_naca0012_2c.cfg_dir   = "rans_uq/naca0012"
    turb_naca0012_2c.cfg_file  = "turb_NACA0012_uq_2c.cfg"
    turb_naca0012_2c.test_iter = 10
    turb_naca0012_2c.test_vals = [-5.230215, 1.262233, 6.052131, 2.395642] #last 4 columns
    turb_naca0012_2c.su2_exec  = "parallel_computation.py -f"
    turb_naca0012_2c.timeout   = 1600
    turb_naca0012_2c.tol       = 0.00001
    test_list.append(turb_naca0012_2c)

    # NACA0012 3c
    turb_naca0012_3c           = TestCase('turb_naca0012_3c')
    turb_naca0012_3c.cfg_dir   = "rans_uq/naca0012"
    turb_naca0012_3c.cfg_file  = "turb_NACA0012_uq_3c.cfg"
    turb_naca0012_3c.test_iter = 10
    turb_naca0012_3c.test_vals = [-5.277132, 1.246269, 6.052474, 2.395999] #last 4 columns
    turb_naca0012_3c.su2_exec  = "parallel_computation.py -f"
    turb_naca0012_3c.timeout   = 1600
    turb_naca0012_3c.tol       = 0.00001
    test_list.append(turb_naca0012_3c)

    # NACA0012 p1c1
    turb_naca0012_p1c1           = TestCase('turb_naca0012_p1c1')
    turb_naca0012_p1c1.cfg_dir   = "rans_uq/naca0012"
    turb_naca0012_p1c1.cfg_file  = "turb_NACA0012_uq_p1c1.cfg"
    turb_naca0012_p1c1.test_iter = 10
    turb_naca0012_p1c1.test_vals = [-5.007909, 1.310935, 6.054741, 2.397402] #last 4 columns
    turb_naca0012_p1c1.su2_exec  = "parallel_computation.py -f"
    turb_naca0012_p1c1.timeout   = 1600
    turb_naca0012_p1c1.tol       = 0.00001
    test_list.append(turb_naca0012_p1c1)

    # NACA0012 p1c2
    turb_naca0012_p1c2           = TestCase('turb_naca0012_p1c2')
    turb_naca0012_p1c2.cfg_dir   = "rans_uq/naca0012"
    turb_naca0012_p1c2.cfg_file  = "turb_NACA0012_uq_p1c2.cfg"
    turb_naca0012_p1c2.test_iter = 10
    turb_naca0012_p1c2.test_vals = [-5.264043, 1.251313, 6.055013, 2.397532] #last 4 columns
    turb_naca0012_p1c2.su2_exec  = "parallel_computation.py -f"
    turb_naca0012_p1c2.timeout   = 1600
    turb_naca0012_p1c2.tol       = 0.00001
    test_list.append(turb_naca0012_p1c2)

    ######################################                                                                                  
    ### Harmonic Balance               ###                                                                                  
    ######################################                                                                                    

    # Description of the regression test 
    harmonic_balance           = TestCase('harmonic_balance')
    harmonic_balance.cfg_dir   = "harmonic_balance"
    harmonic_balance.cfg_file  = "HB.cfg"
    harmonic_balance.test_iter = 25
    harmonic_balance.test_vals = [-1.589862, 3.922099, -0.001443, 0.099456] #last 4 columns
    harmonic_balance.su2_exec  = "parallel_computation.py -f"
    harmonic_balance.timeout   = 1600
    harmonic_balance.tol       = 0.00001
    harmonic_balance.new_output = False
    test_list.append(harmonic_balance)

    # Turbulent pitching NACA 64a010 airfoil
    hb_rans_preconditioning           = TestCase('hb_rans_preconditioning')
    hb_rans_preconditioning.cfg_dir   = "harmonic_balance/hb_rans_preconditioning"
    hb_rans_preconditioning.cfg_file  = "davis.cfg"
    hb_rans_preconditioning.test_iter = 25
    hb_rans_preconditioning.test_vals = [-1.909596, -5.954720, 0.007773, 0.131219] #last 4 columns
    hb_rans_preconditioning.su2_exec  = "parallel_computation.py -f"
    hb_rans_preconditioning.timeout   = 1600
    hb_rans_preconditioning.tol       = 0.00001
    hb_rans_preconditioning.new_output = False
    test_list.append(hb_rans_preconditioning)

    ######################################
    ### Moving Wall                    ###
    ######################################

    # Lid-driven cavity
    cavity           = TestCase('cavity')
    cavity.cfg_dir   = "moving_wall/cavity"
    cavity.cfg_file  = "lam_cavity.cfg"
    cavity.test_iter = 25
    cavity.test_vals = [-5.588455, -0.124966, 0.308126, 0.940895] #last 4 columns
    cavity.su2_exec  = "parallel_computation.py -f"
    cavity.timeout   = 1600
    cavity.tol       = 0.00001
    test_list.append(cavity)

    # Spinning cylinder
    spinning_cylinder           = TestCase('spinning_cylinder')
    spinning_cylinder.cfg_dir   = "moving_wall/spinning_cylinder"
    spinning_cylinder.cfg_file  = "spinning_cylinder.cfg"
    spinning_cylinder.test_iter = 25
    spinning_cylinder.test_vals = [-7.857785, -2.425289, 1.554359, 1.531183] #last 4 columns
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
    square_cylinder.test_vals = [-1.162660, 0.066413, 1.399789, 2.220408] #last 4 columns
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
    sine_gust.test_vals = [-1.977545, 3.481778, -0.001525, -0.007375] #last 4 columns
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
    aeroelastic.test_vals = [0.081326, 0.033214, -0.001666, -0.000155] #last 4 columns
    aeroelastic.su2_exec  = "parallel_computation.py -f"
    aeroelastic.timeout   = 1600
    aeroelastic.tol       = 0.00001
    aeroelastic.unsteady  = True
    test_list.append(aeroelastic)

    # Delayed Detached Eddy Simulation
    ddes_flatplate        = TestCase('ddes_flatplate')
    ddes_flatplate.cfg_dir   = "ddes/flatplate"
    ddes_flatplate.cfg_file  = "ddes_flatplate.cfg"
    ddes_flatplate.test_iter = 10
    ddes_flatplate.test_vals = [-2.714758, -5.883004, -0.215005, 0.023783] #last 4 columns
    ddes_flatplate.su2_exec  = "parallel_computation.py -f"
    ddes_flatplate.timeout   = 1600
    ddes_flatplate.tol       = 0.00001
    ddes_flatplate.unsteady  = True
    test_list.append(ddes_flatplate)    

    # unsteady pitching NACA0015, SA
    unst_inc_turb_naca0015_sa           = TestCase('unst_inc_turb_naca0015_sa')
    unst_inc_turb_naca0015_sa.cfg_dir   = "unsteady/pitching_naca0015_rans_inc"
    unst_inc_turb_naca0015_sa.cfg_file  = "config_incomp_turb_sa.cfg"
    unst_inc_turb_naca0015_sa.test_iter = 1
    unst_inc_turb_naca0015_sa.test_vals = [-2.990702, -6.861616, 1.475736, 0.419771] #last 4 columns
    unst_inc_turb_naca0015_sa.su2_exec  = "parallel_computation.py -f"
    unst_inc_turb_naca0015_sa.timeout   = 1600
    unst_inc_turb_naca0015_sa.tol       = 0.00001
    unst_inc_turb_naca0015_sa.unsteady  = True
    test_list.append(unst_inc_turb_naca0015_sa)

    ######################################
    ### NICFD                          ###
    ######################################	

    # Rarefaction shock wave edge_VW
    edge_VW           = TestCase('edge_VW')
    edge_VW.cfg_dir   = "nicf/edge"
    edge_VW.cfg_file  = "edge_VW.cfg"
    edge_VW.test_iter = 100
    edge_VW.test_vals = [-5.203154, 0.933157, -0.000009, 0.000000] #last 4 columns
    edge_VW.su2_exec  = "parallel_computation.py -f"
    edge_VW.timeout   = 1600
    edge_VW.tol       = 0.00001
    test_list.append(edge_VW)

    # Rarefaction shock wave edge_PPR
    edge_PPR           = TestCase('edge_PPR')
    edge_PPR.cfg_dir   = "nicf/edge"
    edge_PPR.cfg_file  = "edge_PPR.cfg"
    edge_PPR.test_iter = 100
    edge_PPR.test_vals = [-5.385223, 0.755862, -0.000035, 0.000000] #last 4 columns
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
    Jones_tc.test_vals = [-5.316335, 0.355081, 44.772280, 2.269966] #last 4 columns
    Jones_tc.su2_exec  = "parallel_computation.py -f"
    Jones_tc.timeout   = 1600
    Jones_tc.new_output = False
    Jones_tc.tol       = 0.00001
    test_list.append(Jones_tc)

	# Jones APU Turbocharger restart
    Jones_tc_rst           = TestCase('jones_turbocharger_restart')
    Jones_tc_rst.cfg_dir   = "turbomachinery/APU_turbocharger"
    Jones_tc_rst.cfg_file  = "Jones_rst.cfg"
    Jones_tc_rst.test_iter = 5
    Jones_tc_rst.test_vals = [-3.034157, 0.013763, 82.263700, 2.792251] #last 4 columns
    Jones_tc_rst.su2_exec  = "parallel_computation.py -f"
    Jones_tc_rst.timeout   = 1600
    Jones_tc_rst.new_output = False
    Jones_tc_rst.tol       = 0.00001
    test_list.append(Jones_tc_rst)

    # 2D axial stage
    axial_stage2D           = TestCase('axial_stage2D')
    axial_stage2D.cfg_dir   = "turbomachinery/axial_stage_2D"
    axial_stage2D.cfg_file  = "Axial_stage2D.cfg"
    axial_stage2D.test_iter = 20
    axial_stage2D.test_vals = [-1.933241, 5.310060, 73.357650, 1.841414] #last 4 columns
    axial_stage2D.su2_exec  = "parallel_computation.py -f"
    axial_stage2D.timeout   = 1600
    axial_stage2D.new_output = False
    axial_stage2D.tol       = 0.00001
    test_list.append(axial_stage2D)
    
    # 2D transonic stator
    transonic_stator           = TestCase('transonic_stator')
    transonic_stator.cfg_dir   = "turbomachinery/transonic_stator_2D"
    transonic_stator.cfg_file  = "transonic_stator.cfg"
    transonic_stator.test_iter = 20
    transonic_stator.test_vals = [-0.574097, 5.820577, 96.994270, 0.062868] #last 4 columns
    transonic_stator.su2_exec  = "parallel_computation.py -f"
    transonic_stator.timeout   = 1600
    transonic_stator.new_output = False
    transonic_stator.tol       = 0.00001
    test_list.append(transonic_stator)
    
    # 2D transonic stator restart
    transonic_stator_rst           = TestCase('transonic_stator_restart')
    transonic_stator_rst.cfg_dir   = "turbomachinery/transonic_stator_2D"
    transonic_stator_rst.cfg_file  = "transonic_stator_rst.cfg"
    transonic_stator_rst.test_iter = 20
    transonic_stator_rst.test_vals = [-6.797150, -0.755136, 5.007406, 0.0029491] #last 4 columns
    transonic_stator_rst.su2_exec  = "parallel_computation.py -f"
    transonic_stator_rst.timeout   = 1600
    transonic_stator_rst.new_output = False
    transonic_stator_rst.tol       = 0.00001
    test_list.append(transonic_stator_rst)

    ######################################
    ### Sliding Mesh                   ###
    ######################################

    # Uniform flow
    uniform_flow         = TestCase('uniform_flow')
    uniform_flow.cfg_dir   = "sliding_interface/uniform_flow"
    uniform_flow.cfg_file  = "uniform_NN.cfg"
    uniform_flow.test_iter = 5
    uniform_flow.test_vals = [5.000000, 0.000000, -0.188747, -10.631534] #last 4 columns
    uniform_flow.su2_exec  = "parallel_computation.py -f"
    uniform_flow.timeout   = 1600
    uniform_flow.tol       = 0.000001
    uniform_flow.unsteady  = True
    uniform_flow.multizone = True
    test_list.append(uniform_flow) 

    # Channel_2D
    channel_2D           = TestCase('channel_2D')
    channel_2D.cfg_dir   = "sliding_interface/channel_2D"
    channel_2D.cfg_file  = "channel_2D_WA.cfg"
    channel_2D.test_iter = 2
    channel_2D.test_vals = [2.000000, 0.000000, 0.398157, 0.353079, 0.405679] #last 4 columns
    channel_2D.su2_exec  = "parallel_computation.py -f"
    channel_2D.timeout   = 100
    channel_2D.tol       = 0.00001
    channel_2D.unsteady  = True
    channel_2D.multizone = True
    test_list.append(channel_2D)

    # Channel_3D
    channel_3D           = TestCase('channel_3D')
    channel_3D.cfg_dir   = "sliding_interface/channel_3D"
    channel_3D.cfg_file  = "channel_3D_WA.cfg"
    channel_3D.test_iter = 2
    channel_3D.test_vals = [2.000000, 0.000000, 0.620137, 0.505221, 0.415605] #last 4 columns
    channel_3D.su2_exec  = "parallel_computation.py -f"
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
    pipe.test_vals = [0.150025, 0.491954, 0.677756, 0.963980, 1.006936] #last 4 columns
    pipe.su2_exec  = "parallel_computation.py -f"
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
    rotating_cylinders.test_vals = [3.000000, 0.000000, 0.777044, 1.134689, 1.224326] #last 4 columns
    rotating_cylinders.su2_exec  = "parallel_computation.py -f"
    rotating_cylinders.timeout   = 1600
    rotating_cylinders.tol       = 0.00001
    rotating_cylinders.unsteady  = True
    rotating_cylinders.multizone  = True
    test_list.append(rotating_cylinders)

    # Supersonic vortex shedding
    supersonic_vortex_shedding           = TestCase('supersonic_vortex_shedding')
    supersonic_vortex_shedding.cfg_dir   = "sliding_interface/supersonic_vortex_shedding"
    supersonic_vortex_shedding.cfg_file  = "sup_vor_shed_WA.cfg"
    supersonic_vortex_shedding.test_iter = 5
    supersonic_vortex_shedding.test_vals = [5.000000, 0.000000, 1.228298, 1.648277] #last 4 columns
    supersonic_vortex_shedding.su2_exec  = "parallel_computation.py -f"
    supersonic_vortex_shedding.timeout   = 1600
    supersonic_vortex_shedding.tol       = 0.00001
    supersonic_vortex_shedding.unsteady  = True
    supersonic_vortex_shedding.multizone  = True
    test_list.append(supersonic_vortex_shedding)

    # Bars_SST_2D
    bars_SST_2D           = TestCase('bars_SST_2D')
    bars_SST_2D.cfg_dir   = "sliding_interface/bars_SST_2D"
    bars_SST_2D.cfg_file  = "bars.cfg"
    bars_SST_2D.test_iter = 13
    bars_SST_2D.test_vals = [13.000000, -0.590195, -1.955385] #last 4 columns
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
    slinc_steady.test_vals = [19.000000,  -1.766116, -2.206522] #last 4 columns
    slinc_steady.su2_exec  = "SU2_CFD"
    slinc_steady.timeout   = 100
    slinc_steady.tol       = 0.00002
    slinc_steady.multizone = True
    test_list.append(slinc_steady)
    
    # Sliding mesh with incompressible flows (unsteady)
    # slinc_unsteady           = TestCase('slinc_unsteady')
    # slinc_unsteady.cfg_dir   = "sliding_interface/incompressible_unsteady"
    # slinc_unsteady.cfg_file  = "config.cfg"
    # slinc_unsteady.test_iter = 19
    # slinc_unsteady.test_vals = [-3.513701,1.931626,0.000000,0.000000] #last 4 columns
    # slinc_unsteady.su2_exec  = "SU2_CFD"
    # slinc_unsteady.timeout   = 100
    # slinc_unsteady.tol       = 0.00001
    # slinc_unsteady.unsteady  = True
    # test_list.append(slinc_unsteady)    

    ##########################
    ### FEA - FSI          ###
    ##########################   

    # Static beam, 3d
    statbeam3d           = TestCase('statbeam3d')
    statbeam3d.cfg_dir   = "fea_fsi/StatBeam_3d"
    statbeam3d.cfg_file  = "configBeam_3d.cfg"
    statbeam3d.test_iter = 0
    statbeam3d.test_vals = [-8.396797, -8.162206, -8.156102, 64095.0] #last 4 columns
    statbeam3d.su2_exec  = "parallel_computation_fsi.py -f"
    statbeam3d.timeout   = 1600
    statbeam3d.tol       = 0.00001
    test_list.append(statbeam3d)

    # Dynamic beam, 2d
    dynbeam2d           = TestCase('dynbeam2d')
    dynbeam2d.cfg_dir   = "fea_fsi/DynBeam_2d"
    dynbeam2d.cfg_file  = "configBeam_2d.cfg"
    dynbeam2d.test_iter = 6
    dynbeam2d.test_vals = [-3.240015, 2.895057, -0.353146, 6.6127e+04] #last 4 columns
    dynbeam2d.su2_exec  = "parallel_computation_fsi.py -f"
    dynbeam2d.timeout   = 1600
    dynbeam2d.unsteady  = True
    dynbeam2d.tol       = 0.00001
    test_list.append(dynbeam2d)

    # FSI, 2d
    fsi2d           = TestCase('fsi2d')
    fsi2d.cfg_dir   = "fea_fsi/WallChannel_2d"
    fsi2d.cfg_file  = "configFSI.cfg"
    fsi2d.test_iter = 4
    fsi2d.test_vals = [4.000000, 0.000000, -3.801272, -4.123968] #last 4 columns
    fsi2d.su2_exec  = "parallel_computation_fsi.py -f"
    fsi2d.timeout   = 1600
    fsi2d.multizone= True
    fsi2d.unsteady = True
    fsi2d.tol       = 0.00001
    test_list.append(fsi2d)
    
    # FSI, Static, 2D, new mesh solver
    stat_fsi           = TestCase('stat_fsi')
    stat_fsi.cfg_dir   = "fea_fsi/stat_fsi"
    stat_fsi.cfg_file  = "config.cfg"
    stat_fsi.test_iter = 7
    stat_fsi.test_vals = [-3.313322, -4.963786, 0.000000, 46.000000] #last 5 columns
    stat_fsi.su2_exec  = "mpirun -n 2 SU2_CFD"
    stat_fsi.multizone = True
    stat_fsi.timeout   = 1600
    stat_fsi.tol       = 0.00001
    test_list.append(stat_fsi)

    # FSI, Dynamic, 2D, new mesh solver
    dyn_fsi           = TestCase('dyn_fsi')
    dyn_fsi.cfg_dir   = "fea_fsi/dyn_fsi"
    dyn_fsi.cfg_file  = "config.cfg"
    dyn_fsi.test_iter = 4
    dyn_fsi.test_vals = [-4.389734, -4.060117, 0.000000, 64.000000] #last 5 columns
    dyn_fsi.multizone = True
    dyn_fsi.unsteady  = True
    dyn_fsi.su2_exec  = "mpirun -n 2 SU2_CFD"
    dyn_fsi.timeout   = 1600
    dyn_fsi.tol       = 0.00001
    test_list.append(dyn_fsi)

    # FSI, Static, 2D, new mesh solver, restart
    stat_fsi_restart           = TestCase('stat_fsi_restart')
    stat_fsi_restart.cfg_dir   = "fea_fsi/stat_fsi"
    stat_fsi_restart.cfg_file  = "config_restart.cfg"
    stat_fsi_restart.test_iter = 1
    stat_fsi_restart.test_vals = [-3.422307, -4.212725, 0.000000, 46.000000] #last 5 columns
    stat_fsi_restart.su2_exec  = "mpirun -n 2 SU2_CFD"
    stat_fsi_restart.multizone = True
    stat_fsi_restart.timeout   = 1600
    stat_fsi_restart.tol       = 0.00001
    test_list.append(stat_fsi_restart)

    ##########################
    ### Zonal multiphysics ###
    ##########################

    # CHT incompressible
    cht_incompressible           = TestCase('cht_incompressible')
    cht_incompressible.cfg_dir   = "coupled_cht/incomp_2d"
    cht_incompressible.cfg_file  = "cht_2d_3cylinders.cfg"
    cht_incompressible.test_iter = 10
    cht_incompressible.test_vals = [-2.132187, -0.579649, -0.579649, -0.579649] #last 4 columns
    cht_incompressible.su2_exec  = "SU2_CFD"
    cht_incompressible.timeout   = 1600
    cht_incompressible.multizone = True
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
    pywrapper_naca0012.test_vals = [-7.298733, -6.733533, 0.333456, 0.021233] #last 4 columns
    pywrapper_naca0012.su2_exec  = "mpirun -np 2 SU2_CFD.py --parallel -f"
    pywrapper_naca0012.timeout   = 1600
    pywrapper_naca0012.tol       = 0.00001
    test_list.append(pywrapper_naca0012)

    # NACA0012 (SST, FUN3D results for finest grid: CL=1.0840, CD=0.01253)
    pywrapper_turb_naca0012_sst           = TestCase('pywrapper_turb_naca0012_sst')
    pywrapper_turb_naca0012_sst.cfg_dir   = "rans/naca0012"
    pywrapper_turb_naca0012_sst.cfg_file  = "turb_NACA0012_sst.cfg"
    pywrapper_turb_naca0012_sst.test_iter = 10
    pywrapper_turb_naca0012_sst.test_vals = [-15.274205, -6.250384, 1.049989, 0.019166] #last 4 columns
    pywrapper_turb_naca0012_sst.su2_exec  = "mpirun -np 2 SU2_CFD.py --parallel -f"
    pywrapper_turb_naca0012_sst.timeout   = 3200
    pywrapper_turb_naca0012_sst.tol       = 0.00001
    test_list.append(pywrapper_turb_naca0012_sst)

    # Square cylinder
    pywrapper_square_cylinder           = TestCase('pywrapper_square_cylinder')
    pywrapper_square_cylinder.cfg_dir   = "unsteady/square_cylinder"
    pywrapper_square_cylinder.cfg_file  = "turb_square.cfg"
    pywrapper_square_cylinder.test_iter = 3
    pywrapper_square_cylinder.test_vals = [-1.162660, 0.066413, 1.399789, 2.220408] #last 4 columns
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
    pywrapper_aeroelastic.test_vals = [0.081326, 0.033214, -0.001666, -0.000155] #last 4 columns
    pywrapper_aeroelastic.su2_exec  = "mpirun -np 2 SU2_CFD.py --parallel -f"
    pywrapper_aeroelastic.timeout   = 1600
    pywrapper_aeroelastic.tol       = 0.00001
    pywrapper_aeroelastic.unsteady  = True
    test_list.append(pywrapper_aeroelastic)

    # FSI, 2d
    pywrapper_fsi2d           = TestCase('pywrapper_fsi2d')
    pywrapper_fsi2d.cfg_dir   = "fea_fsi/WallChannel_2d"
    pywrapper_fsi2d.cfg_file  = "configFSI.cfg"
    pywrapper_fsi2d.test_iter = 4
    pywrapper_fsi2d.test_vals = [4.000000, 0.000000, -3.801272, -4.123968] #last 4 columns
    pywrapper_fsi2d.su2_exec  = "mpirun -np 2 SU2_CFD.py --nZone 2 --fsi True --parallel -f"
    pywrapper_fsi2d.timeout   = 1600
    pywrapper_fsi2d.unsteady  = True
    pywrapper_fsi2d.multizone = True
    pywrapper_fsi2d.tol       = 0.00001
    test_list.append(pywrapper_fsi2d)

    # Unsteady CHT
    pywrapper_unsteadyCHT               = TestCase('pywrapper_unsteadyCHT')
    pywrapper_unsteadyCHT.cfg_dir       = "py_wrapper/flatPlate_unsteady_CHT"
    pywrapper_unsteadyCHT.cfg_file      = "unsteady_CHT_FlatPlate_Conf.cfg"
    pywrapper_unsteadyCHT.test_iter     = 5
    pywrapper_unsteadyCHT.test_vals     = [-1.614167, 2.245730, -0.000767, 0.175598] #last 4 columns
    pywrapper_unsteadyCHT.su2_exec      = "mpirun -np 2 python launch_unsteady_CHT_FlatPlate.py --parallel -f"
    pywrapper_unsteadyCHT.timeout       = 1600
    pywrapper_unsteadyCHT.tol           = 0.00001
    pywrapper_unsteadyCHT.unsteady      = True
    pywrapper_unsteadyCHT.new_output    = True
    test_list.append(pywrapper_unsteadyCHT)

    # Rigid motion
    pywrapper_rigidMotion               = TestCase('pywrapper_rigidMotion')
    pywrapper_rigidMotion.cfg_dir       = "py_wrapper/flatPlate_rigidMotion"
    pywrapper_rigidMotion.cfg_file      = "flatPlate_rigidMotion_Conf.cfg"
    pywrapper_rigidMotion.test_iter     = 5
    pywrapper_rigidMotion.test_vals     = [-1.614165, 2.242648, -0.037218, 0.173751] #last 4 columns
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
    tutorial_inv_bump.test_vals  = [-1.437425, 4.075857, 0.003000, 0.012720] #last 4 columns
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
    tutorial_inv_wedge.test_vals  = [-0.481460, 5.253008, -0.292159, 0.052922] #last 4 columns
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
    tutorial_inv_onera.test_vals  = [-5.204928, -4.597762, 0.247124, 0.085734] #last 4 columns
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
    tutorial_lam_cylinder.test_vals  = [-6.162141, -0.699617, 0.186570, 69.267308] #last 4 columns
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
    tutorial_lam_flatplate.test_vals  = [-2.821818, 2.657591, -0.400044, 0.029413] #last 4 columns
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
    tutorial_turb_flatplate.test_vals  = [-2.258584, -4.899502, -0.429387, 0.201236] #last 4 columns
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
    tutorial_turb_oneram6.test_vals  = [-4.499497, -11.473637, 0.332666, 0.098280] #last 4 columns
    tutorial_turb_oneram6.su2_exec   = "mpirun -np 2 SU2_CFD"
    tutorial_turb_oneram6.timeout    = 1600
    tutorial_turb_oneram6.tol        = 0.00001
    test_list.append(tutorial_turb_oneram6)

    # Inviscid NACA 0012 Design
    tutorial_design_inv_naca0012            = TestCase('design_inv_naca0012')
    tutorial_design_inv_naca0012.cfg_dir    = "../Tutorials/Inviscid_2D_Unconstrained_NACA0012"
    tutorial_design_inv_naca0012.cfg_file   = "inv_NACA0012_basic.cfg"
    tutorial_design_inv_naca0012.test_iter  = 0
    tutorial_design_inv_naca0012.test_vals  = [-3.585391, -2.989014, 0.134515, 0.208523] #last 4 columns
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
    tutorial_design_turb_rae2822.test_vals  = [-1.700114, -4.941261, 0.218432, 0.190639] #last 4 columns
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

    ##############################################
    ### Method of Manufactured Solutions (MMS) ###
    ##############################################

    # FVM, compressible, laminar N-S
    mms_fvm_ns           = TestCase('mms_fvm_ns')
    mms_fvm_ns.cfg_dir   = "mms/fvm_navierstokes"
    mms_fvm_ns.cfg_file  = "lam_mms_roe.cfg"
    mms_fvm_ns.test_iter = 20
    mms_fvm_ns.test_vals = [-2.851428, 2.192348, 0.000000, 0.000000] #last 4 columns
    mms_fvm_ns.su2_exec  = "mpirun -n 2 SU2_CFD"
    mms_fvm_ns.timeout   = 1600
    mms_fvm_ns.tol       = 0.0001
    test_list.append(mms_fvm_ns)
    
    # FVM, incompressible, euler
    mms_fvm_inc_euler           = TestCase('mms_fvm_inc_euler')
    mms_fvm_inc_euler.cfg_dir   = "mms/fvm_incomp_euler"
    mms_fvm_inc_euler.cfg_file  = "inv_mms_jst.cfg"
    mms_fvm_inc_euler.test_iter = 20
    mms_fvm_inc_euler.test_vals = [-9.128515, -9.441740, 0.000000, 0.000000] #last 4 columns
    mms_fvm_inc_euler.su2_exec  = "mpirun -np 2 SU2_CFD"
    mms_fvm_inc_euler.timeout   = 1600
    mms_fvm_inc_euler.tol       = 0.0001
    test_list.append(mms_fvm_inc_euler)
    
    # FVM, incompressible, laminar N-S
    mms_fvm_inc_ns           = TestCase('mms_fvm_inc_ns')
    mms_fvm_inc_ns.cfg_dir   = "mms/fvm_incomp_navierstokes"
    mms_fvm_inc_ns.cfg_file  = "lam_mms_fds.cfg"
    mms_fvm_inc_ns.test_iter = 20
    mms_fvm_inc_ns.test_vals = [-7.414944, -7.631546, 0.000000, 0.000000] #last 4 columns
    mms_fvm_inc_ns.su2_exec  = "mpirun -np 2 SU2_CFD"
    mms_fvm_inc_ns.timeout   = 1600
    mms_fvm_inc_ns.tol       = 0.0001
    test_list.append(mms_fvm_inc_ns)

    # DG, compressible, euler
    ringleb_dg_euler           = TestCase('ringleb_dg_euler')
    ringleb_dg_euler.cfg_dir   = "mms/dg_ringleb"
    ringleb_dg_euler.cfg_file  = "ringleb_dg.cfg"
    ringleb_dg_euler.test_iter = 100
    ringleb_dg_euler.test_vals = [-5.136652, -4.724941, 0.000000, 0.000000] #last 4 columns
    ringleb_dg_euler.su2_exec  = "SU2_CFD"
    ringleb_dg_euler.timeout   = 1600
    ringleb_dg_euler.tol       = 0.0001
    test_list.append(ringleb_dg_euler)

    # DG, compressible, laminar N-S
    mms_dg_ns           = TestCase('mms_dg_ns')
    mms_dg_ns.cfg_dir   = "mms/dg_navierstokes"
    mms_dg_ns.cfg_file  = "lam_mms_dg.cfg"
    mms_dg_ns.test_iter = 100
    mms_dg_ns.test_vals = [-1.845393, 3.520699, 0.000000, 0.000000] #last 4 columns
    mms_dg_ns.su2_exec  = "SU2_CFD"
    mms_dg_ns.timeout   = 1600
    mms_dg_ns.tol       = 0.0001
    test_list.append(mms_dg_ns)

    # DG, compressible, laminar N-S 3D
    mms_dg_ns_3d           = TestCase('mms_dg_ns_3d')
    mms_dg_ns_3d.cfg_dir   = "mms/dg_navierstokes_3d"
    mms_dg_ns_3d.cfg_file  = "lam_mms_dg_3d.cfg"
    mms_dg_ns_3d.test_iter = 100
    mms_dg_ns_3d.test_vals = [-0.146826, 5.356413, 0.000000, 0.000000] #last 4 columns
    mms_dg_ns_3d.su2_exec  = "SU2_CFD"
    mms_dg_ns_3d.timeout   = 1600
    mms_dg_ns_3d.tol       = 0.0001
    test_list.append(mms_dg_ns_3d)
    
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
    naca0012_def.test_vals = [0.00354532] #residual
    naca0012_def.su2_exec  = "mpirun -n 2 SU2_DEF"
    naca0012_def.timeout   = 1600
    naca0012_def.tol       = 1e-8
    
    pass_list.append(naca0012_def.run_def())
    test_list.append(naca0012_def)
    
    # Inviscid NACA0012 based on SURFACE_FILE input (surface_bump.dat)
    naca0012_def_file            = TestCase('naca0012_def_file')
    naca0012_def_file.cfg_dir   = "deformation/naca0012"
    naca0012_def_file.cfg_file  = "surface_file_NACA0012.cfg"
    naca0012_def_file.test_iter = 10
    naca0012_def_file.test_vals = [0.00354532] #residual
    naca0012_def_file.su2_exec  = "mpirun -n 2 SU2_DEF"
    naca0012_def_file.timeout   = 1600
    naca0012_def_file.tol       = 1e-8
    
    pass_list.append(naca0012_def_file.run_def())
    test_list.append(naca0012_def_file)

    # RAE2822 (mixed tris + quads)
    rae2822_def            = TestCase('rae2822_def')
    rae2822_def.cfg_dir   = "deformation/rae2822"
    rae2822_def.cfg_file  = "def_RAE2822.cfg"
    rae2822_def.test_iter = 10
    rae2822_def.test_vals = [8.2438e-09] #residual
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
    brick_tets_def.test_vals = [0.000954604] #residual
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
    brick_hex_def.test_vals = [0.000199532] #residual
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
    brick_pyra_def.test_vals = [0.00160022] #residual
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
    brick_prism_def.test_vals = [0.00260853] #residual
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
    brick_prism_rans_def.test_vals = [3.10348e-07] #residual
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
    brick_hex_rans_def.test_vals = [3.55635e-06] #residual
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
    cylinder_ffd_def.test_vals = [0.00054847] #residual
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
    sphere_ffd_def.test_vals = [0.00359947] #residual
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
    sphere_ffd_def_bspline.test_vals = [0.00208206] #residual
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
