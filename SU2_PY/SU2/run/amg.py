#!/usr/bin/env python

## \file amg.py
#  \brief python script for running mesh adaptation using the AMG Inria library
#  \author Victorien Menier, Brian Mungu\'ia
#  \version 7.3.0 "Blackbird"
#
# SU2 Project Website: https://su2code.github.io
#
# The SU2 Project is maintained by the SU2 Foundation
# (http://su2foundation.org)
#
# Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
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

import os, shutil, copy, time

from .. import io as su2io
from .. import amginria as su2amg
from .interface import CFD as SU2_CFD

def amg(config):
    """
    Runs the a mesh adaptation loop with the AMG library.

    Inputs:
        config - an SU2 config object
    """

    print('SU2-AMG Anisotropic Mesh Adaptation')

    #--- Check config options related to mesh adaptation

    pyadap_options = [ 'ADAP_SIZES', 'ADAP_SUBITER', 'ADAP_SENSOR', 'ADAP_BACK',
                       'ADAP_HGRAD', 'ADAP_RESIDUAL_REDUCTION', 'ADAP_FLOW_ITER',
                       'ADAP_ADJ_ITER', 'ADAP_CFL', 'ADAP_INV_BACK', 'ADAP_ORTHO',
                       'ADAP_RDG' ]
    required_options = [ 'ADAP_SIZES', 'ADAP_SUBITER', 'ADAP_SENSOR', 'ADAP_HMAX',
                         'ADAP_HMIN', 'MESH_FILENAME', 'RESTART_SOL', 'MESH_OUT_FILENAME' ]

    if not all (opt in config for opt in required_options):
        err = '\n\n## ERROR : Missing options: \n'
        for opt in required_options:
            if not opt in config:
                err += opt + '\n'
        raise AttributeError(err)

    #--- Print adap options

    print(su2amg.print_adap_options(config))

    #--- Target mesh sizes and subiterations at each size

    mesh_sizes = su2amg.get_mesh_sizes(config)
    sub_iter   = su2amg.get_sub_iterations(config)

    #--- Solver iterations/ residual reduction param for each size level

    flow_iter = su2amg.get_flow_iter(config)
    adj_iter  = su2amg.get_adj_iter(config)
    flow_cfl  = su2amg.get_flow_cfl(config)

    adap_sensor = config.ADAP_SENSOR
    sensor_avail = ['MACH', 'PRES', 'MACH_PRES', 'GOAL']

    if adap_sensor not in sensor_avail:
        raise ValueError(f'Unknown adaptation sensor {adap_sensor}. Available options are {sensor_avail}.')

    if len(mesh_sizes) != len(sub_iter):
        raise ValueError(f'Inconsistent number of mesh sizes and sub-iterations. {len(mesh_sizes)} mesh sizes and {len(sub_iter)} sub-iterations provided.')

    #--- Change current directory

    warn = True
    base_dir = os.getcwd()
    adap_dir = './adap'

    if os.path.exists(adap_dir):
        print('./adap exists. Removing old mesh adaptation in 10s.')
        if warn : time.sleep(10)
        shutil.rmtree(adap_dir)
        print(f'The {adap_dir} folder was deleted.')

    dir = f'{adap_dir}/ite0'
    os.makedirs(dir)
    os.chdir(dir)
    os.symlink(os.path.join(base_dir, config.MESH_FILENAME), config.MESH_FILENAME)

    meshfil = config['MESH_FILENAME']

    #--- Format of history file

    history_format = config.TABULAR_FORMAT
    if (history_format == 'TECPLOT'):
        history_filename = os.path.join(base_dir, 'history_adap.dat')
    else:
        history_filename = os.path.join(base_dir, 'history_adap.csv')

    #--- Get mesh dimension

    dim = su2amg.get_su2_dim(meshfil)
    if ( dim != 2 and dim != 3 ):
        raise ValueError('Wrong dimension number.')

    #--- AMG parameters

    config_amg = su2amg.get_amg_config(config, dim)

    #--- Compute initial solution if needed, else link current files

    config_cfd = copy.deepcopy(config)
    config_cfd_ad = copy.deepcopy(config)
    for opt in pyadap_options:
        config_cfd.pop(opt, None)
        config_cfd_ad.pop(opt, None)

    #--- Check config for filenames if restarting
    if config['RESTART_SOL'] == 'YES':
        required_options=['SOLUTION_FILENAME', 'SOLUTION_ADJ_FILENAME']
        if not all (opt in config for opt in required_options):
            err = 'RESTART_SOL is set to YES, but the solution is missing:\n'
            for opt in required_options:
                if not opt in config:
                    err += opt + '\n'
            raise ValueError(err)

        os.symlink(os.path.join(base_dir, config.SOLUTION_FILENAME), config.SOLUTION_FILENAME)

        print('\nInitial CFD solution is provided.')

    else:
        print('\nRunning initial CFD solution.')

    #--- Only allow binary restarts since WRT_BINARY_RESTART is deprecated
    sol_ext = '.dat'

    meshfil = config['MESH_FILENAME']
    solfil  = f'restart_flow{sol_ext}'
    su2amg.set_flow_config_ini(config_cfd, solfil)

    try: # run with redirected outputs
        #--- Run a single iteration of the flow if restarting to get history info
        if config['RESTART_SOL'] == 'YES':
            config_cfd.ITER = 1
            config_cfd.RESTART_CFL = 'YES'

        with su2io.redirect.output('su2.out'): SU2_CFD(config_cfd)

        if config['RESTART_SOL'] == 'YES':
            os.remove(solfil)
            os.symlink(os.path.join(base_dir, config.SOLUTION_FILENAME), solfil)

        #--- Set RESTART_SOL=YES for runs after adaptation
        config_cfd.RESTART_SOL = 'YES'
        config_cfd.RESTART_CFL = 'YES'

        if adap_sensor == 'GOAL':
            adjsolfil = f'restart_adj{sol_ext}'
            su2amg.set_adj_config_ini(config_cfd_ad, solfil, adjsolfil, mesh_sizes[0])

            #--- If restarting, check for the existence of an adjoint restart
            if config['RESTART_SOL'] == 'YES':
                adjsolfil_ini = config_cfd_ad.SOLUTION_ADJ_FILENAME
                func_name          = config.OBJECTIVE_FUNCTION
                suffix             = su2io.get_adjointSuffix(func_name)
                adjsolfil_ini = su2io.add_suffix(adjsolfil_ini, suffix)

                #--- Run an adjoint if the solution file doesn't exist
                if not (os.path.exists(os.path.join(base_dir, adjsolfil_ini))):
                    config_cfd_ad.ITER        = config.ITER
                    config_cfd_ad.RESTART_SOL = 'NO'

                    print('Running initial adjoint CFD solution.')

                #--- Otherwise just compute the metric
                else:
                    os.symlink(os.path.join(base_dir, adjsolfil_ini), adjsolfil_ini)
                    config_cfd_ad.ITER = 0

                    print('Initial adjoint CFD solution is provided.')

            else:
                print('Running initial adjoint CFD solution.')

            with su2io.redirect.output('su2.out'): SU2_CFD(config_cfd_ad)

            func_name      = config.OBJECTIVE_FUNCTION
            suffix         = su2io.get_adjointSuffix(func_name)
            adjsolfil = su2io.add_suffix(adjsolfil, suffix)

            #--- Set RESTART_SOL=YES for runs after adaptation
            config_cfd_ad.RESTART_SOL = 'YES'

    except:
        raise

    #--- Check existence of initial mesh, solution

    required_files = [meshfil, solfil]

    if not all (os.path.exists(fil) for fil in required_files):
        err = "Can't find the following files:\n"
        for fil in required_files:
            if not os.path.exists(fil):
                err += fil + '\n'
        raise Exception(err)

    #--- Start adaptive loop

    global_iter = 0

    #--- Print convergence history

    npoin = su2amg.get_su2_npoin(meshfil)
    su2amg.plot_results(history_format, history_filename, global_iter, npoin)

    print('\nStarting mesh adaptation process.\n')

    nSiz = len(mesh_sizes)
    for iSiz in range(nSiz):
        nSub = int(sub_iter[iSiz])
        for iSub in range(nSub):

            global_iter += 1

            #--- Load su2 mesh

            mesh = su2amg.read_mesh_and_sol(meshfil, solfil)

            #--- Write solution
            su2amg.write_mesh_and_sol('flo.meshb', 'flo.solb', mesh)

            mesh_size = int(mesh_sizes[iSiz])
            if iSub == nSub-1 and iSiz != nSiz-1: mesh_size = int(mesh_sizes[iSiz+1])
            config_amg['size'] = mesh_size

            #--- Use pyAmg interface

            if adap_sensor == 'GOAL':

                #--- Use metric computed from SU2 to drive the adaptation

                metric_wrap = su2amg.create_sensor(mesh, adap_sensor)
                mesh['metric'] = metric_wrap['solution']

                #--- Read and merge adjoint solution to be interpolated

                sol_adj = su2amg.read_sol(adjsolfil, mesh)
                su2amg.merge_sol(mesh, sol_adj)

                del sol_adj

            else:

                #--- Create sensor used to drive the adaptation

                sensor_wrap = su2amg.create_sensor(mesh, adap_sensor)
                mesh['sensor'] = sensor_wrap['solution']

            #--- Adapt mesh with AMG

            mesh_new = su2amg.call_pyamg(mesh, config_amg)

            #--- Remove extra files generated by AMG

            extra_files=['back.meshb','meshp3_smoo.meshb','optim.0.meshb','optim.0.solb','subdom.meshb']
            for file in extra_files:
                try:
                    os.remove(file)
                except OSError:
                    pass

            mesh_new['markers'] = mesh['markers']
            mesh_new['dimension'] = mesh['dimension']
            mesh_new['solution_tag'] = mesh['solution_tag']

            del mesh

            #--- Print mesh sizes
            su2amg.print_adap_table(iSiz, mesh_sizes, iSub, nSub, mesh_new)

            dir = f'./ite{global_iter}'
            os.makedirs(os.path.join('..',dir))
            os.chdir(os.path.join('..',dir))

            meshfil = 'mesh_adap.su2'
            solfil  = f'flo{sol_ext}'

            su2amg.write_mesh_and_sol(meshfil, solfil, mesh_new)

            if adap_sensor == 'GOAL':
                adjsolfil = f'adj{sol_ext}'
                sol_adj = su2amg.split_adj_sol(mesh_new)
                su2amg.write_sol(adjsolfil, sol_adj)

            meshfil_gmf    = 'flo_itp.meshb'
            solfil_gmf     = 'flo_itp.solb'
            su2amg.write_mesh_and_sol(meshfil_gmf, solfil_gmf, mesh_new)

            del mesh_new

            if adap_sensor == 'GOAL':
                solfil_gmf_adj = 'adj_itp.solb'
                su2amg.write_sol(solfil_gmf_adj, sol_adj)
                del sol_adj

            #--- Run su2

            try: # run with redirected outputs

                solfil_ini = f'flo_ini{sol_ext}'
                os.rename(solfil, solfil_ini)

                su2amg.update_flow_config(config_cfd, meshfil, solfil, solfil_ini,
                                          flow_iter[iSiz], flow_cfl[iSiz])

                with su2io.redirect.output('su2.out'): SU2_CFD(config_cfd)

                if not os.path.exists(solfil) :
                    raise RuntimeError('SU2_CFD failed.\n')

                #--- Print convergence history

                npoin = su2amg.get_su2_npoin(meshfil)
                su2amg.plot_results(history_format, history_filename, global_iter, npoin)

                if adap_sensor == 'GOAL':

                    adjsolfil_ini = f'adj_ini{sol_ext}'
                    adjsolfil_ini = su2io.add_suffix(adjsolfil_ini, suffix)
                    os.rename(adjsolfil, adjsolfil_ini)
                    adjsolfil_ini = f'adj_ini{sol_ext}'

                    su2amg.update_adj_config(config_cfd_ad, meshfil, solfil, adjsolfil,
                                             adjsolfil_ini, adj_iter[iSiz], mesh_size)

                    with su2io.redirect.output('su2.out'): SU2_CFD(config_cfd_ad)

                    adjsolfil = su2io.add_suffix(adjsolfil, suffix)

                    if not os.path.exists(adjsolfil) :
                        raise RuntimeError('SU2_CFD_AD failed.\n')

            except:
                raise

    #--- Write final files

    mesh = su2amg.read_mesh_and_sol(meshfil, solfil)
    su2amg.write_mesh_and_sol('flo.meshb', 'flo.solb', mesh)

    os.rename(solfil, os.path.join(base_dir, config.RESTART_FILENAME))
    os.rename(meshfil, os.path.join(base_dir, config.MESH_OUT_FILENAME))

    pad_nul = ' '*15
    print('\nMesh adaptation successfully ended.')
    print(f'Results files: {config.MESH_OUT_FILENAME}\n{pad_nul}{config.RESTART_FILENAME}')
