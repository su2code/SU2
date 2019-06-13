#!/usr/bin/env python

## \file adjoint.py
#  \brief python package for running mesh adaptation using the AMG Inria library
#  \author Victorien Menier
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

import os, sys, shutil, copy, time
import numpy as np

from .. import io   as su2io
from .. import amginria as su2amg
from interface import CFD as SU2_CFD
from interface import MET as SU2_MET
import _amgio as amgio

def amg ( config , kind='' ):
    
    sys.stdout.write("SU2-AMG Anisotropic Mesh Adaptation\n")
        
    #--- Check config options related to mesh adaptation
    
    adap_options = ['ADAP_SIZES', 'ADAP_SUBITE', 'ADAP_SENSOR', \
    'ADAP_BACK', 'ADAP_HMAX', 'ADAP_HMIN', 'ADAP_HGRAD', 'ADAP_RESIDUAL_REDUCTION', 'ADAP_EXT_ITER', 'ADAP_SOURCE','ADAP_PYTHON']
    required_options = ['ADAP_SIZES', 'ADAP_SUBITE', \
    'ADAP_SENSOR', 'MESH_FILENAME', 'RESTART_SOL', 'MESH_OUT_FILENAME']
    
    if not all (opt in config for opt in required_options):
        err = '\n\n## ERROR : Missing options: \n'
        for opt in required_options:
            if not opt in config:
                err += opt + '\n'
        raise RuntimeError , err
    
    # Print adap options
    sys.stdout.write(su2amg.print_adap_options(config, adap_options))
    
    #--- How many iterative loops? Using what prescribed mesh sizes? 
    
    mesh_sizes   = su2amg.get_mesh_sizes(config)
    sub_iter     = su2amg.get_sub_iterations(config)
    
    # solver iterations/ residual reduction param for each size level
    adap_ext_iter = su2amg.get_ext_iter(config)
    adap_res = su2amg.get_residual_reduction(config)

    adap_sensor = config.ADAP_SENSOR
    sensor_avail = ['MACH', 'PRES', 'MACH_PRES', 'GOAL']
    
    if adap_sensor not in sensor_avail:
        raise RuntimeError , 'Unknown adaptation sensor (ADAP_SENSOR option)\n'
        
    if len(mesh_sizes) != len(sub_iter):
        raise RuntimeError , 'Inconsistent number of mesh sizes and sub-iterations'
        
        
    #--- Use the python interface to amg, or the executable?
    
    amg_python = su2amg.get_python_amg(config)
    
    #--- Change current directory
    
    warn = False
    adap_dir = './ADAP'
    cwd = os.getcwd()
        
    if os.path.exists(adap_dir):
        sys.stdout.write('./ADAP exists. Removing old mesh adaptation in 10s.\n')
        sys.stdout.flush()
        if warn : time.sleep(10)
        shutil.rmtree(adap_dir)
    
    os.makedirs(adap_dir)
    os.chdir(adap_dir)
    sys.stdout.write('The %s folder was deleted\n' % adap_dir)
    
    os.symlink(os.path.join(cwd, config.MESH_FILENAME), config.MESH_FILENAME)
    os.symlink(os.path.join(cwd, config.SOLUTION_FLOW_FILENAME), config.SOLUTION_FLOW_FILENAME)
    
    #--- Compute initial solution if needed, else link current files
    
    config_cfd = copy.deepcopy(config)
    for opt in adap_options:
        config_cfd.pop(opt, None)
    config_cfd.LOW_MEMORY_OUTPUT = "NO"
    
    config_cfd.WRT_BINARY_RESTART  = "NO"
    config_cfd.READ_BINARY_RESTART = "NO"
        
    current_mesh     = "Initial_mesh"
    current_solution = "Initial_solution"
        
    if config['RESTART_SOL'] == 'NO':
        
        stdout_hdl = open('ini.out','w') # new targets
        stderr_hdl = open('ini.err','w')
        
        success = False
        val_out = [False]
        
        sys.stdout.write('Running initial CFD solution.\n')
        
        try: # run with redirected outputs
            
            sav_stdout, sys.stdout = sys.stdout, stdout_hdl 
            sav_stderr, sys.stderr = sys.stderr, stderr_hdl
        
            current_mesh     = config['MESH_FILENAME']
            current_solution = "ini_restart_flow.dat"
            
            config_cfd.CONV_FILENAME         = "ini_history"
            config_cfd.RESTART_FLOW_FILENAME = current_solution
            config_cfd.MATH_PROBLEM           = 'DIRECT'
            
            SU2_CFD(config_cfd)
                        
            if adap_sensor == 'GOAL':
                current_solution_adj    = "ini_restart_adj.dat"

                config_cfd.CONV_FILENAME          = "ini_history_adj"
                config_cfd.RESTART_ADJ_FILENAME   = current_solution_adj
                config_cfd.SOLUTION_FLOW_FILENAME = current_solution
                config_cfd.MATH_PROBLEM           = 'DISCRETE_ADJOINT'
                SU2_CFD(config_cfd)

                config_cfd.SOLUTION_ADJ_FILENAME  = current_solution_adj
                config_cfd.ERROR_ESTIMATE         = 'YES'
                config_cfd.MESH_HMAX              = config.ADAP_HMAX
                config_cfd.MESH_HMIN              = config.ADAP_HMIN
                config_cfd.MESH_COMPLEXITY        = int(mesh_sizes[0])
                SU2_MET(config_cfd)

        except:
            sys.stdout = sav_stdout
            sys.stderr = sav_stderr
            raise
        
        sys.stdout = sav_stdout
        sys.stderr = sav_stderr
        
    else:
        required_options=['SOLUTION_FLOW_FILENAME']
        if not all (opt in config for opt in required_options):
            err = '\n\n## ERROR : RESTART_SOL is set to YES, but the solution is missing:\n'
            for opt in required_options:
                if not opt in config:
                    err += opt + '\n'
            raise RuntimeError , err
        
        current_mesh     = config['MESH_FILENAME']
        current_solution = config['SOLUTION_FLOW_FILENAME']
        
        sys.stdout.write('Initial CFD solution is provided.\n')
        
    #--- Check existence of initial mesh, solution
    
    required_files = [current_mesh,current_solution]
    
    if not all (os.path.exists(fil) for fil in required_files):
        err = '\n\n## ERROR : Can\'t find:\n'
        for fil in required_files:
            if not os.path.exists(fil):
                err += fil + '\n'
        raise RuntimeError , err
    
    #--- Start looping
    
    # Get mesh dimension
    dim = su2amg.get_su2_dim(current_mesh)
    if ( dim != 2 and dim != 3 ):
        raise RuntimeError , "Wrong dimension number\n"
    
    #--- AMG parameters
    
    config_amg = dict()
    
    config_amg['hgrad']       = float(config['ADAP_HGRAD'])
    config_amg['hmax']        = float(config['ADAP_HMAX'])
    config_amg['hmin']        = float(config['ADAP_HMIN'])
    config_amg['mesh_in']     = 'current.meshb'
    config_amg['mesh_out']    = 'current.new.meshb'
    config_amg['metric_in']   = ''
    config_amg['sol_in']      = 'current_sensor.solb'
    config_amg['itp_sol_in']  = 'current.solb'
    config_amg['adap_source'] = ''
    
    if 'ADAP_BACK' in config:
        config_amg['adap_back'] = os.path.join(cwd,config['ADAP_BACK'])
        #os.symlink(os.path.join(cwd, config.ADAP_BACK), config.ADAP_BACK)
    else:
        config_amg['adap_back'] = config['MESH_FILENAME']
    
    back_name, back_extension = os.path.splitext(config_amg['adap_back'])
    
    if not os.path.exists(config_amg['adap_back']):
        raise RuntimeError , "\n\n##ERROR : Can't find back mesh: %s.\n\n" % config_amg['adap_back']
    
    if back_extension == ".su2":
        amgio.py_ConvertSU2toInria(config_amg['adap_back'], "", "amg_back")
        config_amg['adap_back'] = "amg_back.meshb"
    
    if 'ADAP_SOURCE' in config:
        config_amg['adap_source'] = os.path.join(cwd,config['ADAP_SOURCE'])
    
    global_iter = 0
    
    sys.stdout.write("\nStarting mesh adaptation process.\n")
    
    for iSiz in range(len(mesh_sizes)):
        
        mesh_size   = int(mesh_sizes[iSiz])
        nSub        = int(sub_iter[iSiz])
                        
        sys.stdout.write("\nIteration %d - Mesh size coefficient %.1lf\n" % (iSiz, mesh_size))
        
        for iSub in range(nSub):
            
            config_amg['size']        = mesh_size
            config_amg['amg_log']     = 'ite%d.amg.out' % (global_iter)
            
            # Prints
            pad_cpt = ("(%d/%d)" % (iSub+1, nSub)).ljust(9)
            pad_nul = "".ljust(9)
            
            #--- Load su2 mesh 
            
            mesh = su2amg.read_mesh(current_mesh, current_solution)
                                    
            if not amg_python : 
                
                #--- If not using the amg python interface, convert the mesh and make system call
                
                su2amg.write_mesh("current.meshb", "current.solb", mesh)
                
                if not os.path.exists("current.solb"):
                    raise RuntimeError , "\n##ERROR : Can't find solution.\n"
                if not os.path.exists("current.meshb"):
                    raise RuntimeError , "\n##ERROR : Can't find mesh.\n"
                
                #--- Get sensor
                                
                sensor = su2amg.create_sensor(mesh, adap_sensor)
                su2amg.write_solution("current_sensor.solb", sensor)
                
                if not os.path.exists("current_sensor.solb"):
                    raise RuntimeError , "\n##ERROR : Can't find adap sensor.\n"
                
                #--- Run amg
                                
                sys.stdout.write("Running amg. Log : %s\n" % config_amg['amg_log'])
                
                if os.path.exists("current.itp.solb"):
                    os.remove("current.itp.solb")
                            
                try :
                    su2amg.amg_call(config_amg)
                except:
                    raise RuntimeError , "\n##ERROR : Call to AMG failed.\n"
                
                if not os.path.exists(config_amg['mesh_out']):
                    raise RuntimeError , "\n##ERROR : Mesh adaptation failed.\n"
                
                if not os.path.exists("current.itp.solb"):
                    raise RuntimeError , "\n##ERROR AMG: Solution interpolation failed.\n"            
                
                #--- Convert output from Inria mesh format to su2
                # Deal with markers
                
                save_markers = mesh['markers']
                del mesh
                
                # Read Inria mesh
                mesh = su2amg.read_mesh(config_amg['mesh_out'], "current.itp.solb")
                mesh['markers'] = save_markers
                
                current_mesh = "ite%d.su2" % global_iter
                current_solution = "ite%d.dat" % global_iter    
                
                su2amg.write_mesh(current_mesh, current_solution, mesh)
                
                if not os.path.exists(current_mesh) or not os.path.exists(current_solution) :
                    raise RuntimeError , "\n##ERROR : Conversion to SU2 failed.\n"

            else :
                
                #--- Use pyAmg interface
                
                try :
                    import pyamg 
                except:
                    sys.stderr.write("## ERROR : Unable to import pyamg module.\n")
                    sys.exit(1)
                
                if adap_sensor == 'GOAL':

                    #--- Use metric computed from SU2 to drive the adaptation

                    config_amg['metric_in']   = 'current.metr.solb'
                    config_amg['sol_in']      = ''
                    config_amg['itp_sol_in']  = 'current.restart.solb'

                    # mesh.pop('sensor',None)

                    sys.stdout.write(' %s Generating adapted mesh using AMG\n' % pad_cpt)
                    
                    # mesh_new = su2amg.amg_call_python(mesh, config_amg)
                    #--- For now, just use executable
                    try :
                        su2amg.amg_call_met(config_amg)
                    except:
                        raise RuntimeError , "\n##ERROR : Call to AMG failed.\n"
                    
                    if not os.path.exists(config_amg['mesh_out']):
                        raise RuntimeError , "\n##ERROR : Mesh adaptation failed.\n"
                    
                    if not os.path.exists("current.restart.itp.solb"):
                        raise RuntimeError , "\n##ERROR AMG: Solution interpolation failed.\n"

                    #--- Convert output from Inria mesh format to su2
                    # Deal with markers
                    
                    save_markers = mesh['markers']
                    del mesh
                    
                    # Read Inria mesh
                    mesh = su2amg.read_mesh(config_amg['mesh_out'], "current.restart.itp.solb")
                    mesh['markers'] = save_markers
                    
                    current_mesh = "ite%d.su2" % global_iter
                    current_solution = "ite%d.dat" % global_iter    
                    
                    su2amg.write_mesh(current_mesh, current_solution, mesh)
                    
                    if not os.path.exists(current_mesh) or not os.path.exists(current_solution) :
                        raise RuntimeError , "\n##ERROR : Conversion to SU2 failed.\n"

                else:
                
                    #--- Create sensor used to drive the adaptation
                    
                    sensor_wrap = su2amg.create_sensor(mesh, adap_sensor)
                    
                    mesh['sensor'] = sensor_wrap['solution']
                    
                    sys.stdout.write(' %s Generating adapted mesh using AMG\n' % pad_cpt)
                    
                    mesh_new = su2amg.amg_call_python(mesh, config_amg)
                                    
                    #--- print mesh size
                    
                    sys.stdout.write(' %s AMG done: %s\n' % (pad_nul, su2amg.return_mesh_size(mesh_new)))
                                    
                    mesh_new['markers'] = mesh['markers']
                    mesh_new['dimension'] = mesh['dimension']
                    
                    current_mesh = "ite%d.su2" % global_iter
                    current_solution = "ite%d.dat" % global_iter
                                    
                    su2amg.write_mesh(current_mesh, current_solution, mesh_new)
                
            #--- Run su2
            
            log = 'ite%d.SU2.'%global_iter
            stdout_hdl = open('%sout'%log,'w') # new targets
            stderr_hdl = open('%serr'%log,'w')
            
            success = False
            val_out = [False]
            
            sys.stdout.write(' %s Running CFD\n' % pad_nul)
        
            try: # run with redirected outputs
            
                sav_stdout, sys.stdout = sys.stdout, stdout_hdl 
                sav_stderr, sys.stderr = sys.stderr, stderr_hdl
                
                current_solution_ini = "ite%d_ini.dat" % global_iter
                os.rename(current_solution, current_solution_ini)
                
                config_cfd.MESH_FILENAME          = current_mesh
                config_cfd.CONV_FILENAME          = "ite%d_history" % global_iter
                config_cfd.SOLUTION_FLOW_FILENAME = current_solution_ini
                config_cfd.RESTART_FLOW_FILENAME  = current_solution
                config_cfd.MATH_PROBLEM           = 'DIRECT'
                config_cfd.RESTART_SOL            = 'YES'
                
                config_cfd.RESIDUAL_REDUCTION = float(adap_res[iSiz])
                config_cfd.EXT_ITER = int(adap_ext_iter[iSiz])
                
                config_cfd.WRT_BINARY_RESTART  = "NO"
                config_cfd.READ_BINARY_RESTART = "NO"
                
                SU2_CFD(config_cfd)
                
                if not os.path.exists(current_solution) :
                    raise RuntimeError , "\n##ERROR : SU2_CFD Failed.\n"
                    
                if adap_sensor == 'GOAL':
                    current_solution_adj = "ite%d_adj.dat" % global_iter

                    config_cfd.CONV_FILENAME          = "ite%d_history_adj" % global_iter
                    config_cfd.RESTART_ADJ_FILENAME   = current_solution_adj
                    config_cfd.SOLUTION_FLOW_FILENAME = current_solution
                    config_cfd.MATH_PROBLEM           = 'DISCRETE_ADJOINT'
                    config_cfd.RESTART_SOL            = 'NO'
                    SU2_CFD(config_cfd)

                    config_cfd.SOLUTION_ADJ_FILENAME  = current_solution_adj
                    config_cfd.ERROR_ESTIMATE         = 'YES'
                    config_cfd.MESH_COMPLEXITY        = int(mesh_sizes[iSiz])
                    SU2_MET(config_cfd)
            
            except:
                sys.stdout = sav_stdout
                sys.stderr = sav_stderr
                raise
            
            sys.stdout = sav_stdout
            sys.stderr = sav_stderr
            
                    
            #--- Print convergence history
            
            plot_format      = config_cfd['OUTPUT_FORMAT']
            plot_extension   = su2io.get_extension(plot_format)
            history_filename = "ite%d_history" % global_iter + plot_extension
            
            history = su2io.read_history(history_filename)
            
            res_flow = history['Res_Flow[0]']
            res_cvg = max(res_flow)-min(res_flow)
            
            del history
            
            sys.stdout.write(' %s CFD done. Residual convergence %.2lf orders of magnitude\n' % (pad_nul, res_cvg))
            
            
            to_remove = ["current.itp.solb", config_amg['mesh_out'], config_amg['sol_in']]
            for fil in to_remove:
                if os.path.exists(fil) : os.remove(fil)
            
            global_iter += 1
    
    os.rename(current_solution,os.path.join(cwd,config.RESTART_FLOW_FILENAME))
    os.rename(current_mesh,os.path.join(cwd,config.MESH_OUT_FILENAME))
    
    sys.stdout.write("\nMesh adaptation successfully ended. Results files:\n")
    sys.stdout.write("%s\n%s\n\n" % (config.MESH_OUT_FILENAME,config.RESTART_FLOW_FILENAME))
    
