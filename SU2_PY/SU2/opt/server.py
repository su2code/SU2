#!/usr/bin/env python

## \file server.py
#  \brief tools for interfacing with scipy
#  \author T. Lukaczyk
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

import os, sys, numpy, time, shutil, glob, traceback

# human readable time stamper
pretty_time = lambda: time.asctime( time.localtime(time.time()))


# todo:
# new project style



# -------------------------------------------------------------------
#  Patiently Sample Design Space 
# -------------------------------------------------------------------
# waits for cfd data request from optimizer

def server( config_filename   , design_filename   ,
            project_filename  , transfer_filename ,
            exchange_location , hot_start = False  ):    

    # -------------------------------------------------------------------
    #  Start Up
    # -------------------------------------------------------------------
    
    # process exchange location
    if exchange_location:
        exchange_location = exchange_location.split(':')
        server_name = exchange_location[0]
        exchange_folder = ':'.join(exchange_location[1:])
    else:
        server_name = ''
        exchange_folder = ''
    
    # check for existing project file
    if os.path.exists(project_filename):
        
        # load project 
        The_Project = libSU2.load_data(project_filename)
        The_Project.folder_self = os.getcwd()        
        
    # or start new project
    else:
        # new design data
        design_init = { 'VARIABLES'  : [] ,
                        'OBJECTIVES' : {} ,
                        'GRADIENTS'  : {}  }
        libSU2.save_data(design_filename,design_init,append=False)
        # start project
        The_Project = Project( config_name = config_filename ,
                               design_name = design_filename  )        
    #: if load/start project   
    
    # make sure to start with waiting
    if not hot_start:
        if server_name:
            os.system('scp -q %s:%s%s ./ ' % (server_name,exchange_folder,transfer_filename) )
            
        Status_set = {'STATUS' : 'WAIT'}
        with libSU2.FileLock(transfer_filename,timeout=100):
            libSU2.save_data(transfer_filename,Status_set,append=True)
            
        if server_name:
            os.system('scp -q %s %s:%s' % (transfer_filename,server_name,exchange_folder) )
    #: if not hot_start    
    
    
    # start log 
    sys.stdout.write( 'Start Server ... \n' )
    sys.stdout.write( pretty_time() + '\n' )
    sys.stdout.write(' \n')
    sys.stdout.flush()
    
    # -------------------------------------------------------------------
    #  Listen for Jobs
    # -------------------------------------------------------------------    
    
    # keep on keepin on
    keepon = True
    while keepon:
        
        # get design data
        if server_name:
            os.system('scp -q %s:%s%s ./ ' % (server_name,exchange_folder,transfer_filename) )
        
        # load design data
        with libSU2.FileLock(transfer_filename,timeout=100):
            Transfer_Data = libSU2.load_data(transfer_filename)
        Status = str( Transfer_Data['STATUS'] )
        
        # -------------------------------------------------------------------
        # Run Case
        if Status == 'RUN':
               
            # log
            sys.stdout.write( 'Run Project ... \n' )
            sys.stdout.write( pretty_time() + '\n' )
            sys.stdout.write( ' \n' )
            sys.stdout.flush()
            
            # setup config deltas
            config_delta = []
            for DV_X in Transfer_Data['VARIABLES']:
                config_delta.append( {'VARIABLES':DV_X} )            
            
            # RUN PROJECT
            try:
                
                # evaluate project
                Transfer_Data,_,_ = The_Project.evaluate(config_delta)
                
                # save project
                libSU2.save_data(project_filename,The_Project)
                
                # save project data
                libSU2.save_data( design_filename, The_Project.design_current, append=False )                
                
            except (KeyboardInterrupt, SystemExit):
                raise
            except Exception,err:
                sys.stdout.write( 'RUN FAILED \n\n' )
                print(traceback.format_exc())
                sys.stdout.write( '\n' )
            #: try SDS.py
            
            # finish up
            sys.stdout.write( pretty_time() + '\n\n' )
            sys.stdout.flush()
            
            # save new transfer data
            Transfer_Data['STATUS'] = 'WAIT'
            with libSU2.FileLock(transfer_filename,timeout=100):
                libSU2.save_data(transfer_filename,Transfer_Data,append=False)
            
            # push design data
            if server_name:
                os.system('scp -q %s %s:%s' % (transfer_filename,server_name,exchange_folder) )
            
        # -------------------------------------------------------------------
        # Stop Case
        elif Status == 'STOP':
            
            sys.stdout.write( 'Caught Stop Signal \n')
            sys.stdout.write( pretty_time() + '\n')
            sys.stdout.write( ' \n')
            sys.stdout.flush()
            
            keepon = False
            
        # -------------------------------------------------------------------
        # Sleep Case
        else:
            
            time.sleep(10)
        
        #: if use case
        
    #: while keepon
        
    sys.stdout.write( 'DONE \n\n' )
    
    return

#: def server()

