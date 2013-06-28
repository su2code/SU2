#!/usr/bin/env python 

## \file patient_designspace.py
#  \brief Python script for running multiple design configurations in multiple sessions
#  \author Trent Lukaczyk, Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
#  \version 2.0.1
#
# Stanford University Unstructured (SU2) Code
# Copyright (C) 2012 Aerospace Design Laboratory
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os, sys, numpy, time, shutil, glob, traceback
from optparse import OptionParser

import libSU2
from tasks_project import Project

# human readable time stamper
pretty_time = lambda: time.asctime( time.localtime(time.time()))



# -------------------------------------------------------------------
#  Main 
# -------------------------------------------------------------------

def main():
    
    # Command Line Options
    parser = OptionParser()
    parser.add_option("-c", "--config",     dest="config_file",
                      help="read config data from FILE", metavar="FILE")
    parser.add_option("-d", "--design",     dest="design_file",
                      help="read design data from FILE", metavar="FILE")    
    parser.add_option("-p", "--project",    dest="project_file",
                      help="read project data from FILE", metavar="FILE")
    parser.add_option("-t", "--transfer",   dest="transfer_file",
                      help="read transfer data from FILE", metavar="FILE")      
    parser.add_option("-e", "--exchange",   dest="exchange_location",  default='',
                      help="optional, SERVER_AND_FOLDER where transfer file can be pulled, example: user@this_server.web:/folder/location/", metavar="SERVER_AND_FOLDER")          
    parser.add_option("-s", "--hot_start",  dest="hot_start",          default="False",
                      help="optional, HOT_START don't initialize design as waiting", metavar="HOT_START")        
    
    (options, args)=parser.parse_args()  
    
    options.hot_start = options.hot_start == 'True'
    
    # Sample Design Space
    patient_designspace ( options.config_file       ,
                          options.design_file       , 
                          options.project_file      ,
                          options.transfer_file     ,
                          options.exchange_location , 
                          options.hot_start          )
    
    return 

#: def main()



# -------------------------------------------------------------------
#  Patiently Sample Design Space 
# -------------------------------------------------------------------
# waits for cfd data request from optimizer

def patient_designspace( config_filename   , design_filename   ,
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
    sys.stdout.write( 'START PATIENT_DESIGNSPACE.PY ... \n' )
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
            sys.stdout.write( 'RUN SAMPLE_DESIGNSPACE.PY ... \n' )
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
                print traceback.format_exc()
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
            
            sys.stdout.write( 'CAUGHT STOP SIGNAL \n')
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

#: def patient_designspace()


# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()
    
    
    
    





