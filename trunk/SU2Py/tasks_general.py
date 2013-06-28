#!/usr/bin/env python 

## \file tasks_general.py
#  \brief Basic class structure for SU2 tasks.
#  \author Trent Lukaczyk, Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
#  \version 2.0.
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

import os, sys, shutil, glob, copy
import numpy
import libSU2


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  GENERAL TASK CLASS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
class General_Task():
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #  ABSTRACTABLE FUNCTIONS
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    # ----------------------------------------------------------------------     
    #  INITIALIZE GENERAL TASK
    # ----------------------------------------------------------------------        
    def __init__( self, config_delta={} ):
        ''' initialize class attributes 
            will start in folder_super
            should not create any new files because 
            the class may be built temporarily 
            as a comparison tester '''        

        self.assets_current = {}
        self.assets_pull    = {}
        self.assets_push    = {}
    
        self.folder_self    = ''
        self.folder_super   = ''
        
        self.config_check   = []
        
        self.config_current = {}
        self.design_current = {}
        
        self.signature      = {}
        
        return
    
    # ----------------------------------------------------------------------     
    #  SET GENERAL TASK SIGNATURE
    # ---------------------------------------------------------------------- 
    def set_signature( self, config_delta={} ):
        ''' assign the properties that represent
            this task's configuration '''
        
        # default signature items
        signature   = {'TYPE':'General_Task'}
        
        # list of keys to copy from a config delta
        config_keys = []
        
        # check in config_delta
        for key in config_keys:
            signature[key] = config_delta[key]
            
        # check in config_current
        for key in self.config_current:
            signature[key] = self.config_current[key]        
        
        return signature
    
    # ----------------------------------------------------------------------     
    #  CHECK GENERAL TASK
    # ----------------------------------------------------------------------     
    def check( self, config_delta, assets_super, design_super ):
        ''' check if job needs to be run or modify super assets ''' 
        
        # modifications to input dictionaries will be by reference
        
        # validate input data, return status True=go
        status = True
        
        ## check if direct solution has been run with requested mesh
        #keys = ['mesh','direct']
        #for key in keys:
            #status = not ( assets_super.has_key(key)        and
                           #self.assets_current.has_key(key) and
                           #assets_super[key] == self.assets_current[key] )        
            #if not status: break
        
        return status   
        
    # ----------------------------------------------------------------------     
    #  CONFIGURE GENERAL TASK
    # ---------------------------------------------------------------------- 
    def configure( self, config_delta, assets_super, design_super ):
        ''' configure any control and data inputs 
            for the analysis task '''
        
        # can also use self.config_current, self.design_current
                
        return
    
    # ----------------------------------------------------------------------     
    #  RUN GENERAL TASK
    # ---------------------------------------------------------------------- 
    def run( self, config_delta, assets_super, design_super ):
        ''' run the analysis task and return an 
            updated design dictionary '''
        
        # can also use self.config_current, self.design_current
        
        # assign design return if needed
        design_return = {}
        
        return design_return

    # ----------------------------------------------------------------------     
    #  PACKUP GENERAL TASK
    # ---------------------------------------------------------------------- 
    def packup( self, config_delta, assets_super, design_super ):
        ''' assign assets to push to super folder 
            and clean uneeded files '''
        
        # can also use self.config_current, self.design_current
        
        # assign assets_return  and config_return if needed
        assets_return = {}
        config_return = {}
        
        return assets_return, config_return
    
    # ----------------------------------------------------------------------     
    #  GENERAL TASK REPRESENTATION
    # ----------------------------------------------------------------------     
    def __repr__( self ):
        return '<Task> General'
    
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #  NON-ABSTRACTABLE FUNCTIONS
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # unless you really want to...
    
    # ----------------------------------------------------------------------     
    #  EVALUATE GENERAL TASK
    # ----------------------------------------------------------------------    
    def evaluate( self, config_delta, assets_super={}, design_super={} ):
        # note: config_delta might be a list
        
        # break pointer references
        config_delta = copy.deepcopy(config_delta)
        assets_super = copy.deepcopy(assets_super)
        design_super = copy.deepcopy(design_super)
            
        # check for self directory
        if not os.path.exists(self.folder_self): 
            os.makedirs(self.folder_self)  # makes intermediate and leaf directories
        
        # move to self directory
        self.folder_super = os.getcwd() # maintains absolute path of super folder at evaluation 
        os.chdir( self.folder_self )    # assumes relative pathing for self folder
        
        # check if needed to run
        if not self.check( config_delta, assets_super, design_super ):
            os.chdir( self.folder_super )
            if not self.config_current['CONSOLE'] == 'QUIET':
                sys.stdout.write('    -skip \n')
            return {}, {}, {}
        
        # pull assets
        self.pull_assets( assets_super )        
                
        # pre configure
        config_delta = self.pre_config( config_delta )
                
        # configure 
        self.configure( config_delta, assets_super, design_super )
        
        # run analysis
        design_return = self.run( config_delta, assets_super, design_super )
        
        # post-run packup
        assets_return, config_return = self.packup( config_delta, assets_super, design_super )
        
        # push assets
        self.push_assets( assets_return )
        
        # return to super directory
        os.chdir( self.folder_super )
        
        # flush output
        sys.stdout.flush()
        
        # done
        return design_return, assets_return, config_return
    
    # ----------------------------------------------------------------------     
    #  PULL GENERAL TASK ASSETS
    # ----------------------------------------------------------------------    
    def pull_assets( self, assets_super ):
        ''' copys files listed in assets_pull to current directory, 
            except for mesh, direct, and adjoint assets, which 
            update the config file with the relative path to the 
            task that holds them''' 
        
        # asset types to point at with filepath (no file copy)
        keys_assets_point = ['mesh','direct','adjoint']
        
        this_assets_pull  = {}
        this_assets_point = {}
        
        # make assets_pull a dictionary for convenience
        if isinstance( self.assets_pull, list ):
            self.assets_pull = dict.fromkeys(self.assets_pull,[])        
        
        # prepare pull asset keys, can only take what super offers
        for key_pull in self.assets_pull.keys():
                        
            # copy asset key and value
            if assets_super.has_key(key_pull):
                # if self.assets_current.has_key(key_pull) :
                    # print 'Warning, %s has %s, pulling %s' % (self.signature['NAME'],self.assets_current[key_pull],assets_super[key_pull])
                # path pointing
                if key_pull in keys_assets_point:
                    this_assets_point[key_pull] = assets_super[key_pull]
                # copying
                else:
                    this_assets_pull[key_pull] = assets_super[key_pull]
                
            # warn if asset not found
            else: 
                # TODO: don't ask for direct if not restarting
                #print 'Super Task did not have asset ' + key_pull  
                pass
                
        #: for each pull key
        
        # --------------------------------------------------------------
        #  ASSETS_PULL
        
        # flatten the asset pull list
        assets_filename = libSU2.flatten_list(this_assets_pull.values())
        assets_filepath = []
                
        # build full file path
        for name in assets_filename:
            assets_filepath.append( os.path.join(self.folder_super,name) )    
        
        # copy files to self folder (should be current folder, could check this)
        for filepath in assets_filepath: shutil.copy( filepath , '.' )
        
        # update current assets dictionary, copied files
        for key,value in this_assets_pull.iteritems():
            # remove all folder paths
            value = os.path.split(value)[-1]
            # store
            self.assets_current[key] = value
            
        # --------------------------------------------------------------
        #  ASSETS_POINT
        
        # build relative back path
        back_path = ""
        this_path = copy.deepcopy( self.folder_self )
        while 1:
            this_path,this_folder = os.path.split(this_path)
            if this_folder:
                if not this_folder == '.':
                    back_path = os.path.join( back_path , '..' )
            else: break
        
        # update current assets dictionary, path pointed files
        for key,value in this_assets_point.iteritems():
            # build relative path
            self.assets_current[key] = os.path.join( back_path , value )
        
        
        # --------------------------------------------------------------
        #  LOAD CONFIG AND DESIGN
        
        # load configuration and design file
        if self.assets_current.has_key('config'):
            self.config_current = libSU2.Get_ConfigParams(self.assets_current['config'])
        if self.assets_current.has_key('design'):
            self.design_current = libSU2.load_data(self.assets_current['design'])   
        
        return


    # ----------------------------------------------------------------------     
    #  PRE-CONFIGURE GENERAL TASK
    # ----------------------------------------------------------------------      
    def pre_config( self, config_delta ):
        ''' assign config updates
            update config filenames, check some conditions 
            mesh asset : check partitions  '''
        
        # -----------------------------------------------------------------
        #  Project Checks
        
        if self.signature['TYPE'] == 'Project':
            # config_delta should be a list
            if not isinstance(config_delta,list):
                config_delta = [config_delta]        
            
            return config_delta
        
        # -----------------------------------------------------------------
        #  Job Checks        
        
        if self.signature['TYPE'] == 'Job':
            # config_delta should be a list
            
            if ( config_delta.has_key('GRADIENTS') and 
                 not isinstance(config_delta['GRADIENTS'],list) ):
                config_delta['GRADIENTS'] = [ config_delta['GRADIENTS']]             
            
            return config_delta        
        
        # else...
        
        # -----------------------------------------------------------------
        #  Task Checks     
            
        config_update = {}
        
        # assign config udpates
        for key in self.config_check:
            if config_delta.has_key(key):
                config_update[key] = config_delta[key]
        
        # check for filenames and assets
        for key in ['mesh','direct','adjoint']:
            if not self.assets_current.has_key(key): continue
            
            asset_name = self.assets_current[key]
            
            # mesh checks
            if key == 'mesh':
                # update mesh filename
                config_update['MESH_FILENAME'] = asset_name
                
                ## check for decomposed mesh        
                #partitions = self.config_current['NUMBER_PART'] 
                #if ( partitions > 1  and 
                     #'DECOMP' in self.config_current['TASKS']
                     #not self.signature['NAME'] in ['Domain Decomposition','Mesh Adaptation'] ):
                    #mesh_supername = libSU2.add_Prefix( asset_name , '%i' )
                    #mesh_partnames = [ mesh_supername % (p+1) for p in range(partitions) ]  
                    
                    ## if partition not found, then no parallel
                    #for mesh_name in mesh_partnames:
                        #if not os.path.exists(mesh_name):
                            #config_update['NUMBER_PART'] = 1
                            #break
                        
            # direct solution checks    
            elif key == 'direct':
                config_update['SOLUTION_FLOW_FILENAME'] = asset_name
        
            # adjoint solution checks
            elif key == 'adjoint':
                config_update['SOLUTION_ADJ_FILENAME'] = asset_name
                
        #: for asset check            
        
        # write any updates
        libSU2.Set_ConfigParams( self.assets_current['config'] , config_update )
        self.config_current.update( config_update )                  
        
        return config_delta
    
    # ----------------------------------------------------------------------     
    #  PUSH GENERAL TASK ASSETS
    # ----------------------------------------------------------------------  
    def push_assets( self, assets_push ):
        ''' moves files listed in assets_push to super directory 
            no relative pathing, mush push file if needed by other tasks '''
        
        #partitions       = self.config_current['NUMBER_PART']
        assets_push = copy.deepcopy(assets_push)
        
        # save design file
        if self.assets_current.has_key('design'):
            libSU2.save_data(self.assets_current['design'],self.design_current)
            
        ## mesh assets, expand if parallel
        #if ( assets_push.has_key('mesh') and 
             #partitions > 1                   and
             #not self.signature['NAME'] in ['Mesh Adaptation'] ):  
            #mesh_supername = libSU2.add_Prefix( assets_push['mesh'] , '%i' )
            #mesh_partnames = [ mesh_supername % (p+1) for p in range(partitions) ]
            #assets_push['mesh'] = mesh_partnames
            
        # adjoint assets, expand by prefix
        if assets_push.has_key('adjoint'):
            adjoint_supername = assets_push['adjoint']
            adjoint_supername = libSU2.add_Prefix(adjoint_supername,'*')
            adjoint_allnames = glob.glob(adjoint_supername)
            assets_push['adjoint'] = adjoint_allnames
        
        # flatten the file list
        assets_filename = libSU2.flatten_list(assets_push.values())
        
        # move files to super folder, overwrite existing files
        for filename in assets_filename: 
            filename_move = os.path.join(self.folder_super,filename)
            if os.path.exists( filename_move ): os.remove( filename_move )
            shutil.move( filename, self.folder_super )
        
        return
    
    # ----------------------------------------------------------------------     
    #  GENERAL TASK BINARY COMPARISON
    # ----------------------------------------------------------------------   
    def __eq__(A,B):
        ''' binary task comparison
            will yield true if two tasks have the same signature
            and false otherwise '''

        A = copy.deepcopy(A.signature)
        B = copy.deepcopy(B.signature)        

        keys_ignore = ['JOB_NUMBER','TASKS','GRADIENTS']
        
        for key in keys_ignore:
            if A.has_key(key):
                del A[key]
            if B.has_key(key):
                del B[key]
        #: for each key_ignore
        
        return A == B


#: class General Task
