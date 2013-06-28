#!/usr/bin/env python 

## \file tasks_SU2.py
#  \brief Python classes and functions for evaluating designs in SU2
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


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  SU2 TASKS MODULE
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import os, sys, shutil, glob, copy
import numpy
import libSU2
from tasks_general      import General_Task
from finite_differences import finite_differences
from mesh_adaptation    import mesh_adaptation
from filter_adjoint     import process_surface_adjoint

SU2_RUN = os.environ['SU2_RUN'] 
sys.path.append( SU2_RUN )

# TODO:
#  mesh decomposition runs at every direct and adjoint solution, fix


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  HIGH LEVEL FUNCTIONS FOR OPTIMIZERS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# -------------------------------------------------------------------
#  OBJECTIVE FUNCTION
# -------------------------------------------------------------------
def eval_f( Variables, The_Project ):
    #print( 'f   ') 
    
    Objective = The_Project.config_current['OPT_OBJFUNC']['OBJECTIVE']
    Scale     = The_Project.config_current['OPT_OBJFUNC']['SCALE']
    Sign      = libSU2.get_ObjFunSign(Objective)
    
    config_delta = {}
    
    config_delta['VARIABLES'] = Variables
    config_delta['TASKS']     = ['DEFORM', 'DIRECT']
    
    design_new,_,_ = The_Project.evaluate(config_delta)
    sys.stdout.flush()

    Value = design_new['OBJECTIVES'][Objective][0][0] * Scale * Sign
    
    return Value

#: def eval_f()

# -------------------------------------------------------------------
#  OBJECTIVE FUNCTION GRADIENT
# -------------------------------------------------------------------
def eval_df( Variables, The_Project ):
    #print( 'df  ')
    
    Objective = The_Project.config_current['OPT_OBJFUNC']['OBJECTIVE']
    Scale     = The_Project.config_current['OPT_OBJFUNC']['SCALE']
    Sign      = libSU2.get_ObjFunSign(Objective)
    
    gradient_type = 'CONT_ADJOINT'
    if 'FINITE_DIFF' in The_Project.config_current['TASKS']:
        gradient_type = 'FINITE_DIFF'
    Tasks = ['DEFORM',gradient_type]
        
    config_delta = {}
    
    config_delta['VARIABLES'] = Variables
    config_delta['TASKS']     = Tasks
    config_delta['GRADIENTS'] = Objective
    
    design_new,_,_ = The_Project.evaluate(config_delta)
    sys.stdout.flush()
    
    Values = numpy.array( design_new['GRADIENTS'][Objective][0] ) * Scale * Sign

    return Values

#: def eval_df()

# -------------------------------------------------------------------
#  EQUALITY CONSTRAINT FUNCTIONS
# -------------------------------------------------------------------
def eval_ceq( Variables, The_Project ):
    #print( 'c   ')
    
    Con_Definition = The_Project.config_current['OPT_CONSTR']['EQUALITY']
    
    if not Con_Definition:
        return numpy.zeros(0)
    
    Constraints = Con_Definition.keys()
    Scales      = [ x['SCALE'] for x in Con_Definition.values() ]
    Limits      = [ x['VALUE'] for x in Con_Definition.values() ]
    
    config_delta = {}
    
    config_delta['VARIABLES'] = Variables
    config_delta['TASKS']     = ['DEFORM', 'DIRECT']
    
    design_new,_,_ = The_Project.evaluate(config_delta)
    sys.stdout.flush()
    
    Values = numpy.zeros([ len(Constraints) , 1 ])
    for i, this_con in enumerate(Constraints):
        Values[i] = ( design_new['OBJECTIVES'][this_con][0][0] - Limits[i] ) * Scales[i]

    return Values

#: def eval_ceq()

# -------------------------------------------------------------------
#  EQUALITY CONSTRAINT FUNCTION GRADIENTS
# -------------------------------------------------------------------
def eval_dceq( Variables, The_Project ):
    #print( 'dc  ')
    
    Con_Definition = The_Project.config_current['OPT_CONSTR']['EQUALITY']
    
    if not Con_Definition:
        return numpy.zeros([0,len(Variables)])    
    
    Constraints = Con_Definition.keys()
    Scales      = [ x['SCALE'] for x in Con_Definition.values() ]
    
    gradient_type = 'CONT_ADJOINT'
    if 'FINITE_DIFF' in The_Project.config_current['TASKS']:
        gradient_type = 'FINITE_DIFF'
        
    Tasks = ['DEFORM',gradient_type]    
    
    config_delta = {}
    
    config_delta['VARIABLES'] = Variables
    config_delta['TASKS']     = Tasks
    config_delta['GRADIENTS'] = Constraints
    
    design_new,_,_ = The_Project.evaluate(config_delta)
    sys.stdout.flush()
    
    Values = numpy.zeros([ len(Constraints) , len(Variables) ])
    for i, this_con in enumerate(Constraints):
        Values[i,:] = numpy.array( design_new['GRADIENTS'][this_con][0] ) * Scales[i]

    return Values

#: def eval_dceq

# -------------------------------------------------------------------
#  INEQUALITY CONSTRAINT FUNCTIONS
# -------------------------------------------------------------------
def eval_cieq( Variables, The_Project ):
    #print( 'c   ')
    
    Con_Definition = The_Project.config_current['OPT_CONSTR']['INEQUALITY']
    
    if not Con_Definition:
        return numpy.zeros(0)    
    
    Constraints = Con_Definition.keys()
    Scales      = [ x['SCALE'] for x in Con_Definition.values() ]
    Signs       = [ libSU2.get_ConSign( x['SIGN'] ) for x in Con_Definition.values() ]
    Limits      = [ x['VALUE'] for x in Con_Definition.values() ]
    
    config_delta = {}
    
    config_delta['VARIABLES'] = Variables
    config_delta['TASKS']     = ['DEFORM', 'DIRECT']
    
    design_new,_,_ = The_Project.evaluate(config_delta)
    sys.stdout.flush()
    
    Values = numpy.zeros([ len(Constraints) , 1 ])
    for i, this_con in enumerate(Constraints):
        Values[i] = ( design_new['OBJECTIVES'][this_con][0][0] - Limits[i] ) * Scales[i] * Signs[i]

    return Values   

#: def eval_cieq()

# -------------------------------------------------------------------
#  INEQUALITY CONSTRAINT FUNCTION GRADIENTS
# -------------------------------------------------------------------
def eval_dcieq( Variables, The_Project ):
    #print( 'dic ')

    Con_Definition = The_Project.config_current['OPT_CONSTR']['INEQUALITY']
    
    if not Con_Definition:
        return numpy.zeros([0,len(Variables)])        
    
    Constraints = Con_Definition.keys()
    Scales      = [ x['SCALE'] for x in Con_Definition.values() ]
    Signs       = [ libSU2.get_ConSign( x['SIGN'] ) for x in Con_Definition.values() ]
    
    gradient_type = 'CONT_ADJOINT'
    if 'FINITE_DIFF' in The_Project.config_current['TASKS']:
        gradient_type = 'FINITE_DIFF'
        
    Tasks = ['DEFORM',gradient_type]
    
    config_delta = {}
    
    config_delta['VARIABLES'] = Variables
    config_delta['TASKS']     = Tasks
    config_delta['GRADIENTS'] = Constraints
    
    design_new,_,_ = The_Project.evaluate(config_delta)
    sys.stdout.flush()
    
    Values = numpy.zeros([ len(Constraints) , len(Variables) ])
    for i, this_con in enumerate(Constraints):
        Values[i,:] = numpy.array( design_new['GRADIENTS'][this_con][0] ) * Scales[i] * Signs[i]

    return Values

#: def eval_dcieq()



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  CLASS MAPPING
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    
# -------------------------------------------------------------------
#  TASK DEPENDENCY
# -------------------------------------------------------------------   
def get_dependency(config_current):
    ''' returns dictionary of list of constructors needed to run task name ''' 
    
    # base dependency
    task_dependency = { 'DECOMP'       : [] ,
                        'DEFORM'       : [] ,
                        'ADAPT'        : [] ,
                        'DIRECT'       : [] ,
                        'CONT_ADJOINT' : [Direct],
                        'FINITE_DIFF'  : [Direct] }

    ## additional cases
    #if config_current['NUMBER_PART'] > 1:
        #task_dependency['DEFORM'].append( SU2_Tasks.Decomp )
        #task_dependency['DIRECT'].append( SU2_Tasks.Decomp )
        #task_dependency['CONT_ADJOINT'].append( SU2_Tasks.Decomp )
        #pass
    
    return task_dependency

#: def get_dependency()
      
# -------------------------------------------------------------------
#  TASK MAPS
# ------------------------------------------------------------------- 
def get_classmap():
    ''' returns dictionary of constructors for task string keys '''
    
    classmap = { 'DECOMP'       : Decomp ,
                 'DEFORM'       : Deform ,
                 'ADAPT'        : Adapt  ,
                 'DIRECT'       : Direct ,
                 'CONT_ADJOINT' : Multiple_Cont_Adjoint ,
                 'FINITE_DIFF'  : Finite_Diff }      
    
    return classmap

#: def get_classmap()
            
            
            
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  SU2 TASKS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
            
            
#   % % %     % % %     % % %     % % %     % % %     % % %     % % %   #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  SU2 MESH DECOMPOSITION TASK
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class Decomp(General_Task): 
    
    # -------------------------------------------------------------------
    #  INITIALIZE DECOMPOSITION TASK
    # -------------------------------------------------------------------
    def __init__( self, assets_super, config_delta ):

        General_Task.__init__(self)

        self.config_current = libSU2.Get_ConfigParams( assets_super['config'] )
        self.config_current.update(config_delta)                      

        self.assets_pull   = ['config','mesh']
        
        self.config_check  = ['NUMBER_PART']
                
        self.folder_self   = 'DECOMP'
        
        self.log_filename  = 'log_Decomp.out'
        
        self.signature     = self.set_signature()
        
        return
        
    # -------------------------------------------------------------------
    #  SET DECOMPOSITION TASK SIGNATURE
    # -------------------------------------------------------------------
    def set_signature( self ):
        ''' assign the properties that represent
            this task's configuration '''
        
        # default signature items
        signature   = {'TYPE' : 'SU2 Task',
                       'NAME' : 'Domain Decomposition'}
        
        # list of keys to copy from a config delta
        config_keys = ['NUMBER_PART']
        
        # assign config signature items
        for key in config_keys:
            signature[key] = self.config_current[key]
        
        return signature
    
    # ----------------------------------------------------------------------     
    #  CHECK DECOMPOSITION TASK
    # ----------------------------------------------------------------------     
    def check( self, config_delta, assets_super, design_super ):
        ''' check if job needs to be run or modify super assets ''' 
        
        # check if the mesh has already been decomposed
        key = 'mesh'
        status = not ( self.assets_current.has_key(key) )
        
        return status    
    
    # -------------------------------------------------------------------
    #  CONFIGURE DECOMPOSITION TASK
    # -------------------------------------------------------------------
    def configure( self, config_delta, assets_super, design_super ):
        ''' Configure SU2_DDC '''
        
        # nothing special to do
        # reqiures NUMBER_PART to be set in project's config file
        
        return
        
    # -------------------------------------------------------------------
    #  RUN DECOMPOSITION TASK
    # -------------------------------------------------------------------
    def run( self, config_delta, assets_super, design_super ):
        ''' Run SU2_DDC '''
        
        # setup run command
        run_Command = "SU2_DDC " + self.assets_current['config']
        
        if self.config_current['CONSOLE'] in ['QUIET','CONCISE']:
            run_Command = run_Command + ' >> ' + self.log_filename
            
        # run command TODO: error checking
        os.system(run_Command)
        
        # no design to return
        design_return = {}
        
        return design_return
    
    # -------------------------------------------------------------------
    #  PACKUP DECOMPOSITION TASK
    # -------------------------------------------------------------------
    def packup( self, config_delta, assets_super, design_super ):
        ''' assign asset push decomposed mesh '''
        
        # mesh is decomposed in place, no assets to push or update
        assets_push = {}

        # remove from pull list for future runs
        del self.assets_pull['mesh']        
        
        return assets_push    
    
    # -------------------------------------------------------------------
    #  DECOMPOSITION TASK REPRESENTATION
    # -------------------------------------------------------------------
    def __repr__ (self):
        return '<Task> SU2 Decomp'
        
#: class Decomp()


#   % % %     % % %     % % %     % % %     % % %     % % %     % % %   #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  SU2 MESH DEFORMATION TASK
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
class Deform(General_Task): 
    
    # -------------------------------------------------------------------
    #  INITIALIZE DEFORMATION TASK
    # -------------------------------------------------------------------
    def __init__( self, assets_super, config_delta ):

        General_Task.__init__(self)
        
        self.config_current = libSU2.Get_ConfigParams( assets_super['config'] )
        self.config_current.update(config_delta)               

        self.assets_pull  = ['config','mesh']
        
        self.config_check = []
                
        self.folder_self  = 'DEFORM'
        
        self.log_filename = 'log_Deform.out'
        
        self.signature    = self.set_signature()
        
        return
        
    # -------------------------------------------------------------------
    #  SET DEFORMATION TASK SIGNATURE
    # -------------------------------------------------------------------
    def set_signature( self ):
        ''' assign the properties that represent
            this task's configuration '''
                    
        # default signature items
        signature = { 'TYPE' : 'SU2 Task',
                      'NAME' : 'Mesh Deformation' }
        
        # list of keys to copy from a config delta
        config_keys = ['DV_VALUE_OLD','DV_VALUE_NEW']
        
        # assign config signature items
        for key in config_keys:
            signature[key] = self.config_current[key]
        
        return signature         

    # ----------------------------------------------------------------------     
    #  CHECK DEFORMATION TASK
    # ----------------------------------------------------------------------     
    def check( self, config_delta, assets_super, design_super ):
        ''' check if job needs to be run or modify super assets ''' 
        
        # check if the mesh has already been deformed
        key = 'mesh'
        status = not ( self.assets_current.has_key(key) )
        
        return status    

    # -------------------------------------------------------------------
    #  CONFIGURE DEFORMATION TASK
    # -------------------------------------------------------------------
    def configure( self, config_delta, assets_super, design_super ):
        ''' Configure SU2_DDC  '''
        
        # nothing special to do
        # needs configs DV_VALUE_OLD and DV_VALUE_NEW
        
        return
        
    # -------------------------------------------------------------------
    #  RUN DEFORMATION TASK
    # -------------------------------------------------------------------
    def run( self, config_delta, assets_super, design_super ):
        ''' Run MDC '''
        
        config_name = self.assets_current['config']
        partitions  = self.config_current['NUMBER_PART']
        
        # prepare command
        #if partitions == 0:
        run_Command = "SU2_MDC " + config_name
        #else:
            #run_Command = os.path.join( SU2_RUN , run_command )
            #run_Command = "mpirun -np %i %s" % ( partitions , run_Command )                
        
        # log file
        if self.config_current['CONSOLE'] in ['QUIET','CONCISE']:
            run_Command = run_Command + ' >> ' + self.log_filename  
            
        # run command 
        os.system(run_Command)
        
        # no design to return
        design_return = {}
                    
        return design_return
    
    # -------------------------------------------------------------------
    #  PACKUP DEFORMATION TASK
    # -------------------------------------------------------------------
    def packup( self, config_delta, assets_super, design_super ):
        ''' assign asset push deformed mesh '''
        
        partitions = self.config_current['NUMBER_PART']
        
        # find deformed output mesh filename
        meshout_filename  = self.config_current['MESH_OUT_FILENAME']
        meshpush_filename = self.config_current['MESH_FILENAME']
        meshpush_filename = libSU2.add_Prefix( meshpush_filename , 'deform')

        # remove pathing
        meshpush_filename = os.path.split(meshpush_filename)[-1]

        # rename meshes with prefixes
        #if partitions > 1:
            #meshout_supername  = libSU2.add_Prefix( meshout_filename,  '%i' )
            #meshpush_supername = libSU2.add_Prefix( meshpush_filename, '%i' )
            #for p in range(partitions): shutil.move( meshout_supername % p , meshpush_supername % p )
        #else:
        shutil.move( meshout_filename , meshpush_filename )

        # mark for push, update current assets
        assets_push = {'mesh' : meshpush_filename }
        self.assets_current.update( assets_push )
        
        config_return = { 'MESH_FILENAME' : meshpush_filename }
        
        # remove from pull list for future runs
        del self.assets_pull['mesh']              
        
        return assets_push, config_return
    
    # -------------------------------------------------------------------
    #  DEFORMATION TASK REPRESENTATION
    # -------------------------------------------------------------------
    def __repr__ (self):
        return '<Task> SU2 Deform'        
        
#: class Deform()


#   % % %     % % %     % % %     % % %     % % %     % % %     % % %   #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  SU2 MESH ADAPTATION TASK
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
class Adapt(General_Task): 
    
    # -------------------------------------------------------------------
    #  INITIALIZE MESH ADAPTATION TASK
    # -------------------------------------------------------------------
    def __init__( self, assets_super, config_delta ):

        General_Task.__init__(self)
        
        self.config_current = libSU2.Get_ConfigParams( assets_super['config'] )
        self.config_current.update(config_delta)               

        self.assets_pull  = ['config','mesh','direct','adjoint']
        
        self.config_check = ['ADAPT_CYCLES','KIND_ADAPT','ADJ_OBJFUNC']
                
        self.folder_self  = 'ADAPT'
        
        self.log_filename = 'log_Adapt.out'
        
        self.signature    = self.set_signature()
        
        return
        
    # -------------------------------------------------------------------
    #  SET MESH ADAPTATION TASK SIGNATURE
    # -------------------------------------------------------------------
    def set_signature( self ):
        ''' assign the properties that represent
            this task's configuration '''
                                
        # default signature items
        signature = { 'TYPE' : 'SU2 Task',
                      'NAME' : 'Mesh Adaptation' }
        
        # list of keys to copy from a config
        config_keys = ['ADAPT_CYCLES','KIND_ADAPT','ADJ_OBJFUNC']
        
        # assign config signature items
        for key in config_keys:
            signature[key] = self.config_current[key]
        
        return signature          
    
    # ----------------------------------------------------------------------     
    #  CHECK MESH ADAPTATION TASK
    # ----------------------------------------------------------------------     
    def check( self, config_delta, assets_super, design_super ):
        ''' check if job needs to be run or modify super assets ''' 
        
        # check if the mesh has already been adapted
        key = 'mesh'
        status = not ( self.assets_current.has_key(key) )
        
        return status    

    # -------------------------------------------------------------------
    #  CONFIGURE MESH ADAPTATION TASK
    # -------------------------------------------------------------------
    def configure( self, config_delta, assets_super, design_super ):
        ''' Configure Mesh Adaptation  '''
             
        kind_adapt       = self.config_current['KIND_ADAPT']
        direct_solution  = self.config_current['SOLUTION_FLOW_FILENAME']
        adjoint_solution = self.config_current['SOLUTION_ADJ_FILENAME']
        adj_objfunc      = self.config_current['ADJ_OBJFUNC']
        direct_solution  = os.path.abspath( direct_solution )
        adjoint_solution = os.path.abspath( adjoint_solution )        
        adj_prefix       = libSU2.get_AdjointPrefix(adj_objfunc)
        adjoint_solution = libSU2.add_Prefix( adjoint_solution , adj_prefix )
        
        assert kind_adapt in ['GRAD_FLOW','GRAD_ADJOINT','GRAD_FLOW_ADJ'] , 'unsupported adaptation kind'
        
        # unlikely that restart will be available for unadapted mesh
        config_delta['RESTART_SOL'] = 'NO'
        
        ## check for restart
        #if self.config_current['RESTART_SOL'] == 'YES':
            #if not os.path.exists(direct_solution):
                #config_delta['RESTART_SOL'] = 'NO'                    
            
            #if kind_adapt in ['GRAD_FLOW','GRAD_FLOW_ADJ']:
                #if not os.path.exists(adjoint_solution):
                    #config_delta['RESTART_SOL'] = 'NO'                     
                    
        #: if check restart
        
        # update config
        libSU2.Set_ConfigParams( self.assets_current['config'], config_delta )
        self.config_current.update(config_delta);
        
        return
        
    # -------------------------------------------------------------------
    #  RUN MESH ADAPTATION TASK
    # -------------------------------------------------------------------
    def run( self, config_delta, assets_super, design_super ):
        ''' Run Mesh Adaptation '''
        
        config_name = self.assets_current['config']
        partitions  = self.config_current['NUMBER_PART']
        cycles      = self.config_current['ADAPT_CYCLES']
        
        # prepare command
        run_Command = 'mesh_adaptation.py -f %s -p %i -c %i' % (config_name,partitions,cycles)
            
        if self.config_current['CONSOLE'] in ['QUIET','CONCISE']:
            run_Command = run_Command + ' >> ' + self.log_filename                  
            
        # run command 
        os.system(run_Command)
                    
        # no design to return
        design_return = {}
        
        return design_return
    
    # -------------------------------------------------------------------
    #  PACKUP MESH ADAPTATION TASK
    # -------------------------------------------------------------------
    def packup( self, config_delta, assets_super, design_super ):
        ''' assign asset push deformed mesh '''
        
        assets_push   = {}
        config_return = {}
        
        # find adapted output mesh filename
        kind_adapt        = self.config_current['KIND_ADAPT']
        meshout_filename  = Mesh_MAC_file = "mesh_MAC_" + self.assets_current['config'].replace(".cfg",".su2")
        meshpush_filename = self.config_current['MESH_FILENAME']
        meshpush_filename = libSU2.add_Prefix( meshpush_filename , 'adapt')
        direct_restart    = self.config_current['RESTART_FLOW_FILENAME']
        direct_solution   = self.config_current['SOLUTION_FLOW_FILENAME']        
        adjoint_restart   = self.config_current['RESTART_ADJ_FILENAME']
        adjoint_solution  = self.config_current['SOLUTION_ADJ_FILENAME']
        cadj_prefix       = libSU2.get_AdjointPrefix( self.config_current['ADJ_OBJFUNC'] )
        
        # remove pathing
        meshpush_filename = os.path.split(meshpush_filename)[-1]
        direct_solution   = os.path.split(direct_solution)[-1]        
        adjoint_solution  = os.path.split(adjoint_solution)[-1]        
        
        # packup mesh
        shutil.move( meshout_filename , meshpush_filename )
        
        # packup direct solution
        shutil.move(direct_restart,direct_solution)
        
        # mark for push
        assets_push['mesh']   = meshpush_filename
        assets_push['direct'] = direct_solution
        
        # config updates
        config_return['MESH_FILENAME']          = meshpush_filename
        config_return['SOLUTION_FLOW_FILENAME'] = direct_solution
        config_return['RESTART_SOL']            = 'YES'
        
        # remove from pull list for future runs
        del self.assets_pull['mesh']              
        del self.assets_pull['direct']              
        
        # adjoint solutions
        if kind_adapt in ['GRAD_ADJOINT','GRAD_FLOW_ADJ']:
            
            # prefix
            restart_prefixed  = libSU2.add_Prefix(adjoint_restart,cadj_prefix)
            solution_prefixed = libSU2.add_Prefix(adjoint_solution,cadj_prefix)
            
            # packup adjoint solutioj
            shutil.move(restart_prefixed,solution_prefixed)
            
            # mark for push
            assets_push['adjoint'] = adjoint_solution
            
            # config updates
            config_return['SOLUTION_ADJ_FILENAME'] = adjoint_solution
            
            # remove from pull list for future runs
            del self.assets_pull['adjoint']
        
        #: if adjoint
            
        # update current assets
        self.assets_current.update(assets_push)
        
        return assets_push, config_return  
    
    # -------------------------------------------------------------------
    #  MESH ADAPTATION TASK REPRESENTATION
    # -------------------------------------------------------------------
    def __repr__ (self):
        return '<Task> SU2 Adapt'        
        
#: class Adapt()


#   % % %     % % %     % % %     % % %     % % %     % % %     % % %   #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  SU2 DIRECT SOLUTION TASK
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
class Direct(General_Task): 
    
    # -------------------------------------------------------------------
    #  INITIALIZE DIRECT SOLUTION TASK
    # -------------------------------------------------------------------
    def __init__( self, assets_super, config_delta ):

        General_Task.__init__(self)
        
        self.config_current = libSU2.Get_ConfigParams( assets_super['config'] )
        self.config_current.update(config_delta)               

        self.assets_pull  = ['config','mesh','direct']
        
        self.config_check = []
                
        self.folder_self  = 'DIRECT'
        
        self.log_filename = 'log_Direct.out'
        
        self.signature    = self.set_signature()
        
        return
        
    # -------------------------------------------------------------------
    #  SET DIRECT SOLUTION TASK SIGNATURE
    # -------------------------------------------------------------------
    def set_signature( self ):
        ''' assign the properties that represent
            this task's configuration '''
                    
        # default signature items
        signature = { 'TYPE' : 'SU2 Task',
                      'NAME' : 'Direct Solution' }
        
        # list of keys to copy from a config delta
        config_keys = []
        
        # assign config signature items
        for key in config_keys:
            signature[key] = self.config_current[key]
        
        return signature          
    
    # ----------------------------------------------------------------------     
    #  CHECK DIRECT SOLUTION  TASK
    # ----------------------------------------------------------------------     
    def check( self, config_delta, assets_super, design_super ):
        ''' check if job needs to be run or modify super assets ''' 
        
        status = False        
        
        # check if direct solution has been run with requested mesh
        keys = ['mesh','direct']
        for key in keys:
            status = not ( self.assets_current.has_key(key) )        
            if status: break
        
        return status    

    # -------------------------------------------------------------------
    #  CONFIGURE DIRECT SOLUTION TASK
    # -------------------------------------------------------------------
    def configure( self, config_delta, assets_super, design_super ):
        ''' Configure SU2_CFD Direct  '''
        
        # assume self.config_current loaded in pull_assets           
        
        direct_solution = self.config_current['SOLUTION_FLOW_FILENAME']
        direct_solution = os.path.abspath( direct_solution )
        
        # set direct math problem
        config_delta['MATH_PROBLEM'] = 'DIRECT'       
             
        # check for restart
        if self.config_current['RESTART_SOL'] == 'YES':
            if not os.path.exists(direct_solution):
                config_delta['RESTART_SOL'] = 'NO'
                                
        # update config
        libSU2.Set_ConfigParams( self.assets_current['config'], config_delta )
        self.config_current.update(config_delta);
        
        # ADAPTATION     
        if 'ADAPT' in self.config_current['TASKS']:
            
            # configure 
            config_adapt = {'KIND_ADAPT':'GRAD_FLOW'}
            This_Adapt = Adapt(self.assets_current,config_adapt)
            
            # run
            print('    %s' % This_Adapt.signature['NAME'])
            (_,assets_new,config_new) = This_Adapt.evaluate(config_adapt,self.assets_current)
            
            # update config and assets
            self.assets_current.update( assets_new )
            config_delta.update( config_new )
        
        # if adaptation
        
        # update config
        libSU2.Set_ConfigParams( self.assets_current['config'], config_delta )
        self.config_current.update(config_delta);        
          
        
        return
        
    # -------------------------------------------------------------------
    #  RUN DIRECT SOLUTION TASK
    # -------------------------------------------------------------------
    def run( self, config_delta, assets_super, design_super ):
        ''' Run Direct '''
        
        config_name      = self.assets_current['config']
        partitions       = self.config_current['NUMBER_PART']
        history_filename = self.config_current['CONV_FILENAME']
        output_format    = self.config_current['OUTPUT_FORMAT']
        plotfile_ext     = libSU2.get_ExtensionName(output_format)
        special_cases    = libSU2.get_SpecialCases(self.config_current)            
        
        # prepare command
        if partitions == 0:
            run_Decomp = ""
            run_Direct = "SU2_CFD %s" % config_name
            run_Merge  = ""
        else:
            run_Decomp = "SU2_DDC %s" % config_name
            run_Direct = os.path.join( SU2_RUN , "SU2_CFD " + config_name )
            run_Direct = "mpirun -np %i %s" % ( partitions , run_Direct )   
            run_Merge  = "merge_solution.py -f %s -p %i" % ( config_name , partitions )
        
        if self.config_current['CONSOLE'] in ['QUIET','CONCISE']:
            run_Decomp = run_Decomp + ' >> ' + self.log_filename  
            run_Direct = run_Direct + ' >> ' + self.log_filename  
            run_Merge  = run_Merge  + ' >> ' + self.log_filename  
            
        # run commands
        if partitions > 1: os.system(run_Decomp)
        os.system(run_Direct)
        if partitions > 1: os.system(run_Merge)
        
        # get objective function values, update design
        ObjFun_Dict = libSU2.get_ObjFunVals( history_filename+plotfile_ext , special_cases )
        
        # make sure outputs are lists for making 2d arrays
        for key,value in ObjFun_Dict.iteritems():
            if not isinstance(value,list):
                ObjFun_Dict[key] = [value]
                
        # update design
        self.design_current['OBJECTIVES'] = ObjFun_Dict            
        
        return self.design_current
    
    # -------------------------------------------------------------------
    #  PACKUP DIRECT SOLUTION TASK
    # -------------------------------------------------------------------
    def packup( self, config_delta, assets_super, design_super ):
        ''' assign asset push deformed mesh '''
        
        assets_push   = {}
        config_return = {}
        
        restart_filename  = self.config_current['RESTART_FLOW_FILENAME']
        solution_filename = self.config_current['SOLUTION_FLOW_FILENAME']
        
        # remove pathing
        solution_filename = os.path.split(solution_filename)[-1]
        
        # rename restart to solution
        shutil.move(restart_filename,solution_filename)
        
        # mark for push
        assets_push['direct'] = solution_filename
        config_return['SOLUTION_FLOW_FILENAME'] = solution_filename
        
        # mark adapted mesh
        if 'ADAPT' in self.config_current['TASKS']:
            assets_push['mesh'] = self.assets_current['mesh']
            config_return['MESH_FILENAME'] = self.assets_current['mesh']
            
        # update current assets
        self.assets_current.update( assets_push )
        
        # remove from pull list for future runs
        del self.assets_pull['mesh']
        del self.assets_pull['direct']
        
        return assets_push, config_return
    
    # -------------------------------------------------------------------
    #  DIRECT SOLUTION TASK REPRESENTATION
    # -------------------------------------------------------------------
    def __repr__ (self):
        return '<Task> SU2 Direct'        
        
#: class Direct()


#   % % %     % % %     % % %     % % %     % % %     % % %     % % %   #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  SU2 MULTIPLE CONTINUOUS ADJOINT SOLUTION TASK
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
class Multiple_Cont_Adjoint(General_Task): 
    # ----------------------------------------------------------------------
    #  INITIALIZE MULTIPLE CONTINUOUS ADJOINTS TASK
    # ----------------------------------------------------------------------    
    def __init__( self, assets_super, config_delta):
        ''' initialize class attributes 
            will start in folder_super
            should not create any new files because
            the class may be built temporarily 
            as a comparison tester '''     
                
        General_Task.__init__(self)
        
        self.config_current = libSU2.Get_ConfigParams( assets_super['config'] )
        self.config_current.update(config_delta)          
        
        self.assets_pull  = ['config','mesh','direct','adjoint'] 
        
        self.config_check = ['GRADIENTS']
        
        self.design_current = {'GRADIENTS':{}}
        
        self.folder_self  = 'CONT_ADJOINT'        
        
        self.signature    = self.set_signature(config_delta)
        
        self.tasks_todo   = []
        self.tasks_done   = []
        
        return
    
    # ----------------------------------------------------------------------
    #  SET MULTIPLE CONTINUOUS ADJOINTS TASK
    # ----------------------------------------------------------------------    
    def set_signature( self, config_delta ):
        ''' assign the properties that represent
            this task's configuration '''
        
        # default signature items
        signature = { 'TYPE' : 'SU2 Task'                              ,
                      'NAME' : 'Multiple Continuous Adjoint Solutions'  } 
        
        # list of keys to copy from a config
        config_keys = []
        
        # assign config signature items
        for key in config_keys:
            signature[key] = self.config_current[key]
        
        return signature        
    
    # ----------------------------------------------------------------------     
    #  CHECK MULTIPLE CONTINUOUS ADJOINTS TASK
    # ----------------------------------------------------------------------     
    def check( self, config_delta, assets_super, design_super ):
        ''' check if job needs to be run or modify super assets ''' 
        
        # modifications to input dictionaries will be by reference
        self.config_current.update(config_delta)
        
        # validate input data, return status True=go
        status = False
        
        # check if solution has been run with requested mesh and direct solution
        keys = ['mesh','direct']
        for key in keys:
            status = not ( self.assets_current.has_key(key) )        
            if status: break
            
        # check if gradients were run
        if not status:
            for gradient in self.config_current['GRADIENTS']:
                status = not ( self.design_current.has_key('GRADIENTS') and
                               gradient in self.design_current['GRADIENTS'] )
                if status: break
        
        return status       
        
    # ----------------------------------------------------------------------
    #  CONFIGURE MULTIPLE CONTINUOUS ADJOINTS TASK
    # ----------------------------------------------------------------------    
    def configure( self, config_delta, assets_super, design_super ):      
        ''' Configure multiple continuous adjoint tasks '''
                
        # tasks to run
        self.tasks_todo  = []
        tasks_todo_names = self.config_current['GRADIENTS']
        
        # process new tasks
        for this_name in tasks_todo_names:
                  
            # set adjoint objective function
            this_config_delta = { 'ADJ_OBJFUNC' : this_name }
                  
            # a potential new task
            Constructor = Cont_Adjoint
            New_Task = Constructor(self.assets_current,this_config_delta)
            
            # list of matching tasks
            Tasks_Check = [ task for task in self.tasks_done if New_Task == task ]
            
            # check for previous runs
            if Tasks_Check:  # task already done
                if New_Task in self.tasks_todo:  # task already scheduled, delete
                    del New_Task; continue
                # else ...
                assert len(Tasks_Check) <= 1 , 'found more than one completed task in multiadjoint'
                del New_Task 
                New_Task = Tasks_Check[0]
            # else ...            
                        
            # add the requested task
            self.tasks_todo.append(New_Task) 
            if not Tasks_Check: self.tasks_done.append(New_Task)
            
        #: for new tasks
        
        return
    
    # ----------------------------------------------------------------------
    #  RUN MULTIPLE CONTINUOUS ADJOINTS TASK
    # ----------------------------------------------------------------------
    def run( self, config_delta, assets_super, design_super ):
        ''' run the analysis task and return an 
            updated design dictionary '''

        # Run the Tasks!!!
        for This_Task in self.tasks_todo:
            
            # run this task
            if not self.config_current['CONSOLE'] == 'QUIET':
                sys.stdout.write( '    %s - %s \n' % (This_Task.signature['NAME'],This_Task.signature['OBJECTIVE']) )
            design_new,assets_new,_ = This_Task.evaluate({},self.assets_current)
            
            # update design and assets
            self.design_current['GRADIENTS'].update( design_new['GRADIENTS'] )
            #self.assets_current.update( assets_new )
        
        #: for each Task
        
        return self.design_current

    # ----------------------------------------------------------------------
    #  PACKUP MULTIPLE CONTINUOUS ADJOINTS TASK
    # ----------------------------------------------------------------------    
    def packup( self, config_delta, assets_super, design_super ):
        ''' assign assets to push to super folder 
            and clean uneeded files '''
        
        gradients        = self.config_current['GRADIENTS'] 
        adjoint_solution = self.config_current['SOLUTION_ADJ_FILENAME']
        
        # remove pathing
        adjoint_solution = os.path.split(adjoint_solution)[-1]
        
        # mark for push
        assets_push = { 'adjoint' : adjoint_solution }        
        config_return = { 'SOLUTION_ADJ_FILENAME' : adjoint_solution }
        
        # update current assets
        self.assets_current.update( assets_push )
        
        # remove from pull list for future runs
        if self.assets_pull.has_key('mesh'): del self.assets_pull['mesh']
        if self.assets_pull.has_key('direct'): del self.assets_pull['direct']
        # allow new adjoints to be pulled
        
        return assets_push, config_return  
    
    # ----------------------------------------------------------------------
    #  MULTIPLE CONTINUOUS ADJOINTS TASK REPRESENATION
    # ----------------------------------------------------------------------    
    def __repr__ (self):
        return '<Task> Multiple Continuous Adjoint Solutions'

#: class Job    
    


#   % % %     % % %     % % %     % % %     % % %     % % %     % % %   #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  SU2 CONTINUOUS ADJOINT SOLUTION TASK
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
class Cont_Adjoint(General_Task): 
    
    # -------------------------------------------------------------------
    #  INITIALIZE CONTINUOUS ADJOINT SOLUTION TASK
    # -------------------------------------------------------------------
    def __init__( self, assets_super, config_delta ):

        General_Task.__init__(self)
        
        self.config_current = libSU2.Get_ConfigParams( assets_super['config'] )
        self.config_current.update(config_delta)               

        objective = self.config_current['ADJ_OBJFUNC']
        prefix    = libSU2.get_AdjointPrefix(objective)
        
        self.assets_pull  = ['config','mesh','direct','adjoint']
        
        self.config_check = ['ADJ_OBJFUNC']
        
        self.folder_self  = objective
        
        self.log_filename = libSU2.add_Prefix( 'log_Adjoint.out' , prefix )
        
        self.FinDiff_Step = 1e-4 # GPC finite difference step
        
        self.signature    = self.set_signature(objective)
        
        return
        
    # -------------------------------------------------------------------
    #  SET CONTINUOUS ADJOINT SOLUTION TASK SIGNATURE
    # -------------------------------------------------------------------
    def set_signature( self , objective ):
        ''' assign the properties that represent
            this task's configuration '''
        
        # default signature items
        signature = { 'TYPE'      : 'SU2 Task'                    ,
                      'NAME'      : 'Continuous Adjoint Solution' ,
                      'OBJECTIVE' : objective                      }
        
        # list of keys to copy from a config
        config_keys = []
        
        # assign config signature items
        for key in config_keys:
            signature[key] = self.config_current[key]
        
        return signature          

    # ----------------------------------------------------------------------     
    #  CHECK CONTINUOUS ADJOINT SOLUTION TASK
    # ----------------------------------------------------------------------     
    def check( self, config_delta, assets_super, design_super ):
        ''' check if job needs to be run or modify super assets ''' 
        
        # validate input data, return status True=go
        status = False
        
        # check if direct solution has been run with requested mesh and adjoint solution
        keys = ['mesh','direct','adjoint']
        for key in keys:
            status = not ( self.assets_current.has_key(key) )        
            if status: break
        
        return status   

    # -------------------------------------------------------------------
    #  CONFIGURE CONTINUOUS ADJOINT SOLUTION TASK
    # -------------------------------------------------------------------
    def configure( self, config_delta, assets_super, design_super ):
        ''' Configure SU2_CFD Continuous Adjoint  '''
        # assume self.config_current loaded in pull_assets           
        
        gradient         = self.signature['OBJECTIVE']
        adjoint_solution = self.config_current['SOLUTION_ADJ_FILENAME']
        cadj_prefix      = libSU2.get_AdjointPrefix( gradient )                
        adjoint_solution = libSU2.add_Prefix(adjoint_solution,cadj_prefix)
        adjoint_solution = os.path.abspath(adjoint_solution)
        
        # set adjoint math problem
        config_delta['MATH_PROBLEM'] = 'ADJOINT'
        config_delta['ADJ_OBJFUNC']  = gradient
        
        # check for restart
        if self.config_current['RESTART_SOL'] == 'YES':
            if not os.path.exists( adjoint_solution ):
                config_delta['RESTART_SOL'] = 'NO'
                                            
        # setup gradient projection
        n_DV = len( self.config_current['DEFINITION_DV']['MARKER'] )
        config_delta['DV_VALUE_OLD'] = numpy.zeros(n_DV)
        config_delta['DV_VALUE_NEW'] = numpy.ones(n_DV) * self.FinDiff_Step    
                                
        # update config
        libSU2.Set_ConfigParams( self.assets_current['config'], config_delta )
        self.config_current.update(config_delta);
        
        # ADAPTATION     
        if 'ADAPT_ADJ' in self.config_current['TASKS']:
            
            # configure
            config_adapt = {'KIND_ADAPT':'GRAD_ADJOINT'}
            This_Adapt = Adapt(self.assets_current,config_adapt)
            
            # run adaptation
            print('      %s' % This_Adapt.signature['NAME'])
            (_,assets_new,config_new) = This_Adapt.evaluate(config_adapt,self.assets_current)
            
            # update assets and configuration
            self.assets_current.update( assets_new )
            config_delta.update( config_new )
            
        # if adaptation
            
        # update config
        libSU2.Set_ConfigParams( self.assets_current['config'], config_delta )
        self.config_current.update(config_delta);            
        
        return
        
    # -------------------------------------------------------------------
    #  RUN CONTINUOUS ADJOINT SOLUTION TASK
    # -------------------------------------------------------------------
    def run( self, config_delta, assets_super, design_super ):
        ''' Run MDC '''
        
        config_name        = self.assets_current['config']
        partitions         = self.config_current['NUMBER_PART']
        gradient_filename  = self.config_current['GRAD_OBJFUNC_FILENAME']
        objective_function = self.config_current['ADJ_OBJFUNC']
        
        # moothing configuration
        smoothing_type     = self.config_current.get('ADJOINT_SMOOTHING','NONE')
        surface_filename   = self.config_current['SURFACE_ADJ_FILENAME']
        filtered_filename  = surface_filename + '_filtered'
        marker_smooth      = 'airfoil'
        chord_length       = 1.0
                
        # prepare command
        if partitions == 0:
            run_Decomp  = ""
            run_Adjoint = "SU2_CFD " + config_name
            run_Project = "SU2_GPC " + config_name
            run_Merge   = ""
        else:
            run_Decomp  = "SU2_DDC %s" % config_name
            run_Adjoint = os.path.join( SU2_RUN , "SU2_CFD " + config_name )
            run_Project = os.path.join( SU2_RUN , "SU2_GPC " + config_name )
            run_Adjoint = "mpirun -np %i %s" % ( partitions , run_Adjoint )
            run_Project = "mpirun -np %i %s" % ( partitions , run_Project )
            run_Merge   = "merge_solution.py -f %s -p %i" % ( config_name , partitions )
        #: if partitions
        run_Filter = 'filter_adjoint.py -f %s -t %s -m %s -c %f' % (config_name,smoothing_type,marker_smooth,chord_length)
        
        if self.config_current['CONSOLE'] in ['QUIET','CONCISE']:
            run_Decomp  = run_Decomp  + ' >> ' + self.log_filename  
            run_Adjoint = run_Adjoint + ' >> ' + self.log_filename  
            run_Project = run_Project + ' >> ' + self.log_filename  
            run_Merge   = run_Merge   + ' >> ' + self.log_filename  
            run_Filter  = run_Filter  + ' >> ' + self.log_filename  
            
        # run solution
        if partitions > 1: os.system( run_Decomp )
        os.system(run_Adjoint)
        if partitions > 1: os.system( run_Merge )
        
        # run filtering
        if not smoothing_type == 'NONE':
            os.system( run_Filter )
            config_update = { 'SURFACE_ADJ_FILENAME' : filtered_filename }
            self.config_current.update(config_update)
            config_delta.update(config_update)
            libSU2.Set_ConfigParams(config_name,config_update)
        
        # run gradient projection
        os.system(run_Project)
                         
        # get gradient values, update design
        gradients     = libSU2.get_GradientVals(gradient_filename)  
        self.design_current['GRADIENTS'] = { objective_function : gradients }

        return self.design_current
    
    # -------------------------------------------------------------------
    #  PACKUP CONTINUOUS ADJOINT SOLUTION TASK
    # -------------------------------------------------------------------
    def packup( self, config_delta, assets_super, design_super ):
        ''' assign asset push deformed mesh '''
        
        assets_push   = {}
        config_return = {}
        
        objective_function = self.config_current['ADJ_OBJFUNC']         
        restart_filename   = self.config_current['RESTART_ADJ_FILENAME']
        solution_filename  = self.config_current['SOLUTION_ADJ_FILENAME'] 
        mesh_filename      = self.config_current['MESH_FILENAME']
        cadj_prefix        = libSU2.get_AdjointPrefix(objective_function)
        mesh_prefixed      = libSU2.add_Prefix(mesh_filename,cadj_prefix)        

        # remove pathing
        solution_filename  = os.path.split(solution_filename)[-1]

        # add prefix
        restart_prefixed   = libSU2.add_Prefix(restart_filename,cadj_prefix)                
        solution_prefixed  = libSU2.add_Prefix(solution_filename,cadj_prefix) 

        # write gradient file
        self.write_gradfile(objective_function)        

        # rename restart to solution
        shutil.move(restart_prefixed,solution_prefixed)
        
        # adapated mesh
        if 'ADAPT_ADJ' in self.config_current['TASKS']:
            shutil.move(mesh_filename,mesh_prefixed)
            assets_push['mesh'] = mesh_prefixed
            config_return['MESH_FILENAME'] = mesh_filename
        
        # mark for push, update current assets
        assets_push['adjoint']                 = solution_filename
        config_return['SOLUTION_ADJ_FILENAME'] = solution_filename
        
        self.assets_current.update( assets_push )
        
        # remove from pull list for future runs
        del self.assets_pull['mesh']
        del self.assets_pull['direct']        
        del self.assets_pull['adjoint']        
                    
        return assets_push, config_return
    
    # -------------------------------------------------------------------
    #  WRITE CONTINUOUS ADJOINT GRADIENT FILE
    # -------------------------------------------------------------------        
    def write_gradfile( self, objective_function ):
        
        Definition_DV      = self.config_current['DEFINITION_DV']
        gradient_filename  = self.config_current['GRAD_OBJFUNC_FILENAME']
        output_format      = self.config_current['OUTPUT_FORMAT']
        cadj_prefix        = libSU2.get_AdjointPrefix(objective_function)
        plotfile_ext       = libSU2.get_ExtensionName(output_format)            
        gradients          = self.design_current['GRADIENTS'][objective_function]
        n_DV               = len(Definition_DV['MARKER'])

        # append extension and plot format
        gradient_filename = os.path.splitext(gradient_filename)[0]
        gradient_filename = gradient_filename+"_"+cadj_prefix+plotfile_ext                        
                    
        # get header and write format information
        header,_ = libSU2.get_GradFileFormat( 'CONTINUOUS_ADJOINT'      ,
                                              output_format            ,
                                              Definition_DV['KIND'][0]  )
        # start gradient file
        Grad_File = open(gradient_filename,'w')
        Grad_File.write(header)
        
        # Write output gradients and dv information 
        for i_DV in range(n_DV):
            _,write_format = libSU2.get_GradFileFormat( 'CONTINUOUS_ADJOINT',
                                                        output_format,
                                                        Definition_DV['KIND'][i_DV] )        
            # merge information for output
            Output_Vector = [i_DV] + [gradients[i_DV]] + [self.FinDiff_Step] + Definition_DV['PARAM'][i_DV]
            Output_Vector = tuple(Output_Vector)
            # write row
            Grad_File.write(write_format % Output_Vector)
        #: for each DV  
        
        Grad_File.close()             
        
        return
    
    # -------------------------------------------------------------------
    #  DIRECT CONTINUOUS ADJOINT TASK REPRESENTATION
    # -------------------------------------------------------------------
    def __repr__ (self):
        return '<Task> SU2 Cont_Adjoint %s' % self.signature['OBJECTIVE']
        
#: class Cont_Adjoint()


#   % % %     % % %     % % %     % % %     % % %     % % %     % % %   #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  SU2 FINITE DIFFERENCING TASK
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
class Finite_Diff(General_Task): 
    
    # -------------------------------------------------------------------
    #  INITIALIZE FINITE DIFFERENCING TASK
    # -------------------------------------------------------------------
    def __init__( self, assets_super, config_delta ):

        General_Task.__init__(self)
        
        self.config_current = libSU2.Get_ConfigParams( assets_super['config'] )
        self.config_current.update(config_delta)               

        self.assets_pull  = ['config','mesh','direct'] # direct solution
        
        self.config_check = []
                
        self.folder_self  = 'FINITE_DIFF'
        
        self.log_filename = 'log_FinDiff.out'
        
        self.signature    = self.set_signature()
        
        return
        
    # -------------------------------------------------------------------
    #  SET FINITE DIFFERENCING TASK SIGNATURE
    # -------------------------------------------------------------------
    def set_signature( self ):
        ''' assign the properties that represent
            this task's configuration '''
        
        # default signature items
        signature = { 'TYPE' : 'SU2 Task',
                      'NAME' : 'Finite Difference' }
        
        # list of keys to copy from a config delta
        config_keys = []
        
        # assign config signature items
        for key in config_keys:
            signature[key] = self.config_current[key]
        
        return signature       
    
    # ----------------------------------------------------------------------     
    #  CHECK FINITE DIFFERENCING TASK
    # ----------------------------------------------------------------------     
    def check( self, config_delta, assets_super, design_super ):
        ''' check if job needs to be run or modify super assets ''' 
        
        # modifications to input dictionaries will be by reference
        
        # validate input data, return status True=go
        status = False
        
        # check if direct solution has been run with requested mesh
        keys = ['mesh','direct']
        for key in keys:
            status = not ( self.assets_current.has_key(key) )        
            if status: break
        
        return status       

    # -------------------------------------------------------------------
    #  CONFIGURE FINITE DIFFERENCING TASK
    # -------------------------------------------------------------------
    def configure( self, config_delta, assets_super, design_super ):
        ''' Configure SU2_CFD Direct  '''
        # assume self.config_current loaded in pull_assets           
                
        # check for restart
        if self.assets_current.has_key('direct'):
            if os.path.exists(self.assets_current['direct']):
                config_delta['RESTART_SOL'] = 'YES'
                                
        # update config
        libSU2.Set_ConfigParams( self.assets_current['config'], config_delta )
        self.config_current.update(config_delta);
        
        return
        
    # -------------------------------------------------------------------
    #  RUN FINITE DIFFERENCING TASK
    # -------------------------------------------------------------------
    def run( self, config_delta, assets_super, design_super ):
        ''' Run Finite_Diff '''
        
        config_name   = self.assets_current['config']
        partitions    = self.config_current['NUMBER_PART']
        fin_diff_step = self.config_current.get('FIN_DIFF_STEP') # defaults 1e-4
        
        # console output
        if self.config_current['CONSOLE'] in ['CONCISE','QUIET']:
            logfile = self.log_filename
        else:
            logfile = None
            
        # run finite differences
        ObjFun_Dict,Gradient_Dict = finite_differences( filename   = config_name   ,
                                                        partitions = partitions    ,
                                                        step       = fin_diff_step ,
                                                        output     = False         ,
                                                        logfile    = logfile        )            
        # save design data
        ### self.design_current['OBJECTIVES'] = ObjFun_Dict
        self.design_current['GRADIENTS']  = Gradient_Dict
        
        return self.design_current
    
    # -------------------------------------------------------------------
    #  PACKUP FINITE DIFFERENCING TASK
    # -------------------------------------------------------------------
    def packup( self, config_delta, assets_super, design_super ):
        ''' assign asset push deformed mesh '''
        
        assets_push   = {}
        config_return = {}
        
        restart_filename  = self.config_current['RESTART_FLOW_FILENAME']
        solution_filename = self.config_current['SOLUTION_FLOW_FILENAME']
        
        # remove pathing
        solution_filename = os.path.split(solution_filename)[-1]
        
        # rename restart to solution
        shutil.move(restart_filename,solution_filename)
        
        # mark for push, update current assets
        assets_push['direct'] = solution_filename
        config_return['SOLUTION_FILENAME'] = solution_filename
        
        # update current assets
        self.assets_current.update(assets_push)

        # remove from pull list for future runs
        del self.assets_pull['mesh']        
        del self.assets_pull['direct']        
        
        return assets_push, config_return
    
    # -------------------------------------------------------------------
    #  FINITE DIFFERENCING TASK REPRESENTATION
    # -------------------------------------------------------------------
    def __repr__ (self):
        return '<Task> SU2 Finite_Diff'        

#: class Finite_Diff()



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  END SU2 TASKS MODULE
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
