#!/usr/bin/env python 

## \file tasks_project.py
#  \brief Python classes for evaluating SU2 projects
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

# http://xkcd.com/353/
# principle of least astonishment

# TODO: 
# data output options: heavy, compress, light
# sorting task list (if given out of order)

# NICE TODO:
# jobs/tasks decide if they will run, no presorting
# restart searching by task, not job

import os, sys, shutil, glob, copy
import numpy
import libSU2
import tasks_su2
from tasks_general      import General_Task
from optparse           import OptionParser
from finite_differences import finite_differences
from mesh_adaptation    import mesh_adaptation
from tasks_su2          import *

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  PROJECT CLASS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
class Project(General_Task):
    
    # ----------------------------------------------------------------------
    #  INITIALIZE PROJECT
    # ----------------------------------------------------------------------    
    def __init__( self, config_name, design_name ):
        ''' initialize class attributes '''      
        
        General_Task.__init__(self)
        
        # read config file
        self.config_current = libSU2.Get_ConfigParams(config_name)
        
        # define assets
        self.assets_current = { 'config' : config_name                          ,   # these are filenames
                                'design' : design_name                          ,
                                'mesh'   : self.config_current['MESH_FILENAME']  }
        
        # job lists
        self.jobs_todo = []
        self.jobs_done = []
        

        self.backup_name = 'backup_project.pkl'        
        
        # restart option
        self.use_restart = self.config_current['RESTART_SOL'] == 'YES'
        
        # check for restart assets in project folder
        if self.use_restart:
            solution_flow = self.config_current['SOLUTION_FLOW_FILENAME']
            solution_adj  = self.config_current['SOLUTION_ADJ_FILENAME']
            
            # filenames
            if os.path.exists(solution_flow): 
                self.assets_current['direct'] = solution_flow
            if len( glob.glob( libSU2.add_Prefix(solution_adj,'*') ) ):
                self.assets_current['adjoint'] = solution_adj
            
            # if the associated files don't exist, then they won't be pulled, 
                
        #: if use_restart
        
        # equivalent area files
        if self.config_current.get('EQUIV_AREA','NO') == 'YES':
            target_filename = 'TargetEA.dat'
            if os.path.exists(target_filename):
                self.assets_current['targetea'] = target_filename
            else:
                sys.stdout.write('Warning: no target equivalent area file \n')        
                
        # working folder is where the project starts
        self.folder_self = '.'
        
        # set signature
        self.signature = self.set_signature()
        
        return
    
    # ----------------------------------------------------------------------
    #  SET PROJECT SIGNATURE
    # ----------------------------------------------------------------------    
    def set_signature( self ):
        ''' assign the properties that represent
            this task's configuration '''
        
        # default signature items
        signature   = { 'TYPE' : 'Project',
                        'NAME' : self.assets_current['config'] }
        
        return signature        
        
    # ----------------------------------------------------------------------
    #  CONFIGURE PROJECT
    # ----------------------------------------------------------------------    
    def configure( self, config_delta, assets_super, design_super ):      
        ''' configure any control and data inputs 
            for the analysis task '''
        
        # build list of jobs to run
        # decide whether to continue an old job
        self.jobs_todo = [] 
        
        # prepare jobs
        for this_config_delta in config_delta:
            
            # check task lists
            if not this_config_delta.has_key('TASKS'):
                this_config_delta['TASKS'] = self.config_current['TASKS']
            else:
                for task in ['ADAPT','ADAPT_ADJ']:
                    if ( task in self.config_current['TASKS'] and 
                         not task in this_config_delta['TASKS'] ):
                        this_config_delta['TASKS'].append(task)
                
            # check gradient lists
            if not this_config_delta.has_key('GRADIENTS'):
                this_config_delta['GRADIENTS'] = self.config_current['GRADIENTS']            
                        
            # new job number
            this_config_delta['JOB_NUMBER'] = len(self.jobs_done)+1    
                    
            # a potential new job
            New_Job = Job(self.assets_current,this_config_delta)
            
            # -----------------------------------------------------------------------
            #  JOB LIST BUILDING            
            
            # look for matching jobs
            Jobs_Check = [ job for job in self.jobs_done if job == New_Job ]
                
            # make new job todo
            if Jobs_Check: 
                # multiple jobs found?
                if len(Jobs_Check) > 1: print 'Warning - found more than one matching job, choosing one'
                del New_Job
                New_Job = Jobs_Check[0]
            #: if Job_Check
            
            # assign job todo
            self.jobs_todo.append( New_Job ) 
            
            # add to done list, avoids a seach durring run
            if not Jobs_Check: self.jobs_done.append( New_Job ) 
                
        #: for each config_delta
        
        return 
    
    # ----------------------------------------------------------------------
    #  RUN PROJECT
    # ----------------------------------------------------------------------    
    def run( self, config_delta, assets_super, design_super ):
        ''' run the analysis task and return an 
            updated design dictionary '''
        
        design_return = { 'VARIABLES'  : [] ,
                          'OBJECTIVES' : {} ,
                          'GRADIENTS'  : {}  }
        
        # Run the Job List!!!
        for This_Job,this_config_delta in zip(self.jobs_todo,config_delta):
            
            # prepare assets for sending
            this_assets_send = copy.deepcopy( self.assets_current )
            
            # search for restart jobs
            if self.use_restart:
                assets_restart = self.find_restart( This_Job, this_config_delta )
                this_assets_send.update( assets_restart )
                # note: if no restart is found, will use project's restart
                
            # --- RUN THIS JOB --- #
            if not self.config_current['CONSOLE'] == 'QUIET':
                sys.stdout.write( 'JOB %i \n' % This_Job.signature['JOB_NUMBER'] )
            this_design,_,_ = This_Job.evaluate(this_config_delta,this_assets_send,{})
            
            # append the nested dictionary
            libSU2.append_nestdict( design_return, this_design )
            
            # save project
            libSU2.save_data( self.backup_name, self )
            
            if not self.config_current['CONSOLE'] == 'QUIET':
                sys.stdout.write('\n')
        
        #: for each Job
        
        
        
        # only return newly updated designs
        return design_return

    # ----------------------------------------------------------------------
    #  SEARCH FOR RESTART
    # ----------------------------------------------------------------------    
    def find_restart( self, The_Job, config_delta ):
        ''' searches for restart job based on signature
            currently only uses 'VARIABLES' '''
        
        # TODO: look for closest restart by each asset instead of by job 
        
        assets_restart = {}
        config_check = self.config_current
        config_check.update(config_delta)        
        restart_direct  = 'DIRECT' in config_check['TASKS']
        restart_adjoint = 'CONT_ADJOINT' in config_check['TASKS']
        
        # minimum distance and current job location
        DV_dist_min = numpy.inf
        DV_thisjob  = numpy.array(The_Job.signature['VARIABLES'])
                
        # search for a restart file
        for This_Restart in self.jobs_done:
            
            if This_Restart == The_Job:
                continue
            
            # make sure this job has direct or adjoint assets
            no_assets = True
            for key in This_Restart.assets_current.keys():
                if ( 'direct' in key or 'adjoint' in key ): 
                    no_assets = False # there may be multiple adjoint_* keys                
            if no_assets: continue
            
            # test DV distance
            DV_restjob = This_Restart.signature['VARIABLES'] 
            DV_dist_this = sum( (numpy.array(DV_restjob) - DV_thisjob)**2 , 0 )
            
            # skip update if design is not closer
            if DV_dist_this > DV_dist_min:
                continue
            
            # new closest job found, store distance and assets
            DV_dist_min = DV_dist_this
            
            # get restart folder
            folder_restart = This_Restart.folder_self
            
            # find available restart assets
            for key,value in This_Restart.assets_current.iteritems():
                if ( 'direct' in key or 'adjoint' in key ): 
                    assets_restart[key] = os.path.join(folder_restart,value)
            
        #: for job in jobs_done    
        
        return assets_restart

    # ----------------------------------------------------------------------
    #  PACKUP PROJECT
    # ----------------------------------------------------------------------    
    def packup( self, config_delta, assets_super, design_super ):
        ''' assign assets to push to super folder 
            and clean uneeded files '''
        
        # clear old design_current
        self.design_current = { 'VARIABLES'  : [] ,
                                'OBJECTIVES' : {} ,
                                'GRADIENTS'  : {}  }
        
        # rewrite project's design_current
        for This_Job in self.jobs_done:
            libSU2.append_nestdict( self.design_current, This_Job.design_current )
        
        # save current design
        libSU2.save_data(self.assets_current['design'],self.design_current)
        
        # save project
        libSU2.save_data( self.backup_name, self )
        
        # this is the root task, no assets to push
        assets_push = {}
        config_return = {}
        
        return assets_push, config_return  
    
    # ----------------------------------------------------------------------
    #  CLEANUP PROJECT
    # ----------------------------------------------------------------------    
    def cleanup( self, config_delta, assets_super, design_super ):
        # nothing to do
        pass
    
    # ----------------------------------------------------------------------
    #  PROJECT REPRESENTATION
    # ----------------------------------------------------------------------    
    def __repr__ (self):
        return '<Task> Project %s' % self.signature['NAME']

#: class Project


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  JOB CLASS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
class Job(General_Task):
    
    # ----------------------------------------------------------------------
    #  INITIALIZE JOB
    # ----------------------------------------------------------------------    
    def __init__( self, assets_super, config_delta):
        ''' initialize class attributes 
            will start in folder_super
            should not create any new files because
            the class may be built temporarily 
            as a comparison tester '''     
                
        General_Task.__init__(self)
        
        # initialize config_current
        self.config_current = libSU2.Get_ConfigParams( assets_super['config'] )
        self.config_current.update(config_delta)          
        
        # define assets keys to pull
        self.assets_pull = ['config','mesh','direct','adjoint'] 
        
        # equivalent area files
        if self.config_current.get('EQUIV_AREA','NO') == 'YES':
            self.assets_pull.append('targetea')
            
        # task lists
        self.tasks_todo = []
        self.tasks_done = []
        
        # set signature
        self.signature = self.set_signature(config_delta)
        
        # set folder name
        self.folder_self = 'JOBS/JOB_%03d' % self.signature['JOB_NUMBER']
        
        return
    
    # ----------------------------------------------------------------------
    #  SET JOB SIGNATURE
    # ----------------------------------------------------------------------    
    def set_signature( self, config_delta ):
        ''' assign the properties that represent
            this task's configuration '''
        
        # default signature items
        signature   = {'TYPE':'Job'}
        
        # list of keys to make sure to copy from config delta
        config_keys = ['JOB_NUMBER','VARIABLES']
        
        # assign config signature items
        for key in config_keys:
            signature[key] = self.config_current[key]
       
        # copy all additional config_delta 
        signature.update(config_delta)
       
        # list needed for job comparison, so we don't have to use np.all()
        if isinstance( signature['VARIABLES'] , numpy.ndarray):
            signature['VARIABLES'] = signature['VARIABLES'].tolist()
        
        return signature        
        
    # ----------------------------------------------------------------------
    #  CONFIGURE JOB
    # ----------------------------------------------------------------------    
    def configure( self, config_delta, assets_super, design_super ):      
        ''' configure any control and data inputs 
            for the analysis task '''
        
        # assume self.config_current, self.design_current loaded in self.pull_assets()        
        
        # interpret variables option
        if config_delta.has_key('VARIABLES'):
            # add to design
            self.design_current['VARIABLES'] = config_delta['VARIABLES']
            # set mesh deformation configs
            config_delta['DV_VALUE_NEW'] = list( config_delta['VARIABLES'] )
            config_delta['DV_VALUE_OLD'] = numpy.zeros( len(config_delta['VARIABLES']) ).tolist()
            # delete from delta
            del config_delta['VARIABLES']
            
        # set deformation parameters
        Definition_DV = self.config_current['DEFINITION_DV']
        config_delta['DV_KIND']   = Definition_DV['KIND']
        config_delta['DV_MARKER'] = Definition_DV['MARKER'][0]
        config_delta['DV_PARAM']  = Definition_DV['PARAM']                    
        
        # update config file and self.
        libSU2.Set_ConfigParams( self.assets_current['config'], config_delta )
        self.config_current.update(config_delta);
        
        # -----------------------------------------------------------------------
        #  TASK LIST BUILDING
        
        # build list of tasks to run
        # check for dependent tasks
        
        # task maps
        task_map = tasks_su2.get_classmap()
        task_dep = tasks_su2.get_dependency(self.config_current)        
        
        # tasks to run
        self.tasks_todo = []
        tasks_todo_names = self.config_current['TASKS']
        
        # process new tasks
        for this_name in tasks_todo_names:
            
            # adaptation done within direct or adjoint
            if this_name in ['ADAPT','ADAPT_ADJ','DECOMP']:
                continue
                  
            # a potential new task
            Constructor = task_map[this_name]
            New_Task = Constructor(self.assets_current,config_delta)
            
            # check for previous runs
            Tasks_Check = [ task for task in self.tasks_done if New_Task == task ]
            if Tasks_Check:  # task already done
                if New_Task in self.tasks_todo:  # task already scheduled, delete
                    del New_Task; continue
                # else ...
                assert len(Tasks_Check) <= 1 , 'found more than one completed task in job'
                del New_Task # incase of pointer leak?
                New_Task = Tasks_Check[0]
            # else ...
            
            # check for dependent tasks
            tasks_dep_constr = task_dep[this_name]
            for Constructor in tasks_dep_constr:
                
                # the constructor
                New_Task_Dep = Constructor(self.assets_current,config_delta)
                
                # check for done tasks
                Tasks_Dep_Check = [ task for task in self.tasks_done if New_Task_Dep == task ]
                if Tasks_Dep_Check:
                    if New_Task_Dep in self.tasks_todo: # dependent task already scheduled
                        del New_Task_Dep; continue
                    # else:
                    assert len(Tasks_Dep_Check) <= 1 , 'found more than one completed task in job'
                    # reassign new task pointer
                    del New_Task_Dep 
                    New_Task_Dep = Tasks_Dep_Check[0]
                # else ...
                
                # add dependent task todo
                self.tasks_todo.append(New_Task_Dep) 
                if not Tasks_Dep_Check: self.tasks_done.append(New_Task_Dep)
                
            #: for dependents
            
            # with all dependent tasks added, now add the requested
            self.tasks_todo.append(New_Task) 
            if not Tasks_Check: self.tasks_done.append(New_Task)
            
        #: for new tasks
        
        return
    
    # ----------------------------------------------------------------------
    #  RUN JOB
    # ----------------------------------------------------------------------
    def run( self,config_delta, assets_super, design_super ):
        ''' run the analysis task and return an 
            updated design dictionary '''
        
        # the updated design
        design_return = {}
        
        # Run the Tasks!!!
        for This_Task in self.tasks_todo:
            
            # run this task
            if not self.config_current['CONSOLE'] == 'QUIET':
                sys.stdout.write( '  %s \n' % This_Task.signature['NAME'] )
            design_new,assets_new,_ = This_Task.evaluate(config_delta,self.assets_current,{})
            
            # update design and assets
            self.design_current.update( design_new )
            self.assets_current.update( assets_new )
            
        #: for each Task
        
        # return newly updated design
        return self.design_current

    # ----------------------------------------------------------------------
    #  PACKUP JOB
    # ----------------------------------------------------------------------    
    def packup( self, config_delta, assets_super, design_super ):
        ''' assign assets to push to super folder 
            and clean uneeded files '''
            
        # no assets to push, dont modify project folder
        assets_push = {}
        config_return = {}
        
        # update pull list
        for key in self.assets_current.keys():
            if self.assets_pull.has_key(key):
                del self.assets_pull[key]        
        
        return assets_push, config_return
    
    # ----------------------------------------------------------------------
    #  CLEANUP PROJECT
    # ----------------------------------------------------------------------    
    def cleanup( self, config_delta, assets_super, design_super ):
        # nothing to do
        pass
    
    # ----------------------------------------------------------------------
    #  JOB REPRESENATION
    # ----------------------------------------------------------------------    
    def __repr__ (self):
        return '<Task> Job %i' % self.signature['JOB_NUMBER']    

#: class Job


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  MAIN
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
def main():
        
    # Command Line Options
    parser = OptionParser()
    parser.add_option("-c", "--config", dest="config",
                      help="read config from FILE", metavar="FILE")
    parser.add_option("-d", "--design", dest="design",
                      help="read design from FILE", metavar="FILE")
    parser.add_option("-p", "--project", dest="project",
                      help="save project to FILE", metavar="FILE")    

    (options, args)=parser.parse_args()

    design_init = { 'VARIABLES'  : [] ,
                    'OBJECTIVES' : {} ,
                    'GRADIENTS'  : {}  }
    libSU2.save_data(options.design,design_init)

    # start project
    The_Project = Project( config_name = options.config ,
                           design_name = options.design  )
    
    # assume no change to config file
    config_baseline = libSU2.Get_ConfigParams(options.config)
    DV_N = len( config_baseline['DEFINITION_DV']['MARKER'] )
    
    config_delta = { 'VARIABLES' : numpy.zeros( DV_N ) }
    
    # evaluate project
    design_new,_,_ = The_Project.evaluate(config_delta)
    
    # save matlab data
    design_name_mat = os.path.splitext( options.design )[0] + '.mat'
    libSU2.save_data(design_name_mat,design_new)
    
    # save project
    libSU2.save_data( options.project , The_Project )
    
    print '\nDONE !'

#: main()


if __name__ == '__main__':
    main()
    
   
