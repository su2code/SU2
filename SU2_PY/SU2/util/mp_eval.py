import os
import multiprocessing as mp
import numpy as np
import sys

if sys.version_info[0] > 2:
    # In Py3, range corresponds to Py2 xrange
    xrange = range

class mp_eval(object):
    
    def __init__(self,function,num_procs=None):
        
        self.__name__ = function.__name__
        
        tasks    = mp.JoinableQueue()
        results  = mp.Queue()
        function = TaskMaster(function)
        
        if num_procs is None:
            num_procs = mp.cpu_count()
            
        procs = [ QueueMaster( tasks, results, function ) 
                  for i in xrange(num_procs) ]
        
        self.tasks    = tasks
        self.results  = results
        self.function = function
        self.procs    = procs
        
        return
    
    def __call__(self,inputs):
        
        tasks   = self.tasks
        results = self.results
        
        if isinstance(inputs,np.ndarray):
            n_inputs = inputs.shape[0]
        elif isinstance(inputs,list):
            n_inputs = len(inputs)
        else:
            raise Exception('unsupported input')
        
        for i_input,this_input in enumerate(inputs):
            this_job = { 'index'  : i_input    ,
                         'input'  : this_input ,
                         'result' : None        }
            tasks.put( this_job )
        #end

        # wait for tasks
        tasks.join() 
        
        # pull results
        result_list = [ [] ]*n_inputs
        for i in xrange(n_inputs):
            result = results.get()
            i_result = result['index']
            result_list[i_result] = result['result']
        
        return result_list

    def __del__(self):
        
        for proc in self.procs:
            self.tasks.put(None)
        self.tasks.join()       
        
        return

class QueueMaster(mp.Process):

    def __init__(self,task_queue,result_queue,task_class=None):
        mp.Process.__init__(self)
        self.task_queue   = task_queue
        self.result_queue = result_queue
        self.task_class   = task_class
        self.daemon       = True
        self.start()

    def run(self):
        proc_name = self.name
        parentPID = os.getppid()
        
        while True:
            
            if os.getppid() != parentPID:
                break # parent died

            this_job = self.task_queue.get()

            if this_job is None:
                self.task_queue.task_done()
                break # kill signal
            
            this_input = this_job['input']
            this_task  = self.task_class

            this_data = this_task(*this_input)
            this_job['result'] = this_data
            self.result_queue.put(this_job)
            
            self.task_queue.task_done()
            
        #: while alive
        
        return

class TaskMaster(object):
    
    def __init__(self, func):
        self.func  = func
    def __call__(self, *arg, **kwarg):  
        # makes object callable
        result = self.func(*arg, **kwarg)
        return result
    def __str__(self):
        return '%s' % self.func

    # pickling
    #def __getstate__(self):
        #dict = self.__dict__.copy()
        #data_dict = cloudpickle.dumps(dict)
        #return data_dict

    #def __setstate__(self,data_dict):
        #self.__dict__ = pickle.loads(data_dict)
        #return    
