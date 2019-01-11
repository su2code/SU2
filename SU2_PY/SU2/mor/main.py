#!/usr/bin/env python

## \file hicks_henne.py
#  \brief Python script for performing sweep of hicks_henne bumps
#  \author J. Lauzon
#  \version 6.0.0 
## python hicks_henne.py -f config_filename.cfg

# imports
import os
import numpy as np
import random as rand
import matplotlib.pyplot as plt
from optparse import OptionParser
import os, sys, shutil, copy, os.path
sys.path.append(os.environ['SU2_RUN'])
import SU2
from bcs import isothermal_bc
from read_flow import Read_Flow
import pickle

def main():
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="filename",
                      help="read config from FILE", metavar="FILE")
    parser.add_option("-n", "--partitions", dest="partitions", default=2,
                      help="number of PARTITIONS", metavar="PARTITIONS")
    parser.add_option("-z", "--zones", dest="nzones", default="1",
                      help="Number of Zones", metavar="ZONES")  
                     
    (options, args)=parser.parse_args()
    options.partitions = int( options.partitions )
    options.nzones = int( options.nzones )
    
    # load config, start state
    config = SU2.io.Config(options.filename)
    state  = SU2.io.State()
    
    # find solution files if they exist
    state.find_files(config)
    
    # prepare config
    config.NUMBER_PART = options.partitions
    config.NZONES      = options.nzones
    
    # generate parameters
    mu_high = 15.0
    mu_low = 1.0
    nparams = 4
    nsamples = 50
    nsnaps = 20
    bc_param = isothermal_bc(mu_high, mu_low, nparams, nsamples)
    samples = bc_param.get_samples_list()
    
    # load precomputed latin hypercube samples (found via matlab)
    #dv_values = np.loadtxt('lhsdesign_6param.csv', dtype=np.float64, delimiter=',')
    
    # update parameter and run solution
    for i, sample in enumerate(samples):
        if i == nsnaps-1:
            break
            
        # Start new folder for each design
        folder = 'DESIGNS/DSN_*'
        _design_number = '%03d'
        folder = folder.rstrip('/')+'/'
        if '*' in folder: folder = SU2.io.next_folder(folder)
        print('New Project: %s' % (folder))
        
        # Update param
        config.MARKER_ISOTHERMAL = sample
        
        # Run solution
        sol = Snapshot(config, config.MESH_FILENAME, folder)
        if i == 0:
            flowdata = np.empty([nsnaps, sol.get_nhighdim])
            
        flowdata[i] = sol.get_data() # nsnaps x nhighdim
        
    pickle.dump(flowdata, open('pickled_data.p','wb'))


class Snapshot:
    """
    Creates a snapshot
    """
    def __init__(self, config, mesh_file, folder):
        # setup snapshot session
        pull = []
        link = [mesh_file]
        config.MESH_FILENAME = mesh_file
        self.flowfile = config.VOLUME_FLOW_FILENAME + '.dat'
        
        self.konfig = config
        self.run(pull, link, folder)
        
        self.results = Read_Flow(self.flowfile)
        
    def run(self, pull, link, folder):
        with SU2.io.redirect_folder(folder,pull,link,force=True) as push:
            info = SU2.run.direct(self.konfig)
        
    def get_data(self):
        return np.asarray(self.results.return_data())
        
    def get_nhighdim(self):
        return results.get_nhighdim()
            
        
    
        


 
if __name__ == "__main__":
    main()