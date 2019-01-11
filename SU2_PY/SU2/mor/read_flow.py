import sys 
import numpy as np
import re
import pickle

class Read_Flow:
    """Read flow.dat file and provide method for saving data
    """
    def __init__(self, data_file):
        """Read a flow.dat file
        """
        
        # load variables and flow data
        with open(data_file) as fp:
            line = fp.readline()
            line = fp.readline() # start with 2nd line
            
            # load variables:
            self.variables = line.strip().split('""')
            self.variables[0] = self.variables[0][-1:]
            print('Found variables: {}'.format(self.variables))
            
            # read how many nodes are in mesh
            line = fp.readline()
            line2 = re.split('= |,',line)
            if line2[0] == 'ZONE NODES':
                self.zone_nodes = int(line2[1])
                print('Found ZONE NODES: {}'.format(self.zone_nodes))
            else:
                print('Note: Possible wrong value entered for zone_nodes')
                print('Looking for ZONE NODES and found: {}'.format(line2[0]))
            
            # initialize data dictionary with variable names
            self.data = {}
            for i in range(len(self.variables)):
                self.data[self.variables[i]] = []
                
            
            node = 1
            while node <= self.zone_nodes:
                line = fp.readline()
                line2 = re.split('\t',line)
                for i in range(len(self.variables)):
                    self.data[self.variables[i]].append(float(line2[i]))
                node += 1
        
        # TODO 
        # write function to return x and y locations since flow data isn't stored in a predictable way
        
    def save_data(self, save_file='pickled_data.p'):
        """Save data, pickled
        """
        pickle.dump(self.data, open(save_file, "wb"))
        
    def save_some_data(self, desired_variables, save_file='pickled_some_data.p'):
        """Only save data with provided variable (string)
        """
        pass
        
    def return_data(self, desired_variable='Conservative_1'):
        """Method to return data specific to desired_variable
        """
        if desired_variable in self.data:
            return self.data[desired_variable]
        else:
            print('Warning: Requested data does not exist')
            return 1
            
    def print_variables(self):
        """Prints out variables, useful if you want to check which ones were found
        """
        print('Found variables: {}'.format(self.variables))
        return 1
            
    def get_nhighdim(self):
        return self.zone_nodes
        
        
        
        
        
 