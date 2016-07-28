#!/usr/bin/env python
# -*-coding:utf-8 -* 

# \file FSI_config.py
#  \brief Python class for handling configuration file for FSI computation.
#  \author THOMAS David, University of Liege, Belgium.  Department of Aerospace and Mechanical Engineering
#  \version BETA

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import os, sys, shutil, copy
from ..util import switch

# ----------------------------------------------------------------------
#  FSI Configuration Class
# ----------------------------------------------------------------------

class FSIConfig:
    """
    Class that contains all the parameters coming from the FSI configuration file.
    Read the file and store all the options into a dictionary.
    """

    def __init__(self,FileName):
        self.ConfigFileName = FileName
        self._ConfigContent = {}
        self.readConfig()

    def __str__(self):
	tempString = str()
	for key, value in self._ConfigContent.items():
	   tempString += "{} = {}\n".format(key,value)
	return tempString

    def __getitem__(self,key):
	return self._ConfigContent[key]

    def __setitem__(self, key, value):
	self._ConfigContent[key] = value

    def readConfig(self):
        input_file = open(self.ConfigFileName)
        while 1:
	    line = input_file.readline()
	    if not line:
	        break
            # remove line returns
            line = line.strip('\r\n')
            # make sure it has useful data
            if (not "=" in line) or (line[0] == '%'):
                continue
            # split across equal sign
            line = line.split("=",1)
            this_param = line[0].strip()
            this_value = line[1].strip()

            for case in switch(this_param):
	        #integer values
		#if case("NDIM")			      : pass
	        #if case("MESH_DEF_LIN_ITER")	      : pass
	        #if case("MESH_DEF_NONLIN_ITER")       : pass
		if case("RESTART_ITER")		      : pass
		if case("NB_EXT_ITER")		      : pass
	        if case("NB_FSI_ITER")		      :
		    self._ConfigContent[this_param] = int(this_value)
		    break

	        #float values
		if case("AITKEN_PARAM")		      : pass
		if case("START_TIME")		      : pass
		if case("UNST_TIMESTEP")	      : pass
		if case("UNST_TIME")		      : pass
	        if case("FSI_TOLERANCE")	      :
		    self._ConfigContent[this_param] = float(this_value)
		    break

	        #string values
		if case("CFD_CONFIG_FILE_NAME")	      : pass
		if case("CSD_SOLVER")		      : pass
		if case("CSD_CONFIG_FILE_NAME")	      : pass
		if case("RESTART_SOL")		      : pass
		if case("MATCHING_MESH")	      : pass
		if case("DISP_PRED")		      : pass
		if case("AITKEN_RELAX")               : pass
	        if case("UNSTEADY_SIMULATION")	      : pass
		if case("INTERNAL_FLOW")	      : 
	        #if case("MESH_DEF_METHOD")	      : pass
		    self._ConfigContent[this_param] = this_value
		    break

 	        if case():
		    print(this_param + " is an invalid option !")
		    break
	    #end for
	


    #def dump()
