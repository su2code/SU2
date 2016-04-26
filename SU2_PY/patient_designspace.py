#!/usr/bin/env python 

## \file patient_designspace.py
#  \brief Python script for running multiple design configurations in multiple sessions
#  \author T. Lukaczyk
#  \version 4.1.2 "Cardinal"
#
# SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
#                      Dr. Thomas D. Economon (economon@stanford.edu).
#
# SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
#                 Prof. Piero Colonna's group at Delft University of Technology.
#                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
#                 Prof. Rafael Palacios' group at Imperial College London.
#
# Copyright (C) 2012-2016 SU2, the open-source CFD code.
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
from optparse import OptionParser
sys.path.append(os.environ['SU2_RUN'])
import SU2

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
    SU2.opt.server ( options.config_file       ,
                     options.design_file       , 
                     options.project_file      ,
                     options.transfer_file     ,
                     options.exchange_location , 
                     options.hot_start          )
    
    return 

#: def main()


# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()
    
    
    
    




