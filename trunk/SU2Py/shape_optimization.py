#!/usr/bin/env python 

## \file shape_optimization.py
#  \brief Python script for performing the shape optimization.
#  \author Francisco Palacios, Tom Economon, Trent Lukaczyk, Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
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

import time, os, sys, shutil, subprocess, numpy, libSU2, copy
from optparse import OptionParser
from scipy.optimize import fmin_bfgs
from scipy.optimize import fmin_slsqp
from parallel_computation import parallel_computation
from continuous_adjoint   import continuous_adjoint
from finite_differences   import finite_differences
from parallel_deformation import parallel_deformation

SU2_RUN = os.environ['SU2_RUN'] 
sys.path.append( SU2_RUN )

counter = 0

parser=OptionParser()
parser.add_option("-f", "--file", dest="filename",
                  help="read config from FILE", metavar="FILE")
parser.add_option("-p", "--partitions", dest="partitions", default=1,
                  help="number of PARTITIONS", metavar="PARTITIONS")
parser.add_option("-g", "--gradient", dest="gradient", default="Adjoint",
                  help="Method for computing the GRADIENT (Adjoint or FinDiff)", metavar="GRADIENT")
parser.add_option("-q", "--quiet", dest="quiet", default="False",
                  help="True/False Quiet all SU2 output (optimizer output only)", metavar="QUIET")
parser.add_option("-c", "--cycle", dest="cycle", default=0,
                  help="number of mesh adaptation CYCLEs", metavar="CYCLE")
parser.add_option("-s", "--step", dest="findiffstep", default=1e-4,
                  help="finite difference STEP", metavar="STEP")

(options, args)=parser.parse_args()

# process inputs
options.partitions = int( options.partitions )
options.cycle      = int( options.cycle )
options.findiffstep = float( options.findiffstep )
parallel = options.partitions > 1
if options.quiet == 'True':
    report = 'quiet'
else: 
    report = 'verbose'

# Read design variable information
params_base = libSU2.Get_ConfigParams( options.filename )

output_format = params_base.get('OUTPUT_FORMAT','TECPLOT')
Definition_DV = params_base.get('DEFINITION_DV')
objfunc       = params_base.get('OBJFUNC'      )
scale         = params_base.get('OBJFUNC_SCALE')
original_grid = params_base.get('MESH_FILENAME')
adaptation    = options.cycle > 0

UnitaryEqConst       = params_base.get('CONST_EQ'      )
UnitaryEqConst_Scale = params_base.get('CONST_EQ_SCALE')
UnitaryEqConst_Value = params_base.get('CONST_EQ_VALUE')

UnitaryIeqConst       = params_base.get('CONST_IEQ'      )
UnitaryIeqConst_Scale = params_base.get('CONST_IEQ_SCALE')
UnitaryIeqConst_Value = params_base.get('CONST_IEQ_VALUE')
UnitaryIeqConst_Sign  = params_base.get('CONST_IEQ_SIGN' )

n_DV = len(Definition_DV['KIND'])

# special physics cases
special_cases = libSU2.get_SpecialCases(params_base)

x_old  = numpy.ones(n_DV) # just to get started
x_new  = numpy.zeros(n_DV)
x_incr = numpy.zeros(n_DV)

Objectives_OLD = {}
Gradients_OLD  = {}

if (output_format == "PARAVIEW"):
    History_file      = "history_opt_"  + options.filename.replace(".cfg",".csv")
    optimization_file = "optimization_" + options.filename.replace(".cfg",".csv")
elif (output_format == "TECPLOT"):
    History_file      = "history_opt_"  + options.filename.replace(".cfg",".plt")
    optimization_file = "optimization_" + options.filename.replace(".cfg",".plt")

Config_MDC_file  = "config_opt_MDC_"  + options.filename
Config_MAC_file  = "config_opt_MAC_"  + options.filename
Config_CFD_file  = "config_opt_CFD_"  + options.filename
Config_DDC_file  = "config_opt_DDC_"  + options.filename
Config_GRAD_file = "config_opt_GRAD_" + options.filename
Mesh_MDC_file    = "mesh_opt_MDC_" + options.filename.replace(".cfg",".su2")

# Prepare Config Files

# MDC
params_set = { 'DV_KIND'           : params_base['DEFINITION_DV']['KIND']      ,
               'DV_MARKER'         : params_base['DEFINITION_DV']['MARKER'][0] , 
               'DV_PARAM'          : params_base['DEFINITION_DV']['PARAM']     ,
               'DV_VALUE_OLD'      : numpy.zeros(n_DV)                         ,  # will always deform from baseline mesh, needed to enable adaptation
               'MESH_OUT_FILENAME' : Mesh_MDC_file                              }
shutil.copy(options.filename,Config_MDC_file)
libSU2.Set_ConfigParams(Config_MDC_file,params_set)   

# MAC
params_set = { 'NUMBER_PART'   : options.partitions ,
               'MESH_FILENAME' : Mesh_MDC_file       }    
shutil.copy(options.filename,Config_MAC_file)   
libSU2.Set_ConfigParams(Config_MAC_file,params_set)    

# CFD
params_set = { 'MATH_PROBLEM'  : 'DIRECT'      ,
               'MESH_FILENAME' : Mesh_MDC_file ,
               'CONV_FILENAME' : History_file.replace(".csv","").replace(".plt","") }
shutil.copy(options.filename,Config_CFD_file)
libSU2.Set_ConfigParams(Config_CFD_file,params_set)   

# GRAD
params_set = { 'MESH_FILENAME' : Mesh_MDC_file ,
               'CONV_FILENAME' : History_file.replace(".csv","").replace(".plt","") }
shutil.copy(options.filename,Config_GRAD_file)
libSU2.Set_ConfigParams(Config_GRAD_file,params_set) 

# DDC
params_set = { 'NUMBER_PART'       : options.partitions ,
               'MESH_FILENAME'     : Mesh_MDC_file       }
shutil.copy(options.filename,Config_DDC_file)
libSU2.Set_ConfigParams( Config_DDC_file, params_set )


# Get header and write format information
header,header_vars,write_format = libSU2.get_OptFileFormat( output_format, special_cases )

# Start optimization file
Opt_file = open(optimization_file,'w')
Opt_file.write('TITLE = "SU2 Optimization"\n')
Opt_file.write(header)
Opt_file.close()



# -------------------------------------------------------------------
# Function for updating the computational grid
# -------------------------------------------------------------------

def update_grid(x):
    #print( 'grid  ') # + str(x) )

    # Copy all the design variables 
    global x_new, x_old, Gradients_OLD, Objectives_OLD
    x_new = [ xx for xx in x ] # avoiding reference copy
    
    delta_x = numpy.array(x_new) - numpy.array(x_old)
    if numpy.dot(delta_x,delta_x) == 0.0:
        #print( 'skip update')
        return False
    
    # Change the parameters of the design variables
    params_set = { 'DV_VALUE_NEW'      : x_new  } #'DV_VALUE_OLD'      : x_old ,
    libSU2.Set_ConfigParams(Config_MDC_file,params_set)       
    
    # Apply grid update
    # print('  SU2_MDC')
    logfile = 'SU2_MDC.out'
    
    if parallel and not adaptation: 
        if options.quiet == 'False':
            logfile = None
            
        parallel_deformation( Config_MDC_file        ,
                              partitions  = options.partitions , 
                              divide_grid = "False"  ,
                              merge_grid  = "True"  ,
                              logfile     = logfile   )        
    else:    
        run_SU2_MDC = "SU2_MDC " + Config_MDC_file
        if options.quiet == 'True':
            run_SU2_MDC = run_SU2_MDC + ' >> ' + logfile
        os.system ( run_SU2_MDC )
        
        
    if adaptation:
        
        #print( '  adapt')
        
        run_adapt = 'mesh_adaptation.py -f %s -p %s -c %s -o True' % ( Config_MAC_file, options.partitions, options.cycle ) 
        run_DDC   = 'SU2_DDC ' + Config_MAC_file 
        if options.quiet == 'True':
            run_adapt = run_adapt + ' >> ' + logfile
            run_DDC   = run_DDC   + ' >> ' + logfile
  
        # ADAPT MESH (will overwrite mesh file)
        os.system( run_adapt )
        
        # decompose new mesh
        os.system( run_DDC )
        

    # Copy solution
    x_old = [ xx for xx in x_new ] # also copy.deepcopy(x_new) works
    
    Objectives_OLD = {}
    Gradients_OLD  = {}    
    
    return True


# -------------------------------------------------------------------
# Function for computing the objective function
# -------------------------------------------------------------------

def f(x):
    #print( 'f     ') # + str(x))

    # Update the computational grid to the design variable state
    updated = update_grid(x)

    # CFD computation to evaluate the objective function
    global counter, Objectives_OLD
    counter = counter + 1

    # run direct problem
    if updated:
        ObjFun_Dict,_ = continuous_adjoint ( filename    = Config_CFD_file    ,
                                             partitions  = options.partitions ,
                                             compute     = 'Direct'           ,
                                             output      = 'True'             ,
                                             divide_grid = 'False'            ,
                                             report      = report             ,
                                             logfile     = 'SU2_CFD.out'       )
        Objectives_OLD = ObjFun_Dict
    else:
        ObjFun_Dict = Objectives_OLD
    
    # copy mesh file
    shutil.copyfile( Mesh_MDC_file, "mesh_ShapeIter_%03d.su2" % counter )

    # To save memory, don't store each iteration for unsteady problems
    if not 'WRT_UNSTEADY' in special_cases: 
        rename_data()
    
    # setup variables for writing
    vars_write = ObjFun_Dict
    vars_write['ITERATION'] = counter  # set to optimization iteration
    
    # values for writing
    write_vals = [ vars_write[key] for key in header_vars ] 
    write_vals = tuple(write_vals)
    
    # Output optimization history
    Opt_file = open(optimization_file,"a")
    Opt_file.write(write_format % write_vals)
    Opt_file.close()    
    
    # return objective function, scaled and sign fixed for maximization problems
    f_val =  ObjFun_Dict[objfunc] * scale * libSU2.get_ObjFunSign(objfunc)

    return f_val



# -------------------------------------------------------------------
# Function for computing the derivative of the objective function
# -------------------------------------------------------------------

def df(x):
    #print( 'df    ') # + str(x))
    
    global Objectives_OLD, Gradients_OLD

    # Update the computational grid to the design variable state
    updated = update_grid(x)
    
    if not Objectives_OLD.has_key(objfunc):
        ObjFun_Dict,_ = continuous_adjoint ( filename    = Config_CFD_file    ,
                                             partitions  = options.partitions ,
                                             compute     = 'Direct'           ,
                                             output      = 'True'            ,
                                             divide_grid = 'False'            ,
                                             report      = report             ,
                                             logfile     = 'SU2_CFD.out'       )
        Objectives_OLD = ObjFun_Dict
        
    if not Gradients_OLD.has_key(objfunc):
        if options.gradient == 'Adjoint':
            # Change grid file for Gradient Computation over the deformed grid
            params_set = { 'ADJ_OBJFUNC'  : objfunc       }
            libSU2.Set_ConfigParams(Config_GRAD_file,params_set) 
        
            # Compute Gradient using continuous adjoint strategy
            _,Gradient_Dict = continuous_adjoint ( filename    = Config_GRAD_file   ,
                                                   partitions  = options.partitions ,
                                                   compute     = 'Adjoint'          , 
                                                   output      = 'False'            ,
                                                   divide_grid = 'False'            ,
                                                   report      = report             ,
                                                   logfile     = 'SU2_ADJ.out'       )             
        elif options.gradient == 'FinDiff':
            _,Gradient_Dict =  finite_differences( filename    = Config_GRAD_file   ,
                                                   partitions  = options.partitions ,
                                                   output      = 'False'            ,
                                                   step        = options.findiffstep,
                                                   logfile     = 'SU_FD.out'       )  
        for key,value in Gradient_Dict.iteritems():
            Gradients_OLD[key] = value
    else:
        Gradient_Dict = Gradients_OLD
    
    # pack gradients, apply scale and sign
    df_vals = numpy.array( Gradient_Dict[ objfunc ] ) * scale * libSU2.get_ObjFunSign(objfunc)
            
    return df_vals
    


# -------------------------------------------------------------------
# Function for computing the equality constraint
# -------------------------------------------------------------------

def eqcons(x):
    #print( 'ceq   ') # + str(x))

    if UnitaryEqConst[0] != "NONE" :
        
        global Objectives_OLD

        # Update the computational grid to the design variable state
        updated = update_grid(x)
        
        # run direct problem
        if updated:
            ObjFun_Dict,_ = continuous_adjoint ( filename    = Config_CFD_file    ,
                                                 partitions  = options.partitions ,
                                                 compute     = 'Direct'           ,
                                                 output      = 'True'            ,
                                                 divide_grid = 'False'            ,
                                                 report      = report             ,
                                                 logfile     = 'SU2_CFD.out'       )
            Objectives_OLD = ObjFun_Dict
        else:
            ObjFun_Dict = Objectives_OLD

        # equality constraint vector
        ceq = numpy.zeros(len(UnitaryEqConst))
        
        # assemble constraint vector
        for i in range(len(UnitaryEqConst)) :
            this_con = UnitaryEqConst[i]
            this_val = UnitaryEqConst_Value[i]
            this_scl = UnitaryEqConst_Scale[i]
            
            # apply scale
            ceq[i] = ( ObjFun_Dict[this_con] - this_val ) * this_scl
            
    else:
        ceq = numpy.zeros(0)

    return ceq       
             
             

# -------------------------------------------------------------------
# Function for computing the derivative of the equality constraint
# -------------------------------------------------------------------

def deqcons(x):
    #print( 'dceq  ') # + str(x))

    # Read the conts function gradient
    dceq = numpy.zeros((len(UnitaryEqConst), n_DV))

    for i in range(len(UnitaryEqConst)) :
        if UnitaryEqConst[i] != "NONE" :
            
            
            global Objectives_OLD, Gradients_OLD
        
            # Update the computational grid to the design variable state
            updated = update_grid(x)
            
            if not Objectives_OLD.has_key(UnitaryEqConst[i]):
                ObjFun_Dict,_ = continuous_adjoint ( filename    = Config_CFD_file    ,
                                                     partitions  = options.partitions ,
                                                     compute     = 'Direct'           ,
                                                     output      = 'True'            ,
                                                     divide_grid = 'False'            ,
                                                     report      = report             ,
                                                     logfile     = 'SU2_CFD.out'       )
                Objectives_OLD = ObjFun_Dict           
        
            if not Gradients_OLD.has_key(UnitaryEqConst[i]):
                if options.gradient == 'Adjoint':
                    # Change grid file for Gradient Computation over the deformed grid
                    params_set = { 'ADJ_OBJFUNC'  : UnitaryEqConst[i] }
                    libSU2.Set_ConfigParams(Config_GRAD_file,params_set) 
                
                    # Compute Gradient using continuous adjoint strategy
                    _,Gradient_Dict = continuous_adjoint ( filename    = Config_GRAD_file   ,
                                                           partitions  = options.partitions ,
                                                           compute     = 'Adjoint'          , 
                                                           output      = 'False'            ,
                                                           divide_grid = 'False'            ,
                                                           report      = report             ,
                                                           logfile     = 'SU2_ADJ.out'       )             
                elif options.gradient == 'FinDiff':
                    _,Gradient_Dict =  finite_differences( filename    = Config_GRAD_file   ,
                                                           partitions  = options.partitions ,
                                                           output      = 'False'            ,
                                                           step        = options.findiffstep,
                                                           logfile     = 'SU2_FD.out'       )  
                for key,value in Gradient_Dict.iteritems():
                    Gradients_OLD[key] = value
            else:
                Gradient_Dict = Gradients_OLD            
                        
            # apply scale 
            dceq[i,:] = numpy.array( Gradient_Dict[ UnitaryEqConst[i] ] ) * UnitaryEqConst_Scale[i]
            
        else:
            dceq = numpy.zeros((0, n_DV))

    return dceq



# -------------------------------------------------------------------
# Function for computing the inequality constraint
# -------------------------------------------------------------------

def ieqcons(x):
    #print( 'cieq  ') # + str(x))

    if UnitaryIeqConst[0].strip() != "NONE" :
        
        global Objectives_OLD

        # Update the computational grid to the design variable state
        updated = update_grid(x)
        
        # run direct problem
        if updated:
            ObjFun_Dict,_ = continuous_adjoint ( filename    = Config_CFD_file    ,
                                                 partitions  = options.partitions ,
                                                 compute     = 'Direct'           ,
                                                 output      = 'True'            ,
                                                 divide_grid = 'False'            ,
                                                 report      = report             ,
                                                 logfile     = 'SU2_CFD.out'       )
            Objectives_OLD = ObjFun_Dict
        else:
            ObjFun_Dict = Objectives_OLD
        
        # equality constraint vector
        cieq = numpy.zeros(len(UnitaryIeqConst))
        
        # build the vector
        for i in range(len(UnitaryIeqConst)) :
            this_con = UnitaryIeqConst[i]
            this_val = UnitaryIeqConst_Value[i]
            this_sgn = UnitaryIeqConst_Sign[i]
            this_scl = UnitaryIeqConst_Scale[i]
            
            if   this_sgn == 'GREATER': this_sgn =  1.0
            elif this_sgn == 'LESS'   : this_sgn = -1.0
            else : raise Exception('unknown inequality constraint sign')
            
            # apply scale and sign
            cieq[i] = ( ObjFun_Dict[this_con] - this_val ) * this_scl * this_sgn
            
    else:
        cieq = numpy.zeros(0)

    return cieq      



# -------------------------------------------------------------------
# Function for computing the derivative of the inequality constraint
# -------------------------------------------------------------------

def dieqcons(x):
    #print( 'dcieq ') # + str(x))

    # Read the conts function gradient
    dcieq = numpy.zeros((len(UnitaryIeqConst), n_DV))

    for i in range(len(UnitaryIeqConst)) :
        if UnitaryIeqConst[i] != "NONE" :
            
            
            global Objectives_OLD, Gradients_OLD
        
            # Update the computational grid to the design variable state
            updated = update_grid(x)
            
            if not Objectives_OLD.has_key(UnitaryIeqConst[i]):
                ObjFun_Dict,_ = continuous_adjoint ( filename    = Config_CFD_file    ,
                                                     partitions  = options.partitions ,
                                                     compute     = 'Direct'           ,
                                                     output      = 'True'            ,
                                                     divide_grid = 'False'            ,
                                                     report      = report             ,
                                                     logfile     = 'SU2_CFD.out'       )
                Objectives_OLD = ObjFun_Dict       
        
            if not Gradients_OLD.has_key(UnitaryIeqConst[i]):
                if options.gradient == 'Adjoint':
                    # Change grid file for Gradient Computation over the deformed grid
                    params_set = { 'ADJ_OBJFUNC'  : UnitaryIeqConst[i] }
                    libSU2.Set_ConfigParams(Config_GRAD_file,params_set) 
                
                    # Compute Gradient using continuous adjoint strategy
                    _,Gradient_Dict = continuous_adjoint ( filename    = Config_GRAD_file   ,
                                                           partitions  = options.partitions ,
                                                           compute     = 'Adjoint'          , 
                                                           output      = 'False'            ,
                                                           divide_grid = 'False'            ,
                                                           report      = report             ,
                                                           logfile     = 'SU2_ADJ.out'       )             
                elif options.gradient == 'FinDiff':
                    _,Gradient_Dict =  finite_differences( filename    = Config_GRAD_file   ,
                                                           partitions  = options.partitions ,
                                                           output      = 'False'            ,
                                                           step        = options.findiffstep,
                                                           logfile     = 'SU2_FD.out'       )  
                for key,value in Gradient_Dict.iteritems():
                    Gradients_OLD[key] = value
            else:
                Gradient_Dict = Gradients_OLD              
              
            if   UnitaryIeqConst_Sign[i] == 'GREATER': this_sign =  1.0
            elif UnitaryIeqConst_Sign[i] == 'LESS'   : this_sign = -1.0
            else : raise Exception('unknown inequality constraint sign')            
            
            # apply scale and sign
            dcieq[i,:] = numpy.array( Gradient_Dict[ UnitaryIeqConst[i] ] ) * UnitaryIeqConst_Scale[i] * this_sign

        else:
            dcieq = numpy.zeros((0, n_DV))

    return dcieq

#: def dieqcons()


    
# -------------------------------------------------------------------
# Function for re-naming data to save
# -------------------------------------------------------------------

def rename_data():
      
    params_get = libSU2.Get_ConfigParams( options.filename )
      
    # list of data names to rename
    append_name = "_ShapeIter_%03d" % counter
    output_names = { 'surface'   : [ params_get['VOLUME_FLOW_FILENAME']  , append_name ] ,
                     'flow'      : [ params_get['SURFACE_FLOW_FILENAME'] , append_name ] ,
                     'chunk'     : [ "deformed_chunk"                    , append_name ] ,
                     'nearfield' : [ "nearfield_flow"                    , append_name ] ,
                     'levelset'  : [ "LevelSet"                          , append_name ]  }
    
    # make second name with appended suffix
    for key in output_names.keys():
        output_names[key][1] = output_names[key][0] + output_names[key][1]
                
    # add extensions
    if output_format == "PARAVIEW" :
        # VTK FORMAT
        this_ext = '.vtk'
        for key in ['surface','flow','chunk','levelset']:
            output_names[key] = [ name + this_ext for name in output_names[key] ]
        # CSV FORMAT
        this_ext = '.csv'
        for key in ['nearfield']:
            output_names[key] = [ name + this_ext for name in output_names[key] ]
    elif output_format == "TECPLOT":
        # PLT FORMAT
        this_ext = '.plt'
        for key in output_names.keys():
            output_names[key] = [ name + this_ext for name in output_names[key] ]        
    #: if output format
    
    # build list of names to rename based on special cases
    keys_to_rename = ['surface','flow']
    if Definition_DV['KIND'][0][0] in ["FFD_CONTROL_POINT","FFD_DIHEDRAL_ANGLE","FFD_TWIST_ANGLE","FFD_ROTATION","FFD_CAMBER","FFD_THICKNESS","FFD_VOLUME"] :
        keys_to_rename.append('chunk')
    if 'EQUIV_AREA' in special_cases:
        keys_to_rename.append('nearfield')
    if 'FREE_SURFACE' in special_cases:
        keys_to_rename.append('levelset')           

    # ok now rename files...
    for key in keys_to_rename:
        os.rename(output_names[key][0],output_names[key][1])
        if output_format == "TECPLOT":
            add_tecplot_iter( output_names[key][1] , counter )
    #: for keys to rename
    
#: def rename_data()



# -------------------------------------------------------------------
# Function for animating optimization files in Tecplot
# -------------------------------------------------------------------

def add_tecplot_iter(filename,this_counter):

    # Open the solution file for this design iteration
    # and find the "ZONE ..." line
    New_file = open("temp_file.plt","w")
    for line in open(filename):
        if "ZONE" in line:
            aug_line = line.replace('ZONE','ZONE STRANDID=' + str(this_counter) + ', SOLUTIONTIME=' + str(this_counter) + ',',1)
            New_file.write(aug_line)
        else:
            New_file.write(line)
    New_file.close()

    # Overwrite old solution with new solution file
    os.rename("temp_file.plt",filename)

#: def add_tecplot_iter()



# -------------------------------------------------------------------
# Main Optimization Call
# -------------------------------------------------------------------

print "\n"
print "-------------------------------------------------------------------------"
print "|                SU2 Suite (Shape Optimization Script)                  |"
print "-------------------------------------------------------------------------"
print "\n"

print objfunc + " objective function with a scale factor of " + str(scale) + "\n"

# decompose initial mesh
startup_DDC_filename = 'config_DDC_startup.cfg'
shutil.copy(Config_DDC_file,startup_DDC_filename)
libSU2.Set_ConfigParams( startup_DDC_filename, {'MESH_FILENAME':original_grid} )

if parallel:
    # print('  SU2_DDC')
    run_DDC = "SU2_DDC " + startup_DDC_filename
    if options.quiet == "True":
        run_DDC = run_DDC + " >> SU2_DDC.out"
    os.system ( run_DDC )

os.remove(startup_DDC_filename)


# Initial guess at the optimum
x0 = numpy.zeros(n_DV)

# Call the SLSPQ optimizer
fmin_slsqp(f, x0, f_eqcons=eqcons, f_ieqcons=ieqcons, bounds=[], fprime=df,
           fprime_eqcons=deqcons, fprime_ieqcons=dieqcons, args=(), iter=100,
           acc=1e-12, iprint=2, full_output=1, epsilon=1.0e-10)



