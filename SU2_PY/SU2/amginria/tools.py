# -------------------------------------------------------------------
#  Imports
# -------------------------------------------------------------------

import numpy as np
from itertools import islice
import sys, os

# --- Prescribed mesh complexities, i.e. desired mesh sizes
def get_mesh_sizes(config):
    return config['PYADAP_COMPLEXITY'].strip('()').split(",")

# --- Return size info from a python mesh structure
def return_mesh_size(mesh):
    
    elt_key  = ['xy', 'xyz',  'Triangles', 'Edges', 'Tetrahedra']
    elt_name = {'xy':'vertices', 'xyz':'vertices',  'Triangles':'triangles', 'Edges':'edges', 'Tetrahedra':'tetrahedra'}
    
    tab_out = []
        
    for elt in elt_key:
        
        if elt not in mesh: continue
        
        nbe = len(mesh[elt])
        
        if nbe > 0:
            tab_out.append("%d %s" % (nbe, elt_name[elt]))
    
    return ', '.join(map(str, tab_out))
    
# --- Use the python interface to amg (YES)? Or the exe (NO)?
def get_python_amg(config):
    
    if 'ADAP_PYTHON' not in config:
        return True
    
    if config['ADAP_PYTHON'] == "YES":
        return True
    elif config['ADAP_PYTHON'] == "NO":
        return False
    else:
        sys.stderr.write("## WARNING : Invalid value for ADAP_PYTHON option. Assuming YES.\n")
        return True

# --- How many sub-iterations per mesh complexity
def get_sub_iterations(config):
    return config['PYADAP_SUBITE'].strip('()').split(",")

# --- What residual reduction for each complexity level
def get_residual_reduction(config):
    if 'PYADAP_RESIDUAL_REDUCTION' in config:
        return config['PYADAP_RESIDUAL_REDUCTION'].strip('()').split(",")
    else:
        nRes = len(config['PYADAP_COMPLEXITY'].strip('()').split(","))
        res = []
        for i in range(nRes):
            res.append(config['RESIDUAL_REDUCTION'])      
        return res 
        
# --- How many SU2 solver iterations for each complexity level
def get_adj_iter(config):
    if 'PYADAP_ADJ_ITER' in config:
        return config['PYADAP_ADJ_ITER'].strip('()').split(",")
    else:
        nExt_iter = len(config['PYADAP_COMPLEXITY'].strip('()').split(","))
        ext_iter = []
        for i in range(nExt_iter):
            ext_iter.append(config['ITER'])        
        return ext_iter

def get_flow_iter(config):
    if 'PYADAP_FLOW_ITER' in config:
        return config['PYADAP_FLOW_ITER'].strip('()').split(",")
    else:
        nExt_iter = len(config['PYADAP_COMPLEXITY'].strip('()').split(","))
        flow_iter = []
        for i in range(nExt_iter):
            flow_iter.append(config['ITER'])        
        return flow_iter

def get_cfl(config):
    if 'PYADAP_CFL' in config:
        return config['PYADAP_CFL'].strip('()').split(",")
    else:
        ncfl = len(config['PYADAP_COMPLEXITY'].strip('()').split(","))
        cfl = []
        for i in range(ncfl):
            cfl.append(config['CFL_NUMBER'])        
        return cfl
    
def print_adap_options(config, kwds):
    prt = '\nMesh adaptation options:\n'
    for kwd in kwds:
        if kwd in config:
            prt += kwd + ' : ' + config[kwd] + '\n'
    prt += '\n'
    return prt
    
def get_su2_dim(filename):
    
    meshfile = open(filename,'r')
    
    def mesh_readlines(n_lines=1):
        fileslice = islice(meshfile,n_lines)
        return list(fileslice)
    
    dim = -1
    
    keepon = True
    while keepon:
        
        line = mesh_readlines()
        
        if not line: 
            keepon = False
            break
        
        # fix white space
        line = line[0]
        line = line.replace('\t',' ')
        line = line.replace('\n',' ')
    
        # skip comments
        if line[0] == "%":
            pass
    
        # number of dimensions
        elif "NDIME=" in line:
            # save to SU2_MESH data
            dim = int( line.split("=")[1].strip() )
            keepon = False
        
    return dim

# --- Merge 2 solutions (e.g. primal and dual)
def merge_sol(mesh0, mesh1):
    mesh0['solution'] = np.hstack((mesh0['solution'], \
                                   mesh1['solution'])).tolist()
    mesh0['solution_tag'] = np.hstack((np.array(mesh0['solution_tag']), \
                                       np.array(mesh1['solution_tag']))).tolist()

# --- Split adjoint solution
def split_adj_sol(mesh):
    nsol = len(mesh['solution_tag'])

    adj_sol = dict()

    for i in range(nsol):
        if "Adjoint" in mesh['solution_tag'][i]:
            iAdj = i

            adj_sol['solution'] = np.delete(np.array(mesh['solution']), np.s_[0:iAdj], axis=1).tolist()
            adj_sol['solution_tag'] = np.delete(np.array(mesh['solution_tag']), np.s_[0:iAdj], axis=0).tolist()

            if 'xyz' in mesh:
                adj_sol['xyz'] = mesh['xyz']
            elif 'xy' in mesh:
                adj_sol['xy'] = mesh['xy']

            adj_sol['dimension'] = mesh['dimension']

            mesh['solution'] = np.delete(np.array(mesh['solution']), np.s_[iAdj:nsol], axis=1).tolist()
            mesh['solution_tag'] = np.delete(np.array(mesh['solution_tag']), np.s_[iAdj:nsol], axis=0).tolist()

            break

    return adj_sol

# ---- Process solution to store desired sensor for adaptation
def create_sensor(solution, sensor):
    
    Dim = solution['dimension']
    Sol = np.array(solution['solution'])
    
    if sensor == "MACH":
        
        iMach = solution['id_solution_tag']['Mach']
        sensor = Sol[:,iMach]
        sensor = np.array(sensor).reshape((len(sensor),1))
        sensor_header = ["Mach"]
        
    elif sensor == "PRES":
        
        iPres = solution['id_solution_tag']['Pressure']
        sensor = Sol[:,iPres]
        sensor = np.array(sensor).reshape((len(sensor),1))        
        sensor_header = ["Pres"]
        
    elif sensor == "MACH_PRES":

        iPres  = solution['id_solution_tag']['Pressure']
        iMach  = solution['id_solution_tag']['Mach']
        mach   = np.array(Sol[:,iMach])
        pres   = np.array(Sol[:,iPres])
        sensor = np.stack((mach, pres), axis=1)    
        sensor_header = ["Mach", "Pres"]

    elif sensor == "GOAL":

        if Dim == 2:
            sensor = Sol[:,-3:]
            sensor = np.array(sensor).reshape((len(sensor),3))
        elif Dim == 3:
            sensor = Sol[:,-6:]
            sensor = np.array(sensor).reshape((len(sensor),6))
        sensor_header = ["Goal"]
                
    else :
        sys.stderr.write("## ERROR : Unknown sensor.\n")
        sys.exit(1)
    
    sensor_wrap = dict()
    
    sensor_wrap['solution_tag'] = sensor_header
    sensor_wrap['xyz'] = solution['xyz']
    
    sensor_wrap['dimension']    = solution['dimension']
    sensor_wrap['solution']     = sensor
    
    return sensor_wrap

def plot_results(history_format, filename, iter):
    """ writes a Tecplot or CSV file for plotting adaptation results
    """

    #--- Format and file name
    if (history_format == 'TECPLOT'):
        solname  = 'history.dat'
    else:
        solname  = 'history.csv'
        
    #--- Write header on first adaptive iteration
    if iter == 0:
        #--- Get header from solution history
        header = ''

        if (history_format == 'TECPLOT'):
            header     = 'VARIABLES='
            headerline = 1
        else:
            headerline = 0

        with open(solname, 'rb') as f:
            for i, line in enumerate(f):
                if i == headerline:
                    break

        header = header + line

        plotfile = open(filename,'w')
        plotfile.write(header)
        plotfile.write('\n')

    # --- Append data on all other iterations
    else:
        plotfile = open(filename,'a')

    #--- Get data from last line of file
    with open(solname, 'rb') as f:
        f.seek(-2, os.SEEK_END). 
        while f.read(1) != b'\n':
            f.seek(-2, os.SEEK_CUR) 
        last_line = f.readline().decode()

    plotfile.write(last_line)
    plotfile.write('\n')

    plotfile.close()
        
