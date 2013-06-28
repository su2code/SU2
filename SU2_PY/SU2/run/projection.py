
import os, sys, shutil, copy

from .. import io   as su2io
from .. import util as su2util
from decompose import decompose as su2decomp
from interface import GPC as SU2_GPC

def projection ( config      ,
                 step = 1e-3  ):
    
    # local copy
    konfig = copy.deepcopy(config)
    
    # decompose
    su2decomp(konfig)
        
    # choose dv values 
    Definition_DV = konfig['DEFINITION_DV']
    n_DV          = len(Definition_DV['KIND'])
    if isinstance(step,list):
        assert len(step) == n_DV , 'unexpected step vector length'
    else:
        step = [step]*n_DV
    dv_old = [0.0]*n_DV # SU2_GPC input requirement, assumes linear superposition of design variables
    dv_new = step
    konfig.unpack_dvs(dv_new,dv_old)
    
    # Run Projection
    SU2_GPC(konfig)
    
    # filenames
    objective      = konfig['ADJ_OBJFUNC']    
    grad_filename  = konfig['GRAD_OBJFUNC_FILENAME']
    output_format  = konfig['OUTPUT_FORMAT']
    plot_extension = su2io.get_extension(output_format)
    adj_suffix     = su2io.get_adjointSuffix(objective)
    grad_plotname  = os.path.splitext(grad_filename)[0] + '_' + adj_suffix + plot_extension    
    
    # read raw gradients
    raw_gradients = su2io.read_gradients(grad_filename)
    os.remove(grad_filename)
    
    # Write Gradients
    data_plot = su2util.ordered_bunch()
    data_plot['VARIABLE']         = range(len(raw_gradients)) 
    data_plot['GRADIENT']     = raw_gradients             
    data_plot['FINDIFF_STEP'] = step                       
    su2util.write_plot(grad_plotname,output_format,data_plot)

    # gradient output dictionary
    gradients = { objective : raw_gradients }
    
    # info out
    info = su2io.State()
    info.GRADIENTS.update( gradients )
    
    return info
