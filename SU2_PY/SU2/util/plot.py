
#from ..io.tools import get_gradFileFormat, get_extension, get_specialCases

def write_plot(filename,plot_format,data_plot,keys_plot=[]):
    
    if not keys_plot:
        keys_plot = data_plot.keys()
    
    header = ('"') + ('","').join(keys_plot) + ('"') + (' \n')
    
    if plot_format == 'TECPLOT':
        header = 'VARIABLES=' + header
    
    n_lines = 0
    
    for i,key in enumerate(keys_plot):
        value = data_plot[key]
        if i == 0:
            n_lines = len(value)
        else:
            assert n_lines == len(value) , 'unequal plot vector lengths'
        
    plotfile = open(filename,'w')
    plotfile.write(header)
        
    for i_line in range(n_lines):
        for j,key in enumerate(keys_plot):
            value = data_plot[key]
            if j > 0: plotfile.write(", ")
            plotfile.write('%s' % value[i_line])
        plotfile.write('\n')
    
    plotfile.close()
    
    return
    
def tecplot(filename,data_plot,keys_plot=[]):
    write_plot(filename,'TECPLOT',data_plot,keys_plot)

def paraview(filename,data_plot,keys_plot=[]):
    write_plot(filename,'PARAVIEW',data_plot,keys_plot)
        
