
import os, sys, shutil, copy
sys.path.append(os.environ['SU2_RUN'])
import SU2

from collections import OrderedDict

# todo:
# verify command line interface
# commenting
# verify optimization, gradients, flow solutions
# verbosity
# plotting
# verbose redirection
# pyopt optimizers
# verify server.py


# needed config options
# DECOMPOSED = YES/NO
# OPT_CONSTRAINT
# OPT_OBJECTIVE
# CONSOLE

# OUTPUT_WEIGHT
# FINDIFF_STEP
# GPC_FINDIFF_STEP
# GRADIENT_METHOD= FINITE_DIFFERENCING, CONTINUOUS_ADJOINT, DISCRETE_ADJOINT
# ADAPTATION= DIRECT, ADJOINT


def main():
    
    #io0()       # working
    #io1()
    level0()    # working
    level1()    # working
    level2()    # working
    level3()    # working
    level4()    # working
    level5()    # working
    
    print 'DONE!'
    
def io0():
    folder='test_io0'; pull='config_NACA0012.cfg'; link='mesh_NACA0012.su2'
    with SU2.io.redirect_folder(folder,pull,link):
        
        config_name = 'config_NACA0012.cfg'
        config = SU2.io.Config(filename=config_name)
        
        print config
        
        config.ADAPT_CYCLES
        config['ADAPT_CYCLES']
        
        config.dump('out.cfg')
        
        konfig = copy.deepcopy(config)
        konfig['TASKS'] = ['TEST']
        konfig['NUMBER_PART'] = 0
        
        config_diff = config.diff(konfig)
        
        print config_diff
        
    
    wait = 0
    
def io1():
    
    option = SU2.io.config.MathProblem()
    
    option = 'DIRECT'
    
    wait = 0

def level0():
    folder='test_level0'; pull='config_NACA0012.cfg'; link='mesh_NACA0012.su2'
    with SU2.io.redirect_folder(folder,pull,link):
        
        # Setup
        config_name = 'config_NACA0012.cfg'
        config = SU2.io.Config(config_name)
        config.EXT_ITER = 9
        config.NUMBER_PART = 2
        
        SU2.run.DDC(config)
        
        SU2.run.CFD(config)
    
def level1():
    folder='test_level1'; pull='config_NACA0012.cfg'; link='mesh_NACA0012.su2'
    with SU2.io.redirect_folder(folder,pull,link):
        
        # Setup
        config_name = 'config_NACA0012.cfg'
        config = SU2.io.Config(config_name)
        config['NUMBER_PART'] = 2
        config['EXT_ITER'] = 9
        
        state = SU2.io.State()
        
        # Deformation
        dv_new = [0.002]*38
        info = SU2.run.deform(config,dv_new)
        state.update(info)
        
        # Direct Solution
        info = SU2.run.direct(config) 
        state.update(info)
        SU2.io.restart2solution(config,state)
        
        # Adjoint Solution
        info = SU2.run.adjoint(config)
        state.update(info)
        SU2.io.restart2solution(config,state)
        
        # Gradient Projection
        info = SU2.run.projection(config)
        state.update(info)
        
        print state
        
        SU2.io.save_data('state.pkl',state)
        data = SU2.io.load_data('state.pkl')
        
        SU2.io.save_data('config.pkl',config)
        data = SU2.io.load_data('config.pkl')
        
    wait = 0
    
def level2():
    folder='test_level2'; pull='config_NACA0012.cfg'; link='mesh_NACA0012.su2'
    with SU2.io.redirect_folder(folder,pull,link):
        
        # Setup
        config_name = 'config_NACA0012.cfg'
        config = SU2.io.Config(config_name)
        config['NUMBER_PART'] = 2
        config['EXT_ITER'] = 9
        dv_new = [0.0]*38
        #dv_new[10] = 0.05
        config.unpack_dvs(dv_new)
        
        state = SU2.io.State()
        
        #with SU2.io.redirect.folder(folder='JOB_001',link='mesh_NACA0012.su2'):
        #    grad = SU2.eval.grad( 'DRAG', 'FINDIFF', config, state )
        
        with SU2.io.redirect_folder(folder='JOB_001',link='mesh_NACA0012.su2'):
            func  = SU2.eval.func( 'LIFT', config, state )
            grads = SU2.eval.grad( 'LIFT', 'ADJOINT', config, state )
        
        with SU2.io.redirect_folder(folder='JOB_001',link='mesh_NACA0012.su2'):
            func  = SU2.eval.func( 'DRAG', config, state ) # will not run direct
            grads = SU2.eval.grad( 'LIFT', 'ADJOINT', config, state ) # will not run adjoint
            grads = SU2.eval.grad( 'DRAG', 'ADJOINT', config, state ) # will run adjoint
    
    wait = 0
    
def level3():
    folder='test_level3'; pull='config_NACA0012.cfg'; link='mesh_NACA0012.su2'
    with SU2.io.redirect_folder(folder,pull,link):
        
        # Setup
        config_name = 'config_NACA0012.cfg'
        config = SU2.io.Config(config_name)
        config['NUMBER_PART'] = 2  
        config['EXT_ITER'] = 9
        
        # initialize design state
        state = SU2.io.State()
        state.find_files(config)
        
        # start design
        design = SU2.eval.Design(config,state)
        
        # run design with dv change
        dv_new = [0.0]*38
        vals = design.obj_f(dv_new)
        vals = design.obj_df(dv_new)
        vals = design.con_ceq(dv_new)
        vals = design.con_dceq(dv_new)
        vals = design.con_cieq(dv_new)
        vals = design.con_dcieq(dv_new)
        vals = design.func('LIFT')
        vals = design.grad('LIFT','ADJOINT')
        
        SU2.io.save_data('design.pkl',design)
        data = SU2.io.load_data('design.pkl')
        
    wait = 0
    
def level4():
    folder='test_level4'; pull='config_NACA0012.cfg'; link='mesh_NACA0012.su2'
    with SU2.io.redirect_folder(folder,pull,link):    
    
        # Setup
        config_name = 'config_NACA0012.cfg'
        config = SU2.io.Config(config_name)
        config['NUMBER_PART'] = 2  
        config['EXT_ITER'] = 9
        config.CONSOLE = 'QUIET'
        
        # initialize design state
        state = SU2.io.State()
        state.find_files(config)
        
        # initialize project
        project = SU2.opt.Project(config,state)
        
        # run project with dv changes
        dv_new = [0.0]*38
        vals = project.obj_f(dv_new)
        vals = project.obj_df(dv_new)

        dv_new = [-0.005]*38
        vals = project.obj_f(dv_new)
        
        dv_new = [0.0]*38
        dv_new[9] = -0.02
        vals = project.obj_f(dv_new)
        
        dv_new = [0.005]*38
        vals = project.obj_f(dv_new) # will not rerun solutions
        
        SU2.io.save_data('project.pkl',project)
        data = SU2.io.load_data('project.pkl')
        
        data = project.data
        
    wait = 0
    print "Done!"
    
def level5():
    folder='test_level5'; pull='config_NACA0012.cfg'; link='mesh_NACA0012.su2'
    with SU2.io.redirect_folder(folder,pull,link):    
        
        # Setup
        config_name = 'config_NACA0012.cfg'
        config = SU2.io.Config(config_name)
        config['NUMBER_PART'] = 2  
        config['EXT_ITER'] = 9
        config['CONSOLE'] = 'CONCISE'
        
        # set optimization problem
        obj = {}
        obj['DRAG'] = {'SCALE':1.e-2}
        
        cons = {}
        cons['EQUALITY'] = {}
        cons['INEQUALITY'] = {}
        cons['INEQUALITY']['LIFT']     = {'SIGN':'>','VALUE':0.328188,'SCALE':1e-1}
        cons['INEQUALITY']['MOMENT_Z'] = {'SIGN':'>','VALUE':0.034068,'SCALE':1e-2}
        
        def_dv = config.DEFINITION_DV
        n_dv   = len(def_dv['KIND'])
        def_dv['SCALE'] = [1.e0]*n_dv
        
        config.OPT_OBJECTIVE  = obj
        config.OPT_CONSTRAINT = cons
                
        # initialize design state
        state = SU2.io.State()
        state.find_files(config)
        
        # initialize project
        project = SU2.opt.Project(config,state)
        
        # optimization setup
        x0 = [0.0]*n_dv
        xb = [] #[[-0.02,0.02]]*n_dv
        its = 20
        
        # optimize
        SU2.opt.SLSQP(project,x0,xb,its)
    
    wait = 0

if __name__ == '__main__':
    main()