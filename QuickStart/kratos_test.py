import SU2


config = SU2.io.Config("inv_NACA0012.cfg")

state = SU2.io.State()
state.find_files(config)
config.NUMBER_PART = 2
config.NZONES= 1
config.GRADIENT_METHOD= 'DISCRETE_ADJOINT'
#config.CONSOLE= 'CONCISE'
project = SU2.opt.Project(config,state)

project.config["MOTION_FILENAME"] = "mesh_motion.dat"
project.state.FILES["MOTION_FILE"] = project.config["MOTION_FILENAME"]

obj = []
dobj = []
#obj.append (project.f( ["DRAG"], True, 1 ))
obj.append  (project.f (["DRAG"], True , 1 ))
dobj.append (project.df(["DRAG"], False, 1 ))
dobj.append (project.df(["DRAG"], True , 2 ))
obj.append  (project.f (["DRAG"], False, 1 ))



#obj_d = { id : [x y z] }
#obj_d = project.df(["LIFT"], False, 2)

print (obj)

print (dobj)
