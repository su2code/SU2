import SU2


config = SU2.io.Config("inv_NACA0012.cfg")

state = SU2.io.State()
state.find_files(config)
config.NUMBER_PART = 2
config.NZONES= 1
project = SU2.opt.Project(config,state)

project.config["MOTION_FILENAME"] = "mesh_motion.dat"
project.state.FILES["MOTION_FILE"] = project.config["MOTION_FILENAME"]
obj = project.f( ["DRAG"], False, 1 )

print(obj)
obj = project.f(["LIFT"], True, 2 )


obj = project.f(["LIFT"], False, 2 )

#obj_d = { id : [x y z] }
#obj_d = project.df(["LIFT"], False, 2)

print (obj)
