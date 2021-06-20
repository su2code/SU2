# FADO script: Finite Differences of unsteady CHT and adjoint run

from FADO import *

# Design variables ----------------------------------------------------- #

nDV = 18
ffd = InputVariable(0.0,PreStringHandler("DV_VALUE="),nDV)

# Parameters ----------------------------------------------------------- #

# The master config `chtMaster.cfg` serves as an SU2 adjoint regression test.
# For a correct gradient validation we need to exchange some options

time_iter_primal = Parameter(["TIME_ITER= 56"],\
                LabelReplacer("TIME_ITER= 54"))
outer_iter_primal = Parameter(["OUTER_ITER= 200"],\
                 LabelReplacer("OUTER_ITER= 100"))
restart_sol_primal = Parameter(["RESTART_SOL= YES"],\
                  LabelReplacer("RESTART_SOL= NO"))

outer_iter_adjoint = Parameter(["OUTER_ITER= 500"],\
                  LabelReplacer("OUTER_ITER= 100"))

# Evaluations ---------------------------------------------------------- #

# Note that correct SU2 version needs to be in PATH

def_command = "SU2_DEF chtMaster.cfg"
cfd_command = "mpirun -n 12 SU2_CFD chtMaster.cfg"

cfd_ad_command = "mpirun -n 12 SU2_CFD_AD chtMaster.cfg"
dot_ad_command = "mpirun -n 12 SU2_DOT_AD chtMaster.cfg"

max_tries = 1

# mesh deformation
deform = ExternalRun("DEFORM",def_command,True) # True means sym links are used for addData
deform.setMaxTries(max_tries)
deform.addConfig("chtMaster.cfg")
deform.addData("fluid.cfg") # zonal cfg's can be symlinked as they are unchanged
deform.addData("solid.cfg")
deform.addData("MeshCHT.su2")
deform.addExpected("mesh_out.su2")

# direct run
direct = ExternalRun("DIRECT",cfd_command,True)
direct.setMaxTries(max_tries)
direct.addConfig("chtMaster.cfg")
direct.addData("fluid.cfg")
direct.addData("solid.cfg")
direct.addData("DEFORM/mesh_out.su2",destination="MeshCHT.su2")
direct.addData("solution_0_00000.dat")
direct.addData("solution_0_00001.dat")
direct.addData("solution_1_00000.dat")
direct.addData("solution_1_00001.dat")
direct.addExpected("solution_0_00055.dat")
direct.addExpected("solution_1_00055.dat")
direct.addParameter(time_iter_primal)
direct.addParameter(outer_iter_primal)
direct.addParameter(restart_sol_primal)

# adjoint run
adjoint = ExternalRun("ADJOINT",cfd_ad_command,True)
adjoint.setMaxTries(max_tries)
adjoint.addConfig("chtMaster.cfg")
adjoint.addData("fluid.cfg") # zonal cfg's can be symlinked
adjoint.addData("solid.cfg")
adjoint.addData("DEFORM/mesh_out.su2", destination="MeshCHT.su2")
# add all primal solution files
for timeIter in range(56): #
    if timeIter < 10:
        timeIter = "0" + str(timeIter)
    adjoint.addData("DIRECT/solution_0_000" + str(timeIter) + ".dat")
    adjoint.addData("DIRECT/solution_1_000" + str(timeIter) + ".dat")
#end
# replace OUTER_ITER= by 500 and TIME_ITER= 56
adjoint.addParameter(outer_iter_adjoint)
adjoint.addExpected("solution_adj_avtp_0_00053.dat")
adjoint.addExpected("solution_adj_avtp_1_00053.dat")

# gradient projection
dot = ExternalRun("DOT",dot_ad_command,True)
dot.setMaxTries(max_tries)
dot.addConfig("chtMaster.cfg")
dot.addData("fluid.cfg") # zonal cfg's can be symlinked
dot.addData("solid.cfg")
dot.addData("DEFORM/mesh_out.su2", destination="MeshCHT.su2")
# add all adjoint solution files
for timeIter in range(54):
    if timeIter < 10:
        timeIter = "0" + str(timeIter)
    dot.addData("ADJOINT/solution_adj_avtp_0_000" + str(timeIter) + ".dat")
    dot.addData("ADJOINT/solution_adj_avtp_1_000" + str(timeIter) + ".dat")
#end
dot.addExpected("of_grad.csv")

# Functions ------------------------------------------------------------ #

tavgT = Function("tavgT", "DIRECT/chtMaster.csv",LabeledTableReader("\"tavg[AvgTemp[1]]\""))
tavgT.addInputVariable(ffd,"DOT/of_grad.csv",TableReader(None,0,(1,0))) # all rows, col 0, don't read the header
tavgT.addValueEvalStep(deform)
tavgT.addValueEvalStep(direct)
tavgT.addGradientEvalStep(adjoint)
tavgT.addGradientEvalStep(dot)

# Driver --------------------------------------------------------------- #

# The input variable is the constraint tolerance which is not used for our purpose of finite differences
driver = ExteriorPenaltyDriver(0.005)
driver.addObjective("min", tavgT)

driver.setWorkingDirectory("DOE")
driver.preprocessVariables()
driver.setStorageMode(True,"DESIGN_")

his = open("doe.his","w",1)
driver.setHistorian(his)

# Simulation Runs ------------------------------------------------------ #

# Primal simulation for each deformed DV
for iLoop in range(0, nDV, 1):
    print("Computing deformed primal ", iLoop, "/", nDV-1)
    x = driver.getInitial()
    x[iLoop] = 1e-4 # DV_VALUE
    driver.fun(x)
#end

# Undeformed/initial primal last in order to have the correct solution in
# the WorkindDirectory for the following adjoint
print("Computing baseline primal")
x = driver.getInitial()
driver.fun(x) # baseline evaluation

# Compute discrete adjoint gradient
print("Computing discrete adjoint gradient")
driver.grad(x)

his.close()

# For results run `python postprocess.py` to get screen output
# of the differences between primal and adjoint simulation.
