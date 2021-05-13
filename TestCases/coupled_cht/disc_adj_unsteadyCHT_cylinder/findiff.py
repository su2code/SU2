# FADO script: Finite Differences of unsteady CHT

from FADO import *

# Design variables ----------------------------------------------------- #

nDV = 18
ffd = InputVariable(0.0,PreStringHandler("DV_VALUE="),nDV)

# Parameters ----------------------------------------------------------- #

# switch input mesh to perform deformation
mesh_in = Parameter(["MESH_FILENAME= mesh_out.su2"],\
       LabelReplacer("MESH_FILENAME= MeshCHT.su2"))

# Evaluations ---------------------------------------------------------- #

def_command = "SU2_DEF chtMaster.cfg" # Note that correct SU2 version needs to be in PATH
cfd_command = "mpirun -n 12 SU2_CFD chtMaster.cfg"

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
direct.addData("DEFORM/mesh_out.su2")
direct.addData("solution_0_00000.dat")
direct.addData("solution_0_00001.dat")
direct.addData("solution_1_00000.dat")
direct.addData("solution_1_00001.dat")
direct.addExpected("solution_0_00055.dat")
direct.addExpected("solution_1_00055.dat")
direct.addParameter(mesh_in)

# Functions ------------------------------------------------------------ #

tavgT = Function("tavgT", "DIRECT/chtMaster.csv",LabeledTableReader("\"tavg[AvgTemp[1]]\""))
tavgT.addInputVariable(ffd,"",None)
tavgT.addValueEvalStep(deform)
tavgT.addValueEvalStep(direct)

# Driver --------------------------------------------------------------- #

driver = ExteriorPenaltyDriver(0.005)
driver.addObjective("min", tavgT)

driver.setWorkingDirectory("DOE")
driver.preprocessVariables()
driver.setStorageMode(True,"DESIGN_")

his = open("doe.his","w",1)
driver.setHistorian(his)

# Finite Differences
x = driver.getInitial()
driver.fun(x) # baseline evaluation

for iLoop in range(0, nDV, 1):
    x = driver.getInitial()
    x[iLoop] = 1e-4 # DV_VALUE
    driver.fun(x)
#end

his.close()
