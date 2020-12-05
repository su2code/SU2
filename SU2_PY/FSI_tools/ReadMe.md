###
Brief explanation on the implemented python interface

for more information please contact nicola.fonzi@polimi.it
###

The current interface is largely based on the work done by David Thomas, but is has
been updated to the latest version of SU2. It is now able to work with the new
driver, singleZoneDriver, and it is able to use the new mesh deformation scheme.

Further, some missing parts have been completed.

The code has been extensively validated and in the folder ../SU2_Nastran/TestCase
a subset of the validation cases are present for reproducibility.

Here, the most important points for the usage of the interface are presented.

The fluid configuration file is set as for a coupled simulation internal to SU2.
Thus, the normal configuration for the fluid zone is provided, plus the
coupling keywords that specify the marker used for the deformation, the marker where
to compute the loads, how to deform the mesh, ecc...

A new marker is introduced, MARKER_MATCH_DEFORM_MESH. This marker is effectively
a symmetry marker for the mesh deformation only. It may be useful in cases where
symmetry in the mesh is required, but not in the fluid simulation. An example may
be the simulation of a plane half-model, in wind tunnel, where the effect of boundary
layer on the tunnel walls must be studied, but a pitch movement of the model is
also allowed. Fluid symmetry cannot be used, but at the same time the mesh should
move on the tunnel wall to match the deformation given by the pitch motion.

Then, other two configuration files must be prepared: one for the solid solver
and one for the interface.

The most important interface configuration keywords are:

NDIM (int): 2 or 3 depending if the model is bidimensional or tridimensional
RESTART_ITER (int): Restart iteration
TIME_TRESHOLD (int): Time iteration after which fluid and structure are coupled
                     in an unsteady simulation
NB_FSI_ITER (int):   Number of max internal iterations to couple fluid and structure
RBF_RADIUS (float):  Radius for the RBF interpolation. It is dimensional (i.e. in meters)
                     and must be set so that at least 5 structural points are always
                     inside a sphere with that radius and centered in any of the
                     structural nodes. The more nodes are included, the better
                     the interpolation. However, with larger radius, the interpolation
                     matrix may become close to singular.
AITKEN_PARAM (float): Under relaxation parameter, between 0 and 1
UNST_TIMESTEP (float): Physical time step size in unsteady simulations, must match
                       the one in the other cfg files.
UNST_TIME (float): Physical simulation time for unsteady problems.
FSI_TOLERANCE (float): Tolerance for inner loop convergence between fluid and structure
CFD_CONFIG_FILE_NAME (string): Path of fluid cfg file
CSD_SOLVER (string): Behaviour of the structural solver to be used. AEROELASTIC if
                     the structural equation of motions must be solved, IMPOSED if
                     a movement of the structure is imposed.
CSD_CONFIG_FILE_NAME (string): Path to solid cfg file
RESTART_SOL (string): YES or NO
MATCHING_MESH (string): YES or NO, the fluid and structural mesh match at the interface
MESH_INTERP_METHOD (string): Interpolation method in case of nonmatching meshes. TPS or RBF
DISP_PRED (string): Displacement predictor order FIRST_ORDER or SECOND_ORDER. To
                    be used in unsteady simulations.
AITKEN_RELAX (string): DYNAMIC or STATIC. It can be automatically changed during
                       the simulation.
TIME_MARCHING (string): YES or NO
