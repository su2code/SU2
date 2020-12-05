###
Brief explanation of the implemented coupling between Nastran and SU2

Please contact nicola.fonzi@polimi.it for more information
###

This structural solver loosely couples the commercial FEM code Nastran with SU2.

First, the FEM model must be prepared, and a modal simulation performed.
In the case control section of the Nastran model (i.e. at the very beginning of
the file), the following lines must be added:

ECHO = SORT
DISPLACEMENT(PRINT,PUNCH)=ALL

This will produce, in the f06 file, an equivalent, ordered, model that will
eventually be read by the python script to create the interface. Further, it will
be created a punch file where all the mode shapes, together with their stiffness,
are stored.

IMPORTANT: The modes should be normalised to unit mass.

Further, in the Nastran model, a SET1 card must be added that contains all the
nodes to be interfaced with aerodynamics. Note that, as of now, only one SET1 card
is allowed. However, this should be sufficiently general not to create issues.

The input and output reference systems of the interface nodes can be defined as
local, but these reference systems must be directly defined with respect to the
global one. Thus, if node x is referring to reference system y, y must be defined
with respect to the global one.

In the structural input file the keyword NMODES must then be defined to select which,
of all the modes in the punch file, to be used.

The solver can work in two ways:

1) It can impose the movement of a mode, with prescribed law, to provide forced
response analysis

2) It can integrate in time the modal equation of motions to study the linearised
structural deformations when the body is surrounded by the flow

Available keyword for the config file:

NMODES (int): number of modes to use in the analysis -> if n modes are available in
             the punch file, but only the first m<n are required, set this to m

IMPOSED_MODE (int): mode with an imposed motion

RESTART_ITER (int): if restart is used, this specifies the iteration to restart

DELTA_T (float): physical time step size to be used in the simulation. Must match
                the one in SU2

MODAL_DAMPING (float): the code is able to add a damping matrix to the system, based
                      on a critical damping. This keyword specifies the amount of damping
                      that can be included: if x% of damping is required, set it
                      to 0.0x

RHO (float): rho parameter for the integrator

TIME_MARCHING (string): YES or NO

MESH_FILE (string): path to the f06 file

PUNCH_FILE (string): path to the pch file

RESTART_SOL (string): YES or NO

IMPOSED_DISP (string): string containing the function for the displacement. Example
                       is sine(2*pi*time)+10

IMPOSED_VEL (string): analytical differentiation of above

IMPOSED_ACC (string): analytical differentiation of above

MOVING_MARKER (string): name for the interface marker

INITIAL_MODES (list): list containing the initial amplitudes of the modes. Example
                      is {0:0.1,1:0.0,3:5.0,...}
