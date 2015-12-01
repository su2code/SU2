#pragma once

#include "../../Common/include/config_structure.hpp"
#include "../../SU2_CFD/include/solver_structure.hpp"
#include "../../SU2_CFD/include/iteration_structure.hpp"
#include "../../SU2_CFD/include/driver_structure.hpp"
#include "../../SU2_CFD/include/integration_structure.hpp"
#include "../../SU2_CFD/include/output_structure.hpp"
#include "../../SU2_CFD/include/numerics_structure.hpp"
#include "../../Common/include/geometry_structure.hpp"
#include "../../Common/include/grid_movement_structure.hpp"
#include <string>

#define ROW_MAJ 0
#define COL_MAJ 1

/*!
 * \class SU2Solver
 * \brief Main API class representing the SU2 solver (fluid part) for staggered external coupling
 *  with other structural solver for FSI computation.
 * \author D. Thomas
 * \version BETA
 */
class SU2Solver{

protected:
    std::string confFile;	/*!< \brief Name of the SU2 (fluid) configuration file (.cfg). */
    ofstream ConvHist_file;     /*!< \brief Convergence history file */
    unsigned short nZone,	/*!< \brief Number of zone. */
    nDim;			/*!< \brief Number of dimension. */
    CDriver *driver;
    CIteration **iteration_container;
    CConfig **config_container;		/*!< \brief Configuration container. */
    CConfig *config;			/*!< \brief Main configuration description. */
    CGeometry ***geometry_container;	/*!< \brief Geometry description container. */
    CNumerics *****numerics_container;	/*!< \brief Numerical method description container. */
    CSolver ****solver_container;	/*!< \brief Solution container. */
    CIntegration ***integration_container;	/*!< \brief Integration description container. */
    COutput *output;		/*!< \brief Output description. */
    CSurfaceMovement **surface_movement;	/*!< \brief Surface movement description. */
    CVolumetricMovement **grid_movement;	/*!< \brief Grid movement description. */
    CFreeFormDefBox*** FFDBox;
    unsigned long nLocalFluidInterfaceVertex, nAllFluidInterfaceVertex;
    unsigned long nSolidInterfaceVertex;    /*!< \brief Number of vertices on the solid interface */
    su2double *interfRigidDispArray;		/*!< \brief Vectors containing the rigid displacement of the FSI interface (6 DOF), used for communication with external solver. */
    su2double* solidInterfaceBuffer;       /*!< \brief Buffer used for communication of solid interface position between MPI processes */
    su2double** solidInterface;            /*!< \brief Array containing coordinates of each node of the solid interface. */
    su2double** fluidSurfaceloads;         /*!< \brief Array containing surface loads at each node of the fluid interface after gathering all the MPI processes. */
    su2double** partFluidSurfaceLoads;     /*!< \brief Array containing surface loads at each node of the partitioned fluid interface. */
    su2double Center[3];                   /*!< \brief Current position of the reference center (motion and moment) */
    su2double Center_n[3];             /*!< \brief Position of the reference center (motion and moment) at previous time step */
    su2double* attitude_i1;                //Used for steady computations, may disappear

public:

    /*!
     * \brief Constructor of the class.
     */
    SU2Solver(std::string str);

    /*!
     * \brief Destructor of the class.
     */
    virtual ~SU2Solver();

    /*!
     * \brief Initialize the SU2 solver (fluid part) for FSI computation.
     * \param[in] FSIComp - Make sure to initialize for FSI computation.
     */
    void initialize( bool FSIComp);

    /*!
     * \brief Exit the SU2 solver, e.g. free the memory after dynamic allocation and closed the
     * opened files.
     */
    void exit();

    /*!
     * \brief Create the communication with external solid solver during FSI pre-processing.
     * \param[in] nVertex - Number of vertices on the solid interface, communicated by external solver.
     */
    void connectToSolidSolver(unsigned long nVertex);

    /*!
     * \brief Broadcast all FSI information and create the FSI interface structure
     */
    void setFSIInterface();

    /*!
     * \brief Used to communicate the rigid body displacement vector to external solver for FSI.
     */
    double* getInterRigidDispArray() const;

    /*!
     * \brief Set the current iteration, update some parameters if required.
     * \param[in] ExtIter - Value of the iteration.
     */
    void setTemporalIteration(unsigned long ExtIter);

    /*!
     * \brief Perform one inner loop (fluid sub-iteration) for dual-time stepping integration.
     * \param[in] ExtIter - Value of the time iteration.
     */
    void dualTimeInnerLoop(unsigned long ExtIter);

    /*!
     * \brief Perform one fluid iteration for steady computation
     * \param[in] ExtIter - Value of the iteration.
     */
    void steadyFluidIteration(unsigned long ExtIter);

    /*!
     * \brief Reset the convergence flag before a new FSI iteration (after a mesh update).
     */
    void resetConvergence();

    /*!
     * \brief Update the dual-time solution.
     */
    void updateDualTime();

    /*!
     * \brief Execute some routines for screen and file output (e.g. write a solution file if required).
     * \param[in] ExtIter - Value of the iteration.
     */
    bool writeSolution(unsigned long ExtIter);

    /*!
     * \brief Apply a rigid body displacement on a surface after communication with the external structural
     * solver (FSI computation).
     */
    void mapRigidDisplacementOnFluidMesh_Old();

    /*!
     * \brief Get the position of the solid interface from the external solid solver.
     * \param[in] solidCoordinate - The coordinates of each node of the solid interface, communicated by external solid solver.
     * \param[in] CenterCoordinate - The coordinate of the center of reference (elastic axis, reference for moment computation, ...), communicated by external solid solver.
     */
    void getSolidDisplacement(double** solidCoordinate, const double* CenterCoordinate);

    /*!
     * \brief Interpolate and map the displacement of the solid interface (solid mesh) on the fluid interface (fluid mesh)
     */
    void mapRigidDisplacementOnFluidMesh();

    /*!
     * \brief Make the one-to-one correspondance between nodes on the solid and fluid interfaces in case of matching meshes.
     */
    unsigned short searchEquivIndex(double** target, unsigned long size_target, unsigned long iPoint);

    /*Search the number of occurence of one parameter in a table (only used for development purpose, may disappear).*/
    int searchNumberOfOccurence(double** target, unsigned long size_target, unsigned long iPoint);

    /*!
     * \brief Dynamically update the mesh (deformation - grid velocity - MG update) for a new unsteady FSI iteration.
     * \param[in] ExtIter - Value of the iteration.
     */
    void dynamicMeshUpdate(unsigned long ExtIter);

    /*!
     * \brief Update the mesh (deformation - MG update) for a new steady FSI iteration or for initial mesh deformation.
     */
    void staticMeshUpdate();

    /*!
     * \brief Set the initial mesh (purely rigid body method) according to the initial condition communicated by the external solid solver.
     * \param[in] restart - Tell weither the solution is restarted or not.
     */
    void setInitialMesh_Old(bool restart);

    /*!
     * \brief Set the initial mesh (general method) according to the initial condition communicated by the external solid solver.
     * \param[in] restart - Tell weither the solution is restarted or not.
     */
    void setInitialMesh(bool restart);

    /*!
     * \brief Ouput a global aerodynamic load - drag.
     * \param[out] Drag - the global drag acting on the monitored surface
     */
    double outputFluidLoads_Drag();

    /*!
     * \brief Ouput a global aerodynamic load - lift.
     * \param[out] Lift - the global drag acting on the monitored surface
     */
    double outputFluidLoads_Lift();

    /*!
     * \brief Ouput global aerodynamic loads and send the data to the external solid solver (used only for purely rigid body motion method, may disappear).
     * \param[in] Array to store the global fluid loads in the solid solver, communicated by the external solid solver.
     */
    void outputGlobalFluidLoads(double* globalFluidLoad);

    /*!
     * \brief Get the surface loads on each node of the fluid interface and gather them to the master process before communication with external solver.
     */
    void mergeSurfaceLoads();

    /*!
     * \brief Ouput the surface loads for each node of the fluid interface and send the data to the external solid solver.
     * \param[in] Array to store the surface fluid loads in the solid solver, communicated by the external solid solver.
     */
    void outputSurfaceLoads(double** solidSurfaceLoads);

    //double computeRigidDisplacementNorm() const;

    /*!
     * \brief Broadcast the rigid displacement of the surface from the root to other processes (used only for purely rigid body motion method, may disappear).
     */
    void broadcastInterfDisp();
};

void MatrixToVec(int order, double** matrix, double* vecteur, int Nrow, int Ncol, int sizeVec);
void VecToMatrix(int order, double** matrix, double* vecteur, int Nrow, int Ncol, int sizeVec);
