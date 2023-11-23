/*!
 * \file docmain.hpp
 * \brief This file contains documentation for Doxygen and does not have any significance with respect to C++.
 * \author F. Palacios
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

/*!
 * \mainpage SU2 version 8.0.0 "Harrier"
 * SU2 suite is an open-source collection of C++ based software tools
 * to perform PDE analysis and PDE constrained optimization.  The toolset is designed with
 * computational fluid dynamics and aerodynamic shape optimization in mind, but is extensible to
 * include other families of governing equations such as potential flow, electrodynamics, chemically reacting
 * flows, and many others.  SU2 is released under an open-source license.
 *
 * The following technical documentation describes the structure and details of the source code for developers.
 * Each module groups classes and functions dedicated to major areas or features of the code base.
 * This documentation helps awnsering questions such as
 * - What are the major classes of the code (architecture)?
 * - How do those classes interact with each other?
 * - What low-level functionality is already implemented that can be re-used?
 *
 * Given the use of polymorphism in SU2, it is usefull to start studying each module from the base class
 * (which usually has a generic name), this will show the inheritance diagram for that family. Moreover,
 * many important functions (defining the interactions between families) are implemented at the base level,
 * the "call graph" of a function will show these interactions (derived classes only specialize certain details).
 * Alternatively, "caller graphs" can be used to navigate to the larger classes that use smaller ones, and
 * thus navigate the architecture from the bottom up (note that this will also be more effective if done
 * for base classes).
 *
 * ### Good practices when documenting classes and functions
 *
 * The groups and subgroups defined in this file form a hierarchy of doxygen modules, when documenting a
 * class or free function it should be added to a group (member functions belong to the group of the class
 * by default), otherwise it will be difficult to navigate to the documentation of the class or function.
 * Therefore, all classes and functions should have a doxygen documentation block with at least \\brief and
 * \\ingroup properties.
 *
 * Note that this is a work in progress and not all classes and functions are currently inserted into groups.
 */

/*!
 * \defgroup Config Description of the Configuration Options
 * \brief Group of variables that can be set using the configuration file.
 */

/*!
 * \defgroup ConvDiscr Discretization of the convective terms
 * \brief Group of classes which define the numerical methods for
 *		  discretizing the convective terms of a Partial Differential Equation.
 *		  There are methods for solving the direct, adjoint and linearized
 *		  systems of equations.
 */

/*!
 * \defgroup ViscDiscr Discretization of the viscous terms
 * \brief Group of classes which define the numerical methods for
 *		  discretizing the viscous terms of a Partial Differential Equation.
 *		  There are methods for solving the direct, adjoint and linearized
 *		  systems of equations.
 */

/*!
 * \defgroup SourceDiscr Discretization of the source terms
 * \brief Group of classes which define the numerical methods for
 *		  discretizing the source terms of a Partial Differential Equation.
 *		  There are methods for solving the direct, adjoint and linearized
 *		  systems of equations.
 */

/*!
 * \defgroup Euler_Equations Solving the Euler equations
 * \brief Group of classes which define the system of Euler equations in
 *		  three formulations: direct, adjoint, and linearized.
 */

/*!
 * \defgroup Navier_Stokes_Equations Solving the Navier-Stokes equations
 * \brief Group of classes which define the system of Navier-Stokes equations in
 *		  three formulations: direct, adjoint, and linearized.
 */

/*!
 * \defgroup Turbulence_Model Solving the turbulence model equations
 * \brief Group of classes which define the turbulence model in
 *		  three formulations: direct, adjoint, and linearized.
 */

/*!
 * \defgroup Scalar_Transport Solving scalar transport equations
 * \brief Classes to solve scalar transport equations.
 */

/*!
 * \defgroup Elasticity_Equations Solving the elasticity equations
 * \brief Group of classes to solve solid deformation problems.
 */

/*!
 * \defgroup Interfaces Multiphysics interfaces
 * \brief Classes for data transfer and interpolation across interfaces between zones.
 */

/*!
 * \defgroup Drivers Iterative solution strategy
 * \brief Group of classes which define the iterative process used to converge
 *        the equations (fluid, turbulence, elasticity, FSI, CHT, etc.).
 *        In general, "Driver" classes use "Iteration" classes to advance one inner
 *        iteration, in turn "Iteration" classes use "Integration" classes to perform
 *        space and time integration. The latter use mostly the "Solvers".
 */

/*!
 * \defgroup PySU2 Python Wrapper functions
 * \brief Functions of the driver classes (\ref Drivers) that can be used to customize SU2
 *        via Python. For example, set custom temperature distributions at boundaries,
 *        deform the mesh, etc.
 */

/*!
 * \defgroup Variable Storing solution variables
 * \brief Classes used to store and access the solution variables of all types of problems.
 */

/*!
 * \defgroup Output Screen, history, and file output
 * \brief Classes used to define the output of the solvers in SU2.
 */

/*!
 * \defgroup DiscAdj Discrete Adjoint
 * \brief Classes and functions used to solve discrete adjoint equations.
 */

/*!
 * \defgroup GradSmooth Gradient Smoothing
 * \brief Classes and functions used to smooth gradients from the discrete adjoint method.
 * \ingroup DiscAdj
 */

/*!
 * \defgroup SpLinSys Sparse linear systems
 * \brief Classes and function to represent and solve large distributed sparse linear systems.
 */

/*!
 * \defgroup FvmAlgos General FVM algorithms
 * \brief Common algorithms used in FVM implementations.
 */

/*!
 * \defgroup FemAlgos General FEM algorithms
 * \brief Common algorithms used in FEM implementations.
 */

/*!
 * \defgroup Toolboxes Utility classes and functions
 * \brief Several classes and functions that implement common operations.
 */

/*!
 * \defgroup GeometryToolbox Geometry toolbox
 * \brief Common geometry operations.
 * \ingroup Toolboxes
 */

/*!
 * \defgroup Containers Data containers
 * \brief Container classes (vectors, matrices, ND-arrays, etc.).
 * \ingroup Toolboxes
 */

/*!
 * \defgroup LookUpInterp Look up and interpolation
 * \brief Data look up and interpolation.
 * \ingroup Toolboxes
 */

/*!
 * \defgroup ADT Alternating Digital Tree
 * \brief Tree-based searches (minimum distance, containment, etc.).
 * \ingroup Toolboxes
 */

/*!
 * \defgroup BLAS Dense linear algebra
 * \brief Linear algebra functions and classes.
 * \ingroup Toolboxes
 */

/*!
 * \defgroup VecExpr Vector math expression templates
 * \brief Expression templates for level-1 BLAS operations.
 * \ingroup BLAS
 */

/*!
 * \defgroup Graph Graph operations
 * \brief Classes to represent graphs and functions to manipulate them (coloring, etc.).
 * \ingroup Toolboxes
 */

/*!
 * \defgroup SIMD Vectorization (SIMD)
 * \brief Classes for explicit (done by the programmer) vectorization (SIMD) of computations.
 * \ingroup Toolboxes
 */

/*!
 * \defgroup Multi-Layer Perceptrons (MLP)
 * \brief Data look up and interpolation via dense, feed-forward multi-layer perceptrons.
 * \ingroup Toolboxes
 */
