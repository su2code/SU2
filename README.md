-----------------------------------------------------------
  SU2 (ver. 6.2.0 "Falcon"): The Open-Source CFD Code
-----------------------------------------------------------

Computational analysis tools have revolutionized the way we design engineering systems, but most established codes are proprietary, unavailable, or prohibitively expensive for many users. The SU2 team is changing this, making multiphysics analysis and design optimization freely available as open-source software and involving everyone in its creation and development. 

For an overview of the technical details in SU2, please see the following AIAA Journal article:

"SU2: An open-source suite for multiphysics simulation and design," AIAA Journal, 54(3):828-846, 2016. http://arc.aiaa.org/doi/10.2514/1.J053813

Please note that this project is released with a [Contributor Code of Conduct](CODE_OF_CONDUCT.md). By participating in this project you agree to abide by its terms.

Continuous Integration:<br/>
[![Regression Testing](https://github.com/su2code/SU2/workflows/Regression%20Testing/badge.svg?branch=develop)](https://github.com/su2code/SU2/actions)


----------------------------------------------------------
  SU2 INTRODUCTION 
----------------------------------------------------------

SU2 is a suite of open-source software tools written in C++ for the numerical solution of partial differential equations (PDE) and performing PDE constrained optimization. 

The primary applications are computational fluid dynamics and aerodynamic shape optimization, but has been extended to treat more general equations such as electrodynamics and chemically reacting flows. 

You will find more information and the latest news in:
   - SU2 Home Page: https://su2code.github.io
   - GitHub repository: https://github.com/su2code
   - CFD Online: http://www.cfd-online.com/Forums/su2/
   - Twitter: https://twitter.com/su2code
   - Facebook: https://www.facebook.com/su2code

---------------------------------------------------
  SU2 INSTALLATION
---------------------------------------------------

To build SU2 from the source code, first open a terminal and execute the ./bootstrap script provided in the root directory of the source distribution in order to set the makefiles for your local system. Then, simply run the './configure', 'make', and 'make install' commands. You can provide an install location using the prefix option to configure. Please note that more detailed instructions on the configure and build processes can be found within the INSTALL file.

----------------------------------------------------------
  SU2 PATH SETUP 
----------------------------------------------------------

SU2 is built using a typical configure/make/make install process. When make install is complete, please be sure to add the $SU2_HOME and $SU2_RUN environment variables, and update your $PATH with $SU2_RUN. 

For example, add these lines to your .bashrc file:

- export SU2_RUN="your_prefix/bin"
- export SU2_HOME="/path/to/SU2vX.X.X/"
- export PATH=$PATH:$SU2_RUN
- export PYTHONPATH=$SU2_RUN:$PYTHONPATH

$SU2_RUN should point to the folder where all binaries and python scripts were installed. This is the prefix you set with the --prefix option to configure. Note that the bin/ directory is automatically added to your prefix path.

$SU2_HOME should point to the root directory of the source code distribution, i.e., /path/to/SU2vX.X.X/.

Thanks for building, and happy optimizing!

- The SU2 Development Team

----------------------------------------------------------
  SU2 DEVELOPERS
----------------------------------------------------------

SU2 is being developed by individuals and organized teams all around the world. 

The current SU2 release has been coordinated by the SU2 International Developers Society with selected contributions from the open-source community.

The main research teams contributing to the current release are:
- Prof. Juan J. Alonso's group at Stanford University.
- Prof. Piero Colonna's group at Delft University of Technology.
- Prof. Nicolas R. Gauger's group at Kaiserslautern U. of Technology.
- Prof. Alberto Guardone's group at Polytechnic University of Milan.
- Prof. Rafael Palacios' group at Imperial College London.
- Prof. Vincent Terrapon's group at the University of Liege.
- Prof. Edwin van der Weide's group at the University of Twente.
- Lab. of New Concepts in Aeronautics at Tech. Inst. of Aeronautics.

Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon, Tim Albring, and the SU2 contributors.
