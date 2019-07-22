Subfolders below this folder contain example applications that use the tecio library.

Building the example applications:

Linux:
------
  cd to the example folder you would like to build.
  Type:
        make <cr>

  Notes:  
    1) Edit the Makefile for the example to change the targets you want to build.
       Not all examples are set up to build using the mpi tecio library.
    2) If you are building against the tecio (or teciompi) library built using the source code
       care package then see instructions in base.make located in this folder.

Windows:
--------
  Developer Studio solution files are provided to help you build the example 
  applications.
