***********************************************
**                   README                  **
***********************************************

To build the TecIO library and/or the pltview utility, 
simply run the Runmake script in this directory.

If customization is needed, it will most likely be done
in GLOBAL.h (to identify machine as 64 bit) and/or in
dataio4.c.  Just look for CRAY in dataio4.c and you
will find most of the critical areas.  Note that the
existing code defined by CRAY is quite old and has
not been in use for some time.

Each example has its own Makefile. You may have to adjust
the variables at the top of the Makefile for your platform.


ReadTec()

The ReadTec() function is included in the tecio library
but is not supported by Tecplot, Inc.  ReadTec is used 
to read Tecplot binary data files (all versions at or 
older than the Tecplot version providing the tecio 
library). See tecsrc/DATAUTIL.h for more information.

The pltview example app gives an example of using ReadTec
to read just the header from a file as well as loading all
field data from a file.

