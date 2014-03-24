/*
 * NOTICE and LICENSE for Tecplot Input/Output Library (TecIO) - OpenFOAM
 *
 * Copyright (C) 1988-2009 Tecplot, Inc.  All rights reserved worldwide.
 *
 * Tecplot hereby grants OpenCFD limited authority to distribute without
 * alteration the source code to the Tecplot Input/Output library, known 
 * as TecIO, as part of its distribution of OpenFOAM and the 
 * OpenFOAM_to_Tecplot converter.  Users of this converter are also hereby
 * granted access to the TecIO source code, and may redistribute it for the
 * purpose of maintaining the converter.  However, no authority is granted
 * to alter the TecIO source code in any form or manner.
 *
 * This limited grant of distribution does not supersede Tecplot, Inc.'s 
 * copyright in TecIO.  Contact Tecplot, Inc. for further information.
 * 
 * Tecplot, Inc.
 * 3535 Factoria Blvd, Ste. 550
 * Bellevue, WA 98006, USA
 * Phone: +1 425 653 1200
 * http://www.tecplot.com/
 *
 */
#if defined EXTERN
#undef EXTERN
#endif
#if defined DATASHRMODULE
#define EXTERN
#else
#define EXTERN extern
#endif

/*
*****************************************************************
*****************************************************************
*******                                                  ********
****** Copyright (C) 1988-2008 Tecplot, Inc.              *******
*******                                                  ********
*****************************************************************
*****************************************************************
*/

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif /* TECPLOTKERNEL */

/*
 * General set of macros for reference count mananagement.
 */
#define IncStructureReference(V)  ((V)->RefCount++)
#define DecStructureReference(V)  ((V)->RefCount--)
#define IsStructureShared(V)      ((V)->RefCount > 1)
#define IsStructureReferenced(V)  ((V)->RefCount > 0)

/*
 * Special set of macros for field data that is having variable sharing between
 * zones tracked. Field data maintains two reference counts: The first,
 * RefCount, is used to keep track of when the field data needs to be
 * deallocated; the second, VarShareRefCount, is used to track variable sharing
 * between zones.
 */
#define IncVarStructureReference(V)  ((V)->VarShareRefCount++)
#define DecVarStructureReference(V)  ((V)->VarShareRefCount--)
#define IsVarStructureShared(V)      ((V)->VarShareRefCount > 1)
#define IsVarStructureReferenced(V)  ((V)->VarShareRefCount > 0)


#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif /* TECPLOTKERNEL */
