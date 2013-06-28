/****************************************************************************
 *
 * Pointwise Plugin utility functions
 *
 * Proprietary software product of Pointwise, Inc.
 * Copyright (c) 1995-2012 Pointwise, Inc.
 * All rights reserved.
 *
 ***************************************************************************/

#ifndef _RTCAEPSUPPORTDATA_H_
#define _RTCAEPSUPPORTDATA_H_

/*! \cond */

/*------------------------------------*/
/* CaeUnsXML format item setup data */
/*------------------------------------*/
CAEP_BCINFO CaeUnsXMLBCInfo[] = {
    { "inflow", 100 },
    { "outflow", 101 },
    { "wall", 102 },
};
/*------------------------------------*/
CAEP_VCINFO CaeUnsXMLVCInfo[] = {
    { "viscous", 200 },
    { "invisid", 201 },
};
/*------------------------------------*/
const char *CaeUnsXMLFileExt[] = {
    "xml"
};

/*! \endcond */



/************************************************************************/
/*! \file
\brief Defines Support Data for the CAEP_RTITEM Array

The file \sf{%rtCaepSupportData.h} defines and initilizes the support data for
the global CAEP_RTITEM \ref caepRtItem[] array. The CAE Plugin SDK uses this data
to implement the functions and behaviors required by the \ref DOXGRP_APICAEP.
If you want to see the SDK implementation details, look in the
\sf{/shared/CAEP/apiCAEP.c} file.

The SDK file \sf{/shared/CAEP/apiCAEP.c} includes \sf{%rtCaepSupportData.h}
prior to the declaration of the \ref caepRtItem[] array as shown below.
\par
\dontinclude apiCAEP.c
\skip "rtCaepSupportData.h"
\until ARRAYSIZE(caepRtItem);

When copied from the \sf{src/plugins/templates/CAEP/} folder to your plugins
project folder, \sf{%rtCaepSupportData.h} will contain the support data needed
for the 3 example CAEP_RTITEM array items. This example support data must be
culled and edited as needed for your plugin's implementation.

The support data used by the CAE Plugin SDK includes the following:

\li Array of valid Boundary Conditions (CAEP_BCINFO) definitions.
\li Array of valid Volume Conditions (CAEP_VCINFO) definitions.
\li Array of valid File Extension definitions.

These support arrays are referenced in rtCaepInitItems.h to initialize the
corresponding data members.

\par Example Support Data Usage

The code segments below show how the example support data is implemented by
the SDK in the \sf{%rtCaepSupportData.h} and rtCaepInitItems.h template files.

The example support data declaration and initialization in
\sf{%rtCaepSupportData.h}:
\par
\dontinclude rtCaepSupportData.h
\skip CaeUnsXMLBCInfo[]
\until "xml"
\skipline };

The support data referenced as initilizers in rtCaepInitItems.h:
\par
\dontinclude rtCaepInitItems.h
\skip == CAEP_BCINFO*
\until ARRAYSIZE(CaeUnsXMLFileExt)
*/

#endif /* _RTCAEPSUPPORTDATA_H_ */
