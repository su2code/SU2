/****************************************************************************
 *
 * Pointwise Plugin utility functions
 *
 * Proprietary software product of Pointwise, Inc.
 * Copyright (c) 1995-2012 Pointwise, Inc.
 * All rights reserved.
 *
 ***************************************************************************/

#ifndef _SITE_H_
#define _SITE_H_

#if !defined(BUILD_PWPLUGIN_DYNLIB)
#   define BUILD_PWPLUGIN_DYNLIB
#endif

#if !defined(PWP_SITE_GROUPID)
#   define PWP_SITE_GROUPID    59012
#endif

#if !defined(PWP_SITE_GROUPNAME)
#   define PWP_SITE_GROUPNAME  "Stanford ADL"
#endif


/*! \file
\brief Site Build Configuration Settings

The plugin SDK uses the values in this file to "brand" the binaries built at
your site.
*/

/*-----------------------------------------------------*/
/*! \def PWP_SITE_GROUPID
\brief Globally unique site identifier 16-bit integer value.

If you plan on releasing a plugin outside your company, you will need to obtain
a PWP_SITE_GROUPID value from Pointwise support. The plugin SDK combines the
PWP_SITE_GROUPID value to construct globally unique ids (GUID).

For example,
the CAEP_FORMATINFO::id values are initialized in rtCaepInitItems.h using
MAKEGUID() as shown in the code below.
\par
\dontinclude rtCaepInitItems.h
\skip == CAEP_FORMATINFO FormatInfo
\until PWP_FILEDEST_FILENAME

\note
Pointwise uses plugin GUIDs at runtime. Pointwise will not be able to load
plugins properly at runtime if any id conflicts are detected.

\note
If you do not define a PWP_SITE_GROUPID value, the SDK will default to
\ref PWP_GRPID_USER.

\sa MAKEGUID(), apiUtils.h, CAEP_FORMATINFO
*/

/*-----------------------------------------------------*/
/*! \def PWP_SITE_GROUPNAME
\brief Your site's identifier string.

This string represents the company or group name you want associated with this
plugin. For instance, all Pointwise distributed plugins use the group name
"Pointwise". The group name is used by the framework to distinguish between
similarly named exporters from different authors.

For example, as shown in the list below, the Pointwise application will display
a CAE exporter's full name as \sf{group/name}.
- Alpha/CGNS
- Alpha/xml
- Beta/CGNS
- Beta/xml

The SDK returns this value in the CAEP_FORMATINFO::group data member as shown
in the code below.
\par
\dontinclude rtCaepInitItems.h
\skip == CAEP_FORMATINFO FormatInfo
\until MAKEGUID

\note
If you do not define a PWP_SITE_GROUPNAME value, the SDK will default to
\ref PWP_GROUPNAME_DEFAULT.

\sa CAEP_FORMATINFO, PwEnumCaeFormat(), PwCaeFormat()
*/
#if DOXYGEN_RUNNING
#   define PWP_SITE_GROUPID    59012
#endif

#if DOXYGEN_RUNNING
#   define PWP_SITE_GROUPNAME  "Stanford ADL"
#endif

#endif
