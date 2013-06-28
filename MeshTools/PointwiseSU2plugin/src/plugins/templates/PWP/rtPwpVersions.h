/****************************************************************************
 *
 * Pointwise Plugin utility functions
 *
 * Proprietary software product of Pointwise, Inc.
 * Copyright (c) 1995-2012 Pointwise, Inc.
 * All rights reserved.
 *
 ***************************************************************************/


#ifndef _RTPWPVERSIONS_H_
#define _RTPWPVERSIONS_H_


/*------------------------------------------------
   PWP api version with which plugin conforms
------------------------------------------------*/
/*! \brief The PWP-API major version value.
*/
#define VERSION_PWP_MAJOR   1

/*! \brief The PWP-API minor version value.
*/
#define VERSION_PWP_MINOR   0

/*! \brief This macro is used for static initialization of a PWP_VERSION
struct to the current \sf{VERSION_PWP_MAJOR} and \sf{VERSION_PWP_MINOR} values.
*/
#define VERSION_PWP_INIT    {VERSION_PWP_MAJOR, VERSION_PWP_MINOR}


/*------------------------------------------------
        plugin software release version
------------------------------------------------*/
/*! \brief The software release major version value.
*/
#define VERSION_LIB_MAJOR   1

/*! \brief The software release minor version value.
*/
#define VERSION_LIB_MINOR   0

/*! \brief This macro is used for static initialization of a PWP_VERSION
struct to the current \sf{VERSION_LIB_MAJOR} and \sf{VERSION_LIB_MINOR} values.
*/
#define VERSION_LIB_INIT    {VERSION_LIB_MAJOR, VERSION_LIB_MINOR}


/************************************************************************/
/*! \file
\brief Defines Implementation Version Information

This file contains 2 version value macro pairs. These values are used
throughout the Plugin SDK. One pair (VERSION_PWP_ prefix) defines the PWP-API
version. The other pair (VERSION_LIB_ prefix) defines the plugin binary release
version.

The plugin developer should edit the VERSION_PWP_MAJOR, VERSION_PWP_MINOR,
VERSION_LIB_MAJOR, and VERSION_LIB_MINOR values as needed.

*/

#endif /* _RTPWPVERSIONS_H_ */
