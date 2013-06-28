/****************************************************************************
 *
 * Plugin Id Registry
 *
 ***************************************************************************/

/*! \file
 \brief Local Site Plugin Registry File

 By design, every plugin must have a unique ID.

 This file is intended to be a central registry of plugin ids already
 in use by your internally developed plugins.

 When adding a new plugin, add a new id to the end of this list and
 then duplicate that value in the .../YourPlugin/rtCaepInitItems.h
 header file for MAKEGUID() to initialize CAEP_FORMATINFO::id.

 \note It is suggested that you do not include this file in any compilation.
 However, if you do use this as an include file, all plugins will recompile
 every time a new ID is added to this file.
*/
/*! \cond sdkINTERNALS */

#define ID_CaeUnsSU2 30

// add new plugin ID above this line

/*! \endcond sdkINTERNALS */
