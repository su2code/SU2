/****************************************************************************
 *
 * Pointwise Plugin utility functions
 *
 * Proprietary software product of Pointwise, Inc.
 * Copyright (c) 1995-2012 Pointwise, Inc.
 * All rights reserved.
 *
 ***************************************************************************/

#ifndef _APICAEGRIDWRITEIMPL_H_
#define _APICAEGRIDWRITEIMPL_H_

#include "apiCAEP.h"
#include "apiCAEPUtils.h"
#include "apiPWP.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! \file
    \brief Starting point for a CAE export.
    
    When a CAE export is initiated, the SDK initializes the plugin and calls
    the function runtimeWrite() in this file.
    
    To implement the exporter logic, The plugin author must edit the copy of
    the file \sf{runtimeWrite.c} located in the
    \sf{src/plugins/MyPlugin} folder.
*/

/*! Called by the SDK to start the export of a grid model.

    \param pRti Pointer to the runtime item instance for this invocation of
                the exporter.
    \param model Handle to the grid model to be exported.
    \param pWriteInfo Pointer to the export settings.

    \return  PWP_TRUE on success.

    \note
        Plugin-specific data members can be added to the CAEP_RTITEM typedef by
        editing the rtCaepInstanceData.h file.

    \sa rtCaepInstanceData.h, rtCaepInitItems.h
*/
PWP_BOOL runtimeWrite(CAEP_RTITEM *pRti, PWGM_HGRIDMODEL model,
                      const CAEP_WRITEINFO *pWriteInfo);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif  /* _APICAEGRIDWRITEIMPL_H_ */
