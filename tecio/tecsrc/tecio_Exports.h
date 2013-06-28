#ifndef tecio_EXPORTS_H
#define tecio_EXPORTS_H

/*
 * This file is a compile shim to allow legacy Tecplot to share common source from the SDK's engine
 * and have a common DLL export interface.
 */

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#else
    #if defined MSWIN
        #if defined MAKEARCHIVE
            #define tecio_API _declspec ( dllexport )
        #else /* !TECPLOTKERNAL && !MAKEARCHIVE */
            #define tecio_API _declspec ( dllimport )
        #endif
    #else
        #define tecio_API
    #endif
#endif

#endif
