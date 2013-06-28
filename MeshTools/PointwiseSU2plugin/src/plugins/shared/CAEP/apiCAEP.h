/****************************************************************************
 *
 * Pointwise CAE Plugin API (CAE-API) v1.0
 *
 * Proprietary software product of Pointwise, Inc.
 * Copyright (c) 1995-2012 Pointwise, Inc.
 * All rights reserved.
 *
 ***************************************************************************/

#ifndef _APICAE_H_
#define _APICAE_H_

#include "apiGridModel.h"
#include "apiPWP.h"


#ifdef __cplusplus
extern "C" {
#endif


/*! \file
    \brief Pointwise CAE Plugin API (CAEP-API)

    This defines the API extension required for all Pointwise CAE plugins.
*/


/***********************************************************/
/***********************************************************/
/*! \defgroup DOXGRP_APICAEP Pointwise CAE Plugin API Specification (CAEP-API)
*
*  In addition to the CAEP-API, all CAE Plugins must also
*  implement the Pointwise Plugin API. This is the base API used by
*  Pointwise to load and manage plugins.
*
*  CAE Plugins also use the Pointwise Grid Model API (PWGM-API) for read-only
*  access to the grid model data being exported.
*
* \sa \link DOXGRP_APIPWP PWP-API specification \endlink,
*     \link DOXGRP_APIPWGM PWGM-API specification \endlink
*
*  @{
*/


/***********************************************************/
/***********************************************************/
/*! \ingroup DOXGRP_APIPWP_BASENAMES
*   The CAEP API specification base name.
*/
#define CAEP_API_EXPORT  "Export-CAE"


/***********************************************************/
/***********************************************************/
/*! \defgroup DOXGRP_APICAEP_TYPES Data Types
*  @{
*/

/*---------------------------------------------------------*/
/*! \brief Boundary condition definition information.
*
* \sa PwCaeEnumBCs()
*/
typedef struct CAEP_BCINFO_t {
    /*! \brief BC physical type name.
    */
    const char *phystype;

    /*! \brief BC physical type id.
    */
    PWP_INT32 id;
}
CAEP_BCINFO;


/*---------------------------------------------------------*/
/*! \brief Output destination types.
*
* \sa CAEP_FORMATINFO
*/
typedef enum CAEP_ENUM_FILEDEST_e {
    PWP_FILEDEST_FILENAME, /*!< \sf{/the/path/to/the/full_filename.ext} */
    PWP_FILEDEST_BASENAME, /*!< \sf{/the/path/to/the/base_filename} (no ext) */
    PWP_FILEDEST_FOLDER,   /*!< \sf{/the/path/to/the/folder_name/} (no filename) */
    /* add new group enum values before this line */

    /*! \cond sdkINTERNALS */
    /*! for Bookkeeping only - NOT valid "id" values
    @{ */
    PWP_FILEDEST_SIZE,
    PWP_FILEDEST_LAST = PWP_FILEDEST_SIZE-1
    /* @} */
    /*! \endcond sdkINTERNALS */
}
CAEP_ENUM_FILEDEST;


/*---------------------------------------------------------*/
/*! \brief Supported CAE dimensionality values.
*
* \sa CAEP_FORMATINFO::allowedDimension2D,
*     CAEP_FORMATINFO::allowedDimension3D,
*     CAEP_WRITEINFO::dimension
*/
typedef enum CAEP_ENUM_DIMENSION_e {
    PWP_DIMENSION_2D, /*!< Plugin supports 2D export */
    PWP_DIMENSION_3D, /*!< Plugin supports 3D export */
    /* add new group enum values before this line */

    /*! \cond sdkINTERNALS */
    /*! for Bookkeeping only - NOT valid "id" values
    @{ */
    PWP_DIMENSION_SIZE,
    PWP_DIMENSION_LAST = PWP_DIMENSION_SIZE-1
    /* @} */
    /*! \endcond sdkINTERNALS */
}
CAEP_ENUM_DIMENSION;


/*---------------------------------------------------------*/
/*! \brief Supported CAE file encoding values.
*
* \sa CAEP_FORMATINFO::allowedFileFormatASCII,
*     CAEP_FORMATINFO::allowedFileFormatBinary,
*     CAEP_FORMATINFO::allowedFileFormatUnformatted,
*     CAEP_WRITEINFO::encoding
*/
typedef enum CAEP_ENUM_ENCODING_e {
    PWP_ENCODING_ASCII,       /*!< Plugin supports ASCII export */
    PWP_ENCODING_BINARY,      /*!< Plugin supports binary export */
    PWP_ENCODING_UNFORMATTED, /*!< Plugin supports Unformatted export */
    /* add new group enum values before this line */

    /*! \cond sdkINTERNALS */
    /*! for Bookkeeping only - NOT valid "id" values
    @{ */
    PWP_ENCODING_SIZE,
    PWP_ENCODING_LAST = PWP_ENCODING_SIZE-1
    /* @} */
    /*! \endcond sdkINTERNALS */
}
CAEP_ENUM_ENCODING;


/*---------------------------------------------------------*/
/*! \brief Supported CAE file precision values.
*
* \sa CAEP_FORMATINFO::allowedDataPrecisionSingle,
*     CAEP_FORMATINFO::allowedDataPrecisionDouble,
*     CAEP_WRITEINFO::precision
*/
typedef enum CAEP_ENUM_PRECISION_e {
    PWP_PRECISION_SINGLE, /*!< Plugin supports single precision */
    PWP_PRECISION_DOUBLE, /*!< Plugin supports double precision */
    /* add new group enum values before this line */

    /*! \cond sdkINTERNALS */
    /*! for Bookkeeping only - NOT valid "id" values
    @{ */
    PWP_PRECISION_SIZE,
    PWP_PRECISION_LAST = PWP_PRECISION_SIZE-1
    /* @} */
    /*! \endcond sdkINTERNALS */
}
CAEP_ENUM_PRECISION;


/*---------------------------------------------------------*/
/*! \brief The information returned for each supported CAEP exporter.
*
*   This information is used to integrate the CAE exporter into
*   the framework.
*
* \sa PwEnumCaeFormat(), PwCaeFormat()
*/
typedef struct CAEP_FORMATINFO_t {
    /*! \brief The plugin's group name.
    *
    * The SDK uses PWP_SITE_GROUPNAME as the default.
    */
    const char *group;

    /*! \brief format Name.
    */
    const char *name;

    /*! \brief format guid.
    */
    PWP_UINT32 id;

    /*! \brief Specifies the desired output destination type.
    *
    * Before an export is requested, the framework is responsible to obtain
    * an output destination in the format dictated by this value.
    *
    * The destination string is passed to PwCaeGridWrite() by the framework
    * through the CAEP_WRITEINFO::fileDest parameter.
    */
    CAEP_ENUM_FILEDEST fileDest;

    /*! \brief Set to PWP_TRUE if separate Condition export supported.
    */
    PWP_BOOL allowedExportConditionsOnly;

    /*! \brief Set to PWP_TRUE if VCs are supported.
    */
    PWP_BOOL allowedVolumeConditions;

    /*! \brief Set to PWP_TRUE if ASCII files supported in export.
    */
    PWP_BOOL allowedFileFormatASCII;

    /*! \brief Set to PWP_TRUE if Binary files supported in export.
    */
    PWP_BOOL allowedFileFormatBinary;

    /*! \brief Set to PWP_TRUE if Unformatted files supported in export.
    */
    PWP_BOOL allowedFileFormatUnformatted;

    /*! \brief Set to PWP_TRUE if single precision data supported in export.
    */
    PWP_BOOL allowedDataPrecisionSingle;

    /*! \brief Set to PWP_TRUE if double precision data supported in export.
    */
    PWP_BOOL allowedDataPrecisionDouble;

    /*! \brief Set to PWP_TRUE if 2D exports supported.
    */
    PWP_BOOL allowedDimension2D;

    /*! \brief Set to PWP_TRUE if 3D exports supported.
    */
    PWP_BOOL allowedDimension3D;
}
CAEP_FORMATINFO;


/*---------------------------------------------------------*/
/*! \brief Volume condition definition information.
*
* \sa PwCaeEnumVCs()
*/
typedef struct CAEP_VCINFO_t {
    /*! \brief VC physical type name.
    */
    const char *phystype;

    /*! \brief VC physical type id.
    */
    PWP_UINT32 id;
}
CAEP_VCINFO;


/*---------------------------------------------------------*/
/*! \brief CAE export write control information.
*
* This data is passed to PwCaeGridWrite() by the framework.
*
* \sa PwCaeGridWrite()
*/
typedef struct CAEP_WRITEINFO_t {
    /*! \brief requested file destination.
    */
    const char *fileDest;

    /*! \brief Set to PWP_TRUE if only Conditions exported.
    */
    PWP_BOOL conditionsOnly;

    /*! \brief export file encoding.
    */
    CAEP_ENUM_ENCODING encoding;

    /*! \brief export precision.
    */
    CAEP_ENUM_PRECISION precision;

    /*! \brief export dimensionality.
    */
    CAEP_ENUM_DIMENSION dimension;

}
CAEP_WRITEINFO;

/*---------------------------------------------------------*/
/*! \brief CAEP exporter instance handle.
*
* Exporter handles are returned by PwCreateCaeById() and PwCreateCaeByName().
*
* This handle is opaque and not intended for direct access by the plugin.
*
* \sa PwCreateCaeById(), PwCreateCaeByName()
*/
PWP_DECLARE_HANDLE(CAEP_EXPORTER);

/*! @} */ /* DOXGRP_APICAEP_TYPES */


/***********************************************************/
/***********************************************************/
/*! \defgroup DOXGRP_APICAEP_FUNCTIONS Functions
*  @{
*/

/*---------------------------------------------------------*/
/*! \brief Create CAE exporter instance with given id.
*
* \param id
*    The plugin-defined CAE exporter id.
*
* \return The CAE exporter handle or NULL if id is invalid.
*
* \sa PwEnumCaeFormat()
*
* \note Valid CAE exporter id's can be obtained from
*         PwEnumCaeFormat(pFormatInfo).
*/
PWP_PROTOTYPE_DECL CAEP_EXPORTER
PwCreateCaeById(PWP_UINT32 id);

/*---------------------------------------------------------*/
/*! \brief Create CAE exporter instance with given name.
*
* \param name
*    The plugin-defined CAE exporter name.
*
* \return The CAE exporter handle or NULL if name is invalid.
* \sa PwEnumCaeFormat()
* \note Valid CAE exporter names can be obtained from
*       PwEnumCaeFormat(pFormatInfo).
*/
PWP_PROTOTYPE_DECL CAEP_EXPORTER
PwCreateCaeByName(const char name[]);

/*---------------------------------------------------------*/
/*! \brief Destroy CAE exporter instance.
*
* \param handle
*    Pointer to a CAEP_EXPORTER handle obtained from PwCreateCaeById() or
*    PwCreateCaeByName().
*
* \sa PwCreateCaeByName(), PwCreateCaeById()
*/
PWP_PROTOTYPE_DECL PWP_VOID
PwDestroyCae(CAEP_EXPORTER *handle);

/*---------------------------------------------------------*/
/*! \brief Enumerate CAEP_FORMATINFO data for all supported CAE exporters.
*
* \param ndx
*    The format index starting with 0.
*
* \param pFormatInfo
*    Pointer to a CAEP_FORMATINFO buffer.
*
* \return The CAE format name (same as info.name). NULL if ndx is not valid.
* \sa PwCaeFormat()
*
* \par "Sample usage:"
* \code
*    CAEP_FORMATINFO info;
*    PWP_UINT32 ndx = 0;
*    while (caep->PwEnumCaeFormat(ndx++, &info)) {
*       ...do something here...
*    }
* \endcode
*/
PWP_PROTOTYPE_DECL const char*
PwEnumCaeFormat(PWP_UINT32 ndx, CAEP_FORMATINFO *pFormatInfo);

/*---------------------------------------------------------*/
/*! \brief Get the number of supported CAE exporters.
*
*   \return The number of supported exporters.
*/
PWP_PROTOTYPE_DECL PWP_UINT32
PwGetCaeFormatCount(); 

/*---------------------------------------------------------*/
/*! \brief Get CAEP_FORMATINFO data for a CAE exporter handle.
*
* \param handle
*    A CAEP_EXPORTER handle obtained from PwCreateCaeById() or
*    PwCreateCaeByName().
*
* \param pFormatInfo
*    Pointer to a CAEP_FORMATINFO buffer.
*
* \return The CAE format name (same as info.name). NULL if ndx is not valid.
*
* \sa PwCaeFormat()
*
* \par "Sample usage:"
* \code
*    CAEP_EXPORTER h = PwCreateCaeByName("name");
*    CAEP_FORMATINFO info;
*    if (caep->PwCaeFormat(h, &info)) {
*       ...do something here...
*    }
* \endcode
*/
PWP_PROTOTYPE_DECL const char*
PwCaeFormat(CAEP_EXPORTER handle, CAEP_FORMATINFO *pFormatInfo);

/*---------------------------------------------------------*/
/*! \brief Test if CAE exporter instance supports the given element type.
*
* \param handle
*    A CAEP_EXPORTER handle obtained from PwCreateCaeById() or
*    PwCreateCaeByName().
*
* \param which
*    The PWGM_ENUM_ELEMTYPE to test.
*
*   \return PWP_TRUE if element type which is supported.
*/
PWP_PROTOTYPE_DECL PWP_BOOL
PwCaeElementType(CAEP_EXPORTER handle, PWGM_ENUM_ELEMTYPE which);

/*---------------------------------------------------------*/
/*! \brief Enumerate CAEP_BCINFO data for a CAE exporter instance.
*
* \param handle
*    A CAEP_EXPORTER handle obtained from PwCreateCaeById() or
*    PwCreateCaeByName().
*
* \param ndx
*    The BC index starting with 0.
*
* \param pBCInfo
*    Pointer to a CAEP_BCINFO buffer.
*
*   \return The BC name (same as pBCInfo->phystype). NULL if
*   handle or ndx is not valid.
*/
PWP_PROTOTYPE_DECL const char*
PwCaeEnumBCs(CAEP_EXPORTER handle, PWP_UINT32 ndx, CAEP_BCINFO *pBCInfo);

/*---------------------------------------------------------*/
/*! \brief Enumerate supported file extensions for a CAE exporter instance.
*
* \param handle
*    A CAEP_EXPORTER handle obtained from PwCreateCaeById() or
*    PwCreateCaeByName().
*
* \param ndx
*    The file extension index starting with 0.
*
*   \return file extension or null if handle or ndx is invalid.
*/
PWP_PROTOTYPE_DECL const char*
PwCaeEnumFileExt(CAEP_EXPORTER handle, PWP_UINT32 ndx);

/*---------------------------------------------------------*/
/*! \brief Enumerate CAEP_VCINFO data for a CAE exporter instance.
*
* \param handle
*    A CAEP_EXPORTER handle obtained from PwCreateCaeById() or
*    PwCreateCaeByName().
*
* \param ndx
*    The VC index starting with 0.
*
* \param pVCInfo
*    Pointer to a CAEP_VCINFO buffer.
*
* \return The VC name (same as pVCInfo->phystype). NULL if
* handle or ndx is not valid.
*/
PWP_PROTOTYPE_DECL const char*
PwCaeEnumVCs(CAEP_EXPORTER handle, PWP_UINT32 ndx, CAEP_VCINFO *pVCInfo);

/*---------------------------------------------------------*/
/*! \brief Get the number of BC's for a CAE exporter instance.
*
* \param handle
*    A CAEP_EXPORTER handle obtained from PwCreateCaeById() or
*    PwCreateCaeByName().
*
*   \return number of BC's
*/
PWP_PROTOTYPE_DECL PWP_UINT32
PwCaeGetBCCount(CAEP_EXPORTER handle);

/*---------------------------------------------------------*/
/*! \brief Get the number of supported file extensions for a CAE exporter
* instance.
*
* \param handle
*    A CAEP_EXPORTER handle obtained from PwCreateCaeById() or
*    PwCreateCaeByName().
*
*   \return number of file extensions
*/
PWP_PROTOTYPE_DECL PWP_UINT32
PwCaeGetFileExtCount(CAEP_EXPORTER handle);

/*---------------------------------------------------------*/
/*! \brief Get the number of VC's for a CAE exporter instance.
*
* \param handle
*    A CAEP_EXPORTER handle obtained from PwCreateCaeById() or
*    PwCreateCaeByName().
*
*   \return number of VC's
*/
PWP_PROTOTYPE_DECL PWP_UINT32
PwCaeGetVCCount(CAEP_EXPORTER handle);

/*---------------------------------------------------------*/
/*! \brief Initiates writing a grid model.
*
* Instructs a plugin to write the given grid model according to the specified
* CAEP_WRITEINFO settings.
*
* \param handle
*    A CAEP_EXPORTER handle obtained from PwCreateCaeById() or
*    PwCreateCaeByName().
*
* \param model
*    A PWGM_HGRIDMODEL handle. The plugin uses this handle to access the grid
*    model data.
*
* \param pWriteInfo
*    Pointer to a CAEP_WRITEINFO settings buffer.
*
* \return PWP_FALSE if write failed.
*
* \sa \link DOXGRP_APIPWGM PWGM-API specification \endlink
*/
PWP_PROTOTYPE_DECL PWP_BOOL
PwCaeGridWrite(CAEP_EXPORTER handle, PWGM_HGRIDMODEL model,
                        const CAEP_WRITEINFO *pWriteInfo);

/*! @} */ /* DOXGRP_APICAEP_FUNCTIONS */

/*! @} */ /* DOXGRP_APICAEP */

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* !_APICAE_H_ */
