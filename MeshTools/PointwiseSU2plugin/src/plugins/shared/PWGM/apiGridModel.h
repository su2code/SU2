#ifndef _APIGRIDMODEL_H_
#define _APIGRIDMODEL_H_
/****************************************************************************
 *
 * Pointwise PWGM API (PWGM-API) v1.0
 *
 * Proprietary software product of Pointwise, Inc.
 * Copyright (c) 1995-2012 Pointwise, Inc.
 * All rights reserved.
 *
 ***************************************************************************/

#include "apiPWP.h"

#if defined(WINDOWS) && defined(SHOW_PWGM_MESSAGES)
#    if !defined(PWGM_HIDE_UNSTRUCTURED_API) && defined(PWGM_HIDE_STRUCTURED_API)
#        pragma message ("FYI: Using Unstructured PWGM API")
#    endif
#    if defined(PWGM_HIDE_UNSTRUCTURED_API) && !defined(PWGM_HIDE_STRUCTURED_API)
#        pragma message ("FYI: Using Structured PWGM API")
#    endif
#    if !defined(PWGM_HIDE_UNSTRUCTURED_API) && !defined(PWGM_HIDE_STRUCTURED_API)
#        pragma message ("WARNING: Using Hybrid PWGM API")
#    endif
#    if defined(PWGM_HIDE_UNSTRUCTURED_API) && defined(PWGM_HIDE_STRUCTURED_API)
#        error ("ERROR: PWGM API is INVISIBLE!")
#    endif
#endif /* WINDOWS && SHOW_PWP_MESSAGES */


/*! \file
    \brief Pointwise Grid Model API Specification (PWGM-API)

    Defines the API used by plugins for read-only access to the grid model
    data.
*/


/*! \cond sdkINTERNALS */
/*! Macro that controls point-of-view for the function signatures and
    public data. For most OS's, this macro does nothing. Used primarily on
    winOS builds.

    For plugin developers, this macro will "export" all functions and data
    from the plugin library module.

    For plugin users, like Pointwise, this macro will "import" all functions
    and data from the plugin library module.
*/
#undef PW_DLL_IMPEXP
#ifdef BUILD_GGPLUGIN
# define PW_DLL_IMPEXP PW_DLLEXPORT
#else
# define PW_DLL_IMPEXP PW_DLLIMPORT
#endif
/*! \endcond sdkINTERNALS */


#ifdef __cplusplus
extern "C" {
#endif

/***********************************************************/
/***********************************************************/
/*! \defgroup DOXGRP_APIPWGM Pointwise Grid Model API Specification (PWGM-API)
*
*  Plugins may use the Pointwise Grid Model API Specification (PWGM-API) to
*  gain read-only access to the grid model data.
*
*  All access to the grid model starts with a \c PWGM_HGRIDMODEL handle. This
*  handle is passed to the plugin by the framework.
*
*  See the \ref DOXGRP_APIPWGM_DATAHANDLES section for more information.
@{
*/

/***********************************************************/
/***********************************************************/
/*! \defgroup DOXGRP_APIPWGM_UNSTRUCTURED Unstructured Grid Model
*
*  An unstructured grid model contains one or more block (PWGM_HBLOCK),
*  domain (PWGM_HDOMAIN), and vertex (PWGM_HVERTEX) instances.
*
*  A block contains one or more 3D, volumetric elements (PWGM_HELEMENT).
*
*  A domain contains one or more 2D, planar elements (PWGM_HELEMENT).
*
*  A vertex defines an XYZ location in space (PWGM_VERTDATA).
*
*  An element is defined by a list of model vertices (PWGM_ELEMDATA).
*
*  Unstructured grid model blocks, domains, and vertices are accessed using the
*  following <tt>PwModXxxx()</tt> functions:
*  \li <tt>PWP_UINT32 PwModBlockCount (PWGM_HGRIDMODEL model)</tt>
*  \li <tt>PWP_UINT32 PwModDomainCount (PWGM_HGRIDMODEL model)</tt>
*  \li <tt>PWGM_HBLOCK PwModEnumBlocks (PWGM_HGRIDMODEL model,
            PWP_UINT32 ndx)</tt>
*  \li <tt>PWGM_HDOMAIN PwModEnumDomains (PWGM_HGRIDMODEL model,
            PWP_UINT32 ndx)</tt>
*  \li <tt>PWGM_HVERTEX PwModEnumVertices (PWGM_HGRIDMODEL model,
            PWP_UINT32 ndx)</tt>
*  \li <tt>PWP_UINT32 PwModVertexCount (PWGM_HGRIDMODEL model)</tt>
*
*  The elements of a block are accessed using the <tt>PwBlkXxxx()</tt>
*  functions:
*  \li <tt>PWP_UINT32 PwBlkElementCount (PWGM_HBLOCK block,
            PWGM_ELEMCOUNTS *pCounts)</tt>
*  \li <tt>PWGM_HELEMENT PwBlkEnumElements (PWGM_HBLOCK block,
            PWP_UINT32 ndx)</tt>
*  \li <tt>PWP_BOOL PwBlkCondition (PWGM_HBLOCK block,
            PWGM_CONDDATA *pCondData)</tt>
*
*  The elements of a domain are accessed using the <tt>PwDomXxxx()</tt>
*  functions:
*  \li <tt>PWP_UINT32 PwDomElementCount (PWGM_HDOMAIN domain,
            PWGM_ELEMCOUNTS *pCounts)</tt>
*  \li <tt>PWGM_HELEMENT PwDomEnumElements (PWGM_HDOMAIN domain,
            PWP_UINT32 ndx)</tt>
*  \li <tt>PWP_BOOL PwDomCondition (PWGM_HDOMAIN domain,
            PWGM_CONDDATA *pCondData)</tt>
*
*  An element's vertex references (topology) are accessed using the
*  <tt>PwElemXxxx()</tt> functions:
*  \li <tt>PWP_BOOL PwElemDataMod (PWGM_HELEMENT element,
            PWGM_ELEMDATA *pElemData)</tt>
*
*  A vertex's data is accessed using the <tt>PwVertXxxx()</tt> functions:
*  \li <tt>PWP_BOOL PwVertDataMod (PWGM_HVERTEX vertex,
            PWGM_VERTDATA *pVertData)</tt>
*  \li <tt>PWP_BOOL PwVertIndexMod (PWGM_HVERTEX vertex,
            PWP_UINT32 *pIndex)</tt>
*  \li <tt>PWP_BOOL PwVertXyzVal (PWGM_HVERTEX vertex, PWGM_ENUM_XYZ which,
            PWGM_XYZVAL *pVal)</tt>
*
* \sa \link DOXGRP_APICAEP CAEP-API specification \endlink
*/


/***********************************************************/
/***********************************************************/
/*! \defgroup DOXGRP_APIPWGM_STRUCTURED Structured Grid Model
*
*  A structured grid model contains one or more block (PWGM_HBLOCK),
*  boundary (PWGM_HBNDRY), and connection (PWGM_CNXNDATA) instances.
*
*  A block is defined by an array of PWGM_VERTDATA instances and a
*  PWGM_BLOCKDATA instance.
*
*  A boundary is defined by a PWGM_BNDRYDATA instance.
*
*  A connection is defined by a PWGM_CNXNDATA instance.
*
*  Structured grid model blocks, and vertices are accessed using the following
*  functions:
*  \li <tt>PWP_UINT32 PwModBlockCount (PWGM_HGRIDMODEL model)</tt>
*  \li <tt>PWGM_HBLOCK PwModEnumBlocks (PWGM_HGRIDMODEL model, PWP_UINT32 ndx)</tt>
*  \li <tt>PWP_BOOL PwBlkSize(PWGM_HBLOCK block, PWGM_STR_SIZE *pSize);</tt>
*  \li <tt>PWP_BOOL PwBlock(PWGM_HBLOCK block, PWGM_BLOCKDATA *pBlockData);</tt>
*  \li <tt>PWP_BOOL PwBlkNdxVertData(PWGM_HBLOCK block, PWGM_INDEX3 ndx3, PWGM_VERTDATA *pVertData);</tt>
*
*  \li <tt>PWP_UINT32 PwModConnectionCount(PWGM_HGRIDMODEL model);</tt>
*  \li <tt>PWGM_HCNXN PwModEnumConnections(PWGM_HGRIDMODEL model, PWP_UINT32 ndx);</tt>
*  \li <tt>PWP_BOOL PwConnection(PWGM_HCNXN connection, PWGM_CNXNDATA *pCnxnData);</tt>
*  \li <tt>PWP_BOOL PwModNdxConnection(PWGM_HGRIDMODEL model, PWP_UINT32 ndx, PWGM_CNXNDATA *pCnxnData);</tt>
*
*  \li <tt>PWP_UINT32 PwModBoundaryCount(PWGM_HGRIDMODEL model);</tt>
*  \li <tt>PWGM_HBNDRY PwModEnumBoundaries(PWGM_HGRIDMODEL model, PWP_UINT32 ndx);</tt>
*  \li <tt>PWP_BOOL PwBoundary(PWGM_HBNDRY boundary, PWGM_BNDRYDATA *pBndryData);</tt>
*  \li <tt>PWP_BOOL PwBndryCondition(PWGM_HBNDRY boundary, PWGM_CONDDATA *pCondData);</tt>
*  \li <tt>PWP_BOOL PwModNdxBoundary(PWGM_HGRIDMODEL model, PWP_UINT32 ndx, PWGM_BNDRYDATA *pBndryData);</tt>
*  \li <tt>PWP_BOOL PwModNdxBoundaryAndCondition(PWGM_HGRIDMODEL model, PWP_UINT32 ndx, PWGM_BNDRYDATA *pBndryData, PWGM_CONDDATA *pCondData);</tt>
*
*  \li <tt>PWP_UINT32 PwBlkBoundaryCount(PWGM_HBLOCK block);</tt>
*  \li <tt>PWGM_HBNDRY PwBlkEnumBoundaries(PWGM_HBLOCK block, PWP_UINT32 ndx);</tt>
*  \li <tt>PWP_BOOL PwBlkNdxBoundary(PWGM_HBLOCK block, PWP_UINT32 ndx, PWGM_BNDRYDATA *pBndryData);</tt>
*  \li <tt>PWP_BOOL PwBlkNdxBoundaryAndCondition(PWGM_HBLOCK block, PWP_UINT32 ndx, PWGM_BNDRYDATA *pBndryData, PWGM_CONDDATA *pCondData);</tt>
*
*  \li <tt>PWP_UINT32 PwBlkConnectionCount(PWGM_HBLOCK block);</tt>
*  \li <tt>PWGM_HCNXN PwBlkEnumConnections(PWGM_HBLOCK block, PWP_UINT32 ndx);</tt>
*  \li <tt>PWP_BOOL PwBlkNdxConnection(PWGM_HBLOCK block, PWP_UINT32 ndx, PWGM_CNXNDATA *pCnxnData);</tt>
*
* \sa \link DOXGRP_APICAEP CAEP-API specification \endlink
*/


/***********************************************************/
/***********************************************************/
/*! \defgroup DOXGRP_APIPWGM_TYPES Data Types
@{
*/

/***********************************************************/
/***********************************************************/
/*! \defgroup DOXGRP_APIPWGM_DATAHANDLES PWGM-API Opaque Data Handle Types
@{
*   Opaque handles used to access grid model data.
*
*   \par "Unstructured Grid Hierarchy:"
\verbatim

    PWGM_HGRIDMODEL
    |
    |---PWGM_HBLOCK[]
    |   |
    |   `---PWGM_HELEMENT[] (3D tet, hex, wedge, or pyramid elements)
    |
    |---PWGM_HDOMAIN[]
    |   |
    |   `---PWGM_HELEMENT[] (2D tri or quad elements)
    |
    `---PWGM_HVERTEX[]

\endverbatim
*
*   \par "Structured Grid Hierarchy:"
\verbatim

    PWGM_HGRIDMODEL
    |
    |---PWGM_HBLOCK[]
    |   |
    |   |---PWGM_VERTDATA[]
    |   |
    |   |---PWGM_HCNXN[]
    |   |   |
    |   |   `---PWGM_CNXNDATA
    |   |
    |   `---PWGM_HBNDRY[]
    |       |
    |       `---PWGM_BNDRYDATA
    |
    |---PWGM_HCNXN[]
    |   |
    |   `---PWGM_CNXNDATA
    |
    `---PWGM_HBNDRY[]
        |
        `---PWGM_BNDRYDATA

\endverbatim
*/


/*! \brief An opaque handle to a grid model.
\sa PWGM_HGRIDMODEL_ISVALID(), PWGM_HGRIDMODEL_INIT
*/
PWP_DECLARE_HANDLE(PWGM_HGRIDMODEL);
/*! \brief returns non-zero value if handle is valid */
#define PWGM_HGRIDMODEL_ISVALID(h) PWP_HANDLE_ISVALID(h)
/*! \cond sdkINTERNALS */
/*! \brief static data init value */
#define PWGM_HGRIDMODEL_INIT       PWP_HANDLE_INIT
/*! \brief assigns internal handle values */
#define PWGM_HGRIDMODEL_SET(h, v)  PWP_HANDLE_SET(h, v)
/*! \brief obtains the underlying PWGM_HGRIDMODEL_IMPL ptr from the handle */
#define PWGM_HGRIDMODEL_GMIMPL(h) ((PWGM_HGRIDMODEL_IMPL*)model)
/*! \endcond sdkINTERNALS */

/*! \brief An opaque handle to a grid block element.
*/
PWP_DECLARE_HELEMGROUP(PWGM_HGRIDMODEL, PWGM_HBLOCK);

/*! \brief Returns non-zero value if handle is valid. */
#define PWGM_HBLOCK_ISVALID(h)    PWP_HEGRP_ISVALID(h)
/*! \cond sdkINTERNALS */
/*! \brief handle init value */
#define PWGM_HBLOCK_INIT          PWP_HEGRP_INIT
/*! \brief assigns internal handle values */
#define PWGM_HBLOCK_SET(h, p, v)  PWP_HEGRP_SET(h, p, v)
/*! \brief obtains the underlying PWGM_HGRIDMODEL_IMPL ptr from the handle */
#define PWGM_HBLOCK_GMIMPL(h) ((PWGM_HGRIDMODEL_IMPL*)PWP_HEGRP_H(h))
/*! \endcond sdkINTERNALS */
/*! \brief obtains the element's parent PWGM_HGRIDMODEL handle */
#define PWGM_HBLOCK_H(h)          PWP_HEGRP_H(h) 
/*! \brief obtains the block's guid from the handle */
#define PWGM_HBLOCK_ID(h)         PWP_HEGRP_ID(h)

/*! \brief An opaque handle to a grid domain element.
*/
PWP_DECLARE_HELEMGROUP(PWGM_HGRIDMODEL, PWGM_HDOMAIN);
/*! \brief returns non-zero value if handle is valid */
#define PWGM_HDOMAIN_ISVALID(h)    PWP_HEGRP_ISVALID(h)
/*! \cond sdkINTERNALS */
/*! \brief handle init value */
#define PWGM_HDOMAIN_INIT          PWP_HEGRP_INIT
/*! \brief assigns internal handle values */
#define PWGM_HDOMAIN_SET(h, p, v)  PWP_HEGRP_SET(h, p, v)
/*! \brief obtains the underlying PWGM_HGRIDMODEL_IMPL ptr from the handle */
#define PWGM_HDOMAIN_GMIMPL(h) ((PWGM_HGRIDMODEL_IMPL*)PWP_HEGRP_H(h))
/*! \endcond sdkINTERNALS */
/*! \brief obtains the element's parent PWGM_HGRIDMODEL handle */
#define PWGM_HDOMAIN_H(h)          PWP_HEGRP_H(h) 
/*! \brief obtains the domain's guid from the handle */
#define PWGM_HDOMAIN_ID(h)         PWP_HEGRP_ID(h)

/*! \brief An opaque handle to a grid vertex element
*/
PWP_DECLARE_HELEMGROUP(PWGM_HGRIDMODEL, PWGM_HVERTEX);
/*! \brief returns non-zero value if handle is valid */
#define PWGM_HVERTEX_ISVALID(h)    PWP_HEGRP_ISVALID(h)
/*! \cond sdkINTERNALS */
/*! \brief handle init value */
#define PWGM_HVERTEX_INIT          PWP_HEGRP_INIT
/*! \brief assigns internal handle values */
#define PWGM_HVERTEX_SET(h, p, v)  PWP_HEGRP_SET(h, p, v)
/*! \brief obtains the underlying PWGM_HGRIDMODEL_IMPL ptr from the handle */
#define PWGM_HVERTEX_GMIMPL(h) ((PWGM_HGRIDMODEL_IMPL*)PWP_HEGRP_H(h))
/*! \endcond sdkINTERNALS */
/*! \brief obtains the element's parent PWGM_HGRIDMODEL handle */
#define PWGM_HVERTEX_H(h)          PWP_HEGRP_H(h) 
/*! \brief obtains the vertex's guid from the handle */
#define PWGM_HVERTEX_ID(h)         PWP_HEGRP_ID(h)

/*! \brief Only used as a generic base handle for PWGM_HELEMENT
*/
PWP_DECLARE_HELEMGROUP(PWGM_HGRIDMODEL, PWGM_HELEMENT_BASE);


/*! \brief Grid element handle declaration */
PWP_DECLARE_HEGRPITEM(PWGM_HELEMENT_BASE, PWGM_HELEMENT);
/*! \brief returns non-zero value if handle is valid */
#define PWGM_HELEMENT_ISVALID(h)        PWP_HEGI_ISVALID(h)
/*! \cond sdkINTERNALS */
/*! \brief handle init value */
#define PWGM_HELEMENT_INIT              PWP_HEGI_INIT
/*! \brief assigns internal handle values */
#define PWGM_HELEMENT_SET(h,p,tt,sv,v)  PWP_HEGI_SET(h,p,tt,sv,v)
/*! \brief obtains the underlying PWGM_HGRIDMODEL_IMPL ptr from the handle */
#define PWGM_HELEMENT_GMIMPL(h) ((PWGM_HGRIDMODEL_IMPL*)PWP_HEGI_H(h))
/*! \endcond sdkINTERNALS */
/*! \brief obtains the element's parent PWGM_HBLOCK or PWGM_HDOMAIN handle */
#define PWGM_HELEMENT_H(h)              PWP_HEGI_H(h)
/*! \brief obtains the element's parent id */
#define PWGM_HELEMENT_PID(h)            PWP_HEGI_PID(h)
/*! \brief obtains the element's parent handle type */
#define PWGM_HELEMENT_PTYPE(h)          PWP_HEGI_PTYPE(h)
/*! \brief obtains the element's guid from the handle */
#define PWGM_HELEMENT_ID(h)             PWP_HEGI_ID(h)


/*! \brief An opaque handle to a structured block boundary.
*/
PWP_DECLARE_HELEMGROUP(PWGM_HGRIDMODEL, PWGM_HBNDRY);
/*! \brief returns non-zero value if handle is valid */
#define PWGM_HBNDRY_ISVALID(h)    PWP_HEGRP_ISVALID(h)
/*! \cond sdkINTERNALS */
    /*! \brief handle init value */
#   define PWGM_HBNDRY_INIT          PWP_HEGRP_INIT
    /*! \brief assigns internal handle values */
#   define PWGM_HBNDRY_SET(h, p, v)  PWP_HEGRP_SET(h, p, v)
    /*! \brief obtains the underlying PWGM_HGRIDMODEL_IMPL ptr from the handle */
#   define PWGM_HBNDRY_GMIMPL(h) ((PWGM_HGRIDMODEL_IMPL*)PWP_HEGRP_H(h))
/*! \endcond sdkINTERNALS */
/*! \brief obtains the element's parent PWGM_HGRIDMODEL handle */
#define PWGM_HBNDRY_H(h)          PWP_HEGRP_H(h) 
/*! \brief obtains the boundary's guid from the handle */
#define PWGM_HBNDRY_ID(h)         PWP_HEGRP_ID(h)


/*! \brief An opaque handle to a structured, inter-block connection.
*/
PWP_DECLARE_HELEMGROUP(PWGM_HGRIDMODEL, PWGM_HCNXN);
/*! \brief returns non-zero value if handle is valid */
#define PWGM_HCNXN_ISVALID(h)    PWP_HEGRP_ISVALID(h)
/*! \cond sdkINTERNALS */
    /*! \brief handle init value */
#   define PWGM_HCNXN_INIT          PWP_HEGRP_INIT
    /*! \brief assigns internal handle values */
#   define PWGM_HCNXN_SET(h, p, v)  PWP_HEGRP_SET(h, p, v)
    /*! \brief obtains the underlying PWGM_HGRIDMODEL_IMPL ptr from the handle */
#   define PWGM_HCNXN_GMIMPL(h) ((PWGM_HGRIDMODEL_IMPL*)PWP_HEGRP_H(h))
/*! \endcond sdkINTERNALS */
/*! \brief obtains the element's parent PWGM_HGRIDMODEL handle */
#define PWGM_HCNXN_H(h)          PWP_HEGRP_H(h) 
/*! \brief obtains the connection's guid from the handle */
#define PWGM_HCNXN_ID(h)         PWP_HEGRP_ID(h)

/*!@} DOXGRP_APIPWGM_DATAHANDLES */


/*---------------------------------------------------------*/
/*       UNSTRUCTURED AND STRUCTURED GRID DATATYPES        */
/*---------------------------------------------------------*/

/*---------------------------------------------------------*/
/*! \brief XYZ component data type.
*
* \sa PwVertXyzVal(), PWGM_VERTDATA
*/
typedef PWP_REAL PWGM_XYZVAL;

/*---------------------------------------------------------*/
/*! \brief Vertex descriptor data type.
*
* \sa PwVertXyzVal(), PwVertDataMod(), PwModEnumVertices()
*/
typedef struct PW_DLL_IMPEXP PWGM_VERTEX_t {
    PWGM_XYZVAL x; /*!< x-component */
    PWGM_XYZVAL y; /*!< y-component */
    PWGM_XYZVAL z; /*!< z-component */
    PWP_UINT32  i; /*!< vertex index in parent's model-space */
}
PWGM_VERTDATA;

/*---------------------------------------------------------*/
/*! \brief Condition descriptor data type.
*
* Block-volume or Domain-boundary condition descriptor data type.
*
* \sa PwBlkCondition(), PwDomCondition()
*/
typedef struct PW_DLL_IMPEXP PWGM_CONDDATA_t {
    const char *name;  /*!< grid-defined condition name */
    PWP_UINT32 id;     /*!< grid-defined condition id */
    const char *type;  /*!< cae-defined condition physical type name */
    PWP_UINT32 tid;    /*!< cae-defined condition id */
}
PWGM_CONDDATA;


/*---------------------------------------------------------*/
/*             UNSTRUCTURED GRID DATATYPES                 */
/*---------------------------------------------------------*/

/*---------------------------------------------------------*/
/*! \brief Element type ids
*
* \sa PWGM_ELEMDATA, PwElemDataMod()
*/
typedef enum PWGM_ENUM_ELEMTYPE_e {
    PWGM_ELEMTYPE_BAR,     /*!< 1D, linear grid element */
    PWGM_ELEMTYPE_HEX,     /*!< 3D, 6-sided (block) grid element */
    PWGM_ELEMTYPE_QUAD,    /*!< 2D, 4-sided grid element */
    PWGM_ELEMTYPE_TRI,     /*!< 2D, 3-sided grid element */
    PWGM_ELEMTYPE_TET,     /*!< 3D, 4-sided (tetrahedral) grid element */
    PWGM_ELEMTYPE_WEDGE,   /*!< 3D, extruded, tri/quad grid element */
    PWGM_ELEMTYPE_PYRAMID, /*!< 3D, 5-sided, quad-based grid element */
    /* add new enum values before this line */

    /*! \cond sdkINTERNALS */
    /*! for Bookkeeping only - NOT valid "element type" values
    @{ */
    PWGM_ELEMTYPE_SIZE,
    PWGM_ELEMTYPE_LAST = PWGM_ELEMTYPE_SIZE-1
    /* @} */
    /*! \endcond sdkINTERNALS */
}
PWGM_ENUM_ELEMTYPE;

/*---------------------------------------------------------*/
/*! \brief Maximum number of verts allowed per element.
*
* \sa PWGM_ELEMDATA
*/
#define PWGM_ELEMDATA_VERT_SIZE  8

/*---------------------------------------------------------*/
/*! \brief Element descriptor data type.
*
* \sa PwElemDataMod()
*/
typedef struct PW_DLL_IMPEXP PWGM_ELEMDATA_t {
    PWGM_ENUM_ELEMTYPE type;
    PWP_UINT32         vertCnt;
    PWGM_HVERTEX       vert[PWGM_ELEMDATA_VERT_SIZE];
    PWP_UINT32         index[PWGM_ELEMDATA_VERT_SIZE];
}
PWGM_ELEMDATA;

/*---------------------------------------------------------*/
/*! \brief Element count information
*
* \sa PwBlkElementCount(), PwDomElementCount()
*/
typedef struct PW_DLL_IMPEXP PWGM_ELEMCOUNTS_t {
    /*! \brief Array of element counts.
    *
    * Index the array using the PWGM_ELEMTYPE_XXX constants or the
    * PWGM_ECNT_Xxx() macros.
    */
    PWP_UINT32 count[PWGM_ELEMTYPE_SIZE];
}
PWGM_ELEMCOUNTS;

/*! \brief Extract the Bar count from a PWGM_ELEMCOUNTS struct.
*
* \param ecs
*    A PWGM_ELEMCOUNTS struct.
*/
#define PWGM_ECNT_Bar(ecs)     (ecs).count[PWGM_ELEMTYPE_BAR]
/*! \brief Extract the Hex count from a PWGM_ELEMCOUNTS struct.
*
* \param ecs
*    A PWGM_ELEMCOUNTS struct.
*/
#define PWGM_ECNT_Hex(ecs)     (ecs).count[PWGM_ELEMTYPE_HEX]
/*! \brief Extract the Quad count from a PWGM_ELEMCOUNTS struct.
*
* \param ecs
*    A PWGM_ELEMCOUNTS struct.
*/
#define PWGM_ECNT_Quad(ecs)    (ecs).count[PWGM_ELEMTYPE_QUAD]
/*! \brief Extract the Tri count from a PWGM_ELEMCOUNTS struct.
*
* \param ecs
*    A PWGM_ELEMCOUNTS struct.
*/
#define PWGM_ECNT_Tri(ecs)     (ecs).count[PWGM_ELEMTYPE_TRI]
/*! \brief Extract the Tet count from a PWGM_ELEMCOUNTS struct.
*
* \param ecs
*    A PWGM_ELEMCOUNTS struct.
*/
#define PWGM_ECNT_Tet(ecs)     (ecs).count[PWGM_ELEMTYPE_TET]
/*! \brief Extract the Wedge count from a PWGM_ELEMCOUNTS struct.
*
* \param ecs
*    A PWGM_ELEMCOUNTS struct.
*/
#define PWGM_ECNT_Wedge(ecs)   (ecs).count[PWGM_ELEMTYPE_WEDGE]
/*! \brief Extract the Pyramid count from a PWGM_ELEMCOUNTS struct.
*
* \param ecs
*    A PWGM_ELEMCOUNTS struct.
*/
#define PWGM_ECNT_Pyramid(ecs) (ecs).count[PWGM_ELEMTYPE_PYRAMID]


/*---------------------------------------------------------*/
/*! \brief XYZ component type ids.
*
* \sa PwVertXyzVal()
*/
typedef enum PWGM_ENUM_XYZ_e {
    PWGM_XYZ_X, /*!< X-component id */
    PWGM_XYZ_Y, /*!< Y-component id */
    PWGM_XYZ_Z, /*!< Z-component id */
    // add new enum values before this line

    /*! \cond sdkINTERNALS */
    /*! for Bookkeeping only - NOT valid "which" values
    @{ */
    PWGM_XYZ_SIZE,
    PWGM_XYZ_LAST = PWGM_XYZ_SIZE-1
    /* @} */
    /*! \endcond sdkINTERNALS */
}
PWGM_ENUM_XYZ;



/*---------------------------------------------------------*/
/*               STRUCTURED GRID DATATYPES                 */
/*---------------------------------------------------------*/

/*! Grid Type IDs */
typedef enum PWGM_ENUM_GRIDTYPE_e {
    PWGM_GRIDTYPE_STRUCTURED,   /*!< Structured grid */
    PWGM_GRIDTYPE_UNSTRUCTURED, /*!< Unstructured grid */

    /*! \cond sdkINTERNALS */
    /*! for Bookkeeping only - NOT valid "which" values
    @{ */
    PWGM_GRIDTYPE_SIZE,
    PWGM_GRIDTYPE_LAST = PWGM_GRIDTYPE_SIZE-1
    /* @} */
    /*! \endcond sdkINTERNALS */
}
PWGM_ENUM_GRIDTYPE;

/*---------------------------------------------------------*/
/*! \brief Structured grid ijk index data type.
*/
typedef struct PW_DLL_IMPEXP PWGM_INDEX3_t {
    PWP_INT32 i; /*!< i-coordinate used for 3D and 2D grids */
    PWP_INT32 j; /*!< j-coordinate used for 3D and 2D grids */
    PWP_INT32 k; /*!< k-coordinate used for 3D grids only */
}
PWGM_INDEX3;

/*---------------------------------------------------------*/
/*! \brief Structured grid ijk size data type.
*/
typedef PWGM_INDEX3 PW_DLL_IMPEXP PWGM_STR_SIZE;

/*---------------------------------------------------------*/
/*! \brief Structured grid ijk range data type.
*/
typedef struct PW_DLL_IMPEXP PWGM_STR_RANGE_t {
    PWGM_INDEX3 beg; /*!< begining index */
    PWGM_INDEX3 end; /*!< ending index */
}
PWGM_STR_RANGE;

/*---------------------------------------------------------*/
/*! \brief The 3D transform matrix data type.
*
* Defines the transform of a structured ijk index across a 3D block connection.
*
* You can obtain the equivalent 2D transform matrix using PwXform3to2().
*/
typedef struct PW_DLL_IMPEXP PWGM_INDEX_XFORM_t {
    PWP_INT32 m[3][4]; /*!< m[row][col] */
}
PWGM_INDEX_XFORM;

/*---------------------------------------------------------*/
/*! \brief The 2D transform matrix data type.
*
* Defines the transform of a structured ij index across a 2D block connection.
*
* You can convert to the equivalent 3D transform matrix using PwXform2to3().
*/
typedef struct PW_DLL_IMPEXP PWGM_INDEX_XFORM2_t {
    PWP_INT32 m[2][3]; /*!< m[row][col] */
}
PWGM_INDEX_XFORM2;

/*---------------------------------------------------------*/
/*! \brief Structured grid block-face ids.
*/
typedef enum PWGM_FACE_ID_e {
    PWGM_FACE_KMIN,  /*!< min K */
    PWGM_FACE_KMAX,  /*!< max K */
    PWGM_FACE_IMIN,  /*!< min I */
    PWGM_FACE_IMAX,  /*!< max I */
    PWGM_FACE_JMIN,  /*!< min J */
    PWGM_FACE_JMAX,  /*!< max J */
    PWGM_FACE_UNSTR, /*!< unstructured */
    // add new enum values before this line

    /*! \cond sdkINTERNALS */
    /*! for Bookkeeping only - NOT valid "which" values
    @{ */
    PWGM_FACE_SIZE,
    PWGM_FACE_LAST = PWGM_FACE_SIZE-1
    /* @} */
    /*! \endcond sdkINTERNALS */
}
PWGM_FACE_ID;

/*---------------------------------------------------------*/
/*! \brief Structured grid block-point association
*/
typedef enum PWGM_CNXNTYPE_e {
    PWGM_CNXNTYPE_ONE_TO_ONE,   /*!< one-to-one grid point match */
    PWGM_CNXNTYPE_MANY_TO_ONE,  /*!< one-to-many grid point match */
    PWGM_CNXNTYPE_MANY_TO_MANY, /*!< many-to-many grid point match */
    PWGM_CNXNTYPE_NONE,         /*!< unspecified */
    // add new enum values before this line

    /*! \cond sdkINTERNALS */
    /*! for Bookkeeping only - NOT valid "which" values
    @{ */
    PWGM_CNXNTYPE_SIZE,
    PWGM_CNXNTYPE_LAST = PWGM_FACE_SIZE-1
    /* @} */
    /*! \endcond sdkINTERNALS */
}
PWGM_CNXNTYPE;

/*---------------------------------------------------------*/
/*! \brief Structured grid boundary data type.
*/
typedef struct PW_DLL_IMPEXP PWGM_BNDRYDATA_t {
    const char     *name; /*!< boundary name */
    PWGM_HBLOCK    block; /*!< boundary block */
    PWGM_FACE_ID   face;  /*!< boundary face id */
    PWGM_STR_RANGE range; /*!< boundary ijk range */
}
PWGM_BNDRYDATA;

/*---------------------------------------------------------*/
/*! \brief Structured grid, inter-block, connection data type.
*/
typedef struct PW_DLL_IMPEXP PWGM_CNXNDATA_t {
    const char       *name;    /*!< connection name */
    PWGM_HBLOCK      block1;   /*!< connection block 1 (source) */
    PWGM_FACE_ID     face1;    /*!< face Id on block 1*/
    PWGM_STR_RANGE   range1;   /*!< ijk connection range in block1 index space */
    PWGM_INDEX_XFORM from1to2; /*!< transforms block1 index to block2 index */
    PWGM_HBLOCK      block2;   /*!< connection block 2 (donor) */
    PWGM_FACE_ID     face2;    /*!< face Id on block2 */
    PWGM_STR_RANGE   range2;   /*!< ijk connection range in block2 index space */
    PWGM_INDEX_XFORM from2to1; /*!< transforms block2 index to block1 index (== inverse(from1to2) */
}
PWGM_CNXNDATA;

/*---------------------------------------------------------*/
/*! \brief Block data type.
*/
typedef struct PW_DLL_IMPEXP PWGM_BLOCKDATA_t {
    const char         *name;    /*!< name */
    PWGM_HBLOCK        block;    /*!< handle */
    PWGM_ENUM_GRIDTYPE gridType; /*!< grid type */
    PWGM_STR_SIZE      size;     /*!< vertex-size */
}
PWGM_BLOCKDATA;

/*---------------------------------------------------------*/
/*! \brief IJK component type ids.
*/
typedef enum PWGM_ENUM_IJK_e {
    PWGM_IJK_I, /*!< I-component id */
    PWGM_IJK_J, /*!< J-component id */
    PWGM_IJK_K, /*!< K-component id */
    // add new enum values before this line

    /*! \cond sdkINTERNALS */
    /*! for Bookkeeping only - NOT valid "which" values
    @{ */
    PWGM_IJK_SIZE,
    PWGM_IJK_LAST = PWGM_IJK_SIZE-1
    /* @} */
    /*! \endcond sdkINTERNALS */
}
PWGM_ENUM_IJK;

/*!@} DOXGRP_APIPWGM_TYPES */


//*********************************************
//   PWGM-API functions
//*********************************************

/***********************************************************/
/***********************************************************/
/*! \defgroup DOXGRP_APIPWGM_FUNCTIONS Functions
*
*  The PWGM-API functions provide read-only access to a grid model.
*
*  Most of the PWGM-API functions are limited to accessing either structured
*  or unstructured grids. However, there are a few that can be used for both
*  grid model types.
@{
*/

/***********************************************************/
/***********************************************************/
/*! \defgroup DOXGRP_APIPWGM_FUNCTIONS_COM Common Grid Functions
*
*  These functions are used to access both structured and unstructured grid
*  models.
@{
*/

/*---------------------------------------------------------*/
/*! \brief Get the number of block elements in the model.
*
* \param model
*    The grid model handle.
*/
PW_DLLEXPORT PWP_UINT32
PwModBlockCount(PWGM_HGRIDMODEL model);

/*---------------------------------------------------------*/
/*! \brief Sequentially enumerate the model block elements
*
* \param model
*    The grid model handle.
*
* \param ndx
*    The block index starting at 0.
*
* \sa PwModBlockCount()
* \par "Sample usage:"
* \code
*    PWP_UINT32 eTotCnt;
*    PWGM_ELEMCOUNTS eCounts;
*    PWP_UINT32 ndx = 0;
*    PWGM_HBLOCK hBlk = PwModEnumBlocks(model, ndx);
*    while (PWGM_HBLOCK_ISVALID(hBlk)) {
*        eTotCnt = PwBlkElementCount(hBlk, &eCounts);
*        // ...etc...
*        hBlk = PwModEnumBlocks(model, ++ndx);
*    }
* \endcode
*/
PW_DLLEXPORT PWGM_HBLOCK
PwModEnumBlocks(PWGM_HGRIDMODEL model, PWP_UINT32 ndx);

/*!@} DOXGRP_APIPWGM_FUNCTIONS_COM */



/***********************************************************/
/***********************************************************/
/*! \defgroup DOXGRP_APIPWGM_FUNCTIONS_UNS Unstructured Grid Functions
*
*  These functions are ONLY used to access unstructured grid models.
@{
*/
#if !defined(PWGM_HIDE_UNSTRUCTURED_API)

/*---------------------------------------------------------*/
/*! \brief Get the number of domain elements in the model.
*
* \param model
*    The grid model handle.
*/
PW_DLLEXPORT PWP_UINT32
PwModDomainCount(PWGM_HGRIDMODEL model);

/*---------------------------------------------------------*/
/*! \brief Sequentially enumerate the model domain elements
*
* \param model
*    The grid model handle.
*
* \param ndx
*    The block index starting at 0.
*
* \sa PwModDomainCount()
* \par "Sample usage:"
* \code
*    PWP_UINT32 eTotCnt;
*    PWGM_ELEMCOUNTS eCounts;
*    PWP_UINT32 ndx = 0;
*    PWGM_HDOMAIN hDom = PwModEnumDomains(model, ndx);
*    while (PWGM_HDOMAIN_ISVALID(hDom)) {
*        eTotCnt = PwDomElementCount(hDom, &eCounts);
*        // ...etc...
*        hDom = PwModEnumDomains(model, ++ndx);
*    }
* \endcode
*
* \note
*   Domain elements are oriented to the interior of the block.
*
* \note
*   Connection elements are oriented to their source connection domain. This
*   is done so that blocks with local coordinate indices have a 1:1 point match
*   across block connections.
*
* \note
*   If a particular CAE plugin uses a different orientation model, the plugin
*   will need to manipulate the data during export.
*
* \par "For 3D Grids"
*   Elements are TRIS and QUADS.
*   The element normals (right-hand rule) point to the block interior.
*
* \par "For 2D Grids"
*   Elements are BARS.
*   The block interior is to the "left" while walking the block's "up" side.
*   Up is defined by the block's local normal (right-hand rule).
*
* \image html elem-orient.png "2D Element Orientation"
*/
PW_DLLEXPORT PWGM_HDOMAIN
PwModEnumDomains(PWGM_HGRIDMODEL model, PWP_UINT32 ndx);

/*---------------------------------------------------------*/
/*! \brief Sequentially enumerate the model vertex elements
*
* \param model
*    The grid model handle.
*
* \param ndx
*    The block index starting at 0.
*
* \sa PwModVertexCount(), PwVertDataMod(), PwVertIndexMod()
* \par "Sample usage:"
* \code
*    PWGM_VERTDATA VertData;
*    PWP_UINT32 ndx = 0;
*    PWGM_HVERTEX hVert = PwModEnumVertices(model, ndx);
*    while (PWGM_HVERTEX_ISVALID(hVert)) {
*        if (PwVertDataMod(hVert, &VertData)) {
*            // ...etc...
*        }
*        hVert = PwModEnumVertices(model, ++ndx);
*    }
* \endcode
*/
PW_DLLEXPORT PWGM_HVERTEX
PwModEnumVertices(PWGM_HGRIDMODEL model, PWP_UINT32 ndx);

/*---------------------------------------------------------*/
/*! \brief Get the number of model vertex elements
*
* \param model
*    The grid model handle.
*
* \return The total number of vertex elements.
*/
PW_DLLEXPORT PWP_UINT32
PwModVertexCount(PWGM_HGRIDMODEL model);

/*---------------------------------------------------------*/
/*! \brief Get the number of block elements
*
* \param block
*    A block handle.
*
* \param pCounts
*    Pointer to a PWGM_ELEMCOUNTS buffer.
*
*   \return The total number of elements
*/
PW_DLLEXPORT PWP_UINT32
PwBlkElementCount(PWGM_HBLOCK block, PWGM_ELEMCOUNTS *pCounts);

/*---------------------------------------------------------*/
/*! \brief Sequentially enumerate the block elements
*
* \param block
*    A block handle.
*
* \param ndx
*    The element index starting at 0.
*
*   \sa PwBlkElementCount(), PwElemDataMod()
*   \par "Sample usage:"
*   \code
*      PWGM_ELEMDATA ElemData;
*      PWP_UINT32 ndx = 0;
*      PWGM_HELEMENT hElem = PwBlkEnumElements(block, ndx);
*      while (PWGM_HELEMENT_ISVALID(hElem)) {
*          if (PwElemDataMod(hElem, &ElemData)) {
*              // ...etc...
*          }
*          hElem = PwBlkEnumElements(model, ++ndx);
*      }
*   \endcode
*/
PW_DLLEXPORT PWGM_HELEMENT
PwBlkEnumElements(PWGM_HBLOCK block, PWP_UINT32 ndx);

/*---------------------------------------------------------*/
/*! \brief Get the block condition data
*
* \param block
*    A block handle.
*
* \param pCondData
*    Pointer to a PWGM_CONDDATA buffer.
*
*   \return PWP_TRUE on success
*/
PW_DLLEXPORT PWP_BOOL
PwBlkCondition(PWGM_HBLOCK block, PWGM_CONDDATA *pCondData);

/*---------------------------------------------------------*/
/*! \brief Get the number of domain elements
*
* \param domain
*    A domain handle.
*
* \param pCounts
*    Pointer to a PWGM_ELEMCOUNTS buffer.
*
*   \return The total number of elements
*/
PW_DLLEXPORT PWP_UINT32
PwDomElementCount(PWGM_HDOMAIN domain, PWGM_ELEMCOUNTS *pCounts);

/*---------------------------------------------------------*/
/*! \brief Sequentially enumerate the domain elements
*
* \param domain
*    A domain handle.
*
* \param ndx
*    The element index starting at 0.
*
* \sa PwDomElementCount()
*/
PW_DLLEXPORT PWGM_HELEMENT
PwDomEnumElements(PWGM_HDOMAIN domain, PWP_UINT32 ndx);

/*---------------------------------------------------------*/
/*! \brief Get the domain condition data
*
* \param domain
*    A domain handle.
*
* \param pCondData
*    Pointer to a PWGM_CONDDATA buffer.
*
* \return PWP_TRUE on success
*/
PW_DLLEXPORT PWP_BOOL
PwDomCondition(PWGM_HDOMAIN domain, PWGM_CONDDATA *pCondData);

/*---------------------------------------------------------*/
/*! \brief Get the vertex data relative to the model's index space
*
* \param vertex
*    A vertex handle.
*
* \param pVertData
*    Pointer to a PWGM_VERTDATA buffer.
*
* \return PWP_TRUE on success
*
* \sa PwModEnumVertices()
*/
PW_DLLEXPORT PWP_BOOL
PwVertDataMod(PWGM_HVERTEX vertex, PWGM_VERTDATA *pVertData);

/*---------------------------------------------------------*/
/*! \brief Get the vertex index relative to the model's index space
*
* \param vertex
*    A vertex handle.
*
* \param pIndex
*    Pointer to a PWP_UINT32 value.
*
* \return PWP_TRUE on success
*
* \sa PwModEnumVertices()
*/
PW_DLLEXPORT PWP_BOOL
PwVertIndexMod(PWGM_HVERTEX vertex, PWP_UINT32 *pIndex);

/*---------------------------------------------------------*/
/*! \brief Get a vertex's x, y, or z component value
*
* \param vertex
*    A vertex handle.
*
* \param which
*    The XYZ component id to retrieve.
*
* \param pVal
*    Pointer to a PWGM_XYZVAL value.
*
* \return PWP_TRUE on success
*
* \sa PwModEnumVertices()
*/
PW_DLLEXPORT PWP_BOOL
PwVertXyzVal(PWGM_HVERTEX vertex, PWGM_ENUM_XYZ which, PWGM_XYZVAL *pVal);


/*---------------------------------------------------------*/
/*! \brief Get the element data relative to the model's index space
*
* \param element
*    A element handle.
*
* \param pElemData
*    Pointer to a PWGM_ELEMDATA buffer.
*
* \return PWP_TRUE on success
*/
PW_DLLEXPORT PWP_BOOL
PwElemDataMod(PWGM_HELEMENT element, PWGM_ELEMDATA *pElemData);

#endif // !defined(PWGM_HIDE_UNSTRUCTURED_API)
/*!@} DOXGRP_APIPWGM_FUNCTIONS_UNS */


/***********************************************************/
/***********************************************************/
/*! \defgroup DOXGRP_APIPWGM_FUNCTIONS_STR Structured Grid Functions
*
*  These functions are ONLY used to access structured grid models.
@{
*/
#if !defined(PWGM_HIDE_STRUCTURED_API)

/*---------------------------------------------------------*/
/*! \brief Get the block's vertex-size.
*
* \param block
*    A block handle.
*
* \param pSize
*    Pointer to a PWGM_STR_SIZE buffer.
*
* \return PWP_FALSE if block is not structured
*/
PW_DLLEXPORT PWP_BOOL
PwBlkSize(PWGM_HBLOCK block, PWGM_STR_SIZE *pSize);

/*---------------------------------------------------------*/
/*! \brief Get the block data.
*
* \param block
*    A block handle.
*
* \param pBlockData
*    Pointer to a PWGM_BLOCKDATA buffer.
*
* \return PWP_TRUE on success.
*/
PW_DLLEXPORT PWP_BOOL
PwBlock(PWGM_HBLOCK block, PWGM_BLOCKDATA *pBlockData);

/*---------------------------------------------------------*/
/*! \brief Get the block's vertex data at the given index location.
*
* \param block
*    A block handle.
*
* \param ndx3
*    The ijk block location.
*
* \param pVertData
*    Pointer to a PWGM_VERTDATA buffer.
*
* \return PWP_FALSE if block is not structured or ndx3 is invalid
*/
PW_DLLEXPORT PWP_BOOL
PwBlkNdxVertData(PWGM_HBLOCK block, PWGM_INDEX3 ndx3, PWGM_VERTDATA *pVertData);

/*---------------------------------------------------------*/
/*! \brief Get the number of structured grid connections in the model.
*
* \param model
*    The grid model handle.
*/
PW_DLLEXPORT PWP_UINT32
PwModConnectionCount(PWGM_HGRIDMODEL model);

/*---------------------------------------------------------*/
/*! \brief Sequentially enumerate the model's connections
*
* \param model
*    The grid model handle.
*
* \param ndx
*    The connection index starting at 0.
*
* \return  valid connection handle on success.
* \par "Sample usage:"
* \code
*    PWP_UINT32 cnxnCnt = PwModConnectionCount(model);
*    PWP_UINT32 ndx = 0;
*    PWGM_HCNXN connection = PwModEnumConnections(model, ndx);
*    while (PWGM_HCNXN_ISVALID(connection)) {
*        // ...etc...
*        connection = PwModEnumConnections(model, ++ndx);
*    }
* \endcode
*/
PW_DLLEXPORT PWGM_HCNXN
PwModEnumConnections(PWGM_HGRIDMODEL model, PWP_UINT32 ndx);

/*---------------------------------------------------------*/
/*! \brief Get the connection data.
*
* \param connection
*    A connection handle.
*
* \param pCnxnData
*    Pointer to a PWGM_CNXNDATA buffer.
*
* \return PWP_TRUE on success.
*/
PW_DLLEXPORT PWP_BOOL
PwConnection(PWGM_HCNXN connection, PWGM_CNXNDATA *pCnxnData);

/*---------------------------------------------------------*/
/*! \brief Get the data for the model's nth connection.
*
* \param model
*    The grid model handle.
*
* \param ndx
*    The model's connection index starting at 0.
*
* \param pCnxnData
*    Pointer to a PWGM_CNXNDATA buffer.
*
* \par "Equivalent to:"
* \code
*    PwConnection(PwModEnumConnections(model, ndx), &cnxnData);
* \endcode
*
* \return PWP_TRUE on success.
*/
PW_DLLEXPORT PWP_BOOL
PwModNdxConnection(PWGM_HGRIDMODEL model, PWP_UINT32 ndx, PWGM_CNXNDATA *pCnxnData);

/*---------------------------------------------------------*/
/*! \brief Get the number of structured grid boundaries in the model.
*
* \param model
*    The grid model handle.
*/
PW_DLLEXPORT PWP_UINT32
PwModBoundaryCount(PWGM_HGRIDMODEL model);

/*---------------------------------------------------------*/
/*! \brief Sequentially enumerate the model's boundaries
*
* \param model
*    The grid model handle.
*
* \param ndx
*    The boundary index starting at 0.
*
* \return  valid boundary handle on success.
* \par "Sample usage:"
* \code
*    PWP_UINT32 bndryCnt = PwModBoundaryCount(model);
*    PWP_UINT32 ndx = 0;
*    PWGM_HBNDRY boundary = PwModEnumBoundaries(model, ndx);
*    while (PWGM_HBNDRY_ISVALID(boundary)) {
*        // ...etc...
*        boundary = PwModEnumBoundaries(model, ++ndx);
*    }
* \endcode
*/
PW_DLLEXPORT PWGM_HBNDRY
PwModEnumBoundaries(PWGM_HGRIDMODEL model, PWP_UINT32 ndx);

/*---------------------------------------------------------*/
/*! \brief Get the boundary data.
*
* \param boundary
*    A boundary handle.
*
* \param pBndryData
*    Pointer to a PWGM_BNDRYDATA buffer.
*
* \return PWP_TRUE on success.
*/
PW_DLLEXPORT PWP_BOOL
PwBoundary(PWGM_HBNDRY boundary, PWGM_BNDRYDATA *pBndryData);

/*---------------------------------------------------------*/
/*! \brief Get the boundary's condition data.
*
* \param boundary
*    A boundary handle.
*
* \param pCondData
*    Pointer to a PWGM_CONDDATA buffer.
*
* \return PWP_TRUE on success.
*/
PW_DLLEXPORT PWP_BOOL
PwBndryCondition(PWGM_HBNDRY boundary, PWGM_CONDDATA *pCondData);

/*---------------------------------------------------------*/
/*! \brief Get the data for the model's nth structured boundary.
*
* \param model
*    The grid model handle.
*
* \param ndx
*    The model's boundary index starting at 0.
*
* \param pBndryData
*    Pointer to a PWGM_BNDRYDATA buffer.
*
* \par "Equivalent to:"
* \code
*    PwBoundary(PwModEnumBoundaries(model, ndx), &bndryData);
* \endcode
*
* \return PWP_TRUE on success.
*/
PW_DLLEXPORT PWP_BOOL
PwModNdxBoundary(PWGM_HGRIDMODEL model, PWP_UINT32 ndx, PWGM_BNDRYDATA *pBndryData);

/*---------------------------------------------------------*/
/*! \brief Get the PWGM_BNDRYDATA and PWGM_CONDDATA for the model's nth
*   structured boundary.
*
* \param model
*    The grid model handle.
*
* \param ndx
*    The model's boundary index starting at 0.
*
* \param pBndryData
*    Pointer to a PWGM_BNDRYDATA buffer or <tt>NULL</tt>.
*
* \param pCondData
*    Pointer to a PWGM_CONDDATA buffer or <tt>NULL</tt>.
*
* \par "Equivalent to:"
* \code
*    hBndry = PwModEnumBoundaries(model, ndx);
*    PwBoundary(hBndry, &bndryData) && PwBndryCondition(hBndry, &condData);
* \endcode
*
* \note
*    It is an error if both <tt>pBndryData</tt> and <tt>pCondData</tt> are <tt>NULL</tt>.
*
* \return PWP_TRUE on success.
*/
PW_DLLEXPORT PWP_BOOL
PwModNdxBoundaryAndCondition(PWGM_HGRIDMODEL model, PWP_UINT32 ndx,
    PWGM_BNDRYDATA *pBndryData, PWGM_CONDDATA *pCondData);

/*---------------------------------------------------------*/
/*! \brief Get the number of boundaries in the block.
*
* \param block
*    A block handle.
*/
PW_DLLEXPORT PWP_UINT32
PwBlkBoundaryCount(PWGM_HBLOCK block);

/*---------------------------------------------------------*/
/*! \brief Sequentially enumerate the block's boundaries
*
* \param block
*    A block handle.
*
* \param ndx
*    The boundary index starting at 0.
*
* \return  valid boundary handle on success.
*
* \par "Sample usage:"
* \code
*    PWP_UINT32 bndryCnt = PwBlkBoundaryCount(block);
*    PWP_UINT32 ndx = 0;
*    PWGM_HBNDRY boundary = PwBlkEnumBoundaries(block, ndx);
*    while (PWGM_HBNDRY_ISVALID(boundary)) {
*        // ...etc...
*        boundary = PwBlkEnumBoundaries(block, ++ndx);
*    }
* \endcode
*/
PW_DLLEXPORT PWGM_HBNDRY
PwBlkEnumBoundaries(PWGM_HBLOCK block, PWP_UINT32 ndx);

/*---------------------------------------------------------*/
/*! \brief Get the data for the block's nth structured boundary.
*
* \param block
*    A block handle.
*
* \param ndx
*    The block's boundary index starting at 0.
*
* \param pBndryData
*    Pointer to a PWGM_BNDRYDATA buffer.
*
* \par "Equivalent to:"
* \code
*    PwBoundary(PwBlkEnumBoundaries(block, ndx), &bndryData);
* \endcode
*
* \return PWP_TRUE on success.
*/
PW_DLLEXPORT PWP_BOOL
PwBlkNdxBoundary(PWGM_HBLOCK block, PWP_UINT32 ndx, PWGM_BNDRYDATA *pBndryData);

/*---------------------------------------------------------*/
/*! \brief Get the PWGM_BNDRYDATA and PWGM_CONDDATA for the block's nth
*          structured boundary.
*
* \param block
*    A block handle.
*
* \param ndx
*    The block's boundary index starting at 0.
*
* \param pBndryData
*    Pointer to a PWGM_BNDRYDATA buffer or <tt>NULL</tt>.
*
* \param pCondData
*    Pointer to a PWGM_CONDDATA buffer or <tt>NULL</tt>.
*
* \par "Equivalent to:"
* \code
*    hBndry = PwBlkEnumBoundaries(block, ndx);
*    PwBoundary(hBndry, &bndryData) && PwBndryCondition(hBndry, &condData);
* \endcode
*
* \note
*    It is an error if both <tt>pBndryData</tt> and <tt>pCondData</tt> are <tt>NULL</tt>.
*
* \return PWP_TRUE on success.
*/
PW_DLLEXPORT PWP_BOOL
PwBlkNdxBoundaryAndCondition(PWGM_HBLOCK block, PWP_UINT32 ndx,
    PWGM_BNDRYDATA *pBndryData, PWGM_CONDDATA *pCondData);

/*---------------------------------------------------------*/
/*! \brief Get the number of connections in the block.
*
* \param block
*    A block handle.
*/
PW_DLLEXPORT PWP_UINT32
PwBlkConnectionCount(PWGM_HBLOCK block);

/*---------------------------------------------------------*/
/*! \brief Sequentially enumerate the block's connections
*
* \param block
*    A block handle.
*
* \param ndx
*    The connection index starting at 0.
*
* \return  valid connection handle on success.
*
* \par "Sample usage:"
* \code
*    PWP_UINT32 cnxnCnt = PwBlkConnectionCount(block);
*    PWP_UINT32 ndx = 0;
*    PWGM_HCNXN connection = PwBlkEnumConnections(block, ndx);
*    while (PWGM_HCNXN_ISVALID(connection)) {
*        // ...etc...
*        connection = PwBlkEnumBoundaries(block, ++ndx);
*    }
* \endcode
*/
PW_DLLEXPORT PWGM_HCNXN
PwBlkEnumConnections(PWGM_HBLOCK block, PWP_UINT32 ndx);

/*---------------------------------------------------------*/
/*! \brief Get the data for the block's nth connection.
*
* \param block
*    A block handle.
*
* \param ndx
*    The block's connection index starting at 0.
*
* \param pCnxnData
*    Pointer to a PWGM_CNXNDATA buffer.
*
* \par "Equivalent to:"
* \code
*    PwConnection(PwBlkEnumConnections(block, ndx), &cnxnData);
* \endcode
*
* \return PWP_TRUE on success.
*/
PW_DLLEXPORT PWP_BOOL
PwBlkNdxConnection(PWGM_HBLOCK block, PWP_UINT32 ndx,
    PWGM_CNXNDATA *pCnxnData);

/***********************************************************/
/***********************************************************/
/*! \defgroup DOXGRP_APIPWGM_FUNCTIONS_STRUTILS Structured Grid Utility Functions
*
*  These utility functions are ONLY used to access structured grid models.
*  These functions are defined locally and do NOT get routed to the framework.
@{
*/
/*---------------------------------------------------------*/
/*! \brief Convert a 2D transform matrix to it's 3D equivalent.
*
* \param pX2
*    Const pointer to a PWGM_INDEX_XFORM2 input buffer.
*
* \param pX3
*    Pointer to a PWGM_INDEX_XFORM output buffer.
*
* \par "The matrix mapping:"
* \code
*    XFORM2       XFORM
*    ------       -------
*    A B C        A B 0 C
*    D E F   ==>  D E 0 F
*                 0 0 1 0
* \endcode
*
* \return PWP_TRUE on success.
*/
PW_DLLEXPORT PWP_BOOL
PwXform2to3(const PWGM_INDEX_XFORM2 *pX2, PWGM_INDEX_XFORM *pX3);

/*---------------------------------------------------------*/
/*! \brief Convert a 3D transform matrix to it's 2D equivalent.
*
* \param pX3
*    Const pointer to a PWGM_INDEX_XFORM input buffer.
*
* \param pX2
*    Pointer to a PWGM_INDEX_XFORM2 output buffer.
*
* \par "The matrix mapping:"
* \code
*    XFORM         XFORM2
*    -------       ------
*    A B C D       A B D
*    E F G H  ==>  E F H
*    I J K L       
* \endcode
*
* \return PWP_TRUE on success.
*/
PW_DLLEXPORT PWP_BOOL
PwXform3to2(const PWGM_INDEX_XFORM *pX3, PWGM_INDEX_XFORM2 *pX2);

/*---------------------------------------------------------*/
/*! \brief Apply a PWGM_INDEX_XFORM transform to a PWGM_INDEX3 value.
*
* \param pX3
*    Const pointer to the PWGM_INDEX_XFORM.
*
* \param ijk
*    The PWGM_INDEX3 value to transform. It will not be changed.
*
* \return The transformed PWGM_INDEX3 value.
*/
PWGM_INDEX3
PwXformApply(const PWGM_INDEX_XFORM *pX3, PWGM_INDEX3 ijk);

/*---------------------------------------------------------*/
/*! \brief For a given localAxis, determine the corresponding axis in
*   the transformed system.
*
* \param pX3
*    Const pointer to the PWGM_INDEX_XFORM.
*
* \param localAxis
*    The local axis of interest.
*
* \param pFlipped
*    Flag set to PWP_TRUE if the returned axis direction is the opposite of
*    localAxis. Pointer may be NULL.
*
* \return The corresponding axis in the transformed system.
*/
PWGM_ENUM_IJK
PwXformFollows(const PWGM_INDEX_XFORM *pX3, PWGM_ENUM_IJK localAxis,
               PWP_BOOL *pFlipped);

/*---------------------------------------------------------*/
/*! \brief Apply a PWGM_INDEX_XFORM2 transform to a PWGM_INDEX3 value.
*
* \param pX2
*    Const pointer to the PWGM_INDEX_XFORM2.
*
* \param ijk
*    The PWGM_INDEX3 value to transform. It will not be changed.
*
* \return The transformed PWGM_INDEX3 value. Only the i and j values will be
*    transformed. The k value will not be changed.
*/
PWGM_INDEX3
PwXform2Apply(const PWGM_INDEX_XFORM2 *pX2, PWGM_INDEX3 ijk);

/*---------------------------------------------------------*/
/*! \brief For a given localAxis, determine the corresponding axis in
*   the transformed system.
*
* \param pX2
*    Const pointer to the PWGM_INDEX_XFORM2.
*
* \param localAxis
*    The local axis of interest.
*
* \param pFlipped
*    Flag set to PWP_TRUE if the returned axis direction is the opposite of
*    localAxis. Pointer may be NULL.
*
* \return The corresponding axis in the transformed system.
*/
PWGM_ENUM_IJK
PwXform2Follows(const PWGM_INDEX_XFORM2 *pX2, PWGM_ENUM_IJK localAxis,
                PWP_BOOL *pFlipped);

/*---------------------------------------------------------*/
/*! \brief Determines if an PWGM_INDEX3 is within a PWGM_STR_RANGE.
*
* \param ijk
*    The PWGM_INDEX3 value to check.
*
* \param pRange
*    Const pointer to the PWGM_STR_RANGE.
*
* \return PWP_TRUE if ijk resides within the limits of pRange.
*/
PWP_BOOL
PwInRange(PWGM_INDEX3 ijk, const PWGM_STR_RANGE *pRange);

/*!@} DOXGRP_APIPWGM_FUNCTIONS_STRUTILS */

#endif // !defined(PWGM_HIDE_STRUCTURED_API)
/*!@} DOXGRP_APIPWGM_FUNCTIONS_STR */



/*! \cond sdkINTERNALS */
/*---------------------------------------------------------*/
/*---------------------------------------------------------*/
/*! PWGM-API function-pointer types (extern C)
*   \sa PWGM_HGRIDMODEL_IMPL
@{*/

/*---------------------------------------------------------*/
/*                  COMMON GRID CALLBACKS                  */
/*---------------------------------------------------------*/

typedef PWP_UINT32
PwModBlockCount_t(PWGM_HGRIDMODEL model);

typedef PWGM_HBLOCK
PwModEnumBlocks_t(PWGM_HGRIDMODEL model, PWP_UINT32 ndx);


/*---------------------------------------------------------*/
/*               UNSTRUCTURED GRID CALLBACKS               */
/*---------------------------------------------------------*/

typedef PWP_UINT32
PwModDomainCount_t(PWGM_HGRIDMODEL model);

typedef PWGM_HDOMAIN
PwModEnumDomains_t(PWGM_HGRIDMODEL model, PWP_UINT32 ndx);

typedef PWGM_HVERTEX
PwModEnumVertices_t(PWGM_HGRIDMODEL model, PWP_UINT32 ndx);

typedef PWP_UINT32
PwModVertexCount_t(PWGM_HGRIDMODEL model);

//---------------------------
typedef PWP_UINT32
PwBlkElementCount_t(PWGM_HBLOCK block, PWGM_ELEMCOUNTS *pCounts);

typedef PWGM_HELEMENT
PwBlkEnumElements_t(PWGM_HBLOCK block, PWP_UINT32 ndx);

typedef PWP_BOOL
PwBlkCondition_t(PWGM_HBLOCK block, PWGM_CONDDATA *pCondData);

//---------------------------
typedef PWP_UINT32
PwDomElementCount_t(PWGM_HDOMAIN domain, PWGM_ELEMCOUNTS *pCounts);

typedef PWGM_HELEMENT
PwDomEnumElements_t(PWGM_HDOMAIN domain, PWP_UINT32 ndx);

typedef PWP_BOOL
PwDomCondition_t(PWGM_HDOMAIN domain, PWGM_CONDDATA *pCondData);

//---------------------------
typedef PWP_BOOL
PwVertDataMod_t(PWGM_HVERTEX vertex, PWGM_VERTDATA *pVertData);

typedef PWP_BOOL
PwVertIndexMod_t(PWGM_HVERTEX vertex, PWP_UINT32 *pIndex);

typedef PWP_BOOL
PwVertXyzVal_t(PWGM_HVERTEX vertex, PWGM_ENUM_XYZ which, PWGM_XYZVAL *pVal);

//---------------------------
typedef PWP_BOOL
PwElemDataMod_t(PWGM_HELEMENT element, PWGM_ELEMDATA *pElemData);



/*---------------------------------------------------------*/
/*               STRUCTURED GRID CALLBACKS                 */
/*---------------------------------------------------------*/

typedef PWP_UINT32
PwModBoundaryCount_t(PWGM_HGRIDMODEL model);

typedef PWGM_HBNDRY
PwModEnumBoundaries_t(PWGM_HGRIDMODEL model, PWP_UINT32 ndx);

typedef PWP_UINT32
PwModConnectionCount_t(PWGM_HGRIDMODEL model);

typedef PWGM_HCNXN
PwModEnumConnections_t(PWGM_HGRIDMODEL model, PWP_UINT32 ndx);

typedef PWP_BOOL
PwModNdxBoundary_t(PWGM_HGRIDMODEL model, PWP_UINT32 ndx,
    PWGM_BNDRYDATA *pBndryData);

typedef PWP_BOOL
PwModNdxBoundaryAndCondition_t(PWGM_HGRIDMODEL model, PWP_UINT32 ndx,
    PWGM_BNDRYDATA *pBndryData, PWGM_CONDDATA *pCondData);

typedef PWP_BOOL
PwModNdxConnection_t(PWGM_HGRIDMODEL model, PWP_UINT32 ndx,
    PWGM_CNXNDATA *pCnxnData);

//---------------------------
typedef PWP_BOOL
PwBlkSize_t(PWGM_HBLOCK block, PWGM_STR_SIZE *pSize);

typedef PWP_BOOL
PwBlkNdxVertData_t(PWGM_HBLOCK block, PWGM_INDEX3 ndx3,
    PWGM_VERTDATA *pVertData);

typedef PWP_UINT32
PwBlkBoundaryCount_t(PWGM_HBLOCK block);

typedef PWGM_HBNDRY
PwBlkEnumBoundaries_t(PWGM_HBLOCK block, PWP_UINT32 ndx);

typedef PWP_UINT32
PwBlkConnectionCount_t(PWGM_HBLOCK block);

typedef PWGM_HCNXN
PwBlkEnumConnections_t(PWGM_HBLOCK block, PWP_UINT32 ndx);

typedef PWP_BOOL
PwBlkNdxBoundary_t(PWGM_HBLOCK block, PWP_UINT32 ndx,
    PWGM_BNDRYDATA *pBndryData);

typedef PWP_BOOL
PwBlkNdxBoundaryAndCondition_t(PWGM_HBLOCK block, PWP_UINT32 ndx,
    PWGM_BNDRYDATA *pBndryData, PWGM_CONDDATA *pCondData);

typedef PWP_BOOL
PwBlkNdxConnection_t(PWGM_HBLOCK block, PWP_UINT32 ndx,
    PWGM_CNXNDATA *pCnxnData);


//---------------------------
typedef PWP_BOOL
PwBoundary_t(PWGM_HBNDRY boundary, PWGM_BNDRYDATA *pBndryData);

typedef PWP_BOOL
PwBndryCondition_t(PWGM_HBNDRY boundary, PWGM_CONDDATA *pCondData);

typedef PWP_BOOL
PwConnection_t(PWGM_HCNXN connection, PWGM_CNXNDATA *pCnxnData);

typedef PWP_BOOL
PwBlock_t(PWGM_HBLOCK block, PWGM_BLOCKDATA *pBlockData);

//---------------------------
typedef PWP_BOOL
PwNOTIMPLEMENTED_t();

/*@}*/


/*---------------------------------------------------------*/
/*! PWGM_HGRIDMODEL implementation
*   \note PWGM_HGRIDMODEL == PWGM_HGRIDMODEL_IMPL*
*/
typedef struct PW_DLL_IMPEXP PWGM_HGRIDMODEL_IMPL_t {
    PwModBlockCount_t *PwModBlockCountCB;
    PwModDomainCount_t *PwModDomainCountCB;
    PwModEnumBlocks_t *PwModEnumBlocksCB;
    PwModEnumDomains_t *PwModEnumDomainsCB;
    PwModEnumVertices_t *PwModEnumVerticesCB;
    PwModVertexCount_t *PwModVertexCountCB;

    PwBlkElementCount_t *PwBlkElementCountCB;
    PwBlkEnumElements_t *PwBlkEnumElementsCB;
    PwNOTIMPLEMENTED_t *PwBlkEnumVerticesCB;
    PwNOTIMPLEMENTED_t *PwBlkVertexCountCB;
    PwBlkCondition_t *PwBlkConditionCB;

    PwDomElementCount_t *PwDomElementCountCB;
    PwDomEnumElements_t *PwDomEnumElementsCB;
    PwNOTIMPLEMENTED_t *PwDomEnumVerticesCB;
    PwNOTIMPLEMENTED_t *PwDomVertexCountCB;
    PwDomCondition_t *PwDomConditionCB;

    PwVertDataMod_t *PwVertDataModCB;
    PwNOTIMPLEMENTED_t *PwVertDataBlkCB;
    PwNOTIMPLEMENTED_t *PwVertDataDomCB;
    PwVertIndexMod_t *PwVertIndexModCB;
    PwNOTIMPLEMENTED_t *PwVertIndexBlkCB;
    PwNOTIMPLEMENTED_t *PwVertIndexDomCB;
    PwVertXyzVal_t *PwVertXyzValCB;

    PwElemDataMod_t *PwElemDataModCB;
    PwNOTIMPLEMENTED_t *PwElemDataBlkCB;
    PwNOTIMPLEMENTED_t *PwElemDataDomCB;

    /*---------------------------------------------------------*/
    /*               STRUCTURED GRID ADDITIONS                 */
    /*---------------------------------------------------------*/

    PwModBoundaryCount_t *PwModBoundaryCountCB;
    PwModEnumBoundaries_t *PwModEnumBoundariesCB;
    PwModConnectionCount_t *PwModConnectionCountCB;
    PwModEnumConnections_t *PwModEnumConnectionsCB;
    PwModNdxBoundary_t *PwModNdxBoundaryCB;
    PwModNdxBoundaryAndCondition_t *PwModNdxBoundaryAndConditionCB;
    PwModNdxConnection_t *PwModNdxConnectionCB;

    PwBlkSize_t *PwBlkSizeCB;
    PwBlkNdxVertData_t *PwBlkNdxVertDataCB;
    PwBlkBoundaryCount_t *PwBlkBoundaryCountCB;
    PwBlkEnumBoundaries_t *PwBlkEnumBoundariesCB;
    PwBlkConnectionCount_t *PwBlkConnectionCountCB;
    PwBlkEnumConnections_t *PwBlkEnumConnectionsCB;
    PwBlkNdxBoundary_t *PwBlkNdxBoundaryCB;
    PwBlkNdxBoundaryAndCondition_t *PwBlkNdxBoundaryAndConditionCB;
    PwBlkNdxConnection_t *PwBlkNdxConnectionCB;

    PwBoundary_t *PwBoundaryCB;
    PwBndryCondition_t *PwBndryConditionCB;
    PwConnection_t *PwConnectionCB;
    PwBlock_t *PwBlockCB;
}
PWGM_HGRIDMODEL_IMPL;

/*! \endcond sdkINTERNALS */

/*!@} DOXGRP_APIPWGM_FUNCTIONS */

/*!@} DOXGRP_APIPWGM */


#ifdef __cplusplus
} // extern "C"
#endif

#endif /* !_APIGRIDMODEL_H_ */
