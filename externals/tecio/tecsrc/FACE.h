/*
******************************************************************
******************************************************************
*******                                                   ********
******  (C) 1988-2010 Tecplot, Inc.                        *******
*******                                                   ********
******************************************************************
******************************************************************
*/
#ifndef _FACE_H_
#define _FACE_H_

#if defined EXTERN
#undef EXTERN
#endif
#if defined FACEMODULE
#define EXTERN
#else
#define EXTERN extern
#endif

namespace tecplot
{
namespace kernel
{
class SubElemValueProducerInterface;
}
}


/**
 */
inline Boolean_t IsCellFaceLogicallyCollapsed(LgIndex_t I1,
                                              LgIndex_t I2,
                                              LgIndex_t I3,
                                              LgIndex_t I4)
{
    return ((I1 == I2 && I3 == I4) ||
            (I1 == I4 && I2 == I3) ||
            (I1 == I3)             ||
            (I2 == I4));
}

/**
 * IMPORTANT NOTE:
 *   A face obscuration of FaceObscuration_LogicallyObscured means that the
 *   face is entirely obscured by either an implicit neighbor for inside faces
 *   of ordered data or an auto generated neighbor for finite element data. In
 *   either case, logical obscuration is not considered if user defined
 *   neighbors have been specified for the face. Therefore, interior faces of
 *   ordered data can have an indication of FaceObscuration_PartiallyObscured.
 */
typedef enum
{
    FaceObscuration_NotObscured,
    FaceObscuration_PartiallyObscured,
    FaceObscuration_EntirelyObscured,
    FaceObscuration_LogicallyObscured,
    END_FaceObscuration_e,
    FaceObscuration_Invalid = BadEnumValue
} FaceObscuration_e;

/**
 */
EXTERN LgIndex_t GetLogicalOrderedNeighbor(LgIndex_t NumIPts,
                                           LgIndex_t NumJPts,
                                           LgIndex_t NumKPts,
                                           LgIndex_t Element,
                                           LgIndex_t Face);

/**
 * Function to determine a cell's neighbor.  It calls FaceNeighborGetSurfaceCellNeighbor()
 * for classic zones.
 */
EXTERN void GetSurfaceCellNeighbor(CZInfo_s const*                                 CZInfo,
                                   CZData_s const*                                 CZData,
                                   LgIndex_t                                       SurfaceCellIndex,
                                   tecplot::kernel::SubElemValueProducerInterface* NodeValueProducer,
                                   ElemFaceOffset_t                                PlaneOrFaceOffset,
                                   ElemFaceOffset_t                                Edge,
                                   LgIndex_t*                                      NeighborSurfaceCellElem,
                                   EntIndex_t*                                     NeighborSurfaceCellZone);
/**
 */
EXTERN FaceObscuration_e GetFaceObscuration(CZInfo_s const* CZInfo,
                                            CZData_s const* CZData,
                                            Set_pa          ActiveRelevantZones,
                                            LgIndex_t       Element,
                                            LgIndex_t       FOffset,
                                            Boolean_t       ConsiderValueBlanking,
                                            Boolean_t       ConsiderIJKBlanking,
                                            Boolean_t       ConsiderDepthBlanking,
                                            Boolean_t       ConsiderClipBlanking);

EXTERN EntIndex_t GetNodesPerElementFace(ZoneType_e ZoneType);

EXTERN EntIndex_t GetFacesPerElement(ZoneType_e ZoneType,
                                     LgIndex_t  IMax,
                                     LgIndex_t  JMax,
                                     LgIndex_t  KMax);

/**
 * @param Points must have MAX_NODES_PER_FACE elements.
 * @param UniquePoints must have MAX_NODES_PER_FACE elements.  
 *   This array is populated with the unique points
 *   from the input Points array.  Values beyond NumUniquePoints-1 may
 *   be uninitialized and should not be accessed.
 * @param NumUniquePoints - This value is populated the number of unique 
 *   points that were placed in the UniquePoints array.
 * @param Is3D is used to help determine if two points are coincident.
 */
EXTERN void GetUniquePoints(XYZ_s const Points[],
                            XYZ_s UniquePoints[],
                            SmInteger_t& NumUniquePoints,
                            Boolean_t Is3D);
EXTERN CollapsedStatus_e GetSurfaceCellCollapsedStatus(CZInfo_s const*                                 CZInfo,
                                                       CZData_s const*                                 CZData,
                                                       tecplot::kernel::SubElemValueProducerInterface* SubElemValueProducer);
EXTERN CollapsedStatus_e GetSurfaceCellCollapsedStatus(CZInfo_s const* CZInfo,
                                                       CZData_s const* CZData,
                                                       LgIndex_t       I1,
                                                       LgIndex_t       I2,
                                                       LgIndex_t       I3,
                                                       LgIndex_t       I4);
EXTERN CollapsedStatus_e GetSurfaceCellLogicalCollapsedStatus(ZoneType_e ZoneType,
                                                              LgIndex_t  I1,
                                                              LgIndex_t  I2,
                                                              LgIndex_t  I3,
                                                              LgIndex_t  I4);
EXTERN CollapsedStatus_e GetSurfEdgeOrVolFaceLogicalCollapsedStatus(NodeMap_pa NodeMap,
                                                                    LgIndex_t  Element,
                                                                    EntIndex_t Face);
#if defined ALLOW_USERDEF_NO_NEIGHBORING_ELEMENT
/**
 */
EXTERN Boolean_t IsUserDefFaceNeighborBoundary(FaceNeighbor_pa FaceNeighbor,
                                               LgIndex_t       Element,
                                               LgIndex_t       Face);
#endif

#endif
