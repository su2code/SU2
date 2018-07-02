/*
*****************************************************************
*****************************************************************
*******                                                  ********
****** Copyright (C) 1988-2010 Tecplot, Inc.              *******
*******                                                  ********
*****************************************************************
*****************************************************************
*/

#if defined EXTERN
#undef EXTERN
#endif
#if defined GEOMMODULE
#define EXTERN
#else
#define EXTERN extern
#endif


/* * macros for checking CoordSys_e * */
#define VALID_RECTANGLE_COORDSYS(sys) \
          (((sys)==CoordSys_Frame) || \
           ((sys)==CoordSys_Grid))
#define VALID_SQUARE_COORDSYS(sys) VALID_RECTANGLE_COORDSYS((sys))
#define VALID_ELLIPSE_COORDSYS(sys) VALID_RECTANGLE_COORDSYS((sys))
#define VALID_CIRCLE_COORDSYS(sys) VALID_ELLIPSE_COORDSYS((sys))
#define VALID_IMAGE_COORDSYS(sys) VALID_RECTANGLE_COORDSYS((sys))
#define VALID_LINESEG_COORDSYS(sys) \
          (((sys)==CoordSys_Frame) || \
           ((sys)==CoordSys_Grid)  || \
           ((sys)==CoordSys_Grid3D))
#define VALID_GEOM_COORDSYS(sys) \
          (((sys)==CoordSys_Frame) || \
           ((sys)==CoordSys_Grid)  || \
           ((sys)==CoordSys_Grid3D))

#define VALID_GEOM_TYPE(geomtype) \
          ( VALID_ENUM((geomtype),GeomType_e) && \
            (geomtype)!=GeomType_LineSegs3D )

#define VALID_GEOM_FIELD_DATA_TYPE(datatype) \
          ( ( (datatype) == FieldDataType_Float ) || \
            ( (datatype) == FieldDataType_Double ) )

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif /* TECPLOTKERNEL */
