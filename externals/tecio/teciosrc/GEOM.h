 #if defined EXTERN
 #undef EXTERN
 #endif
 #if defined ___1618
 #define EXTERN
 #else
 #define EXTERN extern
 #endif
 #define VALID_RECTANGLE_COORDSYS(___3920) \
 (((___3920)==CoordSys_Frame) || \
 ((___3920)==CoordSys_Grid))
 #define VALID_SQUARE_COORDSYS(___3920) VALID_RECTANGLE_COORDSYS((___3920))
 #define VALID_ELLIPSE_COORDSYS(___3920) VALID_RECTANGLE_COORDSYS((___3920))
 #define VALID_CIRCLE_COORDSYS(___3920) VALID_ELLIPSE_COORDSYS((___3920))
 #define VALID_IMAGE_COORDSYS(___3920) (VALID_RECTANGLE_COORDSYS((___3920) || ___3920 == CoordSys_Grid3D))
 #define VALID_LINESEG_COORDSYS(___3920) \
 (((___3920)==CoordSys_Frame) || \
 ((___3920)==CoordSys_Grid)  || \
 ((___3920)==CoordSys_Grid3D))
 #define VALID_GEOM_COORDSYS(___3920) \
 (((___3920)==CoordSys_Frame) || \
 ((___3920)==CoordSys_Grid)  || \
 ((___3920)==CoordSys_Grid3D))
 #define GEOM_USES_V3(___1555) (___1555->___3167 == CoordSys_Grid3D && ___1555->___1652 == GeomType_LineSegs)
 #define VALID_GEOM_TYPE(___1650) \
 ( VALID_ENUM((___1650),GeomType_e) && \
 (___1650)!=GeomType_LineSegs3D )
 #define VALID_GEOM_FIELD_DATA_TYPE(___906) \
 ( ( (___906) == FieldDataType_Float ) || \
 ( (___906) == FieldDataType_Double ) )
