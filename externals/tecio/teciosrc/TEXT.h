 #pragma once
#include "ThirdPartyHeadersBegin.h"
#include <map>
#include <utility>
#include "ThirdPartyHeadersEnd.h"
 #define _TEXT_H_INCLUDED
 #define VALID_TEXT_COORDSYS(___3920)  (((___3920)==CoordSys_Frame)||((___3920)==CoordSys_Grid)||((___3920)==CoordSys_Grid3D))
 #define VALID_TEXT_UNITS(___4266)  (((___4266)==___4269)||((___4266)==___4268)||((___4266)==___4271))
 #define VALID_TEXT_COORDSYS_AND_UNITS(___3171, ___3602) \
 ( VALID_TEXT_COORDSYS((___3171)) && \
 VALID_TEXT_UNITS((___3602)) && \
 ! ((___3171) == CoordSys_Frame && (___3602) == ___4269) )
 #define VALID_FONT_SIZEUNITS(___4266)  (((___4266)==___4269)||((___4266)==___4268)||((___4266)==___4271)||(___4266)==___4267)
