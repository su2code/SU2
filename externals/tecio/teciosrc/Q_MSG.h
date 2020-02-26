 #ifndef Q_MSG_H
 #define Q_MSG_H
 #if defined EXTERN
 #undef EXTERN
 #endif
 #if defined ___3259
 #define EXTERN
 #else
 #define EXTERN extern
 #endif
 #define MAX_STATUS_LINE_MSG_LEN 255
 #define MAX_RUNNING_COORDS_TEXT_LEN 80
#include "TranslatedString.h"
EXTERN ___372 ___4478(const char  *___2873, char       **___2700); EXTERN void Information(tecplot::___4218 format, ...); EXTERN void ___4447(tecplot::___4218 format, ...); EXTERN void ___1177(tecplot::___4218 format, ...);
 #endif 
