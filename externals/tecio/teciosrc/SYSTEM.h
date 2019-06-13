 #ifndef TECPLOT_SYSTEM_H
 #define TECPLOT_SYSTEM_H
 #if defined EXTERN
 #undef EXTERN
 #endif
 #if defined ___3921
 #define EXTERN
 #else
 #define EXTERN extern
 #endif
 #if defined ALLOW_INTERRUPTS_DURING_SLOW_READS_WRITES
 #define MAX_BYTES_PER_CHUNK 4096L 
 #else
 #define MAX_BYTES_PER_CHUNK 131072L 
 #endif
 #if defined TECPLOTKERNEL
size_t ___487(void const* ___3644, size_t      ___2101, size_t      ___2813, FILE*       ___3948, size_t      maxBytesPerChunk = MAX_BYTES_PER_CHUNK);
 #endif
EXTERN ___372 ___4406(const char *___1439); EXTERN ___372 ___2069(const char *___1439); EXTERN ___372 ___1389(const char *F, ___372   ___3570); EXTERN void ___1176(const char *___1439); EXTERN ___372 ___2070(const char *___1395, ___372   ___2054, ___372   ___3572); EXTERN ___372 ___3357(FILE   *File, int64_t ___2224); EXTERN ___372 ___505(FILE     **F, ___372  ___3570); EXTERN ___372 ___2877(FILE       **F, const char *___1439, ___372  ___2054, ___372  ___2001, ___372  ___1473, ___372  ___3570, ___372  ___2003);
 #endif
