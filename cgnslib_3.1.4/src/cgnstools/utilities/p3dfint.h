
#ifndef _P3DFINT_H_
#define _P3DFINT_H_

#ifdef _WIN32

void __stdcall OPENF   (int *,char *,int);
void __stdcall CLOSEF  (void);
void __stdcall READIF  (int *,int *,int *);
void __stdcall READFF  (int *,float *,int *);
void __stdcall READDF  (int *,double *,int *);
void __stdcall READGFF (int *,float *,int *,int *);
void __stdcall READGDF (int *,double *,int *,int *);
void __stdcall WRITEIF (int *,int *,int *);
void __stdcall WRITEFF (int *,float *,int *);
void __stdcall WRITEDF (int *,double *,int *);
void __stdcall WRITEGFF(int *,float *,int *,int *);
void __stdcall WRITEGDF(int *,double *,int *,int *);

#else

#if defined(__linux) || defined(__CYGWIN__)
# ifndef UNDERSCORE
#  define UNDERSCORE
# endif
#endif

#if !defined(CRAY) && !defined(VMS)
# ifdef UNDERSCORE
#  define OPENF    openf_
#  define CLOSEF   closef_
#  define READIF   readif_
#  define READFF   readff_
#  define READDF   readdf_
#  define READGFF  readgff_
#  define READGDF  readgdf_
#  define WRITEIF  writeif_
#  define WRITEFF  writeff_
#  define WRITEDF  writedf_
#  define WRITEGFF writegff_
#  define WRITEGDF writegdf_
# else
#  define OPENF    openf
#  define CLOSEF   closef
#  define READIF   readif
#  define READFF   readff
#  define READDF   readdf
#  define READGFF  readgff
#  define READGDF  readgdf
#  define WRITEIF  writeif
#  define WRITEFF  writeff
#  define WRITEDF  writedf
#  define WRITEGFF writegff
#  define WRITEGDF writegdf
# endif
#endif

#endif
#endif

