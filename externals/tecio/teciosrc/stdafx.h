 #if !defined STDAFX_H_
 # define STDAFX_H_
#include <locale>
 #if defined _WIN32
 #if !defined _DEBUG
 #if !defined NDEBUG
 #define NDEBUG
 #endif
 #endif
 #if !defined MSWIN
 #define MSWIN
 #endif
 #define ___1169 
 #if !defined WINVER
 #define WINVER 0x0501  
 #endif
 #define ___22(___3251,___349) ((___3251)!=NULL)
 #if defined _M_IA64 || defined _M_AMD64
 #define WININT  INT_PTR
 #define WINUINT UINT_PTR
 #else
 #define WININT  int
 #define WINUINT UINT
 #endif
 #define WINCALLBACK ___438
 #endif 
 #endif 
