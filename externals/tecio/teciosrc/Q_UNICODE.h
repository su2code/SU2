 #if !defined Q_UNICODE_H_
 # define Q_UNICODE_H_
 #if defined EXTERN
 #undef EXTERN
 #endif
 #if defined Q_UNICODEMODULE
 #define EXTERN
 #else
 #define EXTERN extern
 #endif
namespace tecplot { EXTERN ___372 ___2073(uint8_t ch); EXTERN ___372 ___2072(uint8_t ch); EXTERN ___372 ___2071(uint8_t ch); EXTERN ___372 ___2051(wchar_t wChar); EXTERN void ___492(); EXTERN ___372 ___2037(const char *S); EXTERN ___372 ___2037(tecplot::___4218 ___4229); EXTERN ___372 ___2016(const char *S); EXTERN ___372 ___2016(tecplot::___4218 S); EXTERN ___372 ___2016(const wchar_t* S);
 #if defined MSWIN
EXTERN std::string  ___2323(std::string& ___3812); EXTERN void ___2624(); EXTERN char *getenv(const char *str);
 #endif
}
 #endif
