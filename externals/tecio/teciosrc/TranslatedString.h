 #ifndef TECPLOT_STRUTIL_TRANSLATEDSTRING
 #define TECPLOT_STRUTIL_TRANSLATEDSTRING
 #if defined MSWIN
 #pragma warning(disable : 4181)
 #endif
namespace tecplot { class ___4218 { public: typedef enum { ___1098, ___1104 } ___2505; ___4218(); static ___4218 ___2765(); static ___4218 ___4217(const char* str, const char* ___4219 = NULL); static ___4218 ___1097(const char* str); virtual ~___4218(); virtual bool ___2035() const; virtual bool ___2036() const; virtual const char* c_str();
 #if defined MSWIN && !defined MAKEARCHIVE
virtual const wchar_t* ___798();
 #endif
virtual operator std::string();
 #if defined MSWIN && !defined MAKEARCHIVE
virtual operator std::wstring();
 #endif
virtual ___4218& operator =(const ___4218& ___2888); ___4218(const ___4218& ___2888);
 #if !defined NO_ASSERTS
virtual bool ___2067() const;
 #endif 
private: ___4218(___4218::___2505 ___2504, const char*            str, const char*            ___4219); void ___1932(___4218::___2505 ___2504, const char*            str, const char*            ___4219); ___4218::___2505  ___2493; bool                    ___2487; std::string             ___2623; std::string            *___2666;
 #if defined MSWIN
std::wstring           *___2675;
 #endif
}; inline ___4218 ___4217(const char* str, const char* ___4219 = NULL) { return ___4218::___4217(str, ___4219); } inline ___4218 ___1097(const char* str) { return ___4218::___1097(str); } inline ___4218 ___1097(const std::string& str) { return ___4218::___1097(str.c_str()); } }
 #endif
