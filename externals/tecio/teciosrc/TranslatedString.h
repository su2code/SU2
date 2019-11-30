 #ifndef TECPLOT_STRUTIL_TRANSLATEDSTRING
 #define TECPLOT_STRUTIL_TRANSLATEDSTRING
 #if defined MSWIN
 #pragma warning(disable : 4181)
 #endif
namespace tecplot { class ___4218 { public: typedef enum { ___1098, ___1104 } ___2505; ___4218(); static ___4218 ___2765(); static ___4218 ___4217(char const* str, char const* ___4219 = NULL); static ___4218 ___1097(char const* str); virtual ~___4218(); virtual bool ___2035() const; virtual bool ___2036() const; virtual char const* c_str() const;
 #if defined MSWIN && !defined MAKEARCHIVE
virtual wchar_t const* ___798() const;
 #endif
virtual operator std::string();
 #if defined MSWIN && !defined MAKEARCHIVE
virtual operator std::wstring();
 #endif
virtual ___4218& operator=(___4218 const& ___2888); ___4218(___4218 const& ___2888);
 #if !defined NO_ASSERTS
virtual bool ___2067() const;
 #endif 
private: ___4218(___4218::___2505 ___2504, char const*            str, char const*            ___4219); void ___1932(___4218::___2505 ___2504, char const*            str, char const*            ___4219); ___4218::___2505  ___2493; bool                    ___2487; std::string             ___2623; mutable std::string*    ___2666;
 #if defined MSWIN
mutable std::wstring*   ___2675;
 #endif
}; inline ___4218 ___4217(char const* str, char const* ___4219 = NULL) { return ___4218::___4217(str, ___4219); } inline ___4218 ___1097(char const* str) { return ___4218::___1097(str); } inline ___4218 ___1097(std::string const& str) { return ___4218::___1097(str.c_str()); } }
 #endif
