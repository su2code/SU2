 #pragma once
 #if defined UNCOPYABLE_CLASS
 #undef UNCOPYABLE_CLASS
 #endif
 #if __cplusplus >= 201103L || (defined _MSC_VER && __cplusplus >= 199711L)
 #define UNCOPYABLE_CLASS(CLASS_NAME) \
 CLASS_NAME(CLASS_NAME const&) = delete;\
 CLASS_NAME& operator=(CLASS_NAME const&) = delete;
 #else
 #define UNCOPYABLE_CLASS(CLASS_NAME) \
 private:\
 CLASS_NAME(CLASS_NAME const&);\
 CLASS_NAME& operator=(CLASS_NAME const&)
 #endif
