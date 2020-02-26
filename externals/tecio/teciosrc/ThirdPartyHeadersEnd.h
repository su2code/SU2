 #if defined(_MSC_VER)
 #pragma warning (pop)
 #if (_MSC_VER > 1600)
 #pragma warning (pop)
 #endif
 #elif defined(__clang__)
 #pragma clang diagnostic pop
 #elif defined(__GNUC__)
 #if ((__GNUC__  >= 4) && (__GNUC_MINOR__ >= 6)) || (__GNUC__ >= 5)
 #pragma GCC diagnostic pop
 #endif
 #endif
