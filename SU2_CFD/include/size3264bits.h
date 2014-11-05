#include <stdint.h>

#ifndef SIZE_TYPE
    #if defined(__LP64__)  
       typedef int64_t SIZE_TYPE;
    #else
       typedef int32_t SIZE_TYPE;
    #endif
#endif

