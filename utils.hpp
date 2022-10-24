#if defined(SUPPORT_TRACELOG)
    #include "tracelog.h"
    #define TRACELOG(level, ...) TraceLog(level, __VA_ARGS__)
#else
    #define TRACELOG(level, ...) (void)0
#endif
