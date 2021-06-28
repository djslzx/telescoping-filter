#ifndef AQF_REMAINDERS_H
#define AQF_REMAINDERS_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include "constants.h"

#if REM_SIZE <= 8
typedef uint8_t rem_t;
#elif REM_SIZE <= 16
typedef uint16_t rem_t;
#elif REM_SIZE <= 32
typedef uint32_t rem_t;
#else
typedef uint64_t rem_t;
#endif

#ifdef __cplusplus
}
#endif

#endif //AQF_REMAINDERS_H
