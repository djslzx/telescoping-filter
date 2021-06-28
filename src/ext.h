// Extension type and constants

#ifndef AQF_EXT_H
#define AQF_EXT_H

typedef struct ext_t {
  uint64_t bits; // store bits of extension
  int len;       // store which bits are valid; if len=0, then ext is empty
} Ext;

// Store 0.875 extension bits per element
#define EXT_CODE_LEN 56
#define EXT_CODE_BYTES (EXT_CODE_LEN >> 3)

#endif //AQF_EXT_H
