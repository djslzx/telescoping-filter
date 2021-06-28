#ifndef TAF_H
#define TAF_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include "constants.h"
#include "remainder.h"

#define SEL_CODE_LEN (56)
#define SEL_CODE_BYTES (SEL_CODE_LEN >> 3)
#define MAX_SELECTOR 6

#define TAF_MODE_NORMAL 0
#define TAF_MODE_ARCD_OVERWRITE 1

typedef struct taf_block_t {
  rem_t remainders[64];
  uint64_t occupieds;
  uint64_t runends;
  size_t offset;
  uint8_t sel_code[SEL_CODE_BYTES];
} TAFBlock;

typedef uint64_t elt_t;

typedef struct remote_elt_t {
  uint64_t elt;
  uint64_t hash;
} Remote_elt ;

typedef struct taf_t {
  size_t p;                     /* fingerprint prefix size = log2(n/E) to get false-pos rate E */
  size_t q;                     /* length of quotient */
  size_t r;                     /* length of remainder */
  size_t nslots;                /* number of slots available (2^q) */
  size_t nblocks;               /* nslots/64 */
  size_t nelts;                 /* number of elements stored  */
  int seed;                     /* seed for Murmurhash */
  TAFBlock* blocks;           /* blocks of 64 remainders with metadata  */
  Remote_elt* remote;           /* array of inserted elements (up to 64 bits) */

  // Extra modes
  int mode;            // mode flag: handle non-adaptive case
} TAF;

void taf_init(TAF *filter, size_t n, int seed);
void taf_destroy(TAF* filter);
int taf_lookup(TAF *filter, elt_t elt);
void taf_insert(TAF *filter, elt_t elt);
void taf_clear(TAF* filter);

// Printing
double taf_load(TAF *filter);
void print_taf(TAF* filter);
void print_taf_metadata(TAF* filter);
void print_taf_block(TAF* filter, size_t block_index);

#ifdef __cplusplus
}
#endif

#endif //TAF_H
