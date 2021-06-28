//
// Created by djsl on 5/2/21.
//

#ifndef AQF_UTAF_H
#define AQF_UTAF_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include "constants.h"
#include "remainder.h"

#define UTAF_MAX_SEL (1 << 8)

typedef struct utaf_block_t {
  rem_t remainders[64];
  uint64_t occupieds;
  uint64_t runends;
  size_t offset;
  uint8_t selectors[64];
} FullTAFBlock;

typedef uint64_t elt_t;

#ifndef TAF_H
typedef struct remote_elt_t {
  uint64_t elt;
  uint64_t hash;
} Remote_elt ;
#endif

typedef struct utaf_t {
  size_t p;                     /* fingerprint prefix size = log2(n/E) to get false-pos rate E */
  size_t q;                     /* length of quotient */
  size_t r;                     /* length of remainder */
  size_t nslots;                /* number of slots available (2^q) */
  size_t nblocks;               /* nslots/64 */
  size_t nelts;                 /* number of elements stored  */
  int seed;                     /* seed for Murmurhash */
  FullTAFBlock* blocks;           /* blocks of 64 remainders with metadata  */
  Remote_elt* remote;           /* array of inserted elements (up to 64 bits) */
} FullTAF;

#define selector(filter, i) ((filter)->blocks[(i)/64].selectors[(i)%64])

void utaf_init(FullTAF *filter, size_t n, int seed);
void utaf_destroy(FullTAF* filter);
int utaf_lookup(FullTAF *filter, elt_t elt);
void utaf_insert(FullTAF *filter, elt_t elt);
void utaf_clear(FullTAF* filter);

// Printing
double utaf_load(FullTAF *filter);
void print_utaf(FullTAF* filter);
void print_utaf_metadata(FullTAF* filter);
void print_utaf_stats(FullTAF* filter);
void print_utaf_block(FullTAF* filter, size_t block_index);

#ifdef __cplusplus
}
#endif

#endif //AQF_UTAF_H
