#ifndef RSQF_H
#define RSQF_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include "constants.h"
#include "remainder.h"

typedef struct rsqf_block_t {
  rem_t remainders[64];
  uint64_t occupieds;
  uint64_t runends;
  size_t offset;
} RSQFBlock;

typedef struct rsqf_t {
  size_t p;                     /* fingerprint prefix size = log2(n/E) to get false-pos rate E */
  size_t q;                     /* length of quotient */
  size_t r;                     /* length of remainder */
  size_t nslots;                /* number of slots available (2^q) */
  size_t nblocks;               /* nslots/64 */
  size_t nelts;                 /* number of elements stored  */
  int seed;                     /* seed for Murmurhash */
  RSQFBlock* blocks;            /* blocks of 64 remainders with metadata  */
} RSQF;

void rsqf_init(RSQF *filter, size_t n, int seed);
void rsqf_destroy(RSQF* filter);
int rsqf_lookup(const RSQF *filter, uint64_t elt);
void rsqf_insert(RSQF *filter, uint64_t elt);
void rsqf_clear(RSQF* filter);

// Printing
double rsqf_load(RSQF* filter);
void print_rsqf(RSQF* filter);
void print_rsqf_metadata(RSQF* filter);
void print_rsqf_block(RSQF* filter, size_t block_index);

#ifdef __cplusplus
}
#endif

#endif // RSQF_H
