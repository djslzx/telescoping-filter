//
// Created by djsl on 5/2/21.
//

#ifndef AQF_FULL_SELQF_H
#define AQF_FULL_SELQF_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include "constants.h"
#include "remainder.h"

#define FULL_SELQF_MAX_SEL (1 << 8)

typedef struct full_selqf_block_t {
  rem_t remainders[64];
  uint64_t occupieds;
  uint64_t runends;
  size_t offset;
  uint8_t selectors[64];
} FullSelQFBlock;

typedef uint64_t elt_t;

#ifndef SELQF_H
typedef struct remote_elt_t {
  uint64_t elt;
  uint64_t hash;
} Remote_elt ;
#endif

typedef struct full_selqf_t {
  size_t p;                     /* fingerprint prefix size = log2(n/E) to get false-pos rate E */
  size_t q;                     /* length of quotient */
  size_t r;                     /* length of remainder */
  size_t nslots;                /* number of slots available (2^q) */
  size_t nblocks;               /* nslots/64 */
  size_t nelts;                 /* number of elements stored  */
  int seed;                     /* seed for Murmurhash */
  FullSelQFBlock* blocks;           /* blocks of 64 remainders with metadata  */
  Remote_elt* remote;           /* array of inserted elements (up to 64 bits) */
} FullSelQF;

#define selector(filter, i) ((filter)->blocks[(i)/64].selectors[(i)%64])

void full_selqf_init(FullSelQF *filter, size_t n, int seed);
void full_selqf_destroy(FullSelQF* filter);
int full_selqf_lookup(FullSelQF *filter, elt_t elt);
void full_selqf_insert(FullSelQF *filter, elt_t elt);
void full_selqf_clear(FullSelQF* filter);

// Printing
double full_selqf_load(FullSelQF *filter);
void print_full_selqf(FullSelQF* filter);
void print_full_selqf_metadata(FullSelQF* filter);
void print_full_selqf_stats(FullSelQF* filter);
void print_full_selqf_block(FullSelQF* filter, size_t block_index);

#ifdef __cplusplus
}
#endif

#endif //AQF_FULL_SELQF_H
