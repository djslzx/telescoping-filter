//
// Created by djsl on 3/26/21.
//

#ifndef EXTQF_H
#define EXTQF_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include "constants.h"
#include "remainder.h"
#include "ext.h"

#ifndef SUMMARIZE_EXTQF
#define SUMMARIZE_EXTQF 1
#endif

typedef struct extqf_block_t {
  rem_t remainders[64];
  uint64_t occupieds;
  uint64_t runends;
  size_t offset;
  uint8_t ext_code[EXT_CODE_BYTES];
} ExtQFBlock;

typedef uint64_t elt_t;

typedef struct extqf_t {
  size_t p;                     /* fingerprint prefix size = log2(n/E) to get false-pos rate E */
  size_t q;                     /* length of quotient */
  size_t r;                     /* length of remainder */
  size_t nslots;                /* number of slots available (2^q) */
  size_t nblocks;               /* nslots/64 */
  size_t nelts;                 /* number of elements stored  */
  int seed;                     /* seed for Murmurhash */
  ExtQFBlock* blocks;           /* blocks of 64 remainders with metadata  */
  elt_t* remote;                /* array of inserted elements (up to 64 bits) */
} ExtQF;

void extqf_init(ExtQF *filter, size_t n, int seed);
void extqf_destroy(ExtQF* filter);
int extqf_lookup(ExtQF *filter, elt_t elt);
void extqf_insert(ExtQF *filter, elt_t elt);
void extqf_clear(ExtQF* filter);

// Printing
double extqf_load(ExtQF *filter);
void print_extqf(ExtQF* filter);
void print_extqf_metadata(ExtQF* filter);
void print_extqf_block(ExtQF* filter, size_t block_index);

#ifdef __cplusplus
}
#endif

#endif //EXTQF_H
