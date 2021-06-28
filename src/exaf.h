//
// Created by djsl on 3/26/21.
//

#ifndef EXAF_H
#define EXAF_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include "constants.h"
#include "remainder.h"
#include "ext.h"

typedef struct exaf_block_t {
  rem_t remainders[64];
  uint64_t occupieds;
  uint64_t runends;
  size_t offset;
  uint8_t ext_code[EXT_CODE_BYTES];
} ExAFBlock;

typedef uint64_t elt_t;

typedef struct exaf_t {
  size_t p;                     /* fingerprint prefix size = log2(n/E) to get false-pos rate E */
  size_t q;                     /* length of quotient */
  size_t r;                     /* length of remainder */
  size_t nslots;                /* number of slots available (2^q) */
  size_t nblocks;               /* nslots/64 */
  size_t nelts;                 /* number of elements stored  */
  int seed;                     /* seed for Murmurhash */
  ExAFBlock* blocks;           /* blocks of 64 remainders with metadata  */
  elt_t* remote;                /* array of inserted elements (up to 64 bits) */
} ExAF;

void exaf_init(ExAF *filter, size_t n, int seed);
void exaf_destroy(ExAF* filter);
int exaf_lookup(ExAF *filter, elt_t elt);
void exaf_insert(ExAF *filter, elt_t elt);
void exaf_clear(ExAF* filter);

// Printing
double exaf_load(ExAF *filter);
void print_exaf(ExAF* filter);
void print_exaf_metadata(ExAF* filter);
void print_exaf_block(ExAF* filter, size_t block_index);

#ifdef __cplusplus
}
#endif

#endif //EXAF_H
