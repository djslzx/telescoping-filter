//
// Created by djsl on 3/26/21.
//

#ifndef SELQF_H
#define SELQF_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include "constants.h"
#include "remainder.h"

#define SEL_CODE_LEN (56)
#define SEL_CODE_BYTES (SEL_CODE_LEN >> 3)
#define MAX_SELECTOR 6
#ifndef SUMMARIZE
#define SUMMARIZE 0
#endif

#define SELQF_MODE_NORMAL 0
#define SELQF_MODE_ARCD_OVERWRITE 1

typedef struct selqf_block_t {
  rem_t remainders[64];
  uint64_t occupieds;
  uint64_t runends;
  size_t offset;
  uint8_t sel_code[SEL_CODE_BYTES];
} SelQFBlock;

typedef uint64_t elt_t;

typedef struct remote_elt_t {
  uint64_t elt;
  uint64_t hash;
} Remote_elt ;

#if SUMMARIZE
typedef struct elt_hist_t {
  int n_fps;       // number of false positives triggered by elt
  int n_rebuilt;  // number of times selector was reset to 0
  int max_sel;     // maximum selector value used for elt
} Elt_hist;

typedef struct selqf_stats_t {
  int n_rebuilds;       // total number of rebuilds
  Elt_hist* elt_hists;  // per-elt history
} SelQF_stats;
#endif

typedef struct selqf_t {
  size_t p;                     /* fingerprint prefix size = log2(n/E) to get false-pos rate E */
  size_t q;                     /* length of quotient */
  size_t r;                     /* length of remainder */
  size_t nslots;                /* number of slots available (2^q) */
  size_t nblocks;               /* nslots/64 */
  size_t nelts;                 /* number of elements stored  */
  int seed;                     /* seed for Murmurhash */
  SelQFBlock* blocks;           /* blocks of 64 remainders with metadata  */
  Remote_elt* remote;           /* array of inserted elements (up to 64 bits) */

  // Extra modes
  int mode;            // mode flag: handle non-adaptive case
#if SUMMARIZE
  SelQF_stats stats;
#endif
} SelQF;

void selqf_init(SelQF *filter, size_t n, int seed);
void selqf_destroy(SelQF* filter);
int selqf_lookup(SelQF *filter, elt_t elt);
void selqf_insert(SelQF *filter, elt_t elt);
void selqf_clear(SelQF* filter);

// Printing
double selqf_load(SelQF *filter);
void print_selqf(SelQF* filter);
void print_selqf_metadata(SelQF* filter);
void print_selqf_block(SelQF* filter, size_t block_index);

#ifdef __cplusplus
}
#endif

#endif //SelQF_H
