//
// Created by djsl on 5/2/21.
//

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <execinfo.h>

#include "murmur3.h"
#include "macros.h"
#include "arcd.h"
#include "utaf.h"
#include "bit_util.h"
#include "set.h"

/**
 * Generate a hash for the input word.
 * Returns the full 128 bit murmurhash.
 * Expects a two-slot array of uint64_t's.
 */
static uint64_t utaf_hash(const FullTAF *filter, elt_t elt) {
  uint64_t buf[2];
  MurmurHash3_x64_128(&elt, 8, filter->seed, buf);
  return buf[0];
}

/**
 * Returns the quotient for a 64-bit fingerprint hash.
 */
static size_t calc_quot(const FullTAF* filter, uint64_t hash) {
  return hash & ONES(filter->q);
}

/**
 * Returns the k-th remainder for h
 */
static rem_t calc_rem(const FullTAF* filter, uint64_t hash, int k) {
  int n_rems = (64 - (int)filter->q)/(int)filter->r;
  if (k >= n_rems) k %= n_rems;
  return (hash >> (filter->q + k * filter->r)) & ONES(filter->r);
}

/* FullTAF Helpers */

/**
 * Print an array of selectors.
 */
static void print_sels(const uint8_t sels[64]) {
  for (int i=0; i<8; i++) {
    printf("   ");
    for (int j=0; j<8; j++) {
      int sel = sels[i*8 + j];
      if (sel == 0) {
        printf(" _");
      } else {
        printf(" %d", sel);
      }
    }
    printf("\n");
  }
}

/**
 * Returns the absolute index of the `rank`-th 1 bit in Q.runends past the start of
 * the block at `block_index`. `rank` indexes from 0.
 *
 * Returns -1 if result is invalid (out of bounds).
 */
static int select_runend(const FullTAF* filter, size_t block_index, size_t rank) {
  assert(block_index < filter->nblocks && "block_index out of bounds");

  size_t step;
  size_t loc = block_index * 64;
  while (1) {
    FullTAFBlock* b = &filter->blocks[loc / 64];
    step = bitselect(b->runends, rank >= 64 ? 63 : (int)rank);
    loc += step;
    if (step != 64 || loc >= filter->nslots) {
      break;
    }
    rank -= popcnt(b->runends);
  }
  if (loc >= filter->nslots) {
    return -1;
  } else {
    return (int)loc;
  }
}

#define RANK_SELECT_EMPTY (-1)
#define RANK_SELECT_OVERFLOW (-2)
/** Performs the blocked equivalent of the unblocked operation
 *    y = select(Q.runends, rank(Q.occupieds, x)).
 *  Note: x indexes from 0.
 *
 *  Return behavior:
 *  - If y <= x, returns Empty
 * - If y > x, returns Full(y)
 * - If y runs off the edge, returns Overflow
 */
static int rank_select(const FullTAF* filter, size_t x) {
  // Exit early if x obviously out of range
  if (x >= filter->nslots) {
    return RANK_SELECT_OVERFLOW;
  }
  size_t block_i = x/64;
  size_t slot_i = x%64;
  FullTAFBlock *b = &filter->blocks[block_i];

  // Compute i + O_i where i = x - (x mod 64)
  if (!GET(b->occupieds, 0) && b->offset == 0 && !GET(b->runends, 0)) {
    // b[0] unoccupied, b.offset = 0, b[0] not a runend =>
    // negative offset
    if (slot_i == 0) {
      return RANK_SELECT_EMPTY;
    }
  } else {
    // non-negative offset
    if (slot_i == 0) {
      return (int)(block_i * 64 + b->offset);
    } else {
      block_i += b->offset/64;
    }
  }

  // Handle case where offset runs off the edge
  if (block_i >= filter->nblocks) {
    return RANK_SELECT_OVERFLOW;
  }

  // Count the number of occupied quotients between i+1 (b.start + i) and j (x)
  uint64_t d = bitrank(b->occupieds, slot_i) - GET(b->occupieds, 0);

  // Advance offset to relevant value for the block that b.offset points to
  size_t offset = b->offset % 64;
  b = &filter->blocks[block_i];

  // Account for the runends in [0, offset] of the new block
  d += bitrank(b->runends, offset);

  // If rank(Q.occupieds, x) == 0, then there's nothing to see here
  if (d == 0) {
    return RANK_SELECT_EMPTY;
  } else {
    // (rank-1) accounts for select's indexing from 0
    int loc = select_runend(filter, block_i, d-1);
    if (loc == -1) {
      return RANK_SELECT_OVERFLOW;
    } else if (loc < x) {
      return RANK_SELECT_EMPTY;
    } else {
      return loc;
    }
  }
}

#define NO_UNUSED (-1)
/**
 * Finds the first unused slot at or after absolute location x.
 */
static int first_unused(const FullTAF* filter, size_t x) {
  while (1) {
    int loc = rank_select(filter, x);
    switch (loc) {
      case RANK_SELECT_EMPTY: return x;
      case RANK_SELECT_OVERFLOW: return NO_UNUSED;
      default:
        if (x <= loc) {
          x = loc + 1;
        } else {
          return x;
        }
    }
  }
}

/**
 * Shift the remainders and runends in [a, b] forward by 1 into [a+1, b+1]
 */
static void shift_rems_and_runends(FullTAF* filter, int a, int b) {
  if (a > b) return;
  for (int i=b; i>=a; i--) {
    remainder(filter, i+1) = remainder(filter, i);
    set_runend_to(filter, i+1, get_runend(filter, i));
  }
  set_runend_to(filter, a, 0);
}

/**
 * Shift the remote elements in [a,b] forward by 1
 */
static void shift_remote_elts(FullTAF* filter, int a, int b) {
  if (a > b) return;
  for (int i=b; i>=a; i--) {
    filter->remote[i+1] = filter->remote[i];
  }
  filter->remote[a].elt = 0;
  filter->remote[a].hash = 0;
}

/**
 * Shift the hash selectors in [a,b] forward by 1
 */
static void shift_sels(FullTAF* filter, int a, int b) {
  if (a > b) return;
  for (int i=b; i>=a; i--) {
    selector(filter, i+1) = selector(filter, i);
  }
  selector(filter, a) = 0;
}

/**
 * Increment all non-negative offsets with targets in [a,b]
 */
static void inc_offsets(FullTAF* filter, size_t a, size_t b) {
  assert(a < filter->nslots && b < filter->nslots);
  // Exit early if invalid range
  if (a > b) {
    return;
  }
  // Start i at the first block after b, clamping it so it doesn't go off the end, and work backwards
  size_t start = min(b/64 + 1, filter->nblocks - 1);
  for (int i = start; i>=0; i--) {
    FullTAFBlock *block = &filter->blocks[i];
    size_t block_start = i * 64;
    // Skip this block if it has a negative offset
    if (!GET(block->occupieds, 0) &&
        block->offset == 0 &&
        !GET(block->runends, 0)) {
      continue;
    }
    // Exit if the target for b.offset is before the interval;
    // if it's within the interval, increment offset
    size_t target = block_start + block->offset;
    if (target < a) {
      break;
    } else if (target <= b) {
      block->offset++;
    }
  }
}

/**
 * Increment non-negative offsets to accommodate insertion of a new run
 * for `quot` at `loc`.
 *
 * Concretely, this function increments unowned offsets in blocks whose
 * first slot `s` is not after `quot`: `s >= quot`.
 */
static void inc_offsets_for_new_run(FullTAF* filter, size_t quot, size_t loc) {
  assert(loc < filter->nslots);
  // Start i at the first block after loc,
  // clamping it so it doesn't go off the end
  size_t start = min(loc/64 + 1, filter->nblocks - 1);
  for (int i=start; i>=0; i--) {
    FullTAFBlock *b = &filter->blocks[i];
    size_t b_start = i*64;
    // Skip this block if it has a negative offset
    if (!GET(b->occupieds, 0) && b->offset == 0 && !GET(b->runends, 0)) {
      continue;
    }
    // Exit if the target for b.offset is before the interval;
    // if the target is within the interval, increment b.offset
    size_t target = b_start + b->offset;
    if (target < loc) {
      break;
    } else if (target == loc && !GET(b->occupieds, 0) && quot <= b_start) {
      b->offset++;
    }
  }
}

static void add_block(FullTAF *filter) {
  // Add block to new_blocks
  FullTAFBlock *new_blocks = realloc(filter->blocks, (filter->nblocks + 1) * sizeof(FullTAFBlock));
  if (new_blocks == NULL) {
    printf("add_block failed to realloc new blocks\n");
    exit(1);
  }
  filter->blocks = new_blocks;
  memset(filter->blocks + filter->nblocks, 0, sizeof(FullTAFBlock));

  // Reallocate remote rep
  Remote_elt *new_remote = realloc(filter->remote,(filter->nslots + 64) * sizeof(Remote_elt));
  if (new_remote == NULL) {
    printf("add_block failed to realloc new remote rep\n");
    exit(1);
  }
  filter->remote = new_remote;
  memset(filter->remote + filter->nslots, 0, 64 * sizeof(Remote_elt));

  // Update counters
  filter->nblocks += 1;
  filter->nslots += 64;
}

/**
 * Adapt a fingerprint at a particular location by incrementing the selector and
 * updating the remainder.
 */
static void adapt_loc(FullTAF *filter, size_t loc) {
  int old_sel = selector(filter, loc);
  int new_sel = (old_sel + 1) % UTAF_MAX_SEL;
  selector(filter, loc) = new_sel;
  remainder(filter, loc) = calc_rem(filter, filter->remote[loc].hash, new_sel);
}

/**
 * Adapt on a query element that collided with a stored fingerprint at loc.
 *
 * Go through the rest of the run and fix any other remaining collisions.
 */
static void adapt(FullTAF *filter, elt_t query, int loc, size_t quot, uint64_t hash) {
  assert(quot <= loc && loc < filter->nslots);
  // Make sure the query elt isn't mapped to an earlier index in the sequence
  for (int i=loc; i>=(int)quot && (i == loc || !get_runend(filter, i)); i--) {
    if (filter->remote[i].elt == query) {
      return;
    }
  }
  // Adapt on all collisions in the run
  for (int i=loc; i>=(int)quot && (i == loc || !get_runend(filter, i)); i--) {
    if (remainder(filter, i) == calc_rem(filter, hash, selector(filter, i))) {
      adapt_loc(filter, i);
    }
  }
}

/* FullTAF */

void utaf_init(FullTAF *filter, size_t n, int seed) {
  filter->seed = seed;
  filter->nelts = 0;
  filter->nblocks = max(1, nearest_pow_of_2(n)/64);
  filter->nslots = filter->nblocks * 64;
  filter->q = (size_t)log2((double)filter->nslots); // nslots = 2^q
  filter->r = REM_SIZE;
  filter->p = filter->q + filter->r;
  filter->blocks = calloc(filter->nblocks, sizeof(FullTAFBlock));
  filter->remote = calloc(filter->nslots, sizeof(Remote_elt));
}

void utaf_destroy(FullTAF* filter) {
  free(filter->blocks);
  free(filter->remote);
  free(filter);
}

void utaf_clear(FullTAF* filter) {
  filter->nelts = 0;
  free(filter->blocks);
  free(filter->remote);
  filter->blocks = calloc(filter->nblocks, sizeof(FullTAFBlock));
  filter->remote = calloc(filter->nslots, sizeof(Remote_elt));
}

static void raw_insert(FullTAF* filter, elt_t elt, uint64_t hash) {
  size_t quot = calc_quot(filter, hash);
  rem_t rem = calc_rem(filter, hash, 0);
  filter->nelts++;

  // Find the appropriate runend
  int r = rank_select(filter, quot);
  switch (r) {
    case RANK_SELECT_EMPTY: {
      set_occupied(filter, quot);
      set_runend(filter, quot);
      remainder(filter, quot) = rem;
      filter->remote[quot].elt = elt;
      filter->remote[quot].hash = hash;
      break;
    }
    case RANK_SELECT_OVERFLOW: {
      printf("FullTAF failed to find runend (nslots=%lu, quot=(block=%lu, slot=%lu))\n",
             filter->nslots, quot/64, quot%64);
      exit(1);
    }
    default: {
      // Find u, the first open slot after r, and
      // shift everything in [r+1, u-1] forward by 1 into [r+2, u],
      // leaving r+1 writable
      size_t u = first_unused(filter, r+1);
      if (u == NO_UNUSED) {
        // Extend filter by one block and use the first empty index
        add_block(filter);
        u = filter->nslots - 64;
      }
      inc_offsets(filter, r+1, u-1);
      shift_rems_and_runends(filter, r + 1, (int)u - 1);
      shift_remote_elts(filter, r + 1, (int)u - 1);
      shift_sels(filter, r + 1, (int)u - 1);

      // Start a new run or extend an existing one
      if (get_occupied(filter, quot)) {
        // quot occupied: extend an existing run
        inc_offsets(filter, r, r);
        unset_runend(filter, r);
      } else {
        // quot unoccupied: start a new run
        inc_offsets_for_new_run(filter, quot, r);
        set_occupied(filter, quot);
      }
      set_runend(filter, r+1);
      remainder(filter, r+1) = rem;
      filter->remote[r+1].elt = elt;
      filter->remote[r+1].hash = hash;
    }
  }
}

static int raw_lookup(FullTAF* filter, elt_t elt, uint64_t hash) {
  size_t quot = calc_quot(filter, hash);

  if (get_occupied(filter, quot)) {
    int loc = rank_select(filter, quot);
    if (loc == RANK_SELECT_EMPTY || loc == RANK_SELECT_OVERFLOW) {
      return 0;
    }
    do {
      int sel = selector(filter, loc);
      rem_t rem = calc_rem(filter, hash, sel);
      if (remainder(filter, loc) == rem) {
        // Check remote
        if (elt != filter->remote[loc].elt) {
          adapt(filter, elt, loc, quot, hash);
        }
        return 1;
      }
      loc--;
    } while (loc >= (int)quot && !get_runend(filter, loc));
  }
  return 0;
}

/**
 * Return 1 if word is in the filter.
 *
 * Exits with 0 immediately if quot(word) is unoccupied.
 * Otherwise, linear probes through the run to see if the
 * run contains rem(word).
 */
int utaf_lookup(FullTAF *filter, elt_t elt) {
  uint64_t hash = utaf_hash(filter, elt);
  return raw_lookup(filter, elt, hash);
}

void utaf_insert(FullTAF *filter, elt_t elt) {
  uint64_t hash = utaf_hash(filter, elt);
  raw_insert(filter, elt, hash);
}

double utaf_load(FullTAF *filter) {
  return (double)filter->nelts/(double)filter->nslots;
}

/* Printing */

void print_utaf_metadata(FullTAF* filter) {
  printf("FILTER METADATA:\n");
  printf("  p=%ld, q=%ld, r=%ld\n",
         filter->p, filter->q, filter->r);
  printf("  nslots=%ld, nblocks=%ld, blocksize=%ld, nelts=%ld\n",
         filter->nslots, filter->nslots/64, sizeof(FullTAFBlock), filter->nelts);
  printf("  seed=%d\n", filter->seed);
  printf("  load factor=%f\n", utaf_load(filter));
}

void print_utaf_block(FullTAF* filter, size_t block_index) {
  assert(0 <= block_index && block_index < filter->nslots/64);
  FullTAFBlock block = filter->blocks[block_index];
  printf("BLOCK %lu:\n", block_index);
  printf("  occupieds=0x%lx\n", block.occupieds);
  printf("  runends=0x%lx\n", block.runends);
  printf("  offset=%ld\n", block.offset);
  printf("  remainders=\n");
  // Print out 8x8
  for (int i=0; i<8; i++) {
    printf("   ");
    for (int j=0; j<8; j++) {
      printf(get_occupied(filter, block_index*64 + i*8 + j) ? "o" : " ");
      printf(get_runend(filter, block_index*64 + i*8 + j) ? "r" : " ");
      printf(" 0x%-*x", (int)(filter->r / 8 + 3), block.remainders[i*8+j]);
    }
    printf("\n");
  }
  printf("  selectors=\n");
  print_sels(filter->blocks[block_index].selectors);
  printf("  remote elts=\n");
  for (int i=0; i<8; i++) {
    printf("   ");
    for (int j=0; j<8; j++) {
      printf(" 0x%-*lx", 8, filter->remote[block_index * 64 + i*8 + j].elt);
    }
    printf("\n");
  }
}

void print_utaf(FullTAF* filter) {
  print_utaf_metadata(filter);
  for (int i=0; i<filter->nblocks; i++) {
    print_utaf_block(filter, i);
  }
}
void print_utaf_stats(FullTAF* filter) {
  printf("FullTAF stats:\n");
  // Hash selector counts
  int max_sel = 0;
  for (int i = 0; i < filter->nslots; i++) {
    int sel = selector(filter, i);
    if (sel > max_sel) {
      max_sel = sel;
    }
  }
  int sel_counts[max_sel];
  for (int i = 0; i < max_sel; i++) {
    sel_counts[i] = 0;
  }
  for (int i = 0; i < filter->nslots; i++) {
    sel_counts[selector(filter, i)]++;
  }
  printf("Hash selector counts:\n");
  for (int i = 0; i < max_sel; i++) {
    printf(" %d: %d (%f%%)\n", i, sel_counts[i],
           100 * (double) sel_counts[i] / (double) filter->nslots);
  }
}

// Tests
//#define TEST_UTAF 1
#ifdef TEST_UTAF

void print_backtrace() {
  void* callstack[128];
  int i, frames = backtrace(callstack, 128);
  char** strs = backtrace_symbols(callstack, frames);
  printf("\n");
  for (i = 0; i < frames; ++i) {
    printf("%s\n", strs[i]);
  }
  free(strs);
}

#define assert_eq(a, b) assert((a) == (b))

#define test_assert_eq(a, b, msg, ...)   \
  if ((a) != (b)) {                 \
    do {                            \
      fprintf(stderr, "Assertion failed: %s != %s: ", #a, #b); \
      fprintf(stderr, msg"\n", __VA_ARGS__); \
      assert_eq(a, b);              \
    } while (0);                    \
  }

#define FullTAF_SEED 32776517

FullTAF *new_utaf(size_t n) {
  FullTAF *filter = malloc(sizeof(FullTAF));
  utaf_init(filter, n, FullTAF_SEED);
  return filter;
}

void test_add_block() {
  printf("Testing %s...", __FUNCTION__);
  FullTAF *filter = new_utaf(64 * 2);
  assert_eq(filter->nslots, 128);
  assert_eq(filter->nblocks, 2);
  add_block(filter);
  // Check metadata
  assert_eq(filter->nslots, 192);
  assert_eq(filter->nblocks, 3);
  // Check new block
  FullTAFBlock b = filter->blocks[2];
  assert_eq(b.occupieds, 0);
  assert_eq(b.runends, 0);
  assert_eq(b.offset, 0);
  for (int i=0; i<64; i++) {
    assert_eq(b.remainders[i], 0);
  }
  // Check remote rep
  for (int i=0; i<64; i++) {
    assert_eq(filter->remote[128 + i].elt, 0);
    assert_eq(filter->remote[128 + i].hash, 0);
  }
  utaf_destroy(filter);
  printf("passed.\n");
}

/// Check that adding a block doesn't overwrite existing data
void test_add_block_no_clobber() {
  printf("Testing %s...", __FUNCTION__);
  FullTAF *filter = new_utaf(128);

  // Setup
  for (int i=0; i<filter->nslots; i++) {
    set_occupied(filter, i);
    set_runend(filter, i);
    remainder(filter, i) = i%16;
    filter->remote[i].elt = i;
    filter->remote[i].hash = i;
  }
  add_block(filter);
  // Check that data in first 2 blocks is preserved
  for (int i=0; i<128; i++) {
    assert(get_occupied(filter, i));
    assert(get_runend(filter, i));
    assert_eq(remainder(filter, i), i%16);
    assert_eq(filter->remote[i].elt, i);
    assert_eq(filter->remote[i].hash, i);
  }
  // Check that 3rd block is empty
  for (int i=128; i<filter->nslots; i++) {
    assert(!get_occupied(filter, i));
    assert(!get_runend(filter, i));
    assert_eq(remainder(filter, i), 0);
    assert_eq(filter->remote[i].elt, 0);
    assert_eq(filter->remote[i].hash, 0);
  }
  // Check filter metadata
  assert_eq(filter->nslots, 192);
  assert_eq(filter->nblocks, 3);

  utaf_destroy(filter);
  printf("passed.\n");
}

void test_adapt_loc_1() {
  printf("Testing %s...", __FUNCTION__);
  FullTAF *filter = new_utaf(128);
  // q=7, r=8
  for (int i=0; i<filter->nslots; i++) {
    assert_eq(selector(filter, i), 0);
  }
  adapt_loc(filter, 0);
  assert_eq(selector(filter, 0), 1);
  for (int i=1; i<64; i++) {
    assert_eq(selector(filter, i), 0);
  }
  utaf_destroy(filter);
  printf("passed.\n");
}

void test_adapt_loc_2() {
  printf("Testing %s...", __FUNCTION__);
  FullTAF *filter = new_utaf(128);
  // q=7, r=8
  for (int i=0; i<64; i++) {
    selector(filter, i) = 0;
  }
  // Adapt a lot
  int limit = 15;
  for (int i=0; i<limit; i++) {
    adapt_loc(filter, i);
  }
  for (int i=0; i<filter->nslots; i++) {
    if (i < limit) {
      assert_eq(selector(filter, i), 1);
    } else {
      assert_eq(selector(filter, i), 0);
    }
  }
  adapt_loc(filter, limit);
  for (int i=0; i<filter->nslots; i++) {
    test_assert_eq(selector(filter, i),
                   i <= limit ? 1 : 0,
                   "i=%d", i);
  }
  utaf_destroy(filter);
  printf("passed.\n");
}

void test_shift_remote_elts() {
  printf("Testing %s...", __FUNCTION__);
  FullTAF *filter = new_utaf(128);
  for (int i=0; i<filter->nslots; i++) {
    assert_eq(filter->remote[i].elt, 0);
    assert_eq(filter->remote[i].hash, 0);
  }
  for (int i=0; i<filter->nslots; i++) {
    filter->remote[i].elt = i;
    filter->remote[i].hash = i;
  }
  // Shift elts in [32, 64+32] to [33, 64+33]
  shift_remote_elts(filter, 32, 64+32);
  for (int i=0; i<=31; i++) {
    assert_eq(filter->remote[i].elt, i);
    assert_eq(filter->remote[i].hash, i);
  }
  assert_eq(filter->remote[32].elt, 0);
  assert_eq(filter->remote[32].hash, 0);
  for (int i=33; i<=64+33; i++) {
    assert_eq(filter->remote[i].elt, i-1);
    assert_eq(filter->remote[i].hash, i-1);
  }
  for (int i=64+34; i<filter->nslots; i++) {
    assert_eq(filter->remote[i].elt, i);
    assert_eq(filter->remote[i].hash, i);
  }
  utaf_destroy(filter);
  printf("passed.\n");
}

/// General integration test: insert and query elts, ensuring that there
/// are no false negatives
void test_insert_and_query() {
  printf("Testing %s...", __FUNCTION__);
  size_t a = 1 << 20;
  double a_s = 100.0; // a/s
  double load = 0.95;
  size_t s = nearest_pow_of_2((size_t)((double)a / a_s));
  s = (size_t)((double)s * load);
  FullTAF* filter = new_utaf(s);

  // Generate query set
  srandom(FullTAF_SEED);
  int nset = (int)(1.5*(double)s);
  Setnode* set = calloc(nset, sizeof(set[0]));
  char str[64];
  for (int i=0; i<s; i++) {
    elt_t elt = random() % a;
    sprintf(str, "%lu", elt);
    set_insert(str, (int)strlen(str), 0, set, nset);
    utaf_insert(filter, elt);
    assert(set_lookup(str, (int)strlen(str), set, nset));
    assert(utaf_lookup(filter, elt));
  }
  // Query [0, a] and ensure that all items in the set return true
  int fps = 0;
  int fns = 0;
  for (int i=0; i<a; i++) {
    elt_t elt = i;
    sprintf(str, "%lu", elt);
    int in_set = set_lookup(str, (int)strlen(str), set, nset);
    int in_utaf = utaf_lookup(filter, elt);
    if (in_set && !in_utaf) {
      fns++;
      uint64_t hash = utaf_hash(filter, elt);
      size_t quot = calc_quot(filter, hash);
      if (get_occupied(filter, quot)) {
        int loc = rank_select(filter, quot);
        if (loc == RANK_SELECT_EMPTY || loc == RANK_SELECT_OVERFLOW) {
          printf("False negative (elt=%lu, 0x%lx): quot=%lu (block=%lu, slot=%lu)"
                 " was occupied but didn't have an associated runend\n",
                 elt, elt, quot, quot/64, quot%64);
          print_utaf_block(filter, quot/64);
          exit(1);
        } else {
          int sel = selector(filter, loc);
          rem_t query_rem = calc_rem(filter, hash, sel);
          rem_t stored_rem = remainder(filter, loc);
          printf("False negative (elt=%lu, 0x%lx): quot=%lu (block=%lu, slot=%lu),"
                 "loc=%d (block=%d, slot=%d); stored rem=0x%hhx doesn't match query rem=0x%hhx\n",
                 elt, elt, quot, quot/64, quot%64, loc, loc/64, loc%64, stored_rem, query_rem);
          print_utaf_metadata(filter);
          print_utaf_block(filter, loc/64);
          exit(1);
        }
      } else {
        printf("False negative (elt=%lu, 0x%lx): quot=%lu (block=%lu, slot=%lu) wasn't occupied\n",
               elt, elt, quot, quot/64, quot%64);
        exit(1);
      }
    } else if (!in_set && in_utaf) {
      fps += in_utaf;
    }
  }
  printf("passed. ");
  printf("FPs: %d (%f%%), FNs: %d (%f%%)\n",
         fps, (double)fps/(double)a * 100, fns, (double)fns/(double)a * 100);
  print_utaf_metadata(filter);
  print_utaf_stats(filter);
  utaf_destroy(filter);
}

void test_insert_and_query_w_repeats() {
  printf("Testing %s...\n", __FUNCTION__);
  int nslots = 1 << 14;
  double load = 0.95;
  double a_s = 100;
  int queries_per_elt = 10;
  printf("nslots=%d, load=%f, a/s=%f, queries_per_elt = %d\n", nslots, load, a_s, queries_per_elt);

  int s = (int)((double)nearest_pow_of_2(nslots) * load);
  int a = (int)((double)s * a_s);
  int n_queries = a * queries_per_elt;

  int fps = 0;  // false positives
  int rfps = 0; // repeated false positives
  int fns = 0;  // false negatives
  int tot_queries = n_queries * queries_per_elt;

  FullTAF *filter = new_utaf(s);
  int nset = (int)(s * 1.5);
  Setnode *set = calloc(nset, sizeof(set[0]));

  srandom(FullTAF_SEED);
  char str[64];
  int len;
  fprintf(stderr, "Initializing membership set and filter...\n");
  for (int i=0; i<s; i++) {
    elt_t elt = random();
    sprintf(str, "%lu", elt);
    len = (int)strlen(str);
    set_insert(str, len, 0, set, nset);
    utaf_insert(filter, elt);
  }
  fprintf(stderr, "Initializing query set...\n");
  elt_t *query_set = calloc(a, sizeof(elt_t));
  for (int i=0; i<a; i++) {
    query_set[i] = random();
  }
  fprintf(stderr, "Querying set and filter...\n");
  int nseen = (int)(s * 1.5);
  Setnode *seen = calloc(nseen, sizeof(seen[0]));
  for (int i=0; i<n_queries; i++) {
    elt_t elt = query_set[random() % a];
    sprintf(str, "%lu", elt);
    len = (int)strlen(str);
    int in_filter = utaf_lookup(filter, elt);
    int in_set = set_lookup(str, len, set, nset);
    if (in_filter && !in_set) {
      fps++;
      if (set_lookup(str, len, seen, nseen)) {
        rfps++;
      } else {
        set_insert(str, len, 0, seen, nseen);
      }
    } else if (!in_filter && in_set) {
      fns++;
      uint64_t hash = utaf_hash(filter, elt);
      size_t quot = calc_quot(filter, hash);
      if (get_occupied(filter, quot)) {
        int loc = rank_select(filter, quot);
        if (loc == RANK_SELECT_EMPTY || loc == RANK_SELECT_OVERFLOW) {
          printf("False negative (elt=%lu, 0x%lx): quot=%lu (block=%lu, slot=%lu)"
                 " was occupied but didn't have an associated runend\n",
                 elt, elt, quot, quot/64, quot%64);
          print_utaf_block(filter, quot/64);
          exit(1);
        } else {
          int sel = selector(filter, loc);
          rem_t query_rem = calc_rem(filter, hash, sel);
          rem_t stored_rem = remainder(filter, loc);
          printf("False negative (elt=%lu, 0x%lx): quot=%lu (block=%lu, slot=%lu),"
                 "loc=%d (block=%d, slot=%d); stored rem=0x%hhx doesn't match query rem=0x%hhx\n",
                 elt, elt, quot, quot/64, quot%64, loc, loc/64, loc%64, stored_rem, query_rem);
          print_utaf_metadata(filter);
          print_utaf_block(filter, loc/64);
          exit(1);
        }
      } else {
        printf("False negative (elt=%lu, 0x%lx): quot=%lu (block=%lu, slot=%lu) wasn't occupied\n",
               elt, elt, quot, quot/64, quot%64);
        exit(1);
      }
    }
  }
  printf("Test results:\n");
  printf("FPs: %d (%f%%), RFPs: %d (%f%%)\n",
         fps, (double)fps/tot_queries, rfps, (double)rfps/tot_queries * 100);
  printf("FNs: %d (%f%%)\n", fns, (double)fns/tot_queries * 100);
  print_utaf_stats(filter);
  utaf_destroy(filter);
  printf("Done testing %s.\n", __FUNCTION__);
}

void test_mixed_insert_and_query_w_repeats() {
  printf("Testing %s...\n", __FUNCTION__);
  int nslots = 1 << 14;
  double load = 0.95;
  double a_s = 100;
  int queries_per_elt = 10;
  printf("nslots=%d, load=%f, a/s=%f, queries_per_elt = %d\n", nslots, load, a_s, queries_per_elt);

  int s = (int)((double)nearest_pow_of_2(nslots) * load);
  int a = (int)((double)s * a_s);
  int n_queries = a * queries_per_elt;

  int fps = 0;  // false positives
  int rfps = 0; // repeated false positives
  int fns = 0;  // false negatives
  int tot_queries = n_queries * queries_per_elt;

  FullTAF *filter = new_utaf(s);
  int nset = (int)(s * 1.5);
  Setnode *set = calloc(nset, sizeof(set[0]));

  srandom(FullTAF_SEED);
  char str[64];
  int len;
  fprintf(stderr, "Initializing query set...\n");
  elt_t *query_set = calloc(a, sizeof(elt_t));
  for (int i=0; i<a; i++) {
    query_set[i] = random();
  }
  fprintf(stderr, "Initializing membership set and filter...\n");
  for (int i=0; i<s/2; i++) {
    elt_t elt = random();
    sprintf(str, "%lu", elt);
    len = (int)strlen(str);
    set_insert(str, len, 0, set, nset);
    utaf_insert(filter, elt);
  }
  fprintf(stderr, "Performing interleaved queries...\n");
  for (int i=0; i<n_queries; i++) {
    elt_t elt = query_set[random() % a];
    utaf_lookup(filter, elt);
  }
  fprintf(stderr, "Finishing initialization of membership set and filter...\n");
  for (int i=s/2; i<s; i++) {
    elt_t elt = random();
    sprintf(str, "%lu", elt);
    len = (int)strlen(str);
    set_insert(str, len, 0, set, nset);
    utaf_insert(filter, elt);
  }
  fprintf(stderr, "Querying set and filter...\n");
  int nseen = (int)(s * 1.5);
  Setnode *seen = calloc(nseen, sizeof(seen[0]));
  for (int i=0; i<n_queries; i++) {
    elt_t elt = query_set[random() % a];
    sprintf(str, "%lu", elt);
    len = (int)strlen(str);
    int in_filter = utaf_lookup(filter, elt);
    int in_set = set_lookup(str, len, set, nset);
    if (in_filter && !in_set) {
      fps++;
      if (set_lookup(str, len, seen, nseen)) {
        rfps++;
      } else {
        set_insert(str, len, 0, seen, nseen);
      }
    } else if (!in_filter && in_set) {
      fns++;
      uint64_t hash = utaf_hash(filter, elt);
      size_t quot = calc_quot(filter, hash);
      if (get_occupied(filter, quot)) {
        int loc = rank_select(filter, quot);
        if (loc == RANK_SELECT_EMPTY || loc == RANK_SELECT_OVERFLOW) {
          printf("False negative (elt=%lu, 0x%lx): quot=%lu (block=%lu, slot=%lu)"
                 " was occupied but didn't have an associated runend\n",
                 elt, elt, quot, quot/64, quot%64);
          print_utaf_block(filter, quot/64);
          exit(1);
        } else {
          int sel = selector(filter, loc);
          rem_t query_rem = calc_rem(filter, hash, sel);
          rem_t stored_rem = remainder(filter, loc);
          printf("False negative (elt=%lu, 0x%lx): quot=%lu (block=%lu, slot=%lu),"
                 "loc=%d (block=%d, slot=%d); stored rem=0x%hhx doesn't match query rem=0x%hhx\n",
                 elt, elt, quot, quot/64, quot%64, loc, loc/64, loc%64, stored_rem, query_rem);
          print_utaf_metadata(filter);
          print_utaf_block(filter, loc/64);
          exit(1);
        }
      } else {
        printf("False negative (elt=%lu, 0x%lx): quot=%lu (block=%lu, slot=%lu) wasn't occupied\n",
               elt, elt, quot, quot/64, quot%64);
        exit(1);
      }
    }
  }
  printf("Test results:\n");
  printf("FPs: %d (%f%%), RFPs: %d (%f%%)\n",
         fps, (double)fps/tot_queries, rfps, (double)rfps/tot_queries * 100);
  printf("FNs: %d (%f%%)\n", fns, (double)fns/tot_queries * 100);
  print_utaf_stats(filter);
  utaf_destroy(filter);
  printf("Done testing %s.\n", __FUNCTION__);
}

// a, b+1 in the same block
void test_shift_sels_single_block() {
  printf("Testing %s...", __FUNCTION__);

  // Setup
  FullTAF *filter = new_utaf(64 * 3);
  for (int i=0; i<64; i++) {
    selector(filter, i) = i % 8 == 0;
  }

  // Shift sels in [0, 62] -> [1, 63] and check
  for (int j=1; j <= 64; j++) {
    shift_sels(filter, 0, 62);
    for (int i=0; i<64; i++) {
      test_assert_eq(selector(filter, i),
                     (i - j) % 8 == 0 && i > j-1,
                     "j=%d, i=%d", j, i);
    }
  }
  utaf_destroy(filter);
  printf("passed.\n");
}

FullTAF* sel_setup() {
  FullTAF* filter = new_utaf(64 * 4);
  for (int i=0; i<filter->nslots; i++) {
    selector(filter, i) = i % 8 == 0;
  }
  return filter;
}

// a, b+1 in different blocks
void test_shift_sels_multi_block() {
  printf("Testing %s...", __FUNCTION__);
  FullTAF *filter = sel_setup();
  // (1) Shift sels in [0, 127] -> [1,128]
  shift_sels(filter, 0, 127);
  for (int i=0; i<filter->nblocks; i++) {
    for (int j=0; j<64; j++) {
      test_assert_eq(selector(filter, i*64 + j),
                     (i < 2) ? ((j-1)%8 == 0) : (i == 2 && j == 0 ? 0 : j%8 == 0),
                     "i=%d, j=%d", i, j);
    }
  }
  utaf_destroy(filter);
  filter = sel_setup();
  // (2) Shift sels in [32, 64+32] -> [33, 64+33]
  shift_sels(filter, 32, 64+32);
  for (int i=0; i<filter->nslots; i++) {
    test_assert_eq(selector(filter, i),
                   i < 32 ? i%8 == 0 :
                   (i == 32 ? 0 :
                    (i <= 64 + 33) ? (i-1)%8 == 0 : i%8 == 0),
                   "i=%d", i);
  }
  utaf_destroy(filter);
  printf("passed.\n");
}

int main() {
  test_add_block();
  test_add_block_no_clobber();
  test_adapt_loc_1();
  test_adapt_loc_2();
  test_shift_remote_elts();
  test_shift_sels_single_block();
  test_shift_sels_multi_block();
  test_insert_and_query();
  test_insert_and_query_w_repeats();
  test_mixed_insert_and_query_w_repeats();
}
#endif // TEST_UTAF
