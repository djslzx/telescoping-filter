#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <execinfo.h>

#include "murmur3.h"
#include "macros.h"
#include "rsqf.h"
#include "bit_util.h"
#include "set.h"

/**
 * Generate a 64-bit hash for the input word.
 */
static uint64_t rsqf_hash(const RSQF *filter, uint64_t elt) {
  uint64_t buf[2];
  MurmurHash3_x64_128(&elt, 8, filter->seed, buf);
  return buf[0];
}

/**
 * Returns the quotient for a 64-bit fingerprint hash.
 */
static size_t calc_quot(const RSQF* filter, uint64_t hash) {
  return hash & ONES(filter->q);
}

/**
 * Returns the remainder for a 64-bit fingerprint hash.
 */
static rem_t calc_rem(const RSQF* filter, const uint64_t hash) {
  return ((hash >> filter->q) & ONES(filter->r));
}

/* RSQF Helpers */

/**
 * Returns the absolute index of the `rank`-th 1 bit in Q.runends past the start of
 * the block at `block_index`. `rank` indexes from 0.
 *
 * Returns -1 if result is invalid (out of bounds).
 */
static int select_runend(const RSQF* filter, size_t block_index, size_t rank) {
  assert(block_index < filter->nblocks && "block_index out of bounds");

  size_t step;
  size_t loc = block_index * 64;
  while (1) {
    RSQFBlock* b = &filter->blocks[loc / 64];
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
static int rank_select(const RSQF* filter, size_t x) {
  // Exit early if x obviously out of range
  if (x >= filter->nslots) {
    return RANK_SELECT_OVERFLOW;
  }
  size_t block_i = x/64;
  size_t slot_i = x%64;
  RSQFBlock *b = &filter->blocks[block_i];

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
static int first_unused(const RSQF* filter, size_t x) {
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
static void shift_rems_and_runends(RSQF* filter, int a, int b) {
  if (a > b) return;
  for (int i=b; i>=a; i--) {
    remainder(filter, i+1) = remainder(filter, i);
    set_runend_to(filter, i+1, get_runend(filter, i));
  }
  set_runend_to(filter, a, 0);
}

/**
 * Increment all non-negative offsets with targets in [a,b]
 */
static void inc_offsets(RSQF* filter, size_t a, size_t b) {
  assert(a < filter->nslots && b < filter->nslots);
  // Exit early if invalid range
  if (a > b) {
    return;
  }
  // Start i at the first block after b, clamping it so it doesn't go off the end, and work backwards
  size_t start = min(b/64 + 1, filter->nblocks - 1);
  for (int i = start; i>=0; i--) {
    RSQFBlock *block = &filter->blocks[i];
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
static void inc_offsets_for_new_run(RSQF* filter, size_t quot, size_t loc) {
  assert(loc < filter->nslots);
  // Start i at the first block after loc,
  // clamping it so it doesn't go off the end
  size_t start = min(loc/64 + 1, filter->nblocks - 1);
  for (int i=start; i>=0; i--) {
    RSQFBlock *b = &filter->blocks[i];
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

static void add_block(RSQF *filter) {
  filter->blocks = realloc(filter->blocks, (filter->nblocks + 1) * sizeof(RSQFBlock));
  memset(filter->blocks + filter->nblocks, 0, sizeof(RSQFBlock));
  filter->nblocks += 1;
  filter->nslots += 64;
}

/* RSQF */

void rsqf_init(RSQF *filter, size_t n, int seed) {
  filter->seed = seed;
  filter->nelts = 0;
  filter->nblocks = max(1, nearest_pow_of_2(n)/64);
  filter->nslots = filter->nblocks * 64;
  filter->q = log2((double)filter->nslots); // nslots = 2^q
  filter->r = REM_SIZE;
  filter->p = filter->q + filter->r;
  filter->blocks = calloc(filter->nblocks, sizeof(RSQFBlock));
}

void rsqf_destroy(RSQF* filter) {
  free(filter->blocks);
  free(filter);
}

static void raw_insert(RSQF* filter, size_t quot, rem_t rem) {
  assert(quot < filter->nslots);
  filter->nelts++;

  // Find the appropriate runend
  int r = rank_select(filter, quot);
  switch (r) {
    case RANK_SELECT_EMPTY: {
      set_occupied_to(filter, quot, 1);
      set_runend_to(filter, quot, 1);
      remainder(filter, quot) = rem;
      break;
    }
    case RANK_SELECT_OVERFLOW: {
      printf("RSQF failed to find runend (nslots=%lu, quot=(block=%lu, slot=%lu))\n",
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
      // Start a new run or extend an existing one
      if (get_occupied(filter, quot)) {
        // quot occupied: extend an existing run
        inc_offsets(filter, r, r);
        set_runend_to(filter, r, 0);
        set_runend_to(filter, r + 1, 1);
        remainder(filter, r+1) = rem;
      } else {
        // quot unoccupied: start a new run
        inc_offsets_for_new_run(filter, quot, r);
        set_occupied_to(filter, quot, 1);
        set_runend_to(filter, r + 1, 1);
        remainder(filter, r+1) = rem;
      }
    }
  }
}

static int raw_lookup(const RSQF* filter, size_t quot, rem_t rem) {
  if (get_occupied(filter, quot)) {
    int loc = rank_select(filter, quot);
    if (loc == RANK_SELECT_EMPTY || loc == RANK_SELECT_OVERFLOW) {
      return 0;
    }
    do {
      if (remainder(filter, loc) == rem) {
        return 1;
      }
      loc--;
    } while (loc >= quot && !get_runend(filter, loc));
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
int rsqf_lookup(const RSQF *filter, uint64_t elt) {
  uint64_t hash = rsqf_hash(filter, elt);
  size_t quot = calc_quot(filter, hash);
  rem_t rem = calc_rem(filter, hash);
  return raw_lookup(filter, quot, rem);
}

void rsqf_insert(RSQF *filter, uint64_t elt) {
  uint64_t hash = rsqf_hash(filter, elt);
  size_t quot = calc_quot(filter, hash);
  rem_t rem = calc_rem(filter, hash);
  raw_insert(filter, quot, rem);
}

double rsqf_load(RSQF *filter) {
  return (double)filter->nelts/(double)filter->nslots;
}

void rsqf_clear(RSQF* filter) {
  filter->nelts = 0;
  free(filter->blocks);
  filter->blocks = calloc(filter->nblocks, sizeof(RSQFBlock));
}

/* Printing */

void print_rsqf_metadata(RSQF* filter) {
  printf("FILTER METADATA:\n");
  printf("  p=%ld, q=%ld, r=%ld\n",
         filter->p, filter->q, filter->r);
  printf("  nslots=%ld, nblocks=%ld, blocksize=%ld, nelts=%ld\n",
         filter->nslots, filter->nslots/64, sizeof(RSQFBlock), filter->nelts);
  printf("  seed=%d\n", filter->seed);
  printf("  load factor=%f\n", rsqf_load(filter));
}

void print_rsqf_block(RSQF* filter, size_t block_index) {
  assert(0 <= block_index && block_index < filter->nslots/64);
  RSQFBlock block = filter->blocks[block_index];
  printf("BLOCK 0x%lx:\n", block_index);
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
}

void print_rsqf(RSQF* filter) {
  print_rsqf_metadata(filter);
  for (int i=0; i<filter->nblocks; i++) {
    print_rsqf_block(filter, i);
  }
}

// Tests
//#define TEST_RSQF 1
#ifdef TEST_RSQF

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

#define RSQF_SEED 32776517

RSQF *new_rsqf(size_t n) {
  RSQF *filter = malloc(sizeof(RSQF));
  rsqf_init(filter, n, RSQF_SEED);
  return filter;
}

void test_calc_quot() {
  printf("Testing %s...", __FUNCTION__);
  RSQF *filter = new_rsqf(128);
  // q = 7, r = 8
  assert_eq(calc_quot(filter, 0), 0);
  assert_eq(calc_quot(filter, 1), 1);
  assert_eq(calc_quot(filter, 0b111110000000), 0);
  assert_eq(calc_quot(filter, 0b000001111111), 0b1111111);
  assert_eq(calc_quot(filter, 0b000001111111), 0b1111111);
  assert_eq(calc_quot(filter, 0b000001010101), 0b1010101);

  rsqf_destroy(filter);
  printf("passed.\n");
}

void test_calc_rem() {
  printf("Testing %s...", __FUNCTION__);
  RSQF *filter = new_rsqf(128);
  // q = 7, r = 8
  assert_eq(calc_rem(filter, 0), 0);
  assert_eq(calc_rem(filter, 0b1111111), 0);
  assert_eq(calc_rem(filter, 0b11111111), 1);
  assert_eq(calc_rem(filter, 0b111111111), 0b11);
  assert_eq(calc_rem(filter, 0b111110000000), 0b11111);
  assert_eq(calc_rem(filter, 0b101010000000), 0b10101);
  assert_eq(calc_rem(filter, 0b1111111110000000), 0b11111111);
  rsqf_destroy(filter);
  printf("passed.\n");
}

void test_select_runend_empty_filter() {
  printf("Testing %s...", __FUNCTION__);
  RSQF *filter = new_rsqf(128);

  assert_eq(select_runend(filter, 0, 0), -1);

  rsqf_destroy(filter);
  printf("passed.\n");
}

void test_select_runend_one_run() {
  printf("Testing %s...", __FUNCTION__);
  RSQF *filter = new_rsqf(128);

  RSQFBlock* b = &filter->blocks[0];
  b->occupieds = 1;
  b->runends = 1;
  b->offset = 0;
  assert_eq(select_runend(filter, 0, 0), 0);

  rsqf_destroy(filter);
  printf("passed.\n");
}

void test_select_runend_mult_runs() {
  printf("Testing %s...", __FUNCTION__);
  RSQF *filter = new_rsqf(128);

  RSQFBlock* b = &filter->blocks[0];
  // note: read bit rep backwards
  b->occupieds = 0b01101;
  b->runends   = 0b11010;
  b->offset    = 1;
  assert_eq(select_runend(filter, 0, 0), 1);
  assert_eq(select_runend(filter, 0, 1), 3);
  assert_eq(select_runend(filter, 0, 2), 4);
  for (int i=3; i<64; i++) {
    assert_eq(select_runend(filter, 0, i), -1);
  }

  rsqf_destroy(filter);
  printf("passed.\n");
}

void test_select_runend_mult_blocks_spanning_run() {
  printf("Testing %s...", __FUNCTION__);
  RSQF *filter = new_rsqf(64 * 3);

  // Filter with run start in block 0 and ending in block 1
  RSQFBlock* b0 = &filter->blocks[0];
  SET(b0->occupieds, 0);

  RSQFBlock* b1 = &filter->blocks[1];
  SET(b1->runends, 0);
  b1->offset = 0;

  assert_eq(select_runend(filter, 0, 0), 64);
  assert_eq(select_runend(filter, 1, 0), 64);
  assert_eq(select_runend(filter, 2, 0), -1);

  rsqf_destroy(filter);
  printf("passed.\n");
}

void test_select_runend_mult_blocks_two_runs() {
  printf("Testing %s...", __FUNCTION__);
  RSQF *filter = new_rsqf(64 * 3);

  // A run starts in b0 and ends in b1, making b1 have nonzero offset
  RSQFBlock* b0 = &filter->blocks[0];
  b0->occupieds = 0b11;
  b0->runends   = 0b01;
  b0->offset    = 0;

  RSQFBlock* b1 = &filter->blocks[1];
  b1->occupieds = 0b10;
  b1->runends   = 0b11;
  b1->offset    = 0;

  assert_eq(select_runend(filter, 0, 0), 0);
  assert_eq(select_runend(filter, 0, 1), 64);
  assert_eq(select_runend(filter, 0, 2), 65);
  assert_eq(select_runend(filter, 0, 3), -1);
  assert_eq(select_runend(filter, 1, 0), 64);
  assert_eq(select_runend(filter, 1, 1), 65);
  assert_eq(select_runend(filter, 1, 2), -1);

  rsqf_destroy(filter);
  printf("passed.\n");
}

void test_rank_select_single_block_empty() {
  printf("Testing %s...", __FUNCTION__);
  RSQF *filter = new_rsqf(64);

  for (int i=0; i<64; i++) {
    assert_eq(rank_select(filter, i), RANK_SELECT_EMPTY);
  }

  rsqf_destroy(filter);
  printf("passed.\n");
}

void test_rank_select_single_block_singleton() {
  printf("Testing %s...", __FUNCTION__);
  RSQF *filter = new_rsqf(64);

  // Filter with one singleton run
  for (int i=0; i<64; i++) {
    // Add singleton run at i
    set_occupied_to(filter, i, 1);
    set_runend_to(filter, i, 1);

    assert_eq(rank_select(filter, i), i);
    for (int j=i+1; j<64; j++) {
      assert_eq(rank_select(filter, j), RANK_SELECT_EMPTY);
    }

    // Clear singleton run at i
    set_occupied_to(filter, i, 0);
    set_runend_to(filter, i, 0);
  }

  rsqf_destroy(filter);
  printf("passed.\n");
}

void test_rank_select_single_block_two_runs() {
  printf("Testing %s...", __FUNCTION__);
  RSQF *filter = new_rsqf(64);

  RSQFBlock* b = &filter->blocks[0];
  b->occupieds = 0b101001;
  b->runends   = 0b110010;
  b->offset    = 1;
  assert_eq(rank_select(filter, 0), 1);
  assert_eq(rank_select(filter, 1), 1);
  assert_eq(rank_select(filter, 2), RANK_SELECT_EMPTY); // empty b/c prev runend is before slot 2
  assert_eq(rank_select(filter, 3), 4);
  assert_eq(rank_select(filter, 4), 4);
  assert_eq(rank_select(filter, 5), 5);
  for (int i=6; i<64; i++) {
    assert_eq(rank_select(filter, i), RANK_SELECT_EMPTY);
  }

  rsqf_destroy(filter);
  printf("passed.\n");
}

void test_rank_select_multi_block_1() {
  printf("Testing %s...", __FUNCTION__);
  RSQF *filter = new_rsqf(64 * 3);

  // Filter with run starting in block 0 and ending in block 1
  set_occupied_to(filter, 0, 1);
  filter->blocks[0].offset = 64;

  set_runend_to(filter, 64, 1);
  filter->blocks[1].offset = 0;

  for (int i=0; i<=64; i++) {
    assert_eq(rank_select(filter, i), 64);
  }
  for (int i=65; i<filter->nslots; i++) {
    assert_eq(rank_select(filter, i), RANK_SELECT_EMPTY);
  }

  rsqf_destroy(filter);
  printf("passed.\n");
}

void test_rank_select_multi_block_2() {
  printf("Testing %s...", __FUNCTION__);
  RSQF *filter = new_rsqf(64 * 3);

  // Run 1: [0, [0]]
  // Run 2: [1, [1,64]]
  // Run 3: [65, [65, 68]]
  // Run 4: [66, [69, 130]]
  set_occupied(filter, 0); // start run 1
  set_runend(filter, 0); // end run 1
  filter->blocks[0].offset = 0;
  set_occupied(filter, 1); // start run 2
  set_runend(filter, 64); // end run 2
  filter->blocks[1].offset = 0;
  set_occupied(filter, 65); // start run 3
  set_runend(filter, 68); // end run 3
  set_occupied(filter, 66); // start run 4
  set_runend(filter, 130); // end run 4
  filter->blocks[2].offset = 2;

  assert_eq(rank_select(filter, 0), 0);
  for (int i=1; i<=64; i++) {
    assert_eq(rank_select(filter, i), 64);
  }
  assert_eq(rank_select(filter, 65), 68);
  for (int i=66; i<=130; i++) {
    assert_eq(rank_select(filter, i), 130);
  }
  for (int i=131; i<filter->nslots; i++) {
    assert_eq(rank_select(filter, i), RANK_SELECT_EMPTY);
  }

  rsqf_destroy(filter);
  printf("passed.\n");
}

void test_first_unused_empty() {
  printf("Testing %s...", __FUNCTION__);
  RSQF *filter = new_rsqf(64 * 3);

  for (int i=0; i<filter->nslots; i++) {
    assert_eq(first_unused(filter, i), i);
  }

  rsqf_destroy(filter);
  printf("passed.\n");
}

void test_first_unused_single() {
  printf("Testing %s...", __FUNCTION__);
  int nslots = 128;
  // Filters with singleton runs at k;
  // first unused should be slot k+1 at x=k, slot x elsewhere
  // unless x = k = nslots-1, in which case there are no unused slots
  for (int k=0; k<nslots; k++) {
    RSQF *filter = new_rsqf(nslots);
    set_occupied(filter, k);
    set_runend(filter, k);
    for (int i=0; i<nslots; i++) {
      int u = first_unused(filter, i);
      if (i == k) {
        if (i == nslots-1) assert_eq(u, NO_UNUSED);
        else assert_eq(u, k+1);
      } else {
        assert_eq(u, i);
      }
    }
    rsqf_destroy(filter);
  }
  printf("passed.\n");
}

// Simulate inserting a single run in [a,b] into a filter of 128 slots
// by modifying metadata bits only
void insert_run(RSQF* filter, size_t a, size_t b) {
  assert(a < 128 && b < 128);
  set_occupied(filter, a);
  set_runend(filter, b);
  if (a==0) {
    filter->blocks[0].offset = b;
  }
  if (a == 64 || (a < 64 && b >= 64)) {
    filter->blocks[1].offset = b-64;
  }
}

void test_first_unused_one_run() {
  printf("Testing %s...", __FUNCTION__);

  // Filter with one run from a to b
  // First unused for run [a,b] at slot x should be
  // b+1 if x in [a,b], else x
  // Except when b=nslots-1, when result should be -1
  RSQF *filter = new_rsqf(128);
  for (int a=0; a<filter->nslots; a++) {
    for (int b=a; b<filter->nslots; b++) {
      // Setup filter with one run in [a,b]
      insert_run(filter, a, b);
      // Test
      for (int i=0; i<filter->nslots; i++) {
        int u = first_unused(filter, i);
        if (i < a || b < i) {
          test_assert_eq(u, i, "run [%d, %d], u=%d, i=%d", a, b, u, i);
        } else if (b == filter->nslots-1){
          test_assert_eq(u, NO_UNUSED, "u=%d, i=%d", u, i);
        } else {
          test_assert_eq(u, b+1, "u=%d, i=%d", u, i);
        }
      }
      // Clear run
      rsqf_clear(filter);
    }
  }
  rsqf_destroy(filter);
  printf("passed.\n");
}

/// Helper for test_first_unused_two_runs
/// Test first_unused for runs [a,b] and [c,d] nonoverlapping
void first_unused_two_runs(size_t a, size_t b, size_t c, size_t d) {
  assert(a < b && b < c && c < d);
  size_t nslots = 128;
  RSQF* filter = new_rsqf(nslots);
  insert_run(filter, a, b);
  insert_run(filter, c, d);
  for (int i=0; i<nslots; i++) {
    int u = first_unused(filter, i);
    if (i < a || i > d || (i > b && i < c)) {
      // i isn't in either interval
      assert_eq(u, i);
    } else {
      if (i >= a && i <= b) {
        // i is in the first interval
        if (c > b+1) {
          // If there's a gap between the two intervals
          assert_eq(u, b+1);
        } else if (d < nslots - 1) {
          // If the two intervals are connected and
          // there's a gap after the second
          assert_eq(u, d+1);
        } else {
          assert_eq(u, NO_UNUSED);
        }
      } else if (i >= c && i <= d && d < nslots - 1) {
        // If i is in the second interval and
        // there's a gap after the second interval
        assert_eq(u, d+1);
      } else {
        assert_eq(u, NO_UNUSED);
      }
    }
  }
  rsqf_destroy(filter);
}

void test_first_unused_two_runs() {
  printf("Testing %s...", __FUNCTION__);
  first_unused_two_runs(0, 63, 65, 127); // center gap
  first_unused_two_runs(0, 64, 66, 127); // center gap
  first_unused_two_runs(1, 63, 64, 127); // left gap
  first_unused_two_runs(1, 64, 65, 127); // left gap
  first_unused_two_runs(0, 63, 64, 126); // right gap
  first_unused_two_runs(0, 64, 65, 126); // right gap
  first_unused_two_runs(1, 63, 65, 127); // left, center gap
  first_unused_two_runs(1, 64, 66, 127); // left, center gap
  first_unused_two_runs(1, 63, 64, 126); // left, right gap
  first_unused_two_runs(1, 64, 65, 126); // left, right gap
  first_unused_two_runs(0, 63, 65, 126); // center, right gap
  first_unused_two_runs(0, 64, 66, 126); // center, right gap
  first_unused_two_runs(1, 63, 65, 126); // left, center, right gap
  first_unused_two_runs(1, 64, 66, 126); // left, center, right gap
  printf("passed.\n");
}

void test_lookup_empty() {
  printf("Testing %s...", __FUNCTION__);
  RSQF *filter = new_rsqf(64 * 3);

  char* str = malloc(20);
  for (int i=0; i<filter->nslots; i++) {
    assert_eq(rsqf_lookup(filter, (uint64_t)i), 0);
  }
  free(str);
  rsqf_destroy(filter);
  printf("passed.\n");
}

void test_lookup_singleton() {
  printf("Testing %s...", __FUNCTION__);
  RSQF *filter = new_rsqf(64 * 3);

  uint64_t elt = 01010101;

  uint64_t hash = rsqf_hash(filter, elt);
  size_t quot = calc_quot(filter, hash);
  rem_t rem = calc_rem(filter, hash);

  RSQFBlock* b = &filter->blocks[quot / 64];
  b->remainders[quot%64] = rem;
  SET(b->occupieds, quot%64);
  SET(b->runends, quot%64);
  assert_eq(rsqf_lookup(filter, elt), 1);

  rsqf_destroy(filter);
  printf("passed.\n");
}

void test_lookup_multi_singletons() {
  printf("Testing %s...", __FUNCTION__);
  RSQF *filter = new_rsqf(64 * 3);
  // I checked that these words' quot/rems don't conflict,
  // so these will have their own singleton runs
  uint64_t elts[] = {0, 1, 2};
  int nelts = 3;

  // Manually insert words
  for (int i=0; i<nelts; i++) {
    // Get quot, rem
    uint64_t elt = elts[i];
    uint64_t hash = rsqf_hash(filter, elt);
    size_t quot = calc_quot(filter, hash);
    rem_t rem = calc_rem(filter, hash);

    RSQFBlock* b = &filter->blocks[quot / 64];
    b->remainders[quot%64] = rem;
    SET(b->occupieds, quot%64);
    SET(b->runends, quot%64);
    assert_eq(rsqf_lookup(filter, elt), 1);
  }

  // Check that words are contained
  for (int i=0; i<nelts; i++) {
    uint64_t elt = elts[i];
    assert_eq(rsqf_lookup(filter, elt), 1);
  }

  rsqf_destroy(filter);
  printf("passed.\n");
}

void test_shift_rems_and_runends() {
  printf("Testing %s...", __FUNCTION__);
  RSQF *filter = new_rsqf(128);
  for (int i=0; i<filter->nslots; i++) {
    remainder(filter, i) = i;
    set_runend_to(filter, i, i%3 == 0);
  }
  // Shift [0, nslots-2] to [1, nslots-1] and clear [0]
  shift_rems_and_runends(filter, 0, (int)(filter->nslots-2));
  assert_eq(remainder(filter, 0), 0);
  assert_eq(get_runend(filter, 0), 0);
  for (int i=1; i<filter->nslots; i++) {
    assert_eq(remainder(filter, i), i-1);
    assert_eq(!!get_runend(filter, i), (i-1)%3 == 0);
  }
  rsqf_destroy(filter);
  printf("passed.\n");
}

RSQF* offset_state_init() {
  RSQF *filter = new_rsqf(64 * 7);
  RSQFBlock* b = filter->blocks;
  // Test cases:
  // - negative offset
  // - zero offset
  //   - elt from prior run at 0
  //   - singleton run at 0
  // - positive offset
  //   - end of prior run
  //   - end of first run in block
  // Run in b0: [0:(0,0)]: zero offset, singleton run
  set_occupied(filter, 0);
  set_runend(filter, 0);
  b[0].offset = 0;
  // Run in b0,b1: [63:(63,64)]: zero offset, end of prior run
  set_occupied(filter, 63);
  set_runend(filter, 64);
  //b[1]->offset = 0;
  b[1].offset = 0;
  // Run in b1: [67: (67,72)]
  set_occupied(filter, 67);
  set_runend(filter, 72);
  // Run in b1: [68: (73,73)]
  set_occupied(filter, 68);
  set_runend(filter, 73);
  // Run in b1, b2: [69: (74,129)]: positive offset, end of prior run
  set_occupied(filter, 69);
  set_runend(filter, 129);
  // Run in b1,b2: [80: (130,130)]: positive offset, end of prior run
  set_occupied(filter, 80);
  set_runend(filter, 130);
  b[2].offset = 2;
  // Run in b3: [192: (192, 194)]: positive offset, end of run at start of block
  set_occupied(filter, 192);
  set_runend(filter, 194);
  b[3].offset = 2;
  // Negative offset for b4
  b[4].offset = 0;
  // Run in b4,b6: [260: (260, 390)]: offset > 64
  set_occupied(filter, 260);
  set_runend(filter, 390);
  b[5].offset = 70;
  b[6].offset = 6;

  return filter;
}

void test_inc_nonneg_offsets_full() {
  printf("Testing %s...", __FUNCTION__);
  RSQF* filter = offset_state_init();
  // Inc all offsets [0, n-2] -> [1, n-1]
  inc_offsets(filter, 0, filter->nslots-1);
  RSQFBlock* b = filter->blocks;
  assert_eq(b[0].offset, 1);
  assert_eq(b[1].offset, 1);
  assert_eq(b[2].offset, 3);
  assert_eq(b[3].offset, 3);
  assert_eq(b[4].offset, 0);
  assert_eq(b[5].offset, 71);
  assert_eq(b[6].offset, 7);
  rsqf_destroy(filter);
  printf("passed.\n");
}

void inc_and_check_offsets_unchanged(RSQF* filter, size_t start, size_t end) {
  inc_offsets(filter, start, end);
  RSQFBlock* b = filter->blocks;
  assert_eq(b[0].offset, 0);
  assert_eq(b[1].offset, 0);
  assert_eq(b[2].offset, 2);
  assert_eq(b[3].offset, 2);
  assert_eq(b[4].offset, 0);
  assert_eq(b[5].offset, 70);
  assert_eq(b[6].offset, 6);
}

void test_inc_nonneg_offsets_untargeted() {
  printf("Testing %s...", __FUNCTION__);
  RSQF* filter = offset_state_init();
  // Inc ranges that nothing is pointing to
  inc_and_check_offsets_unchanged(filter, 1, 63);
  inc_and_check_offsets_unchanged(filter, 65, 129);
  inc_and_check_offsets_unchanged(filter, 131, 193);
  inc_and_check_offsets_unchanged(filter, 195, 389);
  inc_and_check_offsets_unchanged(filter, 391, 447);
  rsqf_destroy(filter);
  printf("passed.\n");
}

void inc_and_check_offsets_match(size_t target, size_t o0, size_t o1, size_t o2,
                                 size_t o3, size_t o4, size_t o5, size_t o6) {
  RSQF* filter = offset_state_init();
  inc_offsets(filter, target, target);
  RSQFBlock* b = filter->blocks;
  test_assert_eq(b[0].offset, o0, "offset=%lu, o=%lu", b[0].offset, o0);
  assert_eq(b[1].offset, o1);
  assert_eq(b[2].offset, o2);
  assert_eq(b[3].offset, o3);
  assert_eq(b[4].offset, o4);
  assert_eq(b[5].offset, o5);
  assert_eq(b[6].offset, o6);
  rsqf_destroy(filter);
}

void test_inc_nonneg_offsets_targeted() {
  printf("Testing %s...", __FUNCTION__);
  // Inc ranges that things are pointing to
  inc_and_check_offsets_match(0, 1, 0, 2, 2, 0, 70, 6);
  inc_and_check_offsets_match(64, 0, 1, 2, 2, 0, 70, 6);
  inc_and_check_offsets_match(130, 0, 0, 3, 2, 0, 70, 6);
  inc_and_check_offsets_match(194, 0, 0, 2, 3, 0, 70, 6);
  inc_and_check_offsets_match(390, 0, 0, 2, 2, 0, 71, 7);
  printf("passed.\n");
}

void test_inc_offsets_negative_target() {
  printf("Testing %s...", __FUNCTION__);
  RSQF *filter = new_rsqf(64);
  RSQFBlock* b = filter->blocks;

  // Check that target doesn't go negative (block_i ==0 and block[0] unoccupied)
  // Negative -> 0
  inc_offsets(filter, 0, 0);
  assert_eq(b[0].offset, 0);

  // 0 -> 1
  set_occupied(filter, 0);
  set_runend(filter, 0);
  inc_offsets(filter, 0, 0);
  assert_eq(b[0].offset, 1);

  // 1 -> 2
  inc_offsets(filter, 1, 1);
  assert_eq(b[0].offset, 2);

  // Check for multiblock case
  rsqf_destroy(filter);
  filter = new_rsqf(64 * 5);
  inc_offsets(filter, 0, filter->nslots-1);
  assert_eq(filter->blocks[0].offset, 0);

  rsqf_destroy(filter);
  printf("passed.\n");
}

void test_inc_offsets_zero_offset() {
  printf("Testing %s...", __FUNCTION__);
  RSQF *filter = new_rsqf(128);
  RSQFBlock* b = filter->blocks;
  set_occupied(filter, 1);
  set_runend(filter, 64);
  b[0].offset = 0;
  b[1].offset = 0;
  inc_offsets(filter, 64, 64);
  assert_eq(b[0].offset, 0); // shouldn't increment, negative
  assert_eq(b[1].offset, 1); // should increment, zero
  rsqf_destroy(filter);
  printf("passed.\n");
}

void test_add_block() {
  printf("Testing %s...", __FUNCTION__);
  RSQF *filter = new_rsqf(64 * 2);
  assert_eq(filter->nslots, 64 * 2);
  assert_eq(filter->nblocks, 2);
  add_block(filter);
  assert_eq(filter->nslots, 64 * 3);
  assert_eq(filter->nblocks, 3);
  RSQFBlock b = filter->blocks[2];
  assert_eq(b.occupieds, 0);
  assert_eq(b.runends, 0);
  assert_eq(b.offset, 0);
  for (int i=0; i<64; i++) {
    assert_eq(b.remainders[i], 0);
  }
  rsqf_destroy(filter);
  printf("passed.\n");
}

int is_pow_of_2(size_t x) {
  size_t mask = x==0 ? 0 : x-1;
  return (x & mask) == 0;
}

/// Insert new runs that don't overlap with anything
void test_raw_insert_new_run() {
  printf("Testing %s...", __FUNCTION__);
  RSQF *filter = new_rsqf(128);
  // Insert runs at powers of 2
  for (int i=0; i<filter->nslots; i++) {
    if (is_pow_of_2(i)) {
      raw_insert(filter, i, i);
    }
  }
  // Check offsets
  assert_eq(filter->blocks[0].offset, 0);
  assert_eq(filter->blocks[1].offset, 0);
  // Check occupieds/runends
  for (int i=0; i<filter->nslots; i++) {
    assert_eq(!!get_occupied(filter, i), is_pow_of_2(i));
    assert_eq(!!get_runend(filter, i), is_pow_of_2(i));
    assert_eq(remainder(filter, i), is_pow_of_2(i) ? i : 0);
  }
  rsqf_destroy(filter);
  printf("passed.\n");
}

/// Initialize a filter of 192 slots with one run 0:[0, 130]
/// where all remainders equal their indices
RSQF* one_long_run() {
  RSQF* filter = new_rsqf(64 * 3);
  RSQFBlock* b = filter->blocks;
  set_occupied(filter, 0);
  set_runend(filter, 130);
  b[0].offset = 130;
  b[1].offset = 130-64;
  b[2].offset = 130-128;
  for (int i=0; i<=130; i++) {
    remainder(filter, i) = i;
  }
  return filter;
}

/// Insert a new run that overlaps with an existing run and shifts
/// multiple offsets
void test_raw_insert_overlapping_run() {
  printf("Testing %s...", __FUNCTION__);
  RSQF* filter;

  // Insert after the run ends, at slot 131
  filter = one_long_run();
  raw_insert(filter, 131, 0xff);
  for (int i=0; i<filter->nslots; i++) {
    assert_eq(!!get_occupied(filter, i), i==0 || i==131);
    assert_eq(!!get_runend(filter, i), i==130 || i==131);
    assert_eq(remainder(filter, i),
              i < 131 ? i : (i == 131 ? 0xff : 0));
  }
  rsqf_destroy(filter);
  // Insert into a new run with quot st 0 < quot < 130
  filter = one_long_run();
  raw_insert(filter, 10, 0xff);
  for (int i=0; i<filter->nslots; i++) {
    assert_eq(!!get_occupied(filter, i), i==0 || i==10);
    assert_eq(!!get_runend(filter, i), i==130 || i==131);
    assert_eq(remainder(filter, i),
              i < 131 ? i : (i == 131 ? 0xff : 0));
  }
  rsqf_destroy(filter);
  // Extend the run (insert with quot=0)
  filter = one_long_run();
  raw_insert(filter, 0, 131);
  for (int i=0; i<filter->nslots; i++) {
    assert_eq(!!get_occupied(filter, i), i==0);
    assert_eq(!!get_runend(filter, i), i==131);
    assert_eq(remainder(filter, i), i <= 131 ? i : 0);
  }
  rsqf_destroy(filter);
  printf("passed.\n");
}

/// Try adding an element into a full filter
void test_raw_insert_extend() {
  printf("Testing %s...", __FUNCTION__);
  RSQF *filter = new_rsqf(128);
  // Place 128 consecutive singleton runs
  for (int i=0; i<filter->nslots; i++) {
    set_occupied(filter, i);
    set_runend(filter, i);
    remainder(filter, i) = (rem_t)i;
  }
  // Inserting another remainder for quot=0 should shift everything else over
  raw_insert(filter, 0, 0xff);
  for (int i=0; i<filter->nslots; i++) {
    assert_eq(!!get_occupied(filter, i), i<128);
    assert_eq(!!get_runend(filter, i), i>0 && i<=128);
    assert_eq(remainder(filter, i),
              i==0 ? 0 : (i==1 ? 0xff : (i <= 128 ? i-1 : 0)));
  }
  rsqf_destroy(filter);
  printf("passed.\n");
}

/// Check that zero offsets are incremented
void test_raw_insert_zero_offset() {
  printf("Testing %s...", __FUNCTION__);
  RSQF *filter = new_rsqf(128);
  RSQFBlock* b = filter->blocks;
  // Run from 1 to 64, insert 2 elements at 0 to shift runends
  set_occupied(filter, 1);
  set_runend(filter, 64);
  for (int i=1; i<=64; i++) {
    remainder(filter, i) = 0xf;
  }
  raw_insert(filter, 0, 0xa);
  raw_insert(filter, 0, 0xb);
  assert(get_occupied(filter, 0)); // insert sets quotient
  assert(get_occupied(filter, 1)); // insert doesn't unset pre-existing quotient
  assert(get_runend(filter, 1)); // new runend for quotient 0
  assert(!get_runend(filter, 64)); // runend at 64 moved
  assert(get_runend(filter, 65)); // runend at 64 moved
  assert_eq(b[1].offset, 1); // offset shifted
  rsqf_destroy(filter);
  printf("passed.\n");
}

void test_insert_repeated() {
  printf("Testing %s...", __FUNCTION__);
  int n = 1 << 10;
  RSQF *filter = new_rsqf(n);
  for (int i=0; i<n; i++) {
    rsqf_insert(filter, 1);
    assert_eq(rsqf_lookup(filter, 1), 1);
  }
  rsqf_destroy(filter);
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
  RSQF* filter = new_rsqf(s);

  // Generate query set
  srand(RSQF_SEED);
  size_t nset = 1.5*s;
  Setnode* set = calloc(nset, sizeof(set[0]));
  char* str = malloc(20);
  for (int i=0; i<s; i++) {
    uint64_t elt = rand() % a;
    sprintf(str, "%lu", elt);
    int len = (int)strlen(str);
    set_insert(str, len, 0, set, nset);
    rsqf_insert(filter, elt);
  }
  // Query [0, a] and ensure that all items in the set return true
  int fps = 0;
  for (int i=0; i<a; i++) {
    uint64_t elt = (uint64_t)i;
    sprintf(str, "%lu", elt);
    int len = (int)strlen(str);
    int in_set = set_lookup(str, len, set, nset);
    int in_rsqf = rsqf_lookup(filter, elt);
    if (in_set) {
      if (!in_rsqf) {
        uint64_t hash = rsqf_hash(filter, elt);
        size_t quot = calc_quot(filter, hash);
        size_t rem = calc_rem(filter, hash);
        printf("False negative: set contains %s, but filter doesn't\n"
               "quot=%lu (block_i=%lu, slot_i=%lu), rem=%zu\n",
               str, quot, quot/64, quot%64, rem);
        print_rsqf_block(filter, quot/64);
        exit(1);
      }
    } else {
      fps += in_rsqf;
    }
  }
  free(str);
  printf("passed. ");
  printf("FPR: %f\n", (double)fps/a);
  print_rsqf_metadata(filter);
  rsqf_destroy(filter);
  set_deallocate(set, nset);
}

void test_template() {
  printf("Testing %s...", __FUNCTION__);
  RSQF *filter = new_rsqf(64 * 3);

  printf("[unimplemented]");

  rsqf_destroy(filter);
  printf("passed.\n");
}

int main() {
  test_calc_quot();
  test_calc_rem();
  test_select_runend_empty_filter();
  test_select_runend_one_run();
  test_select_runend_mult_runs();
  test_select_runend_mult_blocks_spanning_run();
  test_select_runend_mult_blocks_two_runs();
  test_rank_select_single_block_empty();
  test_rank_select_single_block_singleton();
  test_rank_select_single_block_two_runs();
  test_rank_select_multi_block_1();
  test_rank_select_multi_block_2();
  test_first_unused_empty();
  test_first_unused_single();
  test_first_unused_one_run();
  test_first_unused_two_runs();
  test_lookup_empty();
  test_lookup_singleton();
  test_lookup_multi_singletons();
  test_shift_rems_and_runends();
  test_inc_nonneg_offsets_full();
  test_inc_nonneg_offsets_untargeted();
  test_inc_nonneg_offsets_targeted();
  test_inc_offsets_negative_target();
  test_inc_offsets_zero_offset();
  test_add_block();
  test_raw_insert_new_run();
  test_raw_insert_overlapping_run();
  test_raw_insert_extend();
  test_raw_insert_zero_offset();
  test_insert_repeated();
  test_insert_and_query();
}
#endif // TEST_RSQFv
