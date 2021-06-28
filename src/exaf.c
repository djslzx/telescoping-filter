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
#include "exaf.h"
#include "bit_util.h"
#include "set.h"

/**
 * Generate a hash for the input word.
 * Returns the full 128 bit murmurhash.
 * Expects a two-slot array of uint64_t's.
 */
static uint64_t exaf_hash(const ExAF *filter, elt_t elt) {
  uint64_t buf[2];
  MurmurHash3_x64_128(&elt, 8, filter->seed, buf);
  return buf[0];
}

/**
 * Returns the quotient for a 64-bit fingerprint hash.
 */
static size_t calc_quot(const ExAF* filter, uint64_t hash) {
  return hash & ONES(filter->q);
}

/**
 * Returns the remainder for h.
 */
static rem_t calc_rem(const ExAF* filter, const uint64_t hash) {
  return ((hash >> filter->q) & ONES(filter->r));
}

/**
 * Gets `len` bits after the quotient and remainder in `hash`.
 */
static uint64_t calc_ext_bits(const ExAF* filter, uint64_t hash, size_t len) {
  int fp_len = (int)(filter->q + filter->r);
  assert(len > 0 && fp_len + len <= 64);
  return (hash >> fp_len) & ONES(len);
}

/* ExAF Helpers */

/**
 * @return The extension code at block `block_i`, padded with zeros as a uint64_t.
 */
static uint64_t get_ext_code(const ExAF* filter, size_t block_i) {
  uint64_t code = 0;
  memcpy(&code, filter->blocks[block_i].ext_code, EXT_CODE_BYTES);
  return code;
}

/**
 * Set the extension arithmetic code at the `block_i`-th block to the first `CODE_BYTES` bits of `code`.
 */
static void set_ext_code(ExAF* filter, size_t block_i, uint64_t code) {
  memcpy(filter->blocks[block_i].ext_code, &code, EXT_CODE_BYTES);
}

/**
 * Compute the shortest extension from the member's hash that doesn't
 * conflict with the nonmember's hash.
 * @return the length of the extension
 */
static int shortest_diff_ext(ExAF* filter, uint64_t member_hash, uint64_t non_member_hash, Ext* out) {
  uint64_t a = member_hash >> (filter->q + filter->r);
  uint64_t b = non_member_hash >> (filter->q + filter->r);
  if (a == b) {
    return 0;
  } else {
    out->len = tzcnt(a ^ b) + 1;
    out->bits = a & ONES(out->len);
    return out->len;
  }
}

static int ext_matches_hash(ExAF* filter, Ext* ext, uint64_t hash) {
  if (ext->len == 0) {
    return 1;
  } else {
    return calc_ext_bits(filter, hash, ext->len) == ext->bits;
  }
}

/**
 * Returns the absolute index of the `rank`-th 1 bit in Q.runends past the start of
 * the block at `block_index`. `rank` indexes from 0.
 *
 * Returns -1 if result is invalid (out of bounds).
 */
static int select_runend(const ExAF* filter, size_t block_index, size_t rank) {
  assert(block_index < filter->nblocks && "block_index out of bounds");

  size_t step;
  size_t loc = block_index * 64;
  while (1) {
    ExAFBlock* b = &filter->blocks[loc / 64];
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
static int rank_select(const ExAF* filter, size_t x) {
  // Exit early if x obviously out of range
  if (x >= filter->nslots) {
    return RANK_SELECT_OVERFLOW;
  }
  size_t block_i = x/64;
  size_t slot_i = x%64;
  ExAFBlock *b = &filter->blocks[block_i];

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
static int first_unused(const ExAF* filter, size_t x) {
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
static void shift_rems_and_runends(ExAF* filter, int a, int b) {
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
static void shift_remote_elts(ExAF* filter, int a, int b) {
  if (a > b) return;
  for (int i=b; i>=a; i--) {
    filter->remote[i+1] = filter->remote[i];
  }
  filter->remote[a] = 0;
}

/**
 * Helper for `shift_exts`.  Shifts extensions in `[0, b]` in a single block.
 */
static void shift_block_exts(ExAF *filter, int block_i, Ext exts[64], const Ext prev_exts[64], int b) {
  uint64_t code;
  for (int i=b; i > 0; i--) {
    exts[i] = exts[i-1];
  }
  exts[0] = prev_exts[63];
  if (encode_ext(exts, &code) == -1) {
    code = 0;
  }
  set_ext_code(filter, block_i, code);
}

static void inline swap_ptrs(Ext **a, Ext **b) {
  Ext *tmp = *a;
  *a = *b;
  *b = tmp;
}

/**
 * Shift the remainder extensions in [a,b] forward by 1
 */
static void shift_exts(ExAF* filter, int a, int b) {
  if (a > b) return;
  uint64_t code;
  if (a/64 == (b+1)/64) {
    // a and b+1 in the same block
    Ext exts[64];
    decode_ext(get_ext_code(filter, a/64), exts);
    for (int i = (b+1)%64; i > a%64; i--) {
      exts[i] = exts[i-1];
    }
    exts[a%64].bits = 0;
    exts[a%64].len = 0;
    if (encode_ext(exts, &code) == -1) {
      code = 0;
    }
    set_ext_code(filter, a/64, code);
  } else {
    // a and b+1 in different blocks
    Ext* exts = malloc(64 * sizeof(Ext));
    Ext* prev_exts = malloc(64 * sizeof(Ext));
    // (1) last block
    int block_i = (b+1)/64;
    decode_ext(get_ext_code(filter, block_i), exts);
    decode_ext(get_ext_code(filter, block_i - 1), prev_exts);
    shift_block_exts(filter, block_i, exts, prev_exts, (b + 1) % 64);
    swap_ptrs(&exts, &prev_exts);
    // (2) middle blocks
    for (block_i--; block_i > a/64; block_i--) {
      decode_ext(get_ext_code(filter, block_i-1), prev_exts);
      shift_block_exts(filter, block_i, exts, prev_exts, 63);
      swap_ptrs(&exts, &prev_exts);
    }
    // (3) first block
    for (int i=63; i>a%64; i--) {
      exts[i] = exts[i-1];
    }
    exts[a%64].bits = 0;
    exts[a%64].len = 0;
    if (encode_ext(exts, &code) == -1) {
      code = 0;
    }
    set_ext_code(filter, a/64, code);
    free(exts);
    free(prev_exts);
  }
}

/**
 * Increment all non-negative offsets with targets in [a,b]
 */
static void inc_offsets(ExAF* filter, size_t a, size_t b) {
  assert(a < filter->nslots && b < filter->nslots);
  // Exit early if invalid range
  if (a > b) {
    return;
  }
  // Start i at the first block after b, clamping it so it doesn't go off the end, and work backwards
  size_t start = min(b/64 + 1, filter->nblocks - 1);
  for (int i = start; i>=0; i--) {
    ExAFBlock *block = &filter->blocks[i];
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
static void inc_offsets_for_new_run(ExAF* filter, size_t quot, size_t loc) {
  assert(loc < filter->nslots);
  // Start i at the first block after loc,
  // clamping it so it doesn't go off the end
  size_t start = min(loc/64 + 1, filter->nblocks - 1);
  for (int i=start; i>=0; i--) {
    ExAFBlock *b = &filter->blocks[i];
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

static void add_block(ExAF *filter) {
  // Add block to new_blocks
  ExAFBlock *new_blocks = realloc(filter->blocks, (filter->nblocks + 1) * sizeof(ExAFBlock));
  if (new_blocks == NULL) {
    printf("add_block failed to realloc new blocks\n");
    exit(1);
  }
  filter->blocks = new_blocks;
  memset(filter->blocks + filter->nblocks, 0, sizeof(ExAFBlock));

  // Reallocate remote rep
  elt_t *new_remote = realloc(filter->remote,(filter->nslots + 64) * sizeof(elt_t));
  if (new_remote == NULL) {
    printf("add_block failed to realloc new remote rep\n");
    exit(1);
  }
  filter->remote = new_remote;
  memset(filter->remote + filter->nslots, 0, 64 * sizeof(elt_t));

  // Update counters
  filter->nblocks += 1;
  filter->nslots += 64;
}

/**
 * Adapt a fingerprint at a particular location.
 * TODO: consider adding one bit at a time instead of adding the shortest diff ext?
 */
static void adapt_loc(ExAF *filter, size_t loc, uint64_t in_hash, uint64_t out_hash) {
  Ext new_ext;
  if (shortest_diff_ext(filter, in_hash, out_hash, &new_ext) == 0) {
    fprintf(stderr, "Hashes were identical!\n");
    return;
  }
  // Write encoding to the appropriate block
  Ext exts[64];
  decode_ext(get_ext_code(filter, loc/64), exts);
  exts[loc%64] = new_ext;
  uint64_t code;
  if (encode_ext(exts, &code) == -1) {
    // Encoding failed: rebuild
    memset(exts, 0, 64 * sizeof(Ext)); // clear exts
    exts[loc % 64] = new_ext;
    if (encode_ext(exts, &code) == -1) {
      fprintf(stderr, "Encoding failed after rebuild!\n");
      exts[loc % 64].len = 0;
      exts[loc % 64].bits = 0;
      code = 0;
    }
  }
  // Encoding succeeded: update code
  set_ext_code(filter, loc/64, code);
}

/**
 * Adapt on a query element that collided with a stored fingerprint at loc.
 *
 * Go through the rest of the run and fix any other remaining collisions.
 */
static void adapt(ExAF *filter, elt_t query, int loc, size_t quot, rem_t rem, uint64_t hash, Ext exts[64]) {
  assert(quot <= loc && loc < filter->nslots);
  // Make sure the query elt isn't mapped to an earlier index in the sequence
  for (int i=loc; i>=(int)quot && (i == loc || !get_runend(filter, i)); i--) {
    if (filter->remote[i] == query) {
      return;
    }
  }
  // Adapt on all collisions in the run
  for (int i=loc; i>=(int)quot && (i == loc || !get_runend(filter, i)); i--) {
    // Re-decode if at a new block
    if (i != loc && i % 64 == 63) {
      decode_ext(get_ext_code(filter, i/64), exts);
    }
    // Check collision
    Ext ext = exts[i % 64];
    if (remainder(filter, i) == rem && ext_matches_hash(filter, &ext, hash)) {
      // Adapt on hash collision
      uint64_t in_hash = exaf_hash(filter, filter->remote[i]);
      adapt_loc(filter, i, in_hash, hash);
      // TODO: reuse `exts`
    }
  }
}

/* ExAF */

void exaf_init(ExAF *filter, size_t n, int seed) {
  filter->seed = seed;
  filter->nelts = 0;
  filter->nblocks = max(1, nearest_pow_of_2(n)/64);
  filter->nslots = filter->nblocks * 64;
  filter->q = (size_t)log2((double)filter->nslots); // nslots = 2^q
  filter->r = REM_SIZE;
  filter->p = filter->q + filter->r;
  filter->blocks = calloc(filter->nblocks, sizeof(ExAFBlock));
  filter->remote = calloc(filter->nslots, sizeof(elt_t));
}

void exaf_destroy(ExAF* filter) {
  free(filter->blocks);
  free(filter->remote);
  free(filter);
}

void exaf_clear(ExAF* filter) {
  filter->nelts = 0;
  free(filter->blocks);
  free(filter->remote);
  filter->blocks = calloc(filter->nblocks, sizeof(ExAFBlock));
  filter->remote = calloc(filter->nslots, sizeof(elt_t));
}

static void raw_insert(ExAF* filter, elt_t elt, uint64_t hash) {
  size_t quot = calc_quot(filter, hash);
  rem_t rem = calc_rem(filter, hash);
  filter->nelts++;

  // Find the appropriate runend
  int r = rank_select(filter, quot);
  switch (r) {
    case RANK_SELECT_EMPTY: {
      set_occupied(filter, quot);
      set_runend(filter, quot);
      remainder(filter, quot) = rem;
      filter->remote[quot] = elt;
      break;
    }
    case RANK_SELECT_OVERFLOW: {
      printf("ExAF failed to find runend (nslots=%lu, quot=(block=%lu, slot=%lu))\n",
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
       shift_exts(filter, r + 1, (int)u - 1);

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
      filter->remote[r+1] = elt;
    }
  }
}

static int raw_lookup(ExAF* filter, elt_t elt, uint64_t hash) {
  size_t quot = calc_quot(filter, hash);
  rem_t rem = calc_rem(filter, hash);

  if (get_occupied(filter, quot)) {
    int loc = rank_select(filter, quot);
    if (loc == RANK_SELECT_EMPTY || loc == RANK_SELECT_OVERFLOW) {
      return 0;
    }
    // Cache decoded extensions
    Ext decoded[64];
    int decoded_i = -1;
    do {
      if (remainder(filter, loc) == rem) {
        // Refresh cached code
        if (decoded_i != loc/64) {
          decoded_i = loc/64;
          uint64_t code = get_ext_code(filter, loc/64);
          decode_ext(code, decoded);
        }
        // Check if extensions match
        Ext ext = decoded[loc%64];
        if (ext_matches_hash(filter, &ext, hash)) {
          if (elt != filter->remote[loc]) {
            adapt(filter, elt, loc, quot, rem, hash, decoded);
          }
          return 1;
        }
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
int exaf_lookup(ExAF *filter, elt_t elt) {
  uint64_t hash = exaf_hash(filter, elt);
  return raw_lookup(filter, elt, hash);
}

void exaf_insert(ExAF *filter, elt_t elt) {
  uint64_t hash = exaf_hash(filter, elt);
  raw_insert(filter, elt, hash);
}

double exaf_load(ExAF *filter) {
  return (double)filter->nelts/filter->nslots;
}

/* Printing */

void print_exaf_metadata(ExAF* filter) {
  printf("FILTER METADATA:\n");
  printf("  p=%ld, q=%ld, r=%ld\n",
         filter->p, filter->q, filter->r);
  printf("  nslots=%ld, nblocks=%ld, blocksize=%ld, nelts=%ld\n",
         filter->nslots, filter->nslots/64, sizeof(ExAFBlock), filter->nelts);
  printf("  seed=%d\n", filter->seed);
  printf("  load factor=%f\n", exaf_load(filter));
}

void print_exaf_block(ExAF* filter, size_t block_index) {
  assert(0 <= block_index && block_index < filter->nslots/64);
  ExAFBlock block = filter->blocks[block_index];
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
  printf("  extension code=0x%lx\n", get_ext_code(filter, block_index));
  printf("  extensions=\n");
  Ext exts [64];
  decode_ext(get_ext_code(filter, block_index), exts);
  for (int i=0; i<8; i++) {
    printf("   ");
    for (int j=0; j<8; j++) {
      Ext ext = exts[i*8 + j];
      if (ext.len == 0) {
        printf(" _");
      } else {
        printf(" %.*lu", ext.len, ext.bits);
      }
    }
    printf("\n");
  }
  printf("  remote=\n");
  for (int i=0; i<8; i++) {
    printf("   ");
    for (int j=0; j<8; j++) {
      printf(" 0x%-*lx", 8, filter->remote[block_index * 64 + i*8 + j]);
    }
    printf("\n");
  }
}

void print_exaf(ExAF* filter) {
  print_exaf_metadata(filter);
  for (int i=0; i<filter->nblocks; i++) {
    print_exaf_block(filter, i);
  }
}

// Tests
//#define TEST_EXAF 1
#ifdef TEST_EXAF

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

#define EXAF_SEED 32776517

ExAF *new_exaf(size_t n) {
  ExAF *filter = malloc(sizeof(ExAF));
  exaf_init(filter, n, EXAF_SEED);
  return filter;
}

void test_calc_ext() {
  printf("Testing %s...", __FUNCTION__);
  ExAF *filter = new_exaf(128);
  // q = 7, r = 8
  uint64_t hash = exaf_hash(filter, 10101);
  int start = (int)(filter->q + filter->r);
  uint64_t bits;
  for (int i=1; i<64-start; i++) {
    bits = calc_ext_bits(filter, hash, i);
    assert_eq(bits, (hash >> start) & ONES(i));
  }
  exaf_destroy(filter);
  printf("passed.\n");
}

void test_shortest_diff_ext() {
  printf("Testing %s...", __FUNCTION__);
  ExAF *filter = new_exaf(64);
  // q=6, r=8
  Ext ext;
  int out;

  out = shortest_diff_ext(filter, 0, 0, &ext);
  assert_eq(out, 0);
  for (int i=1; i<64 - filter->r - filter->q; i++) {
    // setup h1, h2 to only differ at pos i
    uint64_t h1 = 0;
    uint64_t h2 = 1ull << ((uint64_t)i + filter->p);
    // h1 in, h2 out
    out = shortest_diff_ext(filter, h1, h2, &ext);
    assert_eq(out, i+1);
    assert_eq(ext.len, i+1);
    assert_eq(ext.bits, 0);
    // h1 out, h2 in
    out = shortest_diff_ext(filter, h2, h1, &ext);
    assert_eq(out, i+1);
    assert_eq(ext.len, i+1);
    assert_eq(ext.bits, 1ull << i);
  }
  exaf_destroy(filter);
  printf("passed.\n");
}

void test_ext_matches_hash() {
  printf("Testing %s...", __FUNCTION__);
  ExAF *filter = new_exaf(64);
  // q=6, r=8
  Ext ext;
  ext.bits = 0;
  ext.len = 0;
  assert(ext_matches_hash(filter, &ext, 0b0 << filter->p));
  assert(ext_matches_hash(filter, &ext, 0b1 << filter->p));
  ext.bits = 0;
  ext.len = 1;
  assert(ext_matches_hash(filter, &ext, 0b0 << filter->p));
  assert(ext_matches_hash(filter, &ext, 0b10 << filter->p));
  assert(!ext_matches_hash(filter, &ext, 0b1 << filter->p));
  assert(!ext_matches_hash(filter, &ext, 0b01 << filter->p));
  ext.bits = 1;
  ext.len = 1;
  assert(!ext_matches_hash(filter, &ext, 0b0 << filter->p));
  assert(!ext_matches_hash(filter, &ext, 0b10 << filter->p));
  assert(ext_matches_hash(filter, &ext, 0b1 << filter->p));
  assert(ext_matches_hash(filter, &ext, 0b01 << filter->p));
  ext.bits = 0b11;
  ext.len = 2;
  assert(ext_matches_hash(filter, &ext, 0b1011 << filter->p));
  assert(!ext_matches_hash(filter, &ext, 0b1100 << filter->p));
  exaf_destroy(filter);
  printf("passed.\n");
}

void test_add_block() {
  printf("Testing %s...", __FUNCTION__);
  ExAF *filter = new_exaf(64 * 2);
  assert_eq(filter->nslots, 64 * 2);
  assert_eq(filter->nblocks, 2);
  add_block(filter);
  // Check metadata
  assert_eq(filter->nslots, 64 * 3);
  assert_eq(filter->nblocks, 3);
  // Check new block
  ExAFBlock b = filter->blocks[2];
  assert_eq(b.occupieds, 0);
  assert_eq(b.runends, 0);
  assert_eq(b.offset, 0);
  for (int i=0; i<64; i++) {
    assert_eq(b.remainders[i], 0);
  }
  // Check remote rep
  for (int i=0; i<64; i++) {
    assert_eq(filter->remote[i + 128], 0);
  }
  exaf_destroy(filter);
  printf("passed.\n");
}

/// Check that adding a block doesn't overwrite existing data
void test_add_block_no_clobber() {
  printf("Testing %s...", __FUNCTION__);
  ExAF *filter = new_exaf(64 * 2);
  // Setup
  for (int i=0; i<filter->nslots; i++) {
    set_occupied(filter, i);
    set_runend(filter, i);
    remainder(filter, i) = i%16;
    filter->remote[i] = i;
  }
  add_block(filter);
  // Check that data in first 2 blocks is preserved
  for (int i=0; i<128; i++) {
    assert(get_occupied(filter, i));
    assert(get_runend(filter, i));
    assert_eq(remainder(filter, i), i%16);
    assert_eq(filter->remote[i], i);
  }
  // Check that 3rd block is empty
  for (int i=128; i<filter->nslots; i++) {
    assert(!get_occupied(filter, i));
    assert(!get_runend(filter, i));
    assert_eq(remainder(filter, i), 0);
    assert_eq(filter->remote[i], 0);
  }
  // Check filter metadata
  assert_eq(filter->nslots, 192);
  assert_eq(filter->nblocks, 3);

  exaf_destroy(filter);
  printf("passed.\n");
}

void test_adapt_loc_1() {
  printf("Testing %s...", __FUNCTION__);
  ExAF *filter = new_exaf(128);
  // q=7, r=8
  uint64_t code;
  Ext exts[64];
  code = get_ext_code(filter, 0);
  assert_eq(code, 0);
  decode_ext(get_ext_code(filter, 0), exts);
  for (int i=0; i<64; i++) {
    assert_eq(exts[i].len, 0);
    assert_eq(exts[i].bits, 0);
  }
  adapt_loc(filter, 0, 0, 1 << filter->p);
  decode_ext(get_ext_code(filter, 0), exts);
  assert_eq(exts[0].bits, 0);
  assert_eq(exts[0].len, 1);
  for (int i=1; i<64; i++) {
    assert_eq(exts[i].len, 0);
    assert_eq(exts[i].bits, 0);
  }
  exaf_destroy(filter);
  printf("passed.\n");
}

void test_adapt_loc_2() {
  printf("Testing %s...", __FUNCTION__);
  ExAF *filter = new_exaf(128);
  // q=7, r=8
  Ext exts[64];

  // Push encoding to limit (8 1-bit exts)
  for (int i=0; i<8; i++) {
    adapt_loc(filter, i, 0, 1 << filter->p);
  }
  // Check state
  decode_ext(get_ext_code(filter, 0), exts);
  for (int i=0; i<8; i++) {
    assert_eq(exts[i].len, 1);
    assert_eq(exts[i].bits, 0);
  }
  for (int i=8; i<64; i++) {
    assert_eq(exts[i].len, 0);
    assert_eq(exts[i].bits, 0);
  }

  // Adapt and force encoding to rebuild
  adapt_loc(filter, 8, 0, 0b1000 << filter->p);
  // Check state
  decode_ext(get_ext_code(filter, 0), exts);
  for (int i=0; i<64; i++) {
    if (i == 8) {
      assert_eq(exts[i].len, 4);
      assert_eq(exts[i].bits, 0);
    } else {
      assert_eq(exts[i].len, 0);
      assert_eq(exts[i].bits, 0);
    }
  }

  exaf_destroy(filter);
  printf("passed.\n");
}

void test_adapt_1() {
  printf("Testing %s...", __FUNCTION__);
  ExAF *filter = new_exaf(64);
  Ext exts[64];
  for (int i=0; i<64; i++) {
    exts[i].len = 0;
    exts[i].bits = 0;
  }
  // Try adapting on one elt at start of run, check that elt's ext is updated
  uint64_t in_hash = exaf_hash(filter, 0);
  uint64_t out_hash = in_hash ^ (1 << filter->p); // flip first bit after quot, rem in in_hash
  adapt(filter, 1, 0, 0, 0, out_hash, exts);
  uint64_t code = get_ext_code(filter, 0);
  decode_ext(code, exts);
  assert_eq(exts[0].len, 1);
  assert_eq(exts[0].bits, !!GET(in_hash, filter->q + filter->r));
  for (int i=1; i<filter->nslots; i++) {
    test_assert_eq(exts[i].len, 0, "i=%d", i);
    test_assert_eq(exts[i].bits, 0, "i=%d", i);
  }
  exaf_destroy(filter);
  printf("passed.\n");
}

void test_raw_lookup_1() {
  printf("Testing %s...", __FUNCTION__);
  ExAF *filter = new_exaf(128);

  printf(" [Unimplemented] ");

  exaf_destroy(filter);
  printf("passed.\n");
}

void test_raw_insert_1() {
  printf("Testing %s...", __FUNCTION__);
  ExAF *filter = new_exaf(128);

  printf(" [Unimplemented] ");

  exaf_destroy(filter);
  printf("passed.\n");
}

void test_shift_remote_elts() {
  printf("Testing %s...", __FUNCTION__);
  ExAF *filter = new_exaf(128);
  for (int i=0; i<filter->nslots; i++) {
    assert_eq(filter->remote[i], 0);
  }
  for (int i=0; i<filter->nslots; i++) {
    filter->remote[i] = i;
  }
  // Shift elts in [32, 64+32] to [33, 64+33]
  shift_remote_elts(filter, 32, 64+32);
  for (int i=0; i<=31; i++) {
    assert_eq(filter->remote[i], i);
  }
  assert_eq(filter->remote[32], 0);
  for (int i=33; i<=64+33; i++) {
    assert_eq(filter->remote[i], i-1);
  }
  for (int i=64+34; i<filter->nslots; i++) {
    assert_eq(filter->remote[i], i);
  }
  exaf_destroy(filter);
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
  ExAF* filter = new_exaf(s);

  // Generate query set
  srandom(EXAF_SEED);
  int nset = (int)(1.5*(double)s);
  Setnode* set = calloc(nset, sizeof(set[0]));
  char str[64];
  for (int i=0; i<s; i++) {
    elt_t elt = random() % a;
    sprintf(str, "%lu", elt);
    set_insert(str, (int)strlen(str), 0, set, nset);
    exaf_insert(filter, elt);
    assert(set_lookup(str, (int)strlen(str), set, nset));
    assert(exaf_lookup(filter, elt));
  }
  // Query [0, a] and ensure that all items in the set return true
  int fps = 0;
  int fns = 0;
  for (int i=0; i<a; i++) {
    elt_t elt = i;
    sprintf(str, "%lu", elt);
    int in_set = set_lookup(str, (int)strlen(str), set, nset);
    int in_exaf = exaf_lookup(filter, elt);
    if (in_set) {
      if (!in_exaf) {
        fns++;
        uint64_t hash = exaf_hash(filter, elt);
        size_t quot = calc_quot(filter, hash);
        rem_t rem = calc_rem(filter, hash);
        printf("False negative: set contains %lu (0x%lx), but filter doesn't:"
               " quot=%lu (block_i=%lu, slot_i=%lu), rem=0x%hhx\n",
               elt, elt, quot, quot/64, quot%64, rem);
        print_exaf_metadata(filter);
        print_exaf_block(filter, quot/64);
        exit(1);
      }
    } else {
      fps += in_exaf;
    }
  }
  printf("passed. ");
  printf("FPs: %d (%f%%), FNs: %d (%f%%)\n",
         fps, (double)fps/(double)a * 100, fns, (double)fns/(double)a * 100);
  print_exaf_metadata(filter);
  exaf_destroy(filter);
}

void test_insert_and_query_w_repeats() {
  printf("Testing %s...\n", __FUNCTION__);
  int nslots = 1 << 14;
  double load = 0.95;
  double a_s = 100;
  int queries_per_elt = 10;

  int s = (int)((double)nearest_pow_of_2(nslots) * load);
  int a = (int)((double)s * a_s);
  int n_queries = a * queries_per_elt;

  int fps = 0;  // false positives
  int rfps = 0; // repeated false positives
  int fns = 0;  // false negatives
  int tot_queries = n_queries * queries_per_elt;

  ExAF *filter = new_exaf(s);
  int nset = (int)(s * 1.5);
  Setnode *set = calloc(nset, sizeof(set[0]));

  srandom(EXAF_SEED);
  char str[64];
  int len;
  fprintf(stderr, "Initializing membership set and filter...\n");
  for (int i=0; i<s; i++) {
    elt_t elt = random();
    sprintf(str, "%lu", elt);
    len = (int)strlen(str);
    set_insert(str, len, 0, set, nset);
    exaf_insert(filter, elt);
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
    int in_filter = exaf_lookup(filter, elt);
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
    }
  }
  printf("Test results:\n");
  printf("FPs: %d (%f%%), RFPs: %d (%f%%)\n",
         fps, (double)fps/tot_queries, rfps, (double)rfps/tot_queries * 100);
  printf("FNs: %d (%f%%)\n", fns, (double)fns/tot_queries * 100);

  exaf_destroy(filter);
  printf("Done testing %s.\n", __FUNCTION__);
}

// a, b+1 in the same block
void test_shift_exts_single_block() {
  printf("Testing %s...", __FUNCTION__);
  // Setup
  ExAF *filter = new_exaf(64 * 3);
  Ext exts[64];
  for (int i=0; i<64; i++) {
    exts[i].bits = i % 8 == 0;
    exts[i].len = i % 8 == 0;
  }
  uint64_t code;
  encode_ext(exts, &code);
  set_ext_code(filter, 0, code);

  // Shift exts in [0, 62] -> [1, 63] and check
  for (int j=1; j <= 64; j++) {
    shift_exts(filter, 0, 62);
    decode_ext(get_ext_code(filter, 0), exts);
    for (int i=0; i<64; i++) {
      test_assert_eq(exts[i].bits,
                     (i - j) % 8 == 0 && i > j-1,
                     "j=%d, i=%d", j, i);
      test_assert_eq(exts[i].len,
                     (i - j) % 8 == 0 && i > j-1,
                     "j=%d, i=%d", j, i);
    }
  }
  exaf_destroy(filter);
  printf("passed.\n");
}

void test_swap_exts() {
  printf("Testing %s...", __FUNCTION__);
  Ext* xs = malloc(64 * sizeof(Ext));
  Ext* ys = malloc(64 * sizeof(Ext));
  for (int i=0; i<64; i++) {
    xs[i].bits = i;
    xs[i].len = i;
    ys[i].bits = 64 + i;
    ys[i].len = 64 + i;
  }
  swap_ptrs(&xs, &ys);
  for (int i=0; i<64; i++) {
    assert_eq(ys[i].bits, i);
    assert_eq(ys[i].len, i);
    assert_eq(xs[i].bits, 64+i);
    assert_eq(xs[i].len, 64+i);
  }
  printf("passed.\n");
}

ExAF* ext_setup() {
  ExAF* filter = new_exaf(64 * 4);
  Ext exts[64];
  uint64_t code;
  for (int i=0; i<filter->nblocks; i++) {
    for (int j=0; j<64; j++) {
      exts[j].bits = j % 8 == 0;
      exts[j].len = j % 8 == 0;
    }
    encode_ext(exts, &code);
    set_ext_code(filter, i, code);
  }
  return filter;
}

// a, b+1 in different blocks
void test_shift_exts_multi_block() {
  printf("Testing %s...", __FUNCTION__);
  Ext exts[64];
  ExAF *filter = ext_setup();
  // (1) Shift exts in [0, 127] -> [1,128]
  shift_exts(filter, 0, 127);
  for (int i=0; i<filter->nblocks; i++) {
    decode_ext(get_ext_code(filter, i), exts);
    for (int j=0; j<64; j++) {
      test_assert_eq(exts[j].bits,
                     (i < 2) ? ((j-1)%8 == 0) : (i == 2 && j == 0 ? 0 : j%8 == 0),
                     "i=%d, j=%d", i, j);
      test_assert_eq(exts[j].len,
                     (i < 2) ? ((j-1)%8 == 0) : (i == 2 && j == 0 ? 0 : j%8 == 0),
                     "i=%d, j=%d", i, j);
    }
  }
  exaf_destroy(filter);
  filter = ext_setup();
  // (2) Shift exts in [32, 64+32] -> [33, 64+33]
  shift_exts(filter, 32, 64+32);
  for (int i=0; i<filter->nslots; i++) {
    if (i%64 == 0) {
      decode_ext(get_ext_code(filter, i/64), exts);
    }
    test_assert_eq(exts[i%64].bits,
                   i < 32 ? i%8 == 0 :
                   (i == 32 ? 0 :
                    (i <= 64 + 33) ? (i-1)%8 == 0 : i%8 == 0),
                   "i=%d", i);
    test_assert_eq(exts[i%64].len,
                   i < 32 ? i%8 == 0 :
                   (i == 32 ? 0 :
                    (i <= 64 + 33) ? (i-1)%8 == 0 : i%8 == 0),
                   "i=%d", i);
  }
  exaf_destroy(filter);
  printf("passed.\n");
}

void test_template() {
  printf("Testing %s...", __FUNCTION__);
  ExAF *filter = new_exaf(64 * 3);

  printf("[unimplemented]");

  exaf_destroy(filter);
  printf("passed.\n");
}

int main() {
  test_calc_ext();
  test_shortest_diff_ext();
  test_ext_matches_hash();
  test_add_block();
  test_add_block_no_clobber();
  test_adapt_loc_1();
  test_adapt_loc_2();
  test_adapt_1();
  test_raw_lookup_1();
  test_raw_insert_1();
  test_shift_remote_elts();
  test_shift_exts_single_block();
  test_shift_exts_multi_block();
  test_swap_exts();
  test_insert_and_query();
  test_insert_and_query_w_repeats();
}
#endif // TEST_EXAF