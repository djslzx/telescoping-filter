/* 
   Macros for use in other files
*/
#ifndef MACROS_H
#define MACROS_H

#ifdef __cplusplus
extern "C" {
#endif

/** Set bit at i, indexing from LSB */
#define ONE(i) (1ULL << (i))

/** Get i ones from bit 0 to i-1, equivalent to MASK_CLOSED(0, i-1) */
#define ONES(i) (ONE(i)-1)

/** Bits in the closed interval [a,b] set, all others unset; a <= b */
#define MASK_CLOSED(a,b) (ONES((b)-(a)+1) << (a))

/** Bits in the half-open interval [a,b) set, all others unset; a <= b */
#define MASK_HALF_OPEN(a,b) (ONES((b)-(a)) << (a))

/* Miscellaneous */
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))

/* 
   Access metadata 
   - occupieds, runends index from *LSB*
   - GET returns 0 or a nonzero value, not 0 or 1,
     but this works just fine for conditionals
*/
#define GET(bitarr, i) (ONE(i) & (bitarr))
#define SET(bitarr, i) ((bitarr) |= ONE(i))
#define UNSET(bitarr, i) ((bitarr) &= ~ONE(i))

// Get the block containing absolute index x
#define block_containing(filter, x) ((filter)->blocks[(x)/64])

/* Occupieds (i is an absolute index) */
#define get_occupied(filter, i)                         \
  (GET((filter)->blocks[(i)/64].occupieds, (i)%64))
#define set_occupied_to(filter, i, x)                               \
  ((x) ?                                                            \
   SET((filter)->blocks[(i)/64].occupieds, (i)%64) :                \
   UNSET((filter)->blocks[(i)/64].occupieds, (i)%64))
#define set_occupied(filter, i) SET((filter)->blocks[(i)/64].occupieds, (i)%64)
#define unset_occupied(filter, i) UNSET((filter)->blocks[(i)/64].occupieds, (i)%64)

/* Runends (i is an absolute index) */
#define get_runend(filter, i)                           \
  (GET((filter)->blocks[(i)/64].runends, (i)%64))
#define set_runend_to(filter, i, x)                                 \
  ((x) ?                                                            \
   SET((filter)->blocks[(i)/64].runends, (i)%64) :                  \
   UNSET((filter)->blocks[(i)/64].runends, (i)%64))
#define set_runend(filter, i) SET((filter)->blocks[(i)/64].runends, (i)%64)
#define unset_runend(filter, i) UNSET((filter)->blocks[(i)/64].runends, (i)%64)

/* Remainders */
/* Shorthand to get i-th remainder */
#define remainder(filter, i) ((filter)->blocks[(i)/64].remainders[(i)%64])

// Round v to nearest power of 2
// Pre: v >= 0
// https://graphics.stanford.edu/~seander/bithacks.html
static inline size_t nearest_pow_of_2(size_t v) {
  v--;
  v |= v >> 1UL;
  v |= v >> 2UL;
  v |= v >> 4UL;
  v |= v >> 8UL;
  v |= v >> 16UL;
  v |= v >> 32UL;
  v++;
  return v;
}

#ifdef __cplusplus
}
#endif

#endif
