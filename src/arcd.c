#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "arcd.h"
#include "taf.h"
#include "macros.h"

#define HIGH (~0ull >> (64 - SEL_CODE_LEN))

int encode_ext(const Ext exts[64], uint64_t* code) {
  uint64_t low = 0;
  uint64_t high = HIGH;

  for (int i=0; i<64; i++) {
    Ext ext = exts[i];
    uint64_t range = high - low;
    // Multiply range by ~0.90624 (Pr[ext is empty])
    uint64_t gap = (range >> 1) + (range >> 2) + (range >> 3) + (range >> 5);
    /* printf("range: %f\n", ((double)range)/HIGH); */
    /* printf("big gap: %f\n", ((double)gap)/HIGH); */
    if (ext.len == 0) {
      // If extension is empty, lower top of range
      high = low + gap;
      /* gap = (range >> 5) + (range >> 6); */
      /* /\* gap = range - ((range >> 1) + (range >> 2) + (range >> 3) + (range >> 5)); *\/ */
      /* printf("small gap: %f\n", ((double)gap)/HIGH); */

    } else {
      // If extension is nonempty, raise bottom of range
      low = low + gap;
      // Set gap to range * ~0.04687
      gap = (range >> 5) + (range >> 6);
      /* gap = range - ((range >> 1) + (range >> 2) + (range >> 3) + (range >> 5)); */
      /* printf("small gap: %f\n", ((double)gap)/HIGH); */
      /* exit(1); */
       /* gap = high - low; */
      // Account for probability of extension length:
      // extension length k>0 has probability 2^{-k}
      for (int j=1; j<ext.len; j++) {
        low += gap;
        gap >>= 1;
      }
      // Account for probability of a particular extension of length k:
      // all equally likely -> 1/(2^k)
      gap >>= ext.len; // Divide gap by 2^k
      low += ext.bits * gap; // Take bits-th 1/(2^k)-bit piece
      high = low + gap;
    }
    if (high - low < 2) {
      return -1;
    }
  }
  *code = low;
  return 0;
}

void decode_ext(uint64_t code, Ext exts[64]) {
  uint64_t low = 0;
  uint64_t high = HIGH;
  for (int i=0; i<64; i++) {
    uint64_t range = high - low;
    // Multiply range by ~0.90624 (Pr[ext is empty])
    uint64_t gap = (range >> 1) + (range >> 2) + (range >> 3) + (range >> 5);
    if (low + gap > code) {
      high = low + gap;
      exts[i].len = 0;
      exts[i].bits = 0;
    } else {
      low = low + gap;
      // Multiply range by ~0.04687
      gap = (range >> 5) + (range >> 6);
      /* gap = range - ((range >> 1) + (range >> 2) + (range >> 3) + (range >> 5)); */
      /* gap = range * 0.04701; */

      // Compute k, the length of the extension, by
      // iteratively shrinking the gap in proportion to
      // the probability of each length k=1, 2, ...
      int len = 1;
      while (low + gap <= code) {
        low += gap;
        gap >>= 1;
        len += 1;
      }
      // Get the bits given k, the length of the extension,
      // by dividing the interval into 2^k even pieces and
      // determining which piece the input belongs to
      gap >>= len;
      uint64_t bits = (code - low)/gap;
      low += bits * gap;
      high = low + gap;

      exts[i].bits = bits;
      exts[i].len = len;
    }
  }
}

int encode_sel(const int sels[64], uint64_t *code) {
  uint64_t low = 0;
  uint64_t high = HIGH;

  for(int i=0; i < 64; i++){
    int letter = sels[i];
    if (letter > MAX_SELECTOR) letter %= MAX_SELECTOR;

    uint64_t range = high - low;
    //takes advantage of "fall through" switch behavior
    switch (letter) {
      default:
        return 0;
      case 6:
        low += (range >> 19) + (range >> 20) + (range >> 23);
      case 5:
        low += (range >> 14) + (range >> 16);
      case 4:
        low += (range >> 10) + (range >> 11);
      case 3:
        low += (range >> 6) + (range >> 8);
      case 2:
        low += (range >> 3) + (range >> 4) + (range >> 7) + (range >> 9);
      case 1:
        low += (range >> 1) + (range >> 2) + (range >> 5);
      case 0: ;
    }

    switch (letter) {
      default:
        //high  = low + (range >> 22);
        return 0;
      case 6:
        high  = low + (range >> 24) + (range >> 25) + (range >> 26);
        break;
      case 5:
        high  = low + (range >> 19) + (range >> 20) + (range >> 23);
        break;
      case 4:
        high  = low + (range >> 14) + (range >> 16);
        break;
      case 3:
        high  = low + (range >> 10) + (range >> 11);
        break;
      case 2:
        high  = low + (range >> 6) + (range >> 8);
        break;
      case 1:
        high  = low + (range >> 3) + (range >> 4) + (range >> 7) + (range >> 9);
        break;
      case 0:
        high  = low + (range >> 1) + (range >> 2) + (range >> 5);
        break;
    }
    // Check if ran out of bits
    if(high - low < 2) return -1;
  }
  *code = low;
  return 0;
}

void decode_sel(uint64_t code, int out[64]) {
  uint64_t low = 0;
  uint64_t high = HIGH;

  for (int i=0; i < 64; i++){
    uint64_t range = high - low;
    uint64_t gap = (range >> 1) + (range >> 2) + (range >> 5);

    if(low + gap > code) {
      high = low + gap;
      out[i] = 0;
      continue;
    }
    low += gap;
    gap = (range >> 3) + (range >> 4) + (range >> 7) + (range >> 9);
    if(low + gap > code) {
      high = low + gap;
      out[i] = 1;
      continue;
    }
    low += gap;
    gap = (range >> 6) + (range >> 8);
    if(low + gap > code) {
      high = low + gap;
      out[i] = 2;
      continue;
    }
    low += gap;
    gap = (range >> 10) + (range >> 11);
    if(low + gap > code) {
      high = low + gap;
      out[i] = 3;
      continue;
    }
    low += gap;
    gap = (range >> 14) + (range >> 16);
    if(low + gap > code) {
      high = low + gap;
      out[i] = 4;
      continue;
    }
    low += gap;
    gap = (range >> 19) + (range >> 20) + (range >> 23);
    if(low + gap > code) {
      high = low + gap;
      out[i] = 5;
      continue;
    }
    low += gap;
    gap = (range >> 24) + (range >> 25) + (range >> 26);
    if(low + gap > code) {
      high = low + gap;
      out[i] = 6;
      continue;
    }
    out[i] = 7;
  }
}

/* Tests */
//#define TEST_ARCD 1
#ifdef TEST_ARCD

#define assert_eq(a, b) assert((a) == (b))

#define test_assert_eq(a, b, msg, ...)   \
  if ((a) != (b)) {                 \
    do {                            \
      fprintf(stderr, "Assertion failed: %s != %s: ", #a, #b); \
      fprintf(stderr, msg"\n", __VA_ARGS__); \
      assert_eq(a, b);              \
    } while (0);                    \
  }

#define panic(...) \
  do {             \
    fprintf(stderr, __VA_ARGS__);  \
    exit(1);       \
  } while (0)

/**
 * Convert a string into an ext
 */
void str_to_ext(char* str, Ext* ext) {
  int len = (int)strlen(str);
  if (len == 0) {
    ext->bits = 0;
  } else {
    ext->bits = strtol(str, NULL, 2);
  }
  ext->len = len;
}

/**
 * Convert an array of strs into an array of exts
 */
void strs_to_exts(char* strs[64], Ext exts[64]) {
  for (int i=0; i<64; i++) {
    str_to_ext(strs[i], &exts[i]);
  }
}

int ext_arr_eq(Ext a[64], Ext b[64]) {
  for (int i=0; i<64; i++) {
    if (a[i].len != b[i].len ||
        (a[i].bits & ONES(a[i].len)) != (b[i].bits & ONES(b[i].len))) {
      return 0;
    }
  }
  return 1;
}

void test_strs_to_exts() {
  printf("Testing %s...", __FUNCTION__);
  Ext exts[64];
  char* strs[64] = {
          "", "", "", "", "", "", "", "",
          "0", "0", "0", "0", "0", "0", "0", "0",
          "1", "1", "1", "1", "1", "1", "1", "1",
          "01", "01", "01", "01", "01", "01", "01", "01",
          "10", "10", "10", "10", "10", "10", "10", "10",
          "11", "11", "11", "11", "11", "11", "11", "11",
          "100", "100", "100", "100", "100", "100", "100", "100",
          "1000", "1000", "1000", "1000", "1000", "1000", "1000", "1000",
  };
  strs_to_exts(strs, exts);
  for (int i=0; i<64; i++) {
    switch (i / 8) {
      case 0:
        assert_eq(exts[i].bits, 0);
        assert_eq(exts[i].len, 0);
        break;
      case 1:
        assert_eq(exts[i].bits, 0);
        assert_eq(exts[i].len, 1);
        break;
      case 2:
        assert_eq(exts[i].bits, 1);
        assert_eq(exts[i].len, 1);
        break;
      case 3:
        assert_eq(exts[i].bits, 1);
        assert_eq(exts[i].len, 2);
        break;
      case 4:
        assert_eq(exts[i].bits, 0b10);
        assert_eq(exts[i].len, 2);
        break;
      case 5:
        assert_eq(exts[i].bits, 0b11);
        assert_eq(exts[i].len, 2);
        break;
      case 6:
        assert_eq(exts[i].bits, 0b100);
        assert_eq(exts[i].len, 3);
        break;
      case 7:
        assert_eq(exts[i].bits, 0b1000);
        assert_eq(exts[i].len, 4);
        break;
    }
  }
  printf("passed.\n");
}

/*
 * Check that exts = decode(encode(exts))
 */
void test_encode_decode_w_input(Ext exts[64]) {
  uint64_t code;
  Ext decoded[64];
  int out;
  out = encode_ext(exts, &code);
  if (out < 0) {
    panic("Encoding failed\n");
  } else {
    decode_ext(code, decoded);
    assert(ext_arr_eq(exts, decoded));
  }
}

void test_encode_decode_w_input_expect_fail(Ext exts[64]) {
  uint64_t code;
  int out;
  out = encode_ext(exts, &code);
  assert_eq(out, -1);
}

void test_encode_decode_empty() {
  printf("Testing %s...", __FUNCTION__);
  Ext exts[64];
  char *strs[64] = {
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
  };
  strs_to_exts(strs, exts);
  test_encode_decode_w_input(exts);
  printf("passed.\n");
}

void test_encode_decode_one() {
  printf("Testing %s...", __FUNCTION__);
  Ext exts[64];
  char *strs[64] = {
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "0",
  };
  strs_to_exts(strs, exts);
  test_encode_decode_w_input(exts);
  printf("passed.\n");
}

void test_encode_decode_few() {
  printf("Testing %s...", __FUNCTION__);
  Ext exts[64];
  char *strs[64] = {
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "000", "10", "1", "0",
  };
  strs_to_exts(strs, exts);
  test_encode_decode_w_input(exts);
  printf("passed.\n");
}

void test_encode_decode_many() {
  printf("Testing %s...", __FUNCTION__);
  Ext exts[64];
  char *strs[64] = {
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "0", "0", "0", "0", "0", "0", "0", "0",
  };
  strs_to_exts(strs, exts);
  test_encode_decode_w_input(exts);
  printf("passed.\n");
}

void test_encode_decode_many_rev() {
  printf("Testing %s...", __FUNCTION__);
  Ext exts[64];
  char *strs[64] = {
          "0", "0", "0", "0", "0", "0", "0", "0",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
  };
  strs_to_exts(strs, exts);
  test_encode_decode_w_input(exts);
  printf("passed.\n");
}

void test_encode_decode_long() {
  printf("Testing %s...", __FUNCTION__);
  Ext exts[64];
  char *strs[64] = {
          "1111111111111111111", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
  };
  strs_to_exts(strs, exts);
  test_encode_decode_w_input(exts);
  printf("passed.\n");
}

int will_overflow(Ext exts[64]) {
  uint64_t code;
  Ext decoded[64];
  if (encode_ext(exts, &code) < 0) {
    return 1;
  } else {
    decode_ext(code, decoded);
    assert(ext_arr_eq(exts, decoded));
    return 0;
  }
}

void print_exts(Ext exts[64]) {
  for (int k=0; k<64; k++) {
    printf(" %d", exts[k].len);
  }
  printf("\n");
}


void test_encode_decode_capacity() {
  printf("Testing %s...\n", __FUNCTION__);
  int limit = 20;
  Ext exts[64];
  // Try making ext arrays with different configs to stress test arcd
  // i: length of ext
  // j: number of i-length exts
  // k: iterate over exts
  for (int len=1; len < limit; len++) {
    for (int i=0; i < 64; i++) {
      exts[i].bits = 0;
      exts[i].len = 0;
    }
    int n;
    for (n=1; n < 64; n++) {
      for (int i=0; i < 64; i++) {
        exts[i].len = i < n ? len : 0;
      }
      if (will_overflow(exts)) {
        n--;
        break;
      }
    }
    printf("Can hold %d %d-length exts\n", n, len);
  }
}

void test_encode_decode_too_many() {
  printf("Testing %s...", __FUNCTION__);
  Ext exts[64];
  char *strs[64] = {
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "0",
          "0", "0", "0", "0", "0", "0", "0", "0",
  };
  strs_to_exts(strs, exts);
  test_encode_decode_w_input_expect_fail(exts);
  printf("passed.\n");
}

int main() {
  test_strs_to_exts();
  test_encode_decode_empty();
  test_encode_decode_one();
  test_encode_decode_few();
  test_encode_decode_many();
  test_encode_decode_many_rev();
  test_encode_decode_long();
  test_encode_decode_capacity();
  test_encode_decode_too_many();
}

#endif // TEST_ARCD
