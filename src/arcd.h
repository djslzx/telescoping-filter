#ifndef ARCD_H
#define ARCD_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include "ext.h"

/**
 * Encodes the input extensions, passing the code by reference.
 * @return 0 if encoding succeeds, -1 otherwise.
 */
int encode_ext(const Ext exts[64], uint64_t *code);

/**
 * Decodes the extension code, returning an array of Exts by reference.
  * @return the decoded Ext array
 */
void decode_ext(uint64_t code, Ext exts[64]);

/**
 * Encodes the input selectors, passing code out by ref.
 * @return 0 if encoding succeeds, -1 otherwise.
 */
int encode_sel(const int sels[64], uint64_t *code);

/**
 * Decodes the input, returning an array of selectors by reference.
 * @return the decoded int array
 */
void decode_sel(uint64_t code, int out[64]);

#ifdef __cplusplus
}
#endif

#endif // ARCD_H
